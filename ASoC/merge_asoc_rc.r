#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 27, 2021
# Updated: Jul 20, 2021

options(stringsAsFactors = F, error = traceback)

#' Parse commandline options
#'
#' A function to parse commandline options when using shebang(#!).
#'
getargs <- function() {
    p <- optparse::OptionParser(
        description = "A short R script to estimate allele-specific reads",
        usage = "estimate_asoc.r [options] [FILE1 FILE2 FILE3 ...]"
    )
    p <- optparse::add_option(p, c("-g", "--genotype-file"),
        metavar = "FILE",
        help = "A VCF file including genotypes."
    )
    p <- optparse::add_option(p, c("-m", "--meta-info-file"),
        metavar = "FILE",
        help = "The meta information for each patients. It requires a ID column name used as index."
    )
    p <- optparse::add_option(p, c("-M", "--id-mapping"),
        metavar = "KEY:VAL", default = "DonorID:genoID",
        help = "The idmapping syntax key:value, key is the old name, value is the new value. Default: %default"
    )
    p <- optparse::add_option(p, c("-a", "--asoc-readcounts"),
        metavar = "FILE", default = "readcounts.csv",
        help = "The ASoC estimation results file. Default: %default"
    )
    p <- optparse::add_option(p, c("-p", "--n-cpus"),
        metavar = "INTEGER", type = "integer", default = 4,
        help = "The number of parallel workers to be used. Default: %default"
    )
    p <- optparse::add_option(p, c("--min-sm-rc"),
        metavar = "INTEGER", type = "integer", default = 10,
        help = "The minimum gross read counts of REF and ALT alleles. Default: %default"
    )
    p <- optparse::add_option(p, c("--min-pa-rc"),
        metavar = "INTEGER", type = "integer", default = 2,
        help = "The minimum read counts of REF or ALT alleles. Default: %default"
    )

    argopts <- optparse::parse_args2(p)
    args <- argopts$args
    opts <- argopts$options

    if (!file.exists(opts$asoc_readcounts)) {
      if (is.null(args) || length(args) == 0) {
        stop("When no -a/--asoc-estimation is specified, positional options are required.")
      }

      if (is.null(opts$genotype_file)) {
        stop("When no -a/--asoc-estimation is specified, -g/--genotype is required.")
      }
    }

    return(argopts)
}


#' Read text files
read_txt <- function(fpath) {
    dtfm <- NULL

    if (!is.null(fpath)) {
        fpxt <- tools::file_ext(fpath)
        sep <- ifelse(fpxt == "csv", ",", "")
        dtfm <- read.table(fpath, header = T, sep = sep)
    }

    dtfm
}


#' Main
main <- function() {
    argopts <- getargs()
    ifps <- argopts$args
    opts <- argopts$options

    n_cpus <- opts$n_cpus
    min_sm_rc <- opts$min_sm_rc
    min_pa_rc <- opts$min_pa_rc
    genotype_file <- opts$genotype_file
    asoc_readcounts <- opts$asoc_readcounts

    meta_info_fp <- opts$meta_info_file
    id_mapping <- opts$id_mapping

    # If no options errors, load required libraries.
    library(data.table, verbose = F)
    library(multidplyr, verbose = F)
    library(tidyverse, verbose = F)

    # Some meta information
    metainfo <- fread(meta_info_fp)

    # ID mapping key -> value
    kv_pair <- unlist(str_split(id_mapping, ":"))
    key <- kv_pair[1]
    val <- kv_pair[2]

    # ID mapping file
    idmap <- metainfo %>%
        select(one_of(key, val)) %>%
        deframe()

    # Genotypes
    gntp_dtfm <- fread(genotype_file,
        tmpdir = dirname(genotype_file),
        nThread = n_cpus, data.table = F
    ) %>%
        select(!one_of(c(
            "#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"
        ))) %>%
        filter(!str_detect(ID, "chr")) %>%
        apply(2, FUN = function(x) {
            str_split(x, ":", simplify = T)[, 1]
        }) %>%
        as.data.frame()

    rownames(gntp_dtfm) <- gntp_dtfm[, "ID"]
    gntp_dtfm <- gntp_dtfm %>%
        dplyr::select(-ID) %>%
        as.data.frame()
    cat("Loaded", nrow(gntp_dtfm), "SNPs from the gentype files.\n")
    print(gntp_dtfm[1:5, 1:5])

    if (!file.exists(asoc_readcounts)) {
        # Initialize a cluster
        my_cluster <- new_cluster(n_cpus)
        cluster_library(my_cluster, "dplyr")

        # Load read counts table (from GATK/ASEReadCounter)
        rctabs <- NULL
        for (fp in unique(ifps)) {
            start_time <- Sys.time()
            fp_xt <- tools::file_ext(fp)
            sep <- ifelse(fp_xt == "csv", ",", "")

            # We load read counts results per cell type per smaple
            fp_nm <- str_replace(basename(fp), ".gatkAlleleReadCounts.csv", "")
            donor_id <- str_split(fp_nm, pattern = "-", simplify = T)[, 1]
            celltype <- str_split(fp_nm, pattern = "-", simplify = T)[, 2]

            cat(
                "Merging read counts from", celltype, "cell of", donor_id,
                "into the main table...\n"
            )

            .rctabs <- fread(fp, showProgress = F, nThread = n_cpus) %>%
                dplyr::filter((refCount >= min_pa_rc | altCount >= min_pa_rc) &
                    totalCount >= min_sm_rc) %>%
                (function(d) {
                    cat("Kept", nrow(d), "records.\n")
                    d
                })

            # Pick up information of SNPs in the read counts table
            use_snps <- .rctabs %>%
                dplyr::select(variantID) %>%
                unlist()
            gntp_dtfm_sub <- gntp_dtfm[use_snps, ]
            gntp_rsid <- rownames(gntp_dtfm_sub)

            if (length(gntp_rsid) >= 1e5) {
                rctabs <- .rctabs %>%
                    partition(cluster = my_cluster) %>%
                    mutate(
                        donor_id = donor_id,
                        celltype = celltype,
                        geno_id = mapply(function(i) idmap[i], donor_id),
                        genotype = mapply(
                            function(x, y, gid) {
                                ifelse(x %in% gntp_rsid, gntp_dtfm_sub[x, gid], NA)
                            },
                            variantID, donor_id, geno_id
                        )
                    ) %>%
                    filter(genotype %in% c("1|0", "0|1")) %>%
                    collect() %>%
                    (function(d) {
                        cat("Found heterozygous", nrow(d), "in total.\n")
                        d
                    }) %>%
                    rbind(rctabs)
            } else {
                rctabs <- .rctabs %>%
                    mutate(
                        donor_id = donor_id,
                        celltype = celltype,
                        geno_id = mapply(function(i) idmap[i], donor_id),
                        genotype = mapply(
                            function(x, gid) {
                                ifelse(x %in% gntp_rsid, gntp_dtfm_sub[x, gid], NA)
                            },
                            variantID, geno_id
                        )
                    ) %>%
                    filter(genotype %in% c("1|0", "0|1")) %>%
                    (function(d) {
                        cat("Found heterozygous", nrow(d), "in total.\n")
                        d
                    }) %>%
                    rbind(rctabs)
            }
            # Assign genotype to each SNP
            end_time <- Sys.time()
            interval <- end_time - start_time
            cat("Time eclapsed:", interval, units(interval), "\n\n")
        }
        # Write the estimation result to disk
        if (!is.null(rctabs)) {
            print(rctabs[1:5, 1:5])
            print(dim(rctabs))

            write.csv(rctabs, asoc_readcounts, row.names = F, quote = F)
        } else {
            warning("No allelic open chromatin SNP was found.")
        }
    } else {
        rctabs <- read.csv(asoc_readcounts)
    }

    candiate_snps <- rctabs %>%
        dplyr::group_by(contig, position) %>%
        dplyr::filter(sum(refCount < 3) < 1 & sum(altCount < 3) < 1) %>%
        dplyr::arrange(contig, position, .group = T) %>%
        dplyr::ungroup() %>%
        dplyr::select(variantID) %>%
        dplyr::distinct() %>%
        unlist(use.names = F)

    print(candiate_snps)
}

main()

# vim: set nowrap nu rnu tw=400 ts=4:
