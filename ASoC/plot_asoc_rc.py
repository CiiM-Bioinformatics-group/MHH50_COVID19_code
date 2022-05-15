#!/usr/bin/env python3
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2021
# Updated: Jun 01, 2021

import logging
import argparse
import itertools

from collections import OrderedDict

import matplotlib.pyplot as plt
from pysam import AlignmentFile

# Plasmablast were removed as nearly no ASoC reads can be retrieved.
BLACKLIST_CELLTYPE = ["Plasmablast"]
CONDITION_MAP = {"Severe": "Seve.", "Mild": "Mild", "Post": "Conv.",
                 "active": "Hospitalized", "post": "Convalescent"}


def is_good(read):
    """Wheter a read is a good one.

    We first check the filter used by GATK/ASEReadCounter, then others could
    affect the final read counts.
    """
    # ValidAlignmentStartReadFilter
    if read.reference_start is None: return False
    if read.reference_start <= 0: return False

    # ValidAlignmentEndReadFilter
    if read.reference_end is None : return False
    if read.reference_end - read.reference_start + 1 < 0: return False

    # HasReadGroupReadFilter 
    if all([x[0] != "RG" for x in read.tags]): return False

    # MatchingBasesAndQualsReadFilter
    if len(read.query_qualities) != len(read.query_sequence): return False

    # SeqIsStoredReadFilter
    if read.query_length == 0: return False

    # NotSecondaryAlignmentReadFilter
    if read.is_secondary: return False

    # MappedReadFilter
    if read.is_unmapped: return False

    # # is proper pair
    # if not read.is_proper_pair or not read.is_paired: return False

    # # passed QC
    # if read.is_qcfail: return False

    # # Mate should be mapped
    # if read.mate_is_unmapped: return False

    # # Should be represetative reads
    # if read.is_supplementary: return False

    # Remove duplicates
    if read.is_duplicate: return False

    return True


def shift_by_cigar(cur_pos, read):
    shift = cur_pos - 1 - read.reference_start
    for (oper, length) in read.cigartuples:
        if oper == 1: # I
            shift += length
        elif oper == 2: # D
            shift -= length
        elif oper == 3: #
            shift -= length

    return shift


class FieldMissingError(Exception):
    def __init__(self, msg):
        self.msg = msg
        super(FieldMissingError, self).__init__(msg)

    def __str__(self):
        return f"Missing required fields: {self.msg}!"


class HeaderMissingError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(HeaderMissingError, self).__init__(msg)

    def __str__(self):
        return f"Missing header line! Or you skipped it by {self.msg} lines."


class MultipleSNPIdError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(MultipleSNPIdError, self).__init__(msg)

    def __str__(self):
        return f"Multiple SNP ID given, which is not supported yet! {self.msg}"


class MissingSNPIdError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(MissingSNPIdError, self).__init__(msg)

    def __str__(self):
        return f"Multiple SNP ID given, which is not supported yet! {self.msg}"


def check_file_header(header_fields):
    """Check the header of a file.
    """
    for required_fields in ["contig", "position", "variantID", "refAllele",
                            "altAllele"]:
        if required_fields not in header_fields:
            raise FieldMissingError(required_fields)


def check_rcpool(rcpool):
    """Check whether the read counts pool is good.
    """
    snp_rsid, celltype, condition = [], [], []
    for _snp_rsid, _celltype, _condition, _, _ in rcpool.keys():
        snp_rsid.append(_snp_rsid)
        celltype.append(_celltype)
        condition.append(_condition)

    snp_rsid = list(set(snp_rsid))
    celltype = list(set(celltype))
    condition = list(set(condition))

    n_snp_rsid = len(snp_rsid)
    logging.info(f"Current SNP rsids: {', '.join(snp_rsid)}")

    if n_snp_rsid > 1:
        raise MultipleSNPIdError(", ".join(snp_rsid))
    elif n_snp_rsid < 1:
        raise MissingSNPIdError("Missing SNP rsid!")

    return sorted(celltype), sorted(condition)


def parse_region_file(fpath,
                      sep="\t",
                      skip=0,
                      header=1,
                      add_chr=True,
                      discard_celltype=BLACKLIST_CELLTYPE,
                      bampath="./"):
    """Parse position file."""
    snp_pos_map = {}
    if skip >= header:
        raise HeaderMissingError(skip)

    with open(fpath, "r") as fhandle:
        header_fields = None
        for idx, line in enumerate(fhandle):
            if idx < skip:
                continue

            if idx == header - 1:
                header_fields = line.strip().split(sep)
                check_file_header(header_fields)
                continue

            snp_info = line.strip().split(sep)
            chrom, pos, snpid, ref, alt, donor_id, celltype, condition = snp_info

            if celltype in discard_celltype:
                continue

            chrom = "chr" + chrom if add_chr else chrom
            bamfile = f'{bampath}/{donor_id}-{celltype}_asoc.bam'
            rec_key = (snpid, celltype, condition, chrom, ref, alt, int(pos))
            if rec_key not in snp_pos_map:
                snp_pos_map[rec_key] = [bamfile]
            else:
                snp_pos_map[rec_key].append(bamfile)

    return snp_pos_map


def fetch_flank_rc(htsfile, position, flank=100):
    """Fetch read counts of flank sequence for given ASoC SNP."""
    chrom, pos, *_ = position

    upbound = pos - flank - 1 if pos - flank - 1 > -1 else 0
    dwbound = pos + flank
    logging.debug(f"Position {pos} (1-based) into: {pos-1} (0-based)")
    logging.debug(f"Reads in window({flank*2+1}): {chrom}:{upbound}-{dwbound}")

    flank_readcounts = {}
    for read in htsfile.fetch(chrom, upbound, dwbound):
        read_start = max(read.reference_start, upbound)
        read_end = min(read.reference_end, dwbound)
        for cur_pos in range(read_start, read_end):
            if cur_pos in flank_readcounts:
                flank_readcounts[cur_pos] += 1
            else:
                flank_readcounts[cur_pos] = 1

    return flank_readcounts


def fetch_asoc_rc(htsfile, position):
    """Fetch read counts for each allele of an ASoC SNP.

    Args:
        htsfile
        position

    Returns:
        A dict including two keys: allele_1, allele_2.

    Raises:
        None

    Note:
        The parameter position is 1-based (VCF), but in this function, it's
        substracted one to fit the pysam 0-based scheme.

        The allelic read counts by ASEReadCounter took overlapping paired
        reads and deletions overlapping the position.
        Th CLI seqtech should be used to specify how to deal with overlpping
        paried reads. The --count-overlap-reads-handling options from
        ASEReadCounter is the key to handle the problem.
    """
    chrom, pos, ref, alt = position
    asoc_readcounts = OrderedDict()
    asoc_readcounts[ref] = {}
    asoc_readcounts[alt] = {}

    for read in htsfile.fetch(chrom, pos - 1, pos):
        if not is_good(read):
            continue

        read_start, read_end = read.reference_start, read.reference_end
        het_pos = shift_by_cigar(pos, read)
        het_base = read.query_sequence[het_pos]
        if het_base not in [ref, alt]:
            continue

        for cur_pos in range(read_start, read_end):
            if cur_pos in asoc_readcounts[het_base]:
                asoc_readcounts[het_base][cur_pos] += 1
            else:
                asoc_readcounts[het_base][cur_pos] = 1

    return asoc_readcounts


def fetch_rc(snpid, region_tab, flank_len=200, threads=4, min_parc=1,
             min_ttrc=5):
    """Fetch both ASoC SNP read counts and its flanking read counts.

    Returns:
        A `dict`.
        { 
           ("rs10086", "celltype", "condition", "chr1", 1002): {
               "asoc_rc": {"A": {1000: 9, 1001: 10, 1002: 10},
                           "G": {1000: 3, 1001:  4, 1002:  4}},
               "flank_rc": { 998: 10, 999: 12, 1000: 12, 1001: 14, 1002: 14,
                   1003: 14, 1003: 13 }
            }
        }
    """
    readcounts = {}
    for key, htsfile_pool in region_tab.items():
        if key in readcounts:
            raise KeyError(f"Duplicated key: {key}")

        if snpid == key[0]:
            _, celltype, condition, chrom, ref, alt, pos = key
            snp_pos = (chrom, pos, ref, alt)
            asoc_rc, flank_rc = {ref:{}, alt:{}}, {}
            for htsfile in htsfile_pool:
                mode = "rb" if htsfile.endswith(".bam") else "r"
                with AlignmentFile(htsfile, mode, threads=threads) as samhd:
                    _asoc_rc = fetch_asoc_rc(samhd, snp_pos)
                    asoc_rc = merge_readcounts(_asoc_rc, asoc_rc)

                    _flank_rc = fetch_flank_rc(samhd, snp_pos, flank=flank_len)
                    flank_rc = merge_readcounts(_flank_rc, flank_rc, "flank")

            try:
                bad_rcp = any([asoc_rc[ref][pos]<min_parc,
                               asoc_rc[alt][pos]<min_parc,
                               asoc_rc[ref][pos]+asoc_rc[alt][pos]<min_ttrc])
                if bad_rcp:
                    asoc_rc, flank_rc = {ref:{}, alt:{}}, {}
            except KeyError as e:
                logging.error(f"Bad read counts: {e}")

            # print(celltype, condition, asoc_rc)
            rc_key = (snpid, celltype, condition, chrom, pos)
            readcounts[rc_key] = {"asoc_rc": asoc_rc, "flank_rc": flank_rc}

    return readcounts

def merge_readcounts(rcp_left, rcp_right, ops="asoc"):
    if rcp_right is None:
        return rcp_left

    if ops == "asoc":
        for allele in rcp_left.keys():
            if len(rcp_left[allele]) == len(rcp_right[allele]) ==0:
                pass
            elif len(rcp_left[allele]) == 0:
                for pos in rcp_right[allele]:
                    rcp_left[allele][pos] = rcp_right[allele][pos]
            else:
                for pos in rcp_left[allele]:
                    rcp_left[allele][pos] += rcp_right[allele].get(pos, 0)
    else:
        for pos in rcp_left.keys():
            rcp_left[pos] += rcp_right.get(pos, 0)

        if ops != "flank":
            logging.warning(f"Unknown ops: {ops}, using flank instead.")

    return rcp_left


def plot_rc(rcpool,
            figheight=5,
            figwidth=6,
            cond_order=["active", "post"]):
            #cond_order=["Severe", "Mild", "Post"]):
    """Plot ASoC read counts."""
    celltype, condition = check_rcpool(rcpool)
    condition.sort(key=lambda x: cond_order.index(x))

    n_celltype, n_condition = len(celltype), len(condition)
    if n_celltype == n_condition == 1:
        logging.info("Only one panel will be plotted!")

    axes_idx = list(itertools.product(range(n_condition), range(n_celltype)))

    fig, axes = plt.subplots(ncols=n_celltype, nrows=n_condition, sharex=True,
                             sharey=True, squeeze=False,
                             subplot_kw={"frame_on": True})
    snpid = None
    lines, labels = None, None
    for (snpid, cur_cell, cur_cond, cur_chrom, cur_pos), val in rcpool.items():
        celltype_idx = celltype.index(cur_cell)
        condition_idx = condition.index(cur_cond)
        cur_axe = axes[condition_idx][celltype_idx]

        asoc_rc, flank_rc = val["asoc_rc"], val["flank_rc"]
        cur_axe.fill_between(flank_rc.keys(), 0, flank_rc.values(), lw=0,
                             alpha=0.2, color="#00BA38",
                             label="Gross read counts")

        ref, alt = list(asoc_rc.keys())
        color_dict = {ref: "#00BFC4", alt: "#C77CFF"}
        report_dict = OrderedDict()

        report_dict["chrom"] = cur_chrom
        report_dict["pos"] = cur_pos
        report_dict["snpid"] = snpid
        for allele, asoc_rcpool in asoc_rc.items():
            asoc_pos, asoc_rd = asoc_rcpool.keys(), asoc_rcpool.values()
        
            cur_axe.fill_between(asoc_pos, 0, asoc_rd, lw=0.5,
                                 alpha=0.7, color=color_dict[allele],
                                 label=f"{allele} allele read counts")

            if allele == ref:
                report_dict["ref"] = allele
                report_dict["refc"] = asoc_rcpool.get(cur_pos-1, 0)
            else:
                report_dict["alt"] = allele
                report_dict["altc"] = asoc_rcpool.get(cur_pos-1, 0)

        report_dict["cond"] = cur_cond
        report_dict["celltype"] = cur_cell
        logging.info(",".join([str(x) for x in report_dict.values()]))

        cur_axe.scatter(y=0, x=cur_pos, s=2, c="k", marker="|", label=snpid)

        if lines is None and labels is None:
            lines, labels = cur_axe.get_legend_handles_labels()

    for xi, yi in axes_idx:
        axes[xi][yi].set_xticks([])
        axes[xi][yi].tick_params(length=0.01, labelsize=5)

        cur_cond, cur_cell = condition[xi], celltype[yi]
        if yi == 0:
            axes[xi][yi].set_ylabel(CONDITION_MAP[cur_cond])

        if xi == len(condition) - 1:
            axes[xi][yi].set_xlabel(cur_cell)

        for spine_pos in ["top", "left", "right", "bottom"]:
            axes[xi][yi].spines[spine_pos].set_linewidth(0.0)

    fig.set_tight_layout(True)
    fig.subplots_adjust(wspace=-0.01, hspace=-0.01)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    fig.legend(lines, labels, ncol=1, fontsize="xx-small")

    return fig


def main():
    """The main entry of the script."""
    # CLI arguments
    parser = argparse.ArgumentParser(description="Plot ASoC read counts.")
    parser.add_argument("-r", "--region-file", required=True, metavar="FILE",
                        help="A tab-delimited file including SNP positions.")
    parser.add_argument("-i", "--snp-rsids", required=True, metavar="RSID",
                        nargs="+", help="SNP rsid to be ploted.")
    parser.add_argument("-b", "--bam-path", default="./", metavar="DIR",
                        help="The path to BAM files. Default: %(default)s")
    parser.add_argument("-f", "--flank-len", default=200, type=int,
                        metavar="INT",
                        help="The length of flank sequence in up/down-stream."
                        " of the ASoC SNP." " Default: %(default)s")
    parser.add_argument("-p", "--threads", default=1, type=int,
                        metavar="INT",
                        help="Number of threads used for (de)compression."
                        " Default: %(default)s")
    parser.add_argument("-d", "--discard-celltype", nargs="*", default=None,
                        metavar="STR",
                        help="Cell types will be discarded if region files"
                        " contains cell type. Default: %(default)s")
    parser.add_argument("-a", "--add-chr", action="store_true",
                        help="Add 'chr' the chromosome.")
    parser.add_argument("--seqtech", choices=["se", "pe"], default="pe",
                        metavar="[se, pe]",
                        help="Sequencing technology. Not in usage yet.")
    parser.add_argument("--fig-width", type=float, default=6,
                        metavar="FLOAT",
                        help="The width of the figure. Default: %(default)s")
    parser.add_argument("--fig-height", type=float, default=5,
                        metavar="FLOAT",
                        help="The height of the figure. Default: %(default)s")
    parser.add_argument("-o", "--outdir", default="./", metavar="PATH",
                        help="The output directory. Default: %(default)s")

    options = parser.parse_args()
    region_file = options.region_file
    flank_len = options.flank_len
    bam_path = options.bam_path
    threads = options.threads
    outdir = options.outdir
    snp_rsids = options.snp_rsids
    add_chr = options.add_chr
    discard_celltype = options.discard_celltype
    fig_width = options.fig_width
    fig_height = options.fig_height

    logging.basicConfig(format="{levelname: ^8}| {asctime} | {message}",
                        style="{",
                        datefmt="%Y%m%d,%H:%M:%S",
                        level=logging.INFO)

    if discard_celltype: BLACKLIST_CELLTYPE.extend(discard_celltype)
    logging.info(f"Blacklist: {', '.join(list(set(BLACKLIST_CELLTYPE)))}")

    region_tab = parse_region_file(region_file, sep=",", add_chr=add_chr, bampath=bam_path)
    for snpid in snp_rsids:
        rcpool = fetch_rc(snpid, region_tab, flank_len, threads)
        rcplot = plot_rc(rcpool, figheight=fig_height, figwidth=fig_width)

        rcplot.savefig(f"{outdir}/{snpid}_asoc_read_depth.pdf")


if __name__ == "__main__":
    main()

