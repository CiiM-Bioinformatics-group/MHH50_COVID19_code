#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 14, 2021
# Updated: May 18, 2021

#
## A script to split reads by cell type.
#

import os
import argparse

import pysam

def getargs():
    parser = argparse.ArgumentParser(description="Split reads according to cell type barcode.")
    parser.add_argument("-b", "--bam-files", dest="bam_files", required=True, nargs="+", metavar="FILE",
                        help="The BAM files from which load the reads. Required.")
    parser.add_argument("-c", "--barcode-file", dest="barcode_file", required=True, metavar="FILE",
                        help="A file including barcodes per cell type. Required.")
    parser.add_argument("-t", "--barcode-tag", dest="barcode_tag", default="CB", metavar="STRING",
                        help="The tag where the cell type barcode stored. Default: %(default)s")
    parser.add_argument("--celltype-col", dest="celltype_col", default="celltype", metavar="STRING",
                        help="The columne name of cell type in barcode file. Default: %(default)s")
    parser.add_argument("--barcode-col", dest="barcode_col", default="barcode", metavar="STRING",
                        help="The columne name of barcode in barcode file. Default: %(default)s")
    parser.add_argument("--donor-col", dest="donor_col", default="donor", metavar="STRING",
                        help="The columne name of donor in barcode file. Default: %(default)s")
    parser.add_argument("-n", "--threads", dest="threads", default=10, type=int, metavar="INTEGER",
                        help="Number of threads used for compression and decompression. Default: %(default)s")
    parser.add_argument("-o", "--output-dir", dest="output_dir", default="./", metavar="STRING",
                        help="The directory of the output file. Default: %(default)s")

    return parser.parse_args()


def make_barcode_table(barcode_file, barcode_col, celltype_col, donor_col, sep="\t"):
    """Prepare a barcode to cell type table."""
    barcode_file = os.path.abspath(barcode_file)
    barcode_table = {}
    with open(barcode_file, "r") as fh:
        barcode_idx = celltype_idx = donor_idx = 0
        for idx, line in enumerate(fh):
            line = line.strip()
            if idx == 0:
                header_line = line.split(sep)
                barcode_idx = header_line.index(barcode_col)
                celltype_idx = header_line.index(celltype_col)
                donor_idx = header_line.index(donor_col)
                continue

            if not (barcode_idx != celltype_idx != donor_idx):
                raise ValueError("Barcode, cell type, and donor column should not be the same!")

            record_line = line.split(sep)
            barcode = record_line[barcode_idx]
            celltype = record_line[celltype_idx]
            donor = record_line[donor_idx]

            if barcode in barcode_table and barcode_table[barcode] != celltype:
                raise KeyError("Duplicated barcode! Does the barcode table derive from singlelets?")

            barcode_table[barcode] = (donor, celltype)

    return barcode_table


def make_save_path(output_dir: str, *args):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    fbname = "-".join([a.replace(" ", "") for a in args]) + ".bam"
    return os.path.join(output_dir, fbname)


def split_reads(bam_files, barcode_tag, barcode_table, output_dir, threads):
    """Split reads according to barcode table."""

    threads = max(1, int(threads / 2))
    per_sam_pool = list(set(barcode_table.values()))
    fhpool = {}

    for bam_fp in bam_files:
        bam_fp = os.path.abspath(bam_fp)
        insamfile = pysam.AlignmentFile(bam_fp, "rb", threads=threads)

        for read in insamfile.fetch():
            if read.has_tag(barcode_tag):
                _barcode = read.get_tag(barcode_tag)

                if _barcode not in barcode_table:
                    continue

                _donor, _celltype = barcode_table[_barcode]
                fhkey = (_donor, _celltype)
                if fhkey not in fhpool:
                    optsam_fp = make_save_path(output_dir, _donor, _celltype)
                    fhpool[fhkey] = pysam.AlignmentFile(optsam_fp, "wb",
                                                        template=insamfile,
                                                        threads=threads)

                fhpool[fhkey].write(read)

        for fhkey in per_sam_pool:
            if fhpool[fhkey].is_open:
                fhpool[fhkey].close()

        if insamfile.is_open:
            insamfile.close()


def main():
    args = getargs()
    barcode_file = args.barcode_file
    barcode_tag = args.barcode_tag
    barcode_col = args.barcode_col
    celltype_col = args.celltype_col
    donor_col = args.donor_col
    bam_files = args.bam_files
    output_dir = args.output_dir
    threads = args.threads

    barcode_table = make_barcode_table(barcode_file, barcode_col, celltype_col, donor_col)
    _ = split_reads(bam_files, barcode_tag, barcode_table, output_dir, threads)


if __name__ == "__main__":
    main()

# vim: set nowrap number relativenumber expandtab cindent tw=500 ts=4:
