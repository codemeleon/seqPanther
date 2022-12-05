#!/usr/bin/env python

import pysam
from os import path


def check_sort_and_index_bam(bam_file, tmp_dir):
    """
    Checks file is sorted or not. If not, sorts and indexes.
    """
    stats = pysam.stats(bam_file).splitlines()
    is_sorted = False
    for line in stats:
        if line.startswith("SN"):
            line_split = line.split("\t")
            if line_split[1] == "is sorted":
                if line_split[2] == "1":
                    is_sorted = True
                break
    if is_sorted:
        if not path.exists(bam_file + ".bai"):
            pysam.index(bam_file)
        return bam_file
    else:
        output_prefix = tmp_dir + "/" + path.split(bam_file)[1]
        pysam.sort("-o", output_prefix, bam_file)
        pysam.index(output_prefix)
        return output_prefix
