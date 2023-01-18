#!/usr/bin/env python

import pysam
from os import path


def check_sort_and_index_bam(bam_file, tmp_dir):
    """
    Checks file is sorted or not. If not, sorts and indexes.
    """
    header = pysam.AlignmentFile(bam_file, 'rb').header.as_dict()
    if header['HD']['SO'] == 'coordinate':
        if not path.exists(bam_file + ".bai"):
            pysam.index(bam_file)
        return bam_file
    else:
        output_prefix = tmp_dir + "/" + path.split(bam_file)[1]
        pysam.sort("-o", output_prefix, bam_file)
        pysam.index(output_prefix)
        return output_prefix
