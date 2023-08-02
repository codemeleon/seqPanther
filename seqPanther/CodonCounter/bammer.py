#!/usr/bin/env python
from typing import Union
import pysam
from os import path


def check_sort_and_index_bam(bam_file: str, rid: str,
                             tmp_dir: str) -> Union[str, None]:
    """
    Checks bam file has given reference id.
    It also checks file is sorted or not. If not, sorts and indexes.
    """
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    header = samfile.header.as_dict()
    if rid not in samfile.references:
        print(f"Given reference {rid} not in given bam file {bam_file}")
        print("List of references")
        print(samfile.references)
        print(f"Ignoring {bam_file}")
        return

    if header['HD']['SO'] == 'coordinate':
        if not path.exists(bam_file + ".bai"):
            pysam.index(bam_file)
        return bam_file
    else:
        output_prefix = tmp_dir + "/" + path.split(bam_file)[1]
        pysam.sort("-o", output_prefix, bam_file)
        pysam.index(output_prefix)
        return output_prefix
