#!/usr/bin/env python

# Will fill the missing values in the table

import pysam


def missing(bamfile, ref, pos, offset, type="S"):
    # TODO: add information from user
    """Will resturn codon count for current position

    :bamfile: index bam file
    :ref: reference name
    :pos: 0-based index
    :offset: offset from pos
    :returns: file name, ref, position, offset, codon count

    """
    samfile = pysam.Samfile(bamfile, "rb")
    for pileupcol in samfile.pileup(ref, pos, pos+1, truncate=True):
        if pileupcol.pos != pos:
            continue
        for item in iterable:
            pass
        break
