#!/usr/bin/env python

import pandas as pd
"""Reads gff file buffer"""


def gff2tab(gff_buffer):
    """Convert gff file information to pandas data frame for CDS.

    :gff_buffer: gff file buffer
    :returns: pandas data frame for CDS
    """
    gff_data = gff_buffer.read()
    gff_data = gff_data.split("##FASTA")[0]
    gff_data = gff_data.split("\n")
    gff_data = [x for x in gff_data if (x == "" or x[0] == "#") is False]
    gff_data = [x.split("\t") for x in gff_data]
    gff_data = pd.DataFrame(
        gff_data,
        columns=[
            "seq_id",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    gff_data = gff_data.loc[gff_data["feature"] == "CDS"]
    gff_data["start"] = gff_data["start"].astype(int) - 1
    gff_data["end"] = gff_data["end"].astype(int) - 1
    return gff_data
