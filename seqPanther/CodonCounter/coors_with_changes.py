#!/usr/bin/env python
import pysam
import pandas as pd

import pyfaidx
from os import path

from .subs import sub_table
from .indel_frames import indel_frames


def changed_coordinates(params, bam):
    print(f"Analysing {bam}.")
    rid = params["rid"]
    start = params["start"]
    end = params["end"]
    endlen = params["endlen"]
    sequences = params["sequences"]
    ignore_orphans = params["ignore_orphans"]
    min_mapping_quality = params["min_mapping_quality"]
    min_base_quality = params["min_base_quality"]
    min_seq_depth = params["min_seq_depth"]
    alt_nuc_count = params["alt_nuc_count"]
    ignore_overlaps = params["ignore_overlaps"]

    samfile = pysam.AlignmentFile(bam, "rb")
    if rid not in samfile.references:
        print(f"Given reference {rid} not in given bam file {bam}")
        print("List of references")
        print(samfile.references)
        print(f"Ignoring {bam}")
        return
    iter = samfile.pileup(
        rid,
        start,
        end,
        ignore_orphans=ignore_orphans,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        ignore_overlaps=ignore_overlaps,
        max_depth=1000000,
    )
    coordinates_with_change = {}
    indel_pos_type_size = {"coor": [], "indel": [], "seq": []}
    depth = {"coor": [], "depth": []}
    for pileupcol in iter:
        if pileupcol.n < min_seq_depth:
            continue
        if (pileupcol.pos >= start) & (pileupcol.pos < end):
            # TODO: Include base call quality
            depth["coor"].append(pileupcol.pos)
            depth["depth"].append(pileupcol.n)
            bases = {}
            nuc_indel_count = 0

            for pread in pileupcol.pileups:

                if pread.indel:
                    nuc_indel_count += 1
                    indel_pos_type_size["coor"].append(pileupcol.pos)
                    indel_pos_type_size["indel"].append(pread.indel)
                    if pread.indel > 0:
                        indel_pos_type_size["seq"].append(
                            pread.alignment.
                            query_sequence[pread.query_position +
                                           1:pread.query_position + 1 +
                                           pread.indel])
                    else:
                        indel_pos_type_size["seq"].append("")

                if not pread.is_del and not pread.is_refskip:
                    if (pread.query_position < endlen
                            or (len(pread.alignment.query_sequence) -
                                pread.query_position + 1) < endlen):
                        continue
                    tbase = pread.alignment.query_sequence[
                        pread.query_position]
                    if tbase not in bases:
                        bases[tbase] = {
                            "nuc_count": 0,
                            "codon_count": {},
                        }
                    bases[pread.alignment.query_sequence[
                        pread.query_position]]["nuc_count"] += 1
            # NOTE: Deleting nucleotide which have low frequency

            nucs_to_delete = ""
            for nuc in bases.keys():
                if bases[nuc]["nuc_count"] < alt_nuc_count * pileupcol.n:
                    nucs_to_delete += nuc
            for nuc in nucs_to_delete:
                del bases[nuc]
            if set(bases) - set([sequences[pileupcol.pos].seq]):
                coordinates_with_change[pileupcol.pos] = {
                    "bases": bases,
                    "read_count": pileupcol.n,
                }
    depth = pd.DataFrame(depth)
    indel_pos_type_size = pd.DataFrame(indel_pos_type_size)
    indel_pos_type_size = (indel_pos_type_size.groupby(
        ["coor", "indel",
         "seq"]).size().reset_index().rename(columns={0: "indel_read_count"}))
    indel_pos_type_size = indel_pos_type_size.merge(depth, on="coor")
    indel_pos_type_size = indel_pos_type_size[
        indel_pos_type_size["indel_read_count"] > alt_nuc_count *
        indel_pos_type_size[  # TODO: Replace the alt nuc value with indel alt count
            "depth"]]
    indel_pos_type_size_full = indel_pos_type_size.copy()
    indel_pos_type_size = indel_pos_type_size[indel_pos_type_size["indel"] %
                                              3 == 0]

    return coordinates_with_change, indel_pos_type_size, indel_pos_type_size_full, depth


def coor_with_changes_run(params, bam):
    params["sequences"] = pyfaidx.Fasta(params["ref"])[params["rid"]]
    merged_table = None
    merged_table_nuc = None
    res = changed_coordinates(params, bam)
    print(res[0])
    subs_table = sub_table(res[0], bam, params)
    print(subs_table)
    indelframes = indel_frames(res[1], bam, params)
    merged_table = pd.concat([indelframes[0], indelframes[1], subs_table])
    flb = path.split(bam)[1].split(".")[0]
    res_indel = res[2]
    res_indel.loc[res_indel["indel"] < 0,
                  "seq"] = res_indel.loc[res_indel["indel"] < 0].apply(
                      lambda x: params["sequences"][x["coor"] + 1:x["coor"] + 1
                                                    - x["indel"]].seq,
                      axis=1,
                  )  # NOTE: refefence nucleotide for deletion events
    res_indel["sample"] = flb
    # NOTE: Done till here
    res_sub = res[0]
    res_table = {"pos": [], "read_count": [], "base_count": [], "base_pt": []}
    for pos in res_sub:
        res_table["pos"].append(pos)
        res_table["read_count"].append(res_sub[pos]["read_count"])
        base_count = []
        base_pt = []
        for base in res_sub[pos]["bases"]:
            base_count.append(
                f"{base}: {res_sub[pos]['bases'][base]['nuc_count']}")
            pt_val = "%.f" % (res_sub[pos]["bases"][base]["nuc_count"] *
                              100.0 / res_sub[pos]["read_count"])
            base_pt.append(f"{base}: {pt_val}")
        res_table["base_count"].append(",".join(base_count).replace(" ", ""))
        res_table["base_pt"].append(",".join(base_pt).replace(" ", ""))
    res_table = pd.DataFrame(res_table)
    res_table["sample"] = flb
    res_table["ref_base"] = res_table.apply(
        lambda x: params["sequences"][x["pos"]].seq, axis=1)
    res_indel["indel_read_pt"] = (res_indel["indel_read_count"] * 100.0 /
                                  res_indel["depth"])
    merged_table_nuc = [res_table, res_indel]
    return merged_table, merged_table_nuc, {flb: pd.DataFrame(res[-1])}
