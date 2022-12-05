#!/usr/bin/env python

import pandas as pd
from .codon_table import codon_table
from os import path
import pysam


def sub_table(coordinates_with_change, bam, params):
    sample = path.split(bam)[1]
    alt_nuc_count = params["alt_nuc_count"]
    sequences = params["sequences"]
    rid = params["rid"]
    ignore_orphans = params["ignore_orphans"]
    ignore_overlaps = params["ignore_overlaps"]
    min_base_quality = params["min_base_quality"]

    gff_data = params["gff_data"]
    min_seq_depth = params["min_seq_depth"]
    sample = path.split(bam)[1]
    keys = set(coordinates_with_change)
    samfile = pysam.AlignmentFile(bam, "rb")

    for (
            _,
            row,
    ) in gff_data.iterrows():
        selected_coordinates = keys & set(range(row["start"], row["end"]))
        if selected_coordinates:
            # TODO: Add count of codon. Ignore with ambigious nucleotide and gaps
            # TODO: Check if there is any 3 base indels
            for selected_coordinate in selected_coordinates:
                coordinates_with_change[selected_coordinate]["start"] = row[
                    "start"]
                coordinates_with_change[selected_coordinate]["end"] = row[
                    "end"]
                coordinates_with_change[selected_coordinate]["strand"] = row[
                    "strand"]
                shift = (selected_coordinate - row["start"]) % 3
                iter = samfile.pileup(
                    rid,
                    selected_coordinate,
                    selected_coordinate + 1,
                    ignore_orphans=ignore_orphans,
                    min_mapping_quality=min_base_quality,
                    min_base_quality=min_base_quality,
                    ignore_overlaps=ignore_overlaps,
                )
                # TODO: Discuss with san to use codon count or nuc count to detect proportion
                total_codon_count = 0
                for pileupcol in iter:
                    if pileupcol.pos != selected_coordinate:
                        continue
                    for pread in pileupcol.pileups:
                        if not pread.is_del and not pread.is_refskip:
                            codon = pread.alignment.query_sequence[
                                pread.query_position -
                                shift:pread.query_position - shift + 3]
                            # NOTE: Count codon can keep the same count distribution as the same as bases
                            if (codon in codon_table) and (
                                    pread.alignment.query_sequence[
                                        pread.query_position]
                                    in coordinates_with_change[
                                        pileupcol.pos]["bases"]):

                                if (codon not in coordinates_with_change[
                                        pileupcol.pos]["bases"]
                                    [pread.alignment.query_sequence[
                                        pread.query_position]]["codon_count"]):
                                    coordinates_with_change[
                                        pileupcol.pos]["bases"][
                                            pread.alignment.query_sequence[
                                                pread.query_position]][
                                                    "codon_count"][codon] = 0
                                coordinates_with_change[pileupcol.pos][
                                    "bases"][pread.alignment.query_sequence[
                                        pread.query_position]]["codon_count"][
                                            codon] += 1
                                total_codon_count += 1

                            elif len(codon) > 3:
                                print(codon, selected_coordinate)

                    if pileupcol.pos == selected_coordinate:
                        break
                # NOTE: Removing less less common codons
                coordinates_with_change[selected_coordinate][
                    "total_codon_count"] = total_codon_count

                ref_base = sequences[selected_coordinate].seq
                ref_codon = sequences[
                    selected_coordinate - shift:selected_coordinate - shift +
                    3].seq  # TODO: Integrate the ref codon in dictionary itself
                # NOTE: Reverse complement
                if row["strand"] == "-":
                    ref_codon = str(Seq.Seq(ref_codon).reverse_complement())
                    ref_base = str(Seq.Seq(ref_base).reverse_complement())
                    for k in coordinates_with_change[selected_coordinate][
                            "bases"]:
                        # NOTE: base need to be reverse complemented

                        codon_count = coordinates_with_change[
                            selected_coordinate]["bases"][k]["codon_count"]
                        new_codon_count = {}
                        for codon in codon_count:
                            new_codon_count[str(
                                Seq.Seq(codon).reverse_complement()
                            )] = codon_count[codon]

                        coordinates_with_change[selected_coordinate]["bases"][
                            k]["codon_count"] = new_codon_count
                    new_base = {}
                    for k in coordinates_with_change[selected_coordinate][
                            "bases"]:
                        new_base[str(Seq.Seq(k).reverse_complement(
                        ))] = coordinates_with_change[selected_coordinate][
                            "bases"][k]
                    coordinates_with_change[selected_coordinate][
                        "bases"] = new_base
                coordinates_with_change[selected_coordinate][
                    "ref_codon"] = ref_codon
                coordinates_with_change[selected_coordinate]["codon_pos"] = (
                    selected_coordinate - shift)
                coordinates_with_change[selected_coordinate][
                    "ref_base"] = ref_base

                if row["strand"] == "+":

                    amino_pos = (selected_coordinate - row["start"]) // 3 + 1

                else:
                    amino_pos = (row["end"] - selected_coordinate) // 3 + 1

                coordinates_with_change[selected_coordinate][
                    "amino_pos"] = amino_pos

    ref = {"coor": [], "ref_codon": [], "ref_codon_count": []}
    final_table = {
        "Amino Acid Change": [],
        "Nucleotide Change": [],
        "Codon Change": [],
        "Sample": [],
        "alt_codon": [],
        "alt": [],
        "total": [],
        "coor": [],
        "ref_codon": [],
    }
    # print(coordinates_with_change[21359]['total_codon_count'])
    for coor in coordinates_with_change:
        bases = coordinates_with_change[coor]["bases"]
        read_count = coordinates_with_change[coor]["read_count"]
        ref_codon_bool = False
        for base in bases:
            codon_counts = bases[base]["codon_count"]
            for codon in codon_counts:
                if codon == coordinates_with_change[coor]["ref_codon"]:
                    ref_codon_bool = True
                    ref["coor"].append(coor)
                    ref["ref_codon"].append(
                        coordinates_with_change[coor]["ref_codon"])
                    ref["ref_codon_count"].append(codon_counts[codon])
                    continue
                if (base == coordinates_with_change[coor]["ref_base"]) or (
                        codon_table[coordinates_with_change[coor]["ref_codon"]]
                        == codon_table[codon]):
                    continue
                final_table["Amino Acid Change"].append(
                    f"""{codon_table[coordinates_with_change[coor]['ref_codon']
                        ]}{coordinates_with_change[coor
                            ]['amino_pos']}{codon_table[codon]}""")
                final_table["Nucleotide Change"].append(
                    f'{coor}:{coordinates_with_change[coor]["ref_base"]}>{base}'
                )
                final_table["Codon Change"].append(
                    f"""{coordinates_with_change[coor]["codon_pos"
                        ]}:{coordinates_with_change[coor
                            ]["ref_codon"]}>{codon}""")
                final_table["Sample"].append(sample.split(".bam")[0])
                final_table["alt_codon"].append(codon)
                final_table["alt"].append(codon_counts[codon])
                final_table["total"].append(
                    coordinates_with_change[coor]["total_codon_count"])
                final_table["coor"].append(coor)
                final_table["ref_codon"].append(
                    coordinates_with_change[coor]["ref_codon"])
        # if not ref_codon_bool:
        # ref["coor"].append(coor)
        # print(coor)
        # print(coordinates_with_change[coor])
        # ref["ref_codon"].append(coordinates_with_change[coor]["ref_codon"])
        # ref["ref_codon_count"].append(0)
    final_table = pd.DataFrame(final_table)
    ref = pd.DataFrame(ref)
    final_table = final_table.merge(ref, on=["coor", "ref_codon"], how="inner")
    final_table = final_table[
        (final_table["total"] >= min_seq_depth)
        & (final_table["alt"] > alt_nuc_count * final_table["total"])]
    # print(final_table)
    # print(final_table.columns)
    if not final_table.empty:

        final_table["codon_count"] = final_table[[
            "ref_codon", "ref_codon_count", "alt_codon", "alt"
        ]].apply(
            lambda x: f"""{x['ref_codon']}-{x['ref_codon_count']};{
                x['alt_codon']
         }-{x['alt']}""",
            axis=1,
        )
        final_table["codon_percent"] = final_table[[
            "coor", "ref_codon", "ref_codon_count", "alt_codon", "alt"
        ]].apply(
            lambda x:
            f"""{x['ref_codon']}-{f'%.2f' % (x['ref_codon_count']*100./coordinates_with_change[x['coor']]['total_codon_count'])};{
                x['alt_codon']
         }-{f'%.2f'% (x['alt']*100./coordinates_with_change[x['coor']]['total_codon_count'])}""",
            axis=1,
        )
        # all_codon_count = final_table[]
    final_table = final_table.drop(
        [
            "ref_codon",
            "ref_codon_count",
            "alt_codon",
            "alt",
            "coor",
            "total",
        ],
        axis=1,
    )
    return final_table
