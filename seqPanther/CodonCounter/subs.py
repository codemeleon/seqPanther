#!/usr/bin/env python

import pandas as pd
from .codon_table import codon_table
from os import path
from Bio import Seq


def sub_table(coordinates_with_change, bam, params):
    sample = path.split(bam)[1]
    sequences = params["sequences"]
    alt_codon_frac = params["alt_codon_frac"]

    gff_data = params["gff_data"]
    min_seq_depth = params["min_seq_depth"]
    sample = path.split(bam)[1]
    keys = set(coordinates_with_change)

    for (
            _,
            row,
    ) in gff_data.iterrows():  # TODO: Change it to itertuple for large gff
        selected_coordinates = keys & set(range(row["start"], row["end"] + 1))
        for selected_coordinate in selected_coordinates:
            coordinates_with_change[selected_coordinate]["start"] = row[
                "start"]
            coordinates_with_change[selected_coordinate]["end"] = row["end"]
            coordinates_with_change[selected_coordinate]["strand"] = row[
                "strand"]
            shift = (selected_coordinate - row["start"]) % 3
            ref_base = sequences[selected_coordinate].seq
            ref_codon = sequences[selected_coordinate -
                                  shift:selected_coordinate - shift + 3].seq
            ref_codon_count = 0
            total_codon_count = 0
            codon_counts = {}
            bases = coordinates_with_change[selected_coordinate]["bases"]
            for base in bases.keys():

                codons = bases[base]["codons_count"]
                for extended_codon in codons.keys():
                    codon = extended_codon[2 - shift:5 - shift]
                    if '-' in codon:
                        del codons[extended_codon]
                        continue
                    if codon not in codons:
                        codons[codon] = codons.pop(extended_codon)
                    else:
                        codons[codon] += codons.pop(extended_codon)
                    if codon not in codon_counts:
                        codon_counts[codon] = 0
                    codon_counts[codon] += codons[codon]
                total_codon_count += sum(codons.values())
                try:
                    ref_codon_count += codons[ref_codon]
                except KeyError as _:
                    pass

            # NOTE: Removing less less common codons
            coordinates_with_change[selected_coordinate][
                "total_codon_count"] = total_codon_count

            codons_to_delete = []
            for codon in codon_counts.keys():
                if codon_counts[codon] / total_codon_count < alt_codon_frac:
                    codons_to_delete.append(codon)

            for base in list(bases.key()):
                for codon in codons_to_delete:
                    try:
                        del bases[base]["codons_count"][codon]
                    except KeyError as _:
                        pass

            # NOTE: Reverse complement
            if row["strand"] == "-":
                ref_codon = str(Seq.Seq(ref_codon).reverse_complement())
                ref_base = str(Seq.Seq(ref_base).reverse_complement())
                for k in coordinates_with_change[selected_coordinate]["bases"]:
                    # NOTE: base need to be reverse complemented

                    codon_count = coordinates_with_change[selected_coordinate][
                        "bases"][k]["codon_count"]
                    new_codon_count = {}
                    for codon in codon_count:
                        new_codon_count[str(
                            Seq.Seq(codon).reverse_complement()
                        )] = codon_count[codon]

                    coordinates_with_change[selected_coordinate]["bases"][k][
                        "codon_count"] = new_codon_count
                new_base = {}
                for k in coordinates_with_change[selected_coordinate]["bases"]:
                    new_base[str(Seq.Seq(
                        k).reverse_complement())] = coordinates_with_change[
                            selected_coordinate]["bases"][k]
                coordinates_with_change[selected_coordinate][
                    "bases"] = new_base
            coordinates_with_change[selected_coordinate][
                "ref_codon"] = ref_codon
            coordinates_with_change[selected_coordinate][
                "ref_codon_count"] = ref_codon_count
            coordinates_with_change[selected_coordinate]["codon_pos"] = (
                selected_coordinate - shift)
            coordinates_with_change[selected_coordinate]["ref_base"] = ref_base

            if row["strand"] == "+":

                amino_pos = (selected_coordinate - row["start"]) // 3 + 1

            else:
                amino_pos = (row["end"] - selected_coordinate) // 3 + 1

            coordinates_with_change[selected_coordinate][
                "amino_pos"] = amino_pos

    # NOTE: Dict to Table
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
        "ref_codon_count": [],
    }
    for coor in coordinates_with_change:
        bases = coordinates_with_change[coor]["bases"]
        for base in bases:
            if base == coordinates_with_change[coor]["ref_base"]:
                continue
            codon_counts = bases[base]["codon_count"]
            for codon in codon_counts:
                final_table["Amino Acid Change"].append(
                    f"""{codon_table[coordinates_with_change[coor]['ref_codon']
                        ]}{coordinates_with_change[coor
                            ]['amino_pos']}{codon_table[codon]}"""
                )  # TODO: Correct amino acid position
                final_table["Nucleotide Change"].append(
                    f'{coor+1}:{coordinates_with_change[coor]["ref_base"]}>{base}'
                )
                final_table["Codon Change"].append(
                    f"""{coordinates_with_change[coor]["codon_pos"
                        ]+1}:{coordinates_with_change[coor
                            ]["ref_codon"]}>{codon}"""
                )  # TODO: Correct codon position if it is incorrect
                final_table["Sample"].append(sample.split(".bam")[0])
                final_table["alt_codon"].append(codon)
                final_table["alt"].append(codon_counts[codon])
                final_table["total"].append(
                    coordinates_with_change[coor]["total_codon_count"])
                final_table["coor"].append(
                    coor)  # NOTE: Required later, but also removed
                final_table["ref_codon"].append(
                    coordinates_with_change[coor]["ref_codon"])
                final_table["ref_codon_count"].append(
                    coordinates_with_change[coor]["ref_codon_count"])

    final_table = pd.DataFrame(final_table)
    final_table = final_table[final_table["total"] >= min_seq_depth]
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
