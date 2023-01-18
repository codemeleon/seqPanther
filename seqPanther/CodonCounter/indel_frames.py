#!/usr/bin/env python

import pysam
import pandas as pd
from collections import Counter

from Bio import Seq
import pyfaidx
from os import path


def indel_frames(indel_pos_type_size, bam, params):
    rid = params["rid"]
    gff_data = params["gff_data"]
    sequences = pyfaidx.Fasta(params["ref"])[rid]
    ignore_orphans = params["ignore_orphans"]
    min_mapping_quality = params["min_mapping_quality"]
    min_base_quality = params["min_base_quality"]
    ignore_overlaps = params["ignore_overlaps"]
    alt_nuc_count = params["alt_nuc_count"]
    alt_codon_frac = params["alt_codon_frac"]
    sample = path.split(bam)[1].split(".bam")[0]

    deletion_frame = {}
    insertion_frame = {}
    samfile = pysam.AlignmentFile(bam, "rb")

    for _, row in indel_pos_type_size.iterrows():
        t_gff_data = gff_data[(gff_data["start"] <= row["coor"])
                              & (gff_data["end"] > row["coor"])]
        if not len(t_gff_data):
            print(f'No gff data found for {row["coor"]}')
            continue
        for _, gff_row in t_gff_data.iterrows():
            iter = samfile.pileup(
                rid,
                row["coor"],
                row["coor"] + 1,
                ignore_orphans=ignore_orphans,
                min_mapping_quality=min_mapping_quality,
                min_base_quality=min_base_quality,
                ignore_overlaps=ignore_overlaps,
                max_depth=1000000,
            )

            if row["indel"] < 0:
                for pileupcol in iter:
                    if pileupcol.pos != row["coor"]:
                        continue
                    adjusted_coor = row["coor"] + 1
                    shift = (adjusted_coor - gff_row["start"]) % 3  # + 1
                    r_shift = (3 - shift) % 3
                    ref_sub_seq = sequences[adjusted_coor -
                                            shift:adjusted_coor -
                                            row["indel"] + r_shift].seq
                    amino_pos = (adjusted_coor -
                                 gff_row["start"]) // 3  # - 1 * row["indel"]
                    ref_count = 0
                    depth = 0
                    deleted_codon = []
                    for pread in pileupcol.pileups:
                        if (pread.alignment.reference_start <= pileupcol.pos
                                and pileupcol.pos - row["indel"] <
                                pread.alignment.reference_end):
                            depth += 1
                        if pread.indel:
                            if (pread.indel < 0) and (pread.indel
                                                      == row["indel"]):
                                read_sub_seq = pread.alignment.query_sequence[
                                    pread.query_position - shift +
                                    1:pread.query_position + r_shift + 1]
                                if len(read_sub_seq) % 3 == 0:
                                    deleted_codon.append(read_sub_seq)
                            continue

                        if not pread.is_refskip and not pread.is_del:
                            read_sub_seq = pread.alignment.query_sequence[
                                pread.query_position - shift +
                                1:pread.query_position - row["indel"] +
                                r_shift + 1]
                            # TODO: Check for indels in given surroundings, can be ignored

                            if read_sub_seq == ref_sub_seq:
                                ref_count += 1
                    if pileupcol.pos not in deletion_frame:
                        deletion_frame[pileupcol.pos] = []

                    if gff_row["strand"] == "-":
                        amino_pos = ((gff_row["end"] - gff_row["start"]) // 3 -
                                     (len(ref_sub_seq) -
                                      (shift + r_shift)) // 3 - amino_pos -
                                     (1 if shift else 0))
                        shift, r_shift = r_shift, shift
                        for i, codon in enumerate(deleted_codon):
                            deleted_codon[i] = Seq.Seq(
                                codon).reverse_complement()
                        ref_sub_seq = Seq.Seq(ref_sub_seq).reverse_complement()
                    deleted_codon = dict(Counter(deleted_codon))
                    to_del = []
                    for codon, val in deleted_codon.items():
                        if val * 1. / depth < alt_codon_frac:
                            to_del.append(codon)
                    for codon in to_del:
                        del deleted_codon[codon]

                    deletion_frame[pileupcol.pos].append({  # +1
                        "ref":
                        ref_sub_seq,
                        "shift":
                        shift,
                        "amino_pos":
                        amino_pos,
                        "r_shift":
                        r_shift,
                        "ref_count":
                        ref_count,
                        "alt_count":
                        deleted_codon,
                        "strand":
                        gff_row["strand"],
                        "depth":
                        depth
                    })
                    break

            if row["indel"] > 0:

                for pileupcol in iter:
                    if pileupcol.pos != row["coor"]:
                        continue

                    shift = (row["coor"] - gff_row["start"] + 1) % 3  # + 1
                    amino_pos = (row["coor"] - gff_row["start"] + 1) // 3

                    r_shift = (3 - shift) % 3
                    ref_sub_seq = sequences[row["coor"] - shift +
                                            1:row["coor"] + r_shift + 1].seq
                    inserted_codon = []
                    ref_count = 0
                    depth = 0

                    for pread in pileupcol.pileups:
                        if (pread.alignment.reference_start <= pileupcol.pos
                                and
                                pileupcol.pos < pread.alignment.reference_end):
                            depth += 1

                        if pread.indel > 0:
                            read_sub_seq = pread.alignment.query_sequence[
                                pread.query_position - shift +
                                1:pread.query_position + row["indel"] +
                                r_shift + 1]
                            inserted_codon.append(read_sub_seq)
                        if not pread.is_del and not pread.is_refskip:
                            read_sub_seq = pread.alignment.query_sequence[
                                pread.query_position - shift +
                                1:pread.query_position + r_shift + 1]
                            if read_sub_seq == ref_sub_seq:
                                ref_count += 1

                    if gff_row["strand"] == "-":
                        amino_pos = ((gff_row["end"] - gff_row["start"]) // 3 -
                                     (shift + r_shift) // 3 - amino_pos +
                                     (1 if shift else 0))
                        for i, codon in enumerate(inserted_codon):
                            inserted_codon[i] = Seq.Seq(
                                codon).reverse_complement()
                        ref_sub_seq = Seq.Seq(ref_sub_seq).reverse_complement()
                    inserted_codon = dict(Counter(inserted_codon))
                    to_del = []
                    for codon, val in inserted_codon.items():
                        if val * 1. / depth < alt_codon_frac:
                            to_del.append(codon)
                    for codon in to_del:
                        del inserted_codon[codon]
                    insertion_frame[pileupcol.pos  # - shift
                                    ] = {  # Location is where codon start
                                        "ref": ref_sub_seq,
                                        "amino_pos": amino_pos,
                                        "ref_count": ref_count,
                                        "shift": shift,
                                        "r_shift": r_shift,
                                        "alt_count": inserted_codon,
                                        "strand": gff_row["strand"],
                                        "depth": depth
                                    }
                    break
    delete_final_table = {
        "Amino Acid Change": [],
        "Nucleotide Change": [],
        "Codon Change": [],
        "Sample": [],
        "alt_codon": [],
        "alt_codon_count": [],
        "coor": [],
        "ref_codon": [],
        "ref_codon_count": [],
        "total_codon_count": [],
    }
    print(deletion_frame)
    print(insertion_frame)
    for coor in deletion_frame:
        codon_count = 0
        for x in deletion_frame[coor]:
            codon_count += x["ref_count"] + sum(x["alt_count"].values())

        for info in deletion_frame[coor]:
            for inserted_seq in info["alt_count"]:
                delete_final_table["Amino Acid Change"].append(
                    f"{Seq.Seq(info['ref']).translate()}{info['amino_pos']}{Seq.Seq(inserted_seq).translate()}"
                )
                delete_final_table["Codon Change"].append(
                    f"{coor+(1 if info['shift'] else 2)}:{info['ref']}>{inserted_seq}"
                )
                delete_final_table["Nucleotide Change"].append(
                    f"{coor+(1 if info['shift'] else 2) +(info['shift'] if info['strand']=='+' else info['r_shift'])}:{info['ref'][info['shift']:-1*info['r_shift'] if info['shift'] else len(info['ref'])]}>"
                )
                delete_final_table["coor"].append(coor)
                delete_final_table["ref_codon"].append(info["ref"])
                delete_final_table["ref_codon_count"].append(info["ref_count"])
                delete_final_table["Sample"].append(sample)
                delete_final_table["alt_codon"].append(inserted_seq)
                delete_final_table["alt_codon_count"].append(
                    info["alt_count"][inserted_seq])
                delete_final_table["total_codon_count"].append(codon_count)
    delete_final_table = pd.DataFrame(delete_final_table)
    delete_final_table = delete_final_table[
        delete_final_table["alt_codon_count"] >= alt_nuc_count *
        (delete_final_table["alt_codon_count"] +
         delete_final_table["ref_codon_count"])]
    if not delete_final_table.empty:
        delete_final_table["codon_count"] = delete_final_table.apply(
            lambda x:
            f"{x['ref_codon']}-{x['ref_codon_count']};{x['alt_codon']}-{x['alt_codon_count']}",
            axis=1,
        )
        delete_final_table["codon_percent"] = delete_final_table.apply(
            lambda x:
            f"{x['ref_codon']}-{'%0.2f' % (x['ref_codon_count']*100./x['total_codon_count'])};{x['alt_codon']}-{'%0.2f' % (x['alt_codon_count']*100./x['total_codon_count'])}",
            axis=1,
        )
    delete_final_table = delete_final_table.drop(
        [
            "ref_codon", "ref_codon_count", "alt_codon", "alt_codon_count",
            "coor"
        ],
        axis=1,
    )

    insert_final_table = {
        "Amino Acid Change": [],
        "Nucleotide Change": [],
        "Codon Change": [],
        "Sample": [],
        "alt_codon": [],
        "alt_codon_count": [],
        "coor": [],
        "ref_codon": [],
        "ref_codon_count": [],
    }
    for coor in insertion_frame:
        for inserted_seq in insertion_frame[coor]["alt_count"]:
            insert_final_table["Amino Acid Change"].append(
                f"{Seq.Seq(insertion_frame[coor]['ref']).translate()}{insertion_frame[coor]['amino_pos']}{Seq.Seq(inserted_seq).translate()}"
            )
            insert_final_table["Codon Change"].append(
                f"{coor+1}:{insertion_frame[coor]['ref']}>{inserted_seq}")
            insert_final_table["Nucleotide Change"].append(
                f"{coor+1+insertion_frame[coor]['shift']}:>{inserted_seq[insertion_frame[coor]['shift']:len(inserted_seq)-1*insertion_frame[coor]['r_shift']]}"
            )
            # insert_final_table["Codon Change"].append(f"{}{}>{}")
            insert_final_table["coor"].append(coor)
            insert_final_table["ref_codon"].append(
                insertion_frame[coor]["ref"])
            insert_final_table["ref_codon_count"].append(
                insertion_frame[coor]["ref_count"])
            insert_final_table["Sample"].append(sample)
            insert_final_table["alt_codon"].append(inserted_seq)
            insert_final_table["alt_codon_count"].append(
                insertion_frame[coor]["alt_count"][inserted_seq])
    insert_final_table = pd.DataFrame(insert_final_table)
    insert_final_table = insert_final_table[
        insert_final_table["alt_codon_count"] >= alt_nuc_count *
        (insert_final_table["alt_codon_count"] +
         insert_final_table["ref_codon_count"])]
    if not insert_final_table.empty:
        insert_final_table["codon_count"] = insert_final_table.apply(
            lambda x:
            f"{x['ref_codon']}-{x['ref_codon_count']};{x['alt_codon']}-{x['alt_codon_count']}",
            axis=1,
        )
        insert_final_table["codon_percent"] = insert_final_table.apply(
            lambda x:
            f"{x['ref_codon']}-{'%0.2f' % (x['ref_codon_count']*100./(x['ref_codon_count']+x['alt_codon_count']))};{x['alt_codon']}-{'%0.2f' % (x['alt_codon_count']*100./(x['ref_codon_count']+x['alt_codon_count']))}",
            axis=1,
        )
    insert_final_table = insert_final_table.drop(
        [
            "ref_codon", "ref_codon_count", "alt_codon", "alt_codon_count",
            "coor"
        ],
        axis=1,
    )

    return delete_final_table, insert_final_table
