#!/usr/bin/env python

# Author: Anmol Kiran
# Affiliation: University of Liverpool. UK
# and Malawi-Liverpool-Wellcome Trust, Malawi
# V0.0.1: 20/04/2021

from Bio import Seq, SeqIO
import pandas as pd
import numpy as np
import click
import warnings
from functools import partial

# import re

# from collections import ChainMap
from shutil import copyfile, rmtree

# from functools import partial
import tempfile as tmf
from os import makedirs, path  # , access
from glob import glob

# from multiprocessing import cpu_count, Pool
from collections import defaultdict

# import sys
import subprocess as sb

# import os

from pandas._libs.lib import is_list_like

warnings.filterwarnings("ignore")
# warnings.simplefilter(action="ignore", category=FutureWarning)


__version__ = "0.0.1"

_s_gene_seq = """ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTT
AATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAA
GTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCAT
GCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTAT
TTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCC
CTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTG
GGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGC
ACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAA
TTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTC
CCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTA
CTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGCTTATTAT
GTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGT
GCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCT
AACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTT
TTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCT
GTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGC
TTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGA
AAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTT
GATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGA
GATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCT
TTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAA
CTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTC
AACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGC
AGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCT
TTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAAC
TGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAAT
GTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATT
GGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCC
ATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACA
AATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTAC
ATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCT
TTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACA
CCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCA
TTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTT
GGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACA
GATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGT
GCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTC
TATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACA
GCAAGTGCACTTGGAAAACTTCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTT
AGCTCCAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAGTG
CAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAACAATTAATTAGAGCTGCA
GAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGAT
TTTTGTGGAAAGGGCTATCATCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACT
TATGTCCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTCGT
GAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATGAACCACAAATCATTACT
ACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTG
CAACCTGAATTAGACTCATTCAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTA
GGTGACATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTTGCCAAG
AATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCAGTATATAAAATGGCCATGGTACATT
TGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGT
AGTTGTCTCAAGGGCTGTTGTTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAA
GGAGTCAAATTACATTACACATAA"""

_amb_base = {  # All ambiguous nucleotides
    "N": set(["A", "C", "G", "T"]),
    "R": set(["A", "G"]),
    "Y": set(["T", "C"]),
    "K": set(["G", "T"]),
    "M": set(["A", "C"]),
    "S": set(["G", "C"]),
    "W": set(["A", "T"]),
    "B": set(["C", "G", "T"]),
    "D": set(["A", "G", "T"]),
    "H": set(["A", "C", "T"]),
    "V": set(["A", "C", "G"]),
}

_amb_base_ext = {  # Opposite of above
    tuple(set(["A", "C", "G", "T"])): "N",
    tuple(set(["A", "G"])): "R",
    tuple(set(["T", "C"])): "Y",
    tuple(set(["G", "T"])): "K",
    tuple(set(["A", "C"])): "M",
    tuple(set(["G", "C"])): "S",
    tuple(set(["A", "T"])): "W",
    tuple(set(["C", "G", "T"])): "B",
    tuple(set(["A", "G", "T"])): "D",
    tuple(set(["A", "C", "T"])): "H",
    tuple(set(["A", "C", "G"])): "V",
}


def cmd(command):
    """Runs all the command using Popen communicate.
    :command: in form of list ie ['ls', '-l']
    :return: None
    """
    sb.Popen(command, stdout=sb.DEVNULL, stderr=sb.DEVNULL).communicate()


def min_max(values):
    """Returns minimum and maximum values in the given list/array.

    :values: list/array
    :returns: touple of min and max

    """
    return min(values), max(values)


#           2       20


def trim(mmcount, length, aln_df):
    """Trim the sequence at the end based on mismatched in given length of

    :mmcount: Mismatch count in given range
    :length: Length at both end to explore
    :returns: trimmed alignment dataframe

    """
    # aln_len = length(aln)
    for col in aln_df.columns:
        if col == "ref":
            continue
        min_, max_ = min_max(aln_df[aln_df[col] != "-"].index)
        mismach_locations = aln_df[aln_df[col].isin("-NRYKMSWBDHV")].index
        mismach_locations = mismach_locations[
            (mismach_locations >= min_) & (mismach_locations <= max_)
        ]
        start_mismatch_location = mismach_locations[mismach_locations < length]
        if len(start_mismatch_location) >= mmcount:
            min_ = start_mismatch_location[-1] + 1
        end_mismatch_locations = mismach_locations[
            mismach_locations > (length(aln_df) - length)
        ]
        if len(end_mismatch_locations) >= mmcount:
            max_ = end_mismatch_locations[0] - 1
        aln_df.loc[:min_, col] = "-"
        aln_df.loc[max_:, col] = "-"

    return aln_df


def codon_aln(aln_df):
    """Correct alignment around codon

    :aln_df: normal alignment dataframe
    :returns: colodon alinment dataframe

    """
    df_shape = aln_df.shape[1]
    non_ref_seq = [seq_id for seq_id in df_shape.columns if seq_id != "ref"]
    if df_shape == 2:
        aln_min, aln_max = min_max(aln_df[aln_df[non_ref_seq[0]] != "-"].index)
        # TODO: Check whether is a 0 position if not change the locatio in no indel at next two positions
        #

        pass
    elif df_shape == 3:
        # TODO: Accept nucleotide at overhang and then try to change the codon alignmentt

        pass
    else:
        pass

    return


# TODO: How many codon should be allow to be remove without replacing with N
# TODO: Remove text from last mismatch - I don't think it will make any difference -  But allow to extend the ends for better mapping


def ranges(lst, given_gap=0):
    """
    A generator returns list of range based on given numbers
    [1,2,3,4,5, 10, 11, 12, 13, 14] => [[1,5], [10,14]]
    """
    # print(given_gap, "Anmol")
    lst_sorted = sorted(lst)
    init = 0
    for num in range(1, len(lst_sorted)):
        # reported_gap = lst_sorted[num] - lst_sorted[num - 1] - 1
        if lst_sorted[num] > given_gap + 1 + lst_sorted[num - 1]:
            # or (reported_gap % 3 != 0):
            # gap=0 means overlapping
            yield (lst_sorted[init], lst_sorted[num - 1])
            init = num
    yield (lst_sorted[init], lst_sorted[-1])


def useful_range(lst, gap=10):  # TODO: Add in command
    """Returns first occurance maximum length block range."""
    # TODO: Check with other that if they are happy with it
    # or they want other fragment to be selected

    trange = list(ranges(lst, gap))
    # print(trange, "kiran")
    pre_frag_len = 0
    myrange = trange[0]
    for rng in trange:
        frag_len = rng[1] - rng[0] + 1
        if pre_frag_len < frag_len:
            pre_frag_len = frag_len
            myrange = rng
    return myrange


def df_reverse_complement(dataframe):
    """Reverse Complement of the nucleotide dataframe."""
    dataframe = dataframe.loc[::-1]

    def rvc(nuc):
        """Return complement base."""
        return str(Seq.Seq(nuc).complement())

    def rvc_dict(nuc_dict):
        """Returns complement dictionary."""
        temp_nuc_dict = {}
        for nuc in nuc_dict:
            temp_nuc_dict[rvc(nuc)] = nuc_dict[nuc]
        return temp_nuc_dict

    dataframe["nuc"] = dataframe["nuc"].apply(rvc)
    dataframe["peak"] = dataframe["peak"].apply(rvc_dict)
    return dataframe


def orient(seqfile, ref, tmp_fold):
    """Returns orientation of the sequence"""
    if seqfile.endswith(".ab1"):
        flb = path.split(seqfile)[1].rsplit(".", 1)[0]
        record = SeqIO.read(seqfile, "abi")
        seq = "".join(
            [chr(ascii_val)
             for ascii_val in record.annotations["abif_raw"]["PBAS1"]]
        )
        with open(f"{tmp_fold}/{flb}.fasta", "w") as fout:
            fout.write(f">tmp\n{seq}\n")
        seqfile = f"{tmp_fold}/{flb}.fasta"

    flb = path.split(seqfile)[1].rsplit(".", 1)[0]
    command = [
        "blat",
        "-noHead",
        ref,
        seqfile,
        f"{tmp_fold}/{flb}.psl",
    ]  # Mapping against reference
    cmd(command)
    blat_df = (
        pd.read_table(f"{tmp_fold}/{flb}.psl", header=None)
        .sort_values(0, ascending=False)
        .drop_duplicates(9)
    )
    if len(blat_df) != 1:
        print(f"No match found of {seqfile}. Ignoring")
        return None
    for _, row in blat_df.iterrows():
        if row[8] == "-":
            return "R"
        else:
            return "F"


def ab1_to_fasta_wihout_ref(ab1_list, tmp_fold, res_fold):
    """Convert ab1 to fasta without reference. Not functional at moment"""
    # TODO: If sequence is single just move to results after trimming bad ends

    if len(ab1_list) == 1:
        seqfile = ab1_list[0]
        if seqfile.endswith(".ab1"):
            flb = path.split(seqfile)[1].split(".", 1)[0]
            record = SeqIO.read(seqfile, "abi")
            seq = "".join(
                [
                    chr(ascii_val)
                    for ascii_val in record.annotations["abif_raw"]["PBAS1"]
                ]
            )
            with open(f"{res_fold}/{flb}.fasta", "w") as fout:
                fout.write(f">{flb}\n{seq}\n")
        else:
            copyfile(seqfile, res_fold)
    else:
        flb = path.split(ab1_list[0])[1].split(".", 1)[0]
        # TODO: Check whether both seqyences are same strand
        with open(f"{tmp_fold}/{flb}.fasta", "w") as fout:
            for num, seqfile in enumerate(ab1_list):
                if seqfile.endswith(".ab1"):
                    record = SeqIO.read(seqfile, "abi")
                    seq = "".join(
                        [
                            chr(ascii_val)
                            for ascii_val in record.annotations["abif_raw"]["PBAS1"]
                        ]
                    )
                    if num:
                        seq = Seq.Seq(seq).reverse_complement()
                    fout.write(f">{flb}_{num}\n{seq}\n")
        command = [
            "blat",
            "-noHead",
            f"{tmp_fold}/{flb}.fasta",
            f"{tmp_fold}/{flb}.fasta",
            f"{tmp_fold}/{flb}.psl",
        ]  # Mapping against reference
        cmd(command)
        blat_df = pd.read_table(f"{tmp_fold}/{flb}.psl", header=None).sort_values(
            0, ascending=False
        )
        blat_df = blat_df[blat_df[9] != blat_df[13]].head(2)
        strands = list(set(blat_df[8].values))
        if len(strands) > 1:
            print(f"{flb} has multiple strands. Ignoring")
            return None
        strands = strands[0]
        if strands == "-":
            sequences = {}
            for i, rec in enumerate(SeqIO.parse(f"{tmp_fold}/{flb}.fasta", "fasta")):
                if i:
                    seq = Seq.Seq(rec.seq).reverse_complement()
                else:
                    seq = rec.seq
                sequences[rec.id] = seq
            with open(f"{tmp_fold}/{flb}.fasta", "w") as fout:
                for k in sequences:
                    fout.write(f">{k}\n{sequences[k]}\n")

        command = [
            "muscle",
            "-in",
            f"{tmp_fold}/{flb}.fasta",
            "-out",
            f"{tmp_fold}/{flb}_aln.fasta",
        ]
        cmd(command)
        seq_df = {}
        for rec in SeqIO.parse(f"{tmp_fold}/{flb}_aln.fasta", "fasta"):
            seq_df[rec.id] = list(rec.seq)
        seq_df = pd.DataFrame(seq_df)
        # TODO: add options to extend or not to extend ovelappes
        # TODO: How to select the bases

        # TODO: Look for gaps at the end and trim them
        with open(f"{res_fold}/{flb}.fasta", "w") as fout:
            fout.write(f">{flb}\n{str(seq_df.iloc[0])}\n")  # Change this part


def ab1seq(infile):
    """ab1 to seq trimmed based on reference."""

    bases = {"DATA9": "G", "DATA10": "A", "DATA11": "T", "DATA12": "C"}

    record = SeqIO.read(infile, "abi")
    trace = defaultdict(list)

    for channel in bases:
        trace[bases[channel]] = record.annotations["abif_raw"][channel]

    nuc_df = {"nuc": [], "peak": []}

    for channel in zip(
        record.annotations["abif_raw"]["PBAS1"],
        record.annotations["abif_raw"]["PLOC1"],
    ):
        ambi_base = chr(channel[0])
        nuc_df["nuc"].append(ambi_base)
        if ambi_base in _amb_base:
            # if not peak_selection:
            td = {}  # TODO: Please check what does td reprensts
            for base in _amb_base[ambi_base]:
                td[base] = trace[base][channel[1]]
            nuc_df["peak"].append(td)

        else:
            nuc_df["peak"].append({ambi_base: trace[ambi_base][channel[1]]})
    # NOTE: In case of ambigious nucleotides
    nuc_df = pd.DataFrame(nuc_df)
    peak = nuc_df["peak"].apply(lambda x: np.mean(list(x.values()))).values
    mid_point = int(len(peak) / 2)
    spread = int(mid_point / 2.0)
    peak_mean = np.mean(peak[mid_point - spread: mid_point + spread])
    if peak_mean < 50:
        warnings.warn(
            f"Peak mean is {peak_mean} which is below 50 for {infile}. "
            "This may be due to low quality sequence. "
            "Please check the quality of the sequence."
        )

    correct_peak_pos = np.where(peak > 0.2 * peak_mean)[0]
    # TODO: Autothreshold, check with the sliding window
    # Added 10 extra as towards the end of the sequence quality goes down
    min_pos, max_pos = correct_peak_pos[0], correct_peak_pos[-1]
    nuc_df = nuc_df.loc[min_pos:max_pos]

    if infile.endswith(".R.ab1"):
        nuc_df = df_reverse_complement(nuc_df)
    return nuc_df


def aln_df_with_ref(seq_dict, flb, tmp_fold):
    """Alignment dataframe"""
    inf = f"{tmp_fold}/{flb}.in.fasta"
    otf = f"{tmp_fold}/{flb}.out.fasta"
    with open(inf, "w") as fout:
        for k in seq_dict:
            fout.write(f">{k}\n{seq_dict[k]}\n")
    command = ["muscle", "-in", inf, "-out", otf]
    cmd(command)
    sequences = {}
    for rec in SeqIO.parse(otf, "fasta"):
        sequences[rec.id] = list(str(rec.seq).upper())
    return pd.DataFrame(sequences)


def drop_from_here(arg1):
    """TODO: Docstring for drop_from_here.

    :arg1: TODO
    :returns: TODO

    """
    pass


def merge_base_peak(nuc_df, peak_dict):
    """Merge the peak related information in nucleotide dataframe."""
    # NOTE: Useful only when there is ambigiuity in the base call, rest remove to save memory

    nuc_df["idx"] = list(nuc_df.index)

    # df = enumerate_columns(df)
    # print(df)
    to_drop = ["idx"]
    for nuc in peak_dict:
        # Expecting sequence will not change while aligning, except insertion of gaps
        peak_dict[nuc][f"{nuc}_idx"] = list(
            nuc_df[nuc_df[f"{nuc}"] != "-"].index
        )  # list(range(len(peak_dict[k]))
        nuc_df = nuc_df.merge(
            peak_dict[nuc], left_on="idx", right_on=f"{nuc}_idx", how="outer"
        )
        to_drop += [f"{nuc}_idx", f"{nuc}_nuc"]

    nuc_df = nuc_df.drop(to_drop, axis=1)
    return nuc_df


def rep_paired_base(lst, ambiguous=False):
    """Selecting representative base in presence of forward and
    reverse sanger sequencing data."""
    # NOTE: Select nucleotide based on peak value. Might be incorrect
    if lst["F"] == lst["R"]:
        return lst["F"]
    if lst["F"] == "-":
        if lst["R"] == "-":
            return "-"
        else:
            if ambiguous:
                return lst["R"]
                # amb_base = tuple(set(lst["R_peak"]))
                # return _amb_base_ext[amb_base]
            else:
                bmax = 0
                base = "A"
                peaks = lst["R_peak"]
                for bs in peaks:
                    if peaks[bs] > bmax:
                        bmax = peaks[bs]
                        base = bs
                return base
    else:
        if lst["R"] == "-":
            if ambiguous:
                return lst["F"]
                # amb_base = tuple(set(lst["F_peak"]))
                # return _amb_base_ext[amb_base]
            else:
                # TODO: Correct issues because of gap
                bmax = 0
                base = "A"
                peaks = lst["F_peak"]
                for bs in peaks:
                    if peaks[bs] > bmax:
                        bmax = peaks[bs]
                        base = bs
                return base
        else:
            # if ambiguous:
            # print(lst["F_peak"], lst["R_peak"], "Kkiran", lst)
            amb_base = tuple(list(set(lst["F_peak"]) & set(lst["R_peak"])))
            if not amb_base:
                amb_base = tuple(list(set(lst["F_peak"]) | set(lst["R_peak"])))
            if len(amb_base) == 1:
                return amb_base[0]
            else:
                if ambiguous:
                    return _amb_base_ext[amb_base]
                else:
                    tbase = "A"
                    tval = 0
                    for bs in amb_base:
                        if bs in lst["F_peak"]:
                            if lst["F_peak"][bs] > tval:
                                tval = lst["F_peak"][bs]
                                tbase = bs
                        if bs in lst["R_peak"]:
                            if lst["R_peak"][bs] > tval:
                                tval = lst["R_peak"][bs]
                                tbase = bs

                    return tbase


def aln_clean(aln_df, gap=15, ambiguous=False):  # , at, res_fold):
    """Clearning alihnment dataframe to remove unneccessary gaps and
    alignements."""
    sang_type = None
    rev = {"R": "F", "F": "R"}
    if "F" in aln_df and "R" in aln_df:
        idx = aln_df[~((aln_df["F"] == "-") & (aln_df["R"] == "-"))].index
        idx_f = aln_df[aln_df["F"] != "-"].index
        idx_r = aln_df[aln_df["R"] != "-"].index
        u_range_f = list(useful_range(list(idx_f), gap))
        aln_df.loc[: u_range_f[0] - 1, "F"] = "-"
        aln_df.loc[u_range_f[-1] + 1:, "F"] = "-"
        u_range_r = list(useful_range(list(idx_r), gap))
        aln_df.loc[: u_range_r[0] - 1, "R"] = "-"
        aln_df.loc[u_range_r[-1] + 1:, "R"] = "-"
        locations_to_check = [
            u_range_f[0],
            u_range_f[0] + 1,
            u_range_r[0],
            u_range_r[0] + 1,
        ]
        for ltc in locations_to_check:
            for dir in ["F", "R"]:
                if aln_df.loc[ltc, dir] == "-" and aln_df.loc[ltc, rev[dir]] != "-":
                    if aln_df.loc[ltc, rev[dir]] != aln_df.loc[ltc, "ref"]:
                        if (aln_df.loc[ltc, rev[dir]] in _amb_base) and (
                            aln_df.loc[ltc, "ref"]
                            in _amb_base[aln_df.loc[ltc, rev[dir]]]
                        ):
                            pass
                        else:
                            aln_df.loc[:ltc, rev[dir]] = "-"
        locations_to_check = [
            u_range_f[-1],
            u_range_f[-1] - 1,
            u_range_r[-1],
            u_range_r[-1] - 1,
        ]
        for ltc in locations_to_check:
            for dir in ["F", "R"]:
                if aln_df.loc[ltc, dir] == "-" and aln_df.loc[ltc, rev[dir]] != "-":
                    if aln_df.loc[ltc, rev[dir]] != aln_df.loc[ltc, "ref"]:
                        if (aln_df.loc[ltc, rev[dir]] in _amb_base) and (
                            aln_df.loc[ltc, "ref"]
                            in _amb_base[aln_df.loc[ltc, rev[dir]]]
                        ):
                            pass
                        else:
                            aln_df.loc[ltc:, rev[dir]] = "-"

        # print(aln_df.loc[u_range_r[-1] + 1 : u_range_r[-1] + 30, "R"])
        # print(u_range_f, u_range_r, "dead zone")
        u_range = [
            np.min([u_range_f[0], u_range_r[0]]),
            np.max([u_range_f[-1], u_range_r[-1]]),
        ]
        aln_df.loc[: u_range[0] - 1, ["F", "R"]] = "-"
        aln_df.loc[u_range[1] + 1:, ["F", "R"]] = "-"
        sang_type = "P"  # For paired
        # print(aln_df.head())
        # TODO: Must check in case of mismatch and ambigious it remove it carefully
        # print(idx)
    else:
        if "F" in aln_df:
            idx = aln_df[aln_df["F"] != "-"].index
            u_range = list(useful_range(list(idx), gap))
            sang_type = "F"
            aln_df.loc[: u_range[0] - 1, "F"] = "-"
            aln_df.loc[u_range[1] + 1:, "F"] = "-"
        else:
            idx = aln_df[aln_df["R"] != "-"].index
            u_range = list(useful_range(list(idx), gap))
            sang_type = "R"
            aln_df.loc[: u_range[0] - 1, "R"] = "-"
            aln_df.loc[u_range[1] + 1:, "R"] = "-"

    # TODO: This part need to be corrected
    for col in aln_df.columns:

        if col not in ["F", "R"]:
            continue
        ambi_indexes = aln_df[aln_df[col].isin(list("NRYKMSWBDHV"))].index
        ambi_indexes = ambi_indexes[
            (ambi_indexes >= u_range[0]) & (ambi_indexes <= u_range[-1])
        ]
        # print(u_range, ambi_indexes, "anmol", u_range)

        if len(ambi_indexes):
            print(ambi_indexes)
            ambi_index_left = ambi_indexes[ambi_indexes < u_range[0] + 20]
            ambi_index_right = ambi_indexes[ambi_indexes > u_range[-1] - 20]
            if ambi_index_left.any():
                u_range[0] = ambi_index_left.max()
            if ambi_index_right.any():
                u_range[-1] = ambi_index_right.min()
            # print(u_range, "testing")
            # exit()
            aln_df.loc[: u_range[0] - 1, col] = "-"
            aln_df.loc[u_range[1] + 1:, col] = "-"
    # print(aln_df.loc[u_range_r[-1]: u_range_r[-1] + 30, "R"])
    # TODO: Avoid gap

    # # TODO: Perform average base call for gaps

    if sang_type in "FR":
        # TODO: Select if both shows the same base
        aln_df["consensus"] = aln_df[sang_type].values
        aln_df.loc[: u_range[0] - 1, "consensus"] = aln_df.loc[
            : u_range[0] - 1, "ref"
        ].values
        aln_df.loc[u_range[1] + 1:, "consensus"] = aln_df.loc[
            u_range[1] + 1:, "ref"
        ].values
        # TODO: Create a file for ambigious nucleotides
        if not ambiguous:
            ambi_indexes = aln_df[aln_df[sang_type].isin(
                list("NRYKMSWBDHV"))].index
            # Delete 3 in range of 20

            for ambi_index in ambi_indexes:
                bmax = 0
                base = "A"
                peaks = aln_df.loc[ambi_index, f"{sang_type}_peak"]
                for bs in peaks:
                    if peaks[bs] > bmax:
                        bmax = peaks[bs]
                        base = bs
                aln_df.loc[ambi_index, "consensus"] = base

        insert_ranges = aln_df.loc[u_range[0]: u_range[1]]
        insert_ranges = insert_ranges[insert_ranges["ref"] == "-"].index
        if insert_ranges.any():
            insert_ranges = ranges(insert_ranges)
            for insert_range in insert_ranges:
                if (insert_range[1] - insert_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (insert_range[1] - insert_range[0] + 1) < 3:
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "-"
                    else:
                        # TODO: Talk to san one more time
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "N"
        del_ranges = aln_df.loc[u_range[0]: u_range[1]]
        del_ranges = del_ranges[del_ranges[sang_type] == "-"].index
        if del_ranges.any():
            del_ranges = ranges(del_ranges)
            for del_range in del_ranges:
                if (del_range[1] - del_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (del_range[1] - del_range[0] + 1) < 3:
                        aln_df.loc[
                            del_range[0]: del_range[1], "consensus"
                        ] = aln_df.loc[del_range[0]: del_range[1], "ref"].values
                    else:
                        aln_df.loc[del_range[0]
                            : del_range[1], "consensus"] = "N"

    else:
        mistmatched_indexes = aln_df[aln_df["F"] != aln_df["R"]].index
        # TODO: Reduce the size of the index
        f_range = aln_df[aln_df["F"] != "-"].index
        f_range = f_range.min(), f_range.max()
        r_range = aln_df[aln_df["R"] != "-"].index
        r_range = r_range.min(), r_range.max()
        f_end = np.max([f_range[0], r_range[0]])
        r_end = np.min([f_range[1], r_range[1]])
        mistmatched_indexes = mistmatched_indexes[
            (mistmatched_indexes >= f_end) & (mistmatched_indexes <= r_end)
        ]
        for mismatched_index in mistmatched_indexes:
            # print(
            # "test",
            # mistmatched_indexes,
            # aln_df.loc[mismatched_index, "F"],
            # aln_df.loc[mismatched_index, "R"],
            # )

            if (
                aln_df.loc[mismatched_index, "F"] == "-"
                or aln_df.loc[mismatched_index, "R"] == "-"
            ):
                continue
            if (
                aln_df.loc[mismatched_index - 1, "F"] == "-"
                and aln_df.loc[mismatched_index - 1, "R"]
                == aln_df.loc[mismatched_index, "F"]
            ):
                aln_df.loc[mismatched_index, "F"] = "-"
                aln_df.loc[mismatched_index - 1, "F"] = aln_df.loc[
                    mismatched_index - 1, "R"
                ]
                aln_df.loc[mismatched_index - 1, "F_peak"] = [
                    aln_df.loc[mismatched_index, "F_peak"]
                ]

                aln_df.loc[mismatched_index, "F_peak"] = np.nan

            elif (
                aln_df.loc[mismatched_index + 1, "F"] == "-"
                and aln_df.loc[mismatched_index + 1, "R"]
                == aln_df.loc[mismatched_index, "F"]
            ):
                aln_df.loc[mismatched_index, "F"] = "-"
                aln_df.loc[mismatched_index + 1, "F"] = aln_df.loc[
                    mismatched_index + 1, "R"
                ]

                aln_df.loc[mismatched_index + 1, "F_peak"] = [
                    aln_df.loc[mismatched_index, "F_peak"]
                ]
                aln_df.loc[mismatched_index, "F_peak"] = np.nan

            if (
                aln_df.loc[mismatched_index - 1, "R"] == "-"
                and aln_df.loc[mismatched_index - 1, "F"]
                == aln_df.loc[mismatched_index, "R"]
            ):
                aln_df.loc[mismatched_index, "R"] = "-"
                aln_df.loc[mismatched_index - 1, "R"] = aln_df.loc[
                    mismatched_index - 1, "F"
                ]
                aln_df.loc[mismatched_index - 1, "R_peak"] = [
                    aln_df.loc[mismatched_index, "R_peak"]
                ]
                aln_df.loc[mismatched_index, "R_peak"] = np.nan
            elif (
                aln_df.loc[mismatched_index + 1, "R"] == "-"
                and aln_df.loc[mismatched_index + 1, "F"]
                == aln_df.loc[mismatched_index, "R"]
            ):
                # print(
                # mismatched_index,
                # aln_df.loc[mismatched_index + 1, "R"],
                # aln_df.loc[mismatched_index + 1, "F"],
                # aln_df.loc[mismatched_index, "R"],
                # "TEsting",
                # )
                aln_df.loc[mismatched_index, "R"] = "-"
                aln_df.loc[mismatched_index + 1, "R"] = aln_df.loc[
                    mismatched_index + 1, "F"
                ]
                aln_df.loc[mismatched_index + 1, "R_peak"] = [
                    aln_df.loc[mismatched_index, "R_peak"]
                ]
                aln_df.loc[mismatched_index, "R_peak"] = np.nan
        # aln_df.to_csv("test.csv")

        mistmatched_indexes = aln_df[aln_df["F"] != aln_df["R"]].index
        mistmatched_indexes = mistmatched_indexes[
            (mistmatched_indexes >= f_end) & (mistmatched_indexes <= r_end)
        ]

        # Corrections around single nucleotide indels
        # print(f_range, "abd", r_range, f_end, r_end, "KK")
        for mismatched_index in mistmatched_indexes:
            for rv in rev:

                if (
                    (aln_df.loc[mismatched_index, f"{rv}"] == "-")
                    & (aln_df.loc[mismatched_index - 1, f"{rv}"] != "-")
                    & (aln_df.loc[mismatched_index + 1, f"{rv}"] != "-")
                ):
                    if aln_df.loc[mismatched_index, f"{rev[rv]}"] not in _amb_base:
                        aln_df.loc[mismatched_index, f"{rv}"] = aln_df.loc[
                            mismatched_index, f"{rev[rv]}"
                        ]
                    else:
                        if (
                            aln_df.loc[mismatched_index + 1, f"{rev[rv]}"]
                            == aln_df.loc[mismatched_index - 1, f"{rev[rv]}"]
                        ) and (
                            aln_df.loc[mismatched_index + 1, f"{rev[rv]}"]
                            in _amb_base[aln_df.loc[mismatched_index, f"{rev[rv]}"]]
                        ):
                            aln_df.loc[mismatched_index, f"{rv}"] = aln_df.loc[
                                mismatched_index + 1, f"{rev[rv]}"
                            ]
                            aln_df.loc[mismatched_index, f"{rev[rv]}"] = aln_df.loc[
                                mismatched_index + 1, f"{rev[rv]}"
                            ]
                        else:
                            aln_df.loc[mismatched_index, f"{rv}"] = aln_df.loc[
                                mismatched_index, f"{rev[rv]}"
                            ]
                            # TODO: insert a variable inticating that only on was present

        # TODO: Mismach near ends

        # aln_df.to_csv("test2.csv")

        # exit(1)
        rep_paired_base_fn = partial(rep_paired_base, ambiguous=ambiguous)
        aln_df["consensus"] = aln_df.apply(rep_paired_base_fn, axis=1)

        insert_ranges = aln_df.loc[u_range[0]: u_range[1]]
        insert_ranges = insert_ranges[
            (insert_ranges["ref"] == "-")
            & ((insert_ranges["F"] != "-") & (insert_ranges["R"] != "-"))
        ].index
        if insert_ranges.any():
            insert_ranges = ranges(insert_ranges)
            for insert_range in insert_ranges:
                if (insert_range[1] - insert_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (insert_range[1] - insert_range[0] + 1) < 3:
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "-"
                    else:
                        # TODO: Talk to san one more time
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "N"

        del_ranges = aln_df.loc[u_range[0]: u_range[1]]
        del_ranges = del_ranges[
            (del_ranges["ref"] != "-")
            & ((del_ranges["F"] == "-") & (del_ranges["R"] == "-"))
        ].index
        if del_ranges.any():
            del_ranges = ranges(del_ranges)
            for del_range in del_ranges:
                if (del_range[1] - del_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (del_range[1] - del_range[0] + 1) < 3:
                        aln_df.loc[
                            del_range[0]: del_range[1], "consensus"
                        ] = aln_df.loc[del_range[0]: del_range[1], "ref"].values
                    else:
                        aln_df.loc[del_range[0]
                            : del_range[1], "consensus"] = "N"
    consesus_range = aln_df[aln_df["consensus"] != "-"].index
    consensus_range = np.min(consesus_range), np.max(consesus_range)

    aln_df.loc[: consensus_range[0] - 1, "consensus"] = aln_df.loc[
        : consensus_range[0] - 1, "ref"
    ]
    aln_df.loc[consensus_range[1] + 1:, "consensus"] = aln_df.loc[
        consensus_range[1] + 1:, "ref"
    ]
    gaps_index = aln_df[aln_df["consensus"] == "-"].index
    gaps_index = gaps_index[(gaps_index >= u_range[0])
                            & (gaps_index <= u_range[1])]
    u_range[-1] -= len(gaps_index)
    aln_df = aln_df[aln_df["consensus"] != "-"]

    return aln_df, u_range


def fasta_map2ref(infile, gap, tmp_fold, n3, idb):
    # TODO: Use some part for arranging the sequence
    """Integrates Sanger fasta to refgene

    :args: infile, outfile, tmp_fold, idb, gap


    """
    sequences = {}
    for rec in SeqIO.parse(infile, "fasta"):
        if infile.endswith(".R.fasta"):  # Generates revese complement
            sequences[rec.id] = rec.seq.reverse_complement()
        else:
            sequences[rec.id] = rec.seq
    cds = True  # TODO: get this information from reference sequence file
    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
        sequences["ref"] = rec.seq
        # if rec.description.split()[1] == "CDS":
        # cds = True

    flb = path.split(infile)[1].split(".")[0]

    aln_df = aln_df_with_ref(sequences, flb, tmp_fold)
    # print(aln_df)
    mapped_index = aln_df[aln_df[flb] != "-"].index
    u_range = useful_range(mapped_index, gap)
    aln_df["concensus"] = "-"
    aln_df.loc[: u_range[0], "concensus"] = aln_df.loc[: u_range[0], "ref"]
    aln_df.loc[u_range[0]: u_range[1], "concensus"] = aln_df.loc[
        u_range[0]: u_range[1], flb
    ]
    aln_df.loc[u_range[1]:, "concensus"] = aln_df.loc[u_range[1]:, "ref"]
    # TODO: Use cds value here
    if idb in ["del", "both"]:
        del_sites = aln_df[aln_df[flb] == "-"].index.values
        if len(del_sites):
            del_ranges = ranges(del_sites)
            del_ranges = [  # Selecting internal deletions
                rng for rng in del_ranges if (rng[0] > u_range[0] & rng[1] < u_range[1])
            ]
            # Del codon is accepted when codon deletion is allowed else deletions are filled with gaps
            for rng in del_ranges:
                if n3:
                    if (rng[1] - rng[0] + 1) % 3 != 0:
                        aln_df.loc[rng[0]: rng[1], "concensus"] = "N"
                        # aln_df.loc[
                        # rng[0]: rng[1], "ref"
                        # ]
                # else: # TODO: Talk to SAN, does he want insert ref, of Ns or leave as gaps
                # aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                # rng[0]: rng[1], "ref"
                # ]

    if idb in ["ins", "both"]:
        ins_sites = aln_df[aln_df["ref"] == "-"].index.values
        if len(ins_sites):
            ins_ranges = ranges(ins_sites)
            ins_ranges = [  # Internal inserts
                rng for rng in ins_ranges if (rng[0] > u_range[0] & rng[1] < u_range[1])
            ]
            for rng in ins_ranges:
                if n3:
                    if (rng[1] - rng[0] + 1) % 3 != 0:
                        aln_df.loc[
                            rng[0]: rng[1], "concensus"
                        ] = aln_df.loc[  # TODO: ask san whether he wants to insert Ns or leave reported nucletides
                            rng[0]: rng[1], "ref"
                        ]
                # else:
                # aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                # rng[0]: rng[1], "ref"
                # ]

    seq = "".join(aln_df["concensus"])  # .replace("-", "N")
    outfile = f"{tmp_fold}/sanger_converted_fasta/{flb}.fasta"

    with open(outfile, "w") as fout:
        fout.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")


def ab1_2seq_map2ref(infiles, gap, tmp_fold):  # , amb, at, res_fold):
    """TODO: Docstring for ab1_2seq.

    :infiles: TODO
    :returns: TODO

    """
    ab1seq_dfs = {}
    tsequences = {}
    flb = path.split(infiles[0])[1].split(".")[0]
    for fl in infiles:
        if fl.endswith(".F.ab1"):
            ab1seq_dfs["F"] = ab1seq(fl)  # TODO: add the condition
            tsequences["F"] = "".join(ab1seq_dfs["F"]["nuc"].values)
        else:
            ab1seq_dfs["R"] = ab1seq(fl)
            tsequences["R"] = "".join(ab1seq_dfs["R"]["nuc"].values)
    # TODO: Keep the ref name as ref
    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
        tsequences[rec.id] = str(rec.seq)

    for k in ab1seq_dfs:
        ab1seq_dfs[k].columns = [f"{k}_{col}" for col in ab1seq_dfs[k].columns]

    aln_with_peak = merge_base_peak(
        aln_df_with_ref(tsequences, flb, tmp_fold), ab1seq_dfs
    )
    # exit()
    # add file names as well, amb, at, res_fold)
    aln_with_peak, u_range = aln_clean(aln_with_peak, gap)
    # aln_with_peak.to_csv(f"{flb}.csv")
    # aln_with_peak.to_csv("testxx.csv")
    # exit(0)
    seq = "".join(list(aln_with_peak["consensus"].values))
    # TODO: Generate sequence and exit
    outfile = path.split(infiles[0])[1].split(".")[0]
    outfile = f"{tmp_fold}/sanger_converted_fasta/{outfile}.fasta"

    output_file = open(outfile, "w")
    output_file.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")
    output_file.close()


def ab2fasta(
    sang_dict, tmp_fold, gap, key, n3, idb  # , bc="neigh"
):  # Base criteria, max, neighbors, mixed # Inputfiles paired and none paired
    # sanger_outputs, tmp_fold, gap
    """Retains fasta and converts ab1 to fasta"""
    # print(key, sang_dict)
    infiles = sang_dict[key]

    if len(infiles) == 1 and infiles[0].endswith(".fasta"):
        fasta_map2ref(infiles[0], gap, tmp_fold, n3, idb)

    else:
        # TODO: Fin a way to usilise only this part to generate fasta
        ab1_2seq_map2ref(infiles, gap, tmp_fold)


def files_and_groups(sanger_files):
    """List fasta and ab1 sequences as dictionary.

    :sanger_files: List of files in the folder
    :returns: dictionary of files as paired/single and fasta/ab1

    """
    file_groups = {}
    for file_ in sanger_files:
        flx = path.split(file_)[1]
        flb = flx.split(".")[0]
        if flb not in file_groups:
            file_groups[flb] = []
        file_groups[flb].append(file_)
    return file_groups


# TODO: generate sequence in absence of reference.
# TODO: Integration if sequence is not from coding region.


def non_overlapping_ids(asseblies, ab1s):
    """Check for ovelapping and non-ovelapping ids and generates csv table

    :asseblies: Fasta assembly containing folder
    :ab1s: Sanger generated ab1 or fasta files
    :returns: Pandas dataframe

    """
    # Assembly IDS
    assembly_ids = []
    for fl in glob(f"{asseblies}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            assembly_ids.append(rec.id)

    # Sanger sequences IDs
    sanger_fasta = []
    for fl in glob(f"{ab1s}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            sanger_fasta.append(rec.id)
    sanger_fasta_missing_assembly = ",".join(
        set(sanger_fasta) - set(assembly_ids))

    # Sanger forward ab1 IDs
    sanger_ab1_f = []
    for fl in glob(f"{ab1s}/*.F.ab1"):
        sanger_ab1_f.append(path.split(fl)[1].split(".F.ab1")[0])
    sanger_ab1_f_missing_assembly = ",".join(
        set(sanger_fasta) - set(sanger_ab1_f))

    # Sanger Reverse ab1 IDs
    sanger_ab1_r = []
    for fl in glob(f"{ab1s}/*.R.ab1"):
        sanger_ab1_r.append(path.split(fl)[1].split(".R.ab1")[0])
    sanger_ab1_r_missing_assembly = ",".join(
        set(sanger_fasta) - set(sanger_ab1_r))

    data_frame = {"assembly": [], "ab1_Forward": [],
                  "ab1_Reverse": [], "fasta": []}
    for assembly_id in assembly_ids:
        data_frame["assembly"].append(assembly_id)

        if assembly_id in sanger_ab1_f:
            data_frame["ab1_Forward"].append(1)
        else:
            data_frame["ab1_Forward"].append(0)

        if assembly_id in sanger_ab1_r:
            data_frame["ab1_Reverse"].append(1)
        else:
            data_frame["ab1_Reverse"].append(0)

        if assembly_id in sanger_fasta:
            data_frame["fasta"].append(1)
        else:
            data_frame["fasta"].append(0)

    deduct = False

    if (
        sanger_ab1_f_missing_assembly
        or sanger_ab1_r_missing_assembly
        or sanger_fasta_missing_assembly
    ):
        deduct = True
        data_frame["assembly"].append("No Assembly")
        data_frame["ab1_Forward"].append(sanger_ab1_f_missing_assembly)
        data_frame["ab1_Reverse"].append(sanger_ab1_r_missing_assembly)
        data_frame["fasta"].append(sanger_fasta_missing_assembly)

    data_frame = pd.DataFrame(data_frame)

    # Check for overlap
    if deduct:
        is_overlap = (
            data_frame.iloc[:-1][["ab1_Forward",
                                  "ab1_Reverse", "fasta"]].sum().sum()
        )
    else:
        is_overlap = data_frame[["ab1_Forward",
                                 "ab1_Reverse", "fasta"]].sum().sum()

    if not is_overlap:
        return pd.DataFrame()
    return data_frame


# TODO: Generate min-max stats for each sequence and provide warnings


def integrate_in_assembly(outputfold, tmp_fold, sample_id):
    """Mergre sange sequences in NGS assemblies
    :outputfold: Final output folder
    :tempfold: Intermediate files generated by other part of the script
    :sample_id: Sample with NGS assembly and sange sequencing

    """
    # TODO: create psl folder in temp folder

    assembly = f"{tmp_fold}/assemblies/{sample_id}.fasta"
    if not path.exists(assembly):
        return
    sanger = f"{tmp_fold}/sanger_converted_fasta/{sample_id}.fasta"
    psl_file = f"{tmp_fold}/tmp/{sample_id}.psl"
    command = ["blat", "-noHead", assembly, sanger, psl_file]
    cmd(command)

    sanger_seq = {}
    sanger_seq_desc = {}
    for rec in SeqIO.parse(sanger, "fasta"):
        sanger_seq_desc[rec.id] = rec.description
        sanger_seq[rec.id] = str(rec.seq)
    org_seq = {}
    for rec in SeqIO.parse(assembly, "fasta"):
        org_seq[rec.id] = str(rec.seq)
    blat_df = pd.read_table(psl_file, header=None)

    blat_df = (
        blat_df[blat_df[9] == blat_df[13]]
        .sort_values(0, ascending=False)
        .drop_duplicates(9)
    )

    for _, row in blat_df.iterrows():
        start, end = list(map(int, sanger_seq_desc[row[9]].split()[1:]))
        qstarts = np.array(list(map(int, row[19][:-1].split(","))))
        tstarts = np.array(list(map(int, row[20][:-1].split(","))))
        block = np.array(list(map(int, row[18][:-1].split(","))))
        qends = qstarts + block
        # tends = tstarts + block
        my_start, my_end = None, None
        for ps in range(row[17]):
            if start >= qstarts[ps] and start <= qends[ps]:
                my_start = tstarts[ps] + start - qstarts[ps]
            if end >= qstarts[ps] and end <= qends[ps]:
                my_end = tstarts[ps] + end - qstarts[ps]
        # tseq = org_seq[row[13]][my_start:my_end]
        # range_gen = re.finditer("N+", tseq)
        # n_range = [match.span() for match in range_gen]
        # print(my_start, my_end, start, end, n_range, "anmol")

        org_seq[row[13]] = (
            org_seq[row[13]][:my_start]
            + sanger_seq[row[9]][start:end]
            + org_seq[row[13]][my_end:]
        )

        # print("Please report at a bug at")
        # print("https://github.com/krisp-kwazulu-natal/" "seqPatcher/issues")
    with open(f"{outputfold}/{sample_id}.fasta", "w") as fout:
        for k in org_seq:
            fout.write(f">{k}\n{org_seq[k]}\n")


@click.command()
@click.option(
    "-s",
    "--sanger-ab1",
    "sa_ab1",
    help="Folder containing Sanger sequencing trace files or Fasta files"
    " generated from trace files.",
    type=str,
    default="ab1",
    show_default=True,
)  # Convert this to folder
# "/home/devil/Documents/San/Corona/Merging/Sanger/12April2021"
# @click.option("-fa", help="Fasta output file.
# If not given, only sequences will be printed in terminal",
#               type=str, default=None, show_default=True)
@click.option(
    "-a",
    "-assemblies-foder",
    "asf",
    help="Folder containing HTS generate incomplete assembly Fasta files",
    type=str,
    default="assemblies",
    show_default=True,
)
@click.option(
    "-o",
    "--out-dir",
    "outd",
    help="Result output Folder",
    type=str,
    default="Results",
    show_default=True,
)
@click.option(
    "-t",
    "--tab",
    "tab",
    help="CSV file for overlapping assemblies and Sanger IDs."
    " If not specified, stdout.",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-O",
    "--output-fasta",
    "ss",
    help="Ouput file name for Fasta from sanger ab1",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-R",
    "--ref-gene-fasta-file",
    "rf",
    help="Refence gene in fasta format",
    type=str,
    default=None,
    show_default=True,
)
# @click.option(
# "-p",
# "--peak-value",
# "pv",
# help="Minmum value for peak. if sequence with peak not covering length of minimum of 50 nucleode length, will return error message",
# default=50,
# type=int,
# show_default=True,
# )
# TODO: ask user if they want accoring to proxymity to end
# TODO: Check the stats of peak decrease over length and generate the math equation for that in case of ambious nucleotide is not utilised
# TODO: Whats happend ins case of ambigious nucleotide is not utilised
# @click.option(
# "-n",
# "--cpu",
# "cpu",
# help="Number of CPU to use",
# type=int,
# default=1,
# show_default=True,
# )
@click.option(
    "-c",
    "--clean-intermediate",
    "ci",
    help="Remove intermediate files",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-g",
    "--gap-allowed",
    "gap",
    help="Minimum gap length between aligned fragment to consider the alignment continuous",
    type=int,
    default=10,
    show_default=True,
)
@click.option(
    "-3",
    "--only-3-nuc",
    "n3",
    help="Allow multiple 3 nucleotide InDels else replace with reference nucleotides or Ns  ",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-x",
    "--indel-selection",
    "idb",
    help="Replace Insertion, Deletion or Both",
    type=click.Choice(["del", "ins", "both"]),
    default="del",
    show_default=True,
    multiple=False,
)
@click.option(
    "-m",
    "--allow-ambi-nuc",
    "amb",
    help="Allow ambigious nucleotide integration, if present in both forward and reverse sequences. Else nucleotide will be calculated.",
    type=bool,
    default=False,
    show_default=True,
)
# @click.option(
# "-M",
# "--ambigious-base-table",
# "at",
# help="Generate table of ambigious nucletide in reads and their replacement"
# " in concestion alongth with the position in consesus",
# type=bool,
# default=False,
# show_default=True,
# )
# TODO: Suggest option to integrate ambigious nucleotide
@click.version_option(__version__)
def run(
    sa_ab1, asf, outd, tab, ss, rf, ci, gap, n3, idb, amb  # , at
):  # , fa, asb, al, bscpu,
    # print(sa_ab1, asf, outd, tab, ss, rf, cpu, ci, gap, n3, idb)
    """
    Reports nucleotide sequence from Sanger chromatogram data based on user
    provided parameters and integrate that in assembly generated using NGS
    data"""
    # TODO: Integrate multi core system
    # if cpu < 1:
    # print("Number of CPU use given is < 1.")
    # print("Using default CPU value of 1")
    # cpu = 1
    # elif cpu > cpu_count() - 1:
    # print("Given cpu usage is more or equal to cpus avalable on system")
    # print(f"Setting CPU usage to {cpu_count() - 1 }")
    # cpu = cpu_count() - 1
    # else:
    # pass

    # pool = Pool(cpu)

    if not sa_ab1:
        exit("Sanger data folder is not given. Exiting . . . .")
    if not path.exists(sa_ab1) or not path.isdir(sa_ab1):
        exit(
            f"Given sanger data folder {sa_ab1} doesn't exist or path is not a folder."
            " Exiting . . . . ."
        )

    if not asf:
        exit("Assembly folder not given. Exiting . . . . . . .")
    if not path.exists(asf) or not path.isdir(asf):
        exit(
            f"Given assembly folder {asf} doesn't exist or path is not a folder."
            " Exiting . . . . ."
        )

    if not rf:
        print(
            "Reference sequence file is not given."
            " Considering sars-cov-2 spike protein sequence as reference"
        )
    elif not path.exists(rf) or not path.isfile(rf):
        print(
            f"Given reference file {rf} doesn't exist or path is not a file.")

        print("Considering sars-cov-2 spike protein sequence as reference")
        rf = None

    # tmp_fold = "tmp"
    tmp_fold = tmf.mkdtemp()

    # ----------Housekeeping----------------
    # Creating temporary  files and folders

    ref_path = f"{tmp_fold}/ref.fasta"
    makedirs(outd, exist_ok=True)
    for folder in [
        "assemblies",
        "sanger_raw",
        "sanger_converted_fasta",
        "sanger_final_fasta",
        "tmp",
    ]:
        makedirs(f"{tmp_fold}/{folder}", exist_ok=True)

    # Copying ref fasta
    # TODO: If ref fasta not given, convert ab1 with them and provide region of overlaps
    with open(ref_path, "w") as fout:
        if not rf:
            fout.write(f">ref CDS 0\n{_s_gene_seq}\n")
        else:
            seq_count = 0
            for rec in SeqIO.parse(rf, "fasta"):
                seq_count += 1
                if seq_count > 1:
                    if ci:
                        rmtree(tmp_fold)
                    exit(
                        f"{rf} contains more than 1 sequence. "
                        "Expect only one."
                        " Exiting."
                    )
                seq_desc = rec.description.split()
                # TODO: Include in future documentation
                seq = rec.seq
                if len(seq_desc) == 3 and seq_desc[1] == "CDS":
                    seq = seq[int(seq_desc[2]):]
                    fout.write(f">{rec.id} CDS 0\n{seq}\n")
                else:
                    fout.write(f">{rec.id}\n{seq}\n")

            if not seq_count:
                if ci:
                    rmtree(tmp_fold)
                exit(
                    f"{rf} contains 0 (zero) sequence. " "Expect only one." " Exiting."
                )

    sanger_files = glob(f"{sa_ab1}/*")
    if not sanger_files:
        if ci:
            rmtree(tmp_fold)
        exit(f"No file found in {sa_ab1} folder. Exiting . . . .")

    sanger_names = []
    for fl in sanger_files:
        if fl.endswith(".fasta"):
            for rec in SeqIO.parse(fl, "fasta"):
                if rec.id not in sanger_names:
                    sanger_names.append(rec.id)
        if fl.endswith(".ab1"):
            flb = path.split(fl)[1].split(".")[0]
            if flb not in sanger_names:
                sanger_names.append(flb)
    # print(sanger_names)
    # exit()

    ss = "Anmol.fasta"
    for fl in glob(f"{sa_ab1}/*"):
        if fl.endswith(".fasta"):  # TODO: Allow fa, faa, fna, and other fasta formats
            for rec in SeqIO.parse(fl, "fasta"):
                with open(f"{tmp_fold}/tmp/{rec.id}.fasta", "w") as fout:
                    fout.write(f">{rec.id}\n{rec.seq}\n")

                l_r = orient(
                    f"{tmp_fold}/tmp/{rec.id}.fasta",
                    ref_path,
                    f"{tmp_fold}/tmp",
                )

                with open(f"{tmp_fold}/sanger_raw/{rec.id}.{l_r}.fasta", "w") as fout:
                    fout.write(f">{rec.id}\n{rec.seq}\n")
        if fl.endswith(".ab1"):
            fl_e = path.split(fl)[1]
            flb = fl_e.split(".")[0]
            l_r = orient(fl, ref_path, f"{tmp_fold}/tmp")
            copyfile(fl, f"{tmp_fold}/sanger_raw/{flb}.{l_r}.ab1")
    sanger_outputs = files_and_groups(glob(f"{tmp_fold}/sanger_raw/*"))
    # print(sanger_outputs)
    for k in sanger_outputs:
        # print(k)
        ab2fasta(sanger_outputs, tmp_fold, gap, k, n3, idb)
        # exit()

    if ss:  # WARNING: Does it matter? If this doesn't work while code be not functional
        with open(ss, "w") as fout:
            for fl in glob(f"{tmp_fold}/sanger_converted_fasta/*"):
                for rec in SeqIO.parse(fl, "fasta"):
                    coors = rec.description.split()[1:]
                    coors = int(coors[0]), int(coors[1])
                    fout.write(f">{rec.id}\n{rec.seq[coors[0]:coors[1]]}\n")

    assembly_files = glob(f"{asf}/*.fasta")
    if not assembly_files:
        if ci:
            rmtree(tmp_fold)
        exit(f"No file found in {asf} folder. Exiting . . . .")
    assembly_names = []

    for fl in assembly_files:
        for rec in SeqIO.parse(fl, "fasta"):
            assembly_names.append(rec.id)
    # else:
    # exit(f"No file is assembly folder {asf}. Exiting . . . .")
    if not assembly_names:
        # TODO: Should run even assemblies are not give and produce fasta from sanger seq
        if ci:
            rmtree(tmp_fold)
        exit("No assembly sequence found. Exiting . . . . .")

    common_ids = set(assembly_names) & set(sanger_names)
    if not common_ids:
        if ci:
            rmtree(tmp_fold)
        exit(
            "Genome assembly and sanger sequencing data doesn't have common"
            " id(s). Exiting..."
        )

    # Copying assembly to tmp folder
    for fl in glob(f"{asf}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            if rec.id in common_ids:
                with open(f"{tmp_fold}/assemblies/{rec.id}.fasta", "w") as fout:
                    fout.write(f">{rec.id}\n{rec.seq}\n")

    seq_id_df = non_overlapping_ids(
        f"{tmp_fold}/assemblies", f"{tmp_fold}/sanger_raw")

    seq_id_df = seq_id_df[
        ~(
            (seq_id_df["assembly"] == "No Assembly")
            | (
                (seq_id_df["ab1_Forward"] == 0)
                & (seq_id_df["ab1_Reverse"] == 0)
                & (seq_id_df["fasta"] == 0)
            )
        )
    ]

    if tab:
        seq_id_df.to_csv(tab, index=False)
    else:
        print(seq_id_df.to_csv(index=False))

    seq_id_df = seq_id_df[
        (
            ((seq_id_df["ab1_Forward"] == 1) | (seq_id_df["ab1_Reverse"] == 1))
            & (seq_id_df["fasta"] == 0)
        )
        | (
            ((seq_id_df["ab1_Forward"] == 0) & (seq_id_df["ab1_Reverse"] == 0))
            & (seq_id_df["fasta"] == 1)
        )
    ]
    print("The patcher executed for ..")
    print(",".join(seq_id_df["assembly"]))

    # assemblies = glob(f"{tmp_fold}/assemblies/*.fasta")

    # sanger_outputs = files_and_groups(glob(f"{tmp_fold}/sanger_raw/*"))
    # for k in sanger_outputs:
    # ab2fasta(sanger_outputs, tmp_fold, gap, k, n3, idb)

    for id_ in sanger_outputs:
        integrate_in_assembly(outd, tmp_fold, id_)

    if ci:
        rmtree(tmp_fold)


if __name__ == "__main__":
    run()
