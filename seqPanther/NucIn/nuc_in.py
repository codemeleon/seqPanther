#!/usr/bin/env python3

import click
from Bio import SeqIO
import pandas as pd
from glob import glob
from os import path, makedirs, system
import itertools
from pandas.core.indexes.range import RangeIndex
from typing import Generator, Tuple

from tempfile import mkdtemp


def ranges(i: RangeIndex) -> Generator[Tuple[int, int], None, None]:
    for _, b in itertools.groupby(enumerate(i),
                                  lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


@click.command(context_settings={'help_option_names': ["-h", "--help"]},
               no_args_is_help=True)
@click.option('-r',
              '--ref',
              "ref",
              help="Reference fasta file",
              default=None,
              show_default=True,
              type=str,
              required=True)
@click.option(
    '-i',
    '--rid',
    "rid",
    help="Reference sequence id. If not give, first will be selected",
    default='ref',
    show_default=True,
    type=str)
@click.option("-c",
              "--con",
              "consensus",
              help="consenus fasta file or folder containing fasta files",
              default=None,
              show_default=True,
              required=True)
@click.option('--tabd',
              '-t',
              "tabd",
              help="Nucletide substitution table folder",
              default=None,
              show_default=True,
              type=str,
              required=True)
@click.option('--outd',
              '-o',
              "outd",
              help="Output directory",
              default=".",
              show_default=True,
              type=str,
              required=True)
@click.option(
    '--merge',
    '-m',
    "merge",
    help="Output directory",
    default=False,
    show_default=True,
    type=bool,
)
def run(ref: str, rid: str, tabd: str, outd: str, consensus: str,
        merge: bool) -> None:
    """
    Integrate changes in nucleotide sequences.\n
    fasta: file for folder containing fasta files name ending with .fasta\n
    tab: file containing changes in format of\n
        SampleA\tC:21301:A,C:23063:C,T:25312:G,G:22188:T\n
    outfile: Output fasta file. All the sequences will be in one file

    """

    # NOTE: Expecting the /tmp folder is writable
    temp_dir = mkdtemp(prefix='nucin')
    makedirs(outd, exist_ok=True)
    if not path.exists(consensus):
        exit(f"{consensus} doesn't exists.")

    if ref:
        if not path.exists(ref):
            exit(f"{ref} path doesn't exist.")
        if not path.isfile(ref):
            exit(f"{ref} is not a file")
        try:
            ref_seqs = SeqIO.parse(ref, "fasta")
        except Exception as e:
            print(e)
            exit(f"{ref} is not fasta file.")
        else:
            ref_fasta = None
            for rec in ref_seqs:
                if rec.id == rid:
                    ref_fasta = rec.seq
                    break
            if not ref_fasta:
                exit(f"Given reference '{rid}' is not in fasta file '{ref}'")

    else:
        exit("Reference file path not given.")
    sample_seq = {}
    is_fold = False
    if path.isfile(consensus):
        flb = path.split(consensus)[-1].split(".")[0]
        try:
            for rec in SeqIO.parse(consensus, "fasta"):
                sample_seq[flb] = rec.seq
        except Exception as e:
            print(e)
            exit(f"Consenus file {consensus} is not fasta.")
    else:
        is_fold = True
        for ss in glob(f"{consensus}/*.fasta"):
            flb = ss.split("/")[-1].split(".")[0]
            try:
                for rec in SeqIO.parse(ss, "fasta"):
                    sample_seq[flb] = rec.seq
                    break
            except Exception as e:
                print(e)
                print(f"Consenus file {ss} is not fasta. ignoring...")

    for fl in glob(f"{tabd}/*.tsv"):
        df = pd.read_table(fl)
        df['coor'] -= 1  # NOTE: Cor 0-based index
        for samp in df["Sample"].unique():
            with open(f"{temp_dir}/test.fasta", "w") as fout:
                fout.write(f">ref\n{ref_fasta}\n>query\n{sample_seq[samp]}\n")
            if path.isfile(f"{temp_dir}/test_mafft.fasta"):
                system(f"rm {temp_dir}/test_mafft.fasta")
            system(
                f"mafft --auto {temp_dir}/test.fasta > {temp_dir}/test_mafft.fasta 2>/dev/null"
            )
            seq_df = {}
            for rec in SeqIO.parse(f"{temp_dir}/test_mafft.fasta", "fasta"):
                seq_df[rec.id] = list(rec.seq)
            seq_df = pd.DataFrame(seq_df)
            seq_df["pos"] = seq_df.index
            for rng in ranges(seq_df[seq_df["ref"] == "-"].index):
                seq_df.loc[rng[0]:rng[1],
                           "pos"] = -1  # NOTE:  to remember gaps
                seq_df.loc[seq_df["pos"] > rng[1],
                           "pos"] -= (rng[1] - rng[0] + 1)

            samp_df = df[df["Sample"] == samp]

            tdf = samp_df[samp_df["type"] == "sub"]
            for row in tdf.itertuples():
                seq_df.loc[seq_df["pos"] == row.coor,
                           "query"] = row.sub.split(":")[1]

            tdf = samp_df[samp_df["type"] == "del"]
            for row in tdf.itertuples():
                frag = row.sub.split(':')[1]
                seq_df.loc[(seq_df["pos"] >= row.coor) &
                           (seq_df["pos"] < row.coor + len(frag)),
                           "query"] = '-'

            tdf = samp_df[samp_df["type"] == "ins"].sort_values(
                "coor", ascending=False)
            for row in tdf.itertuples():
                frag = row.sub.split(':')[1]
                seq_df.loc[seq_df["pos"] == row.coor,
                           "query"] += frag  # TODO: Merge the frag here
            seq = "".join(seq_df["query"].values).replace("-", "")
            sample_seq[samp] = seq
    if merge:

        if not is_fold:
            with open(f"{outd}/{path.split(consensus)[1]}", "w") as fout:
                for sid in sample_seq:
                    fout.write(f">{sid}\n{sample_seq[sid]}\n")
        else:
            with open(f"{outd}/{path.split(consensus)[1]}_merged.fasta",
                      "w") as fout:
                for sid in sample_seq:
                    fout.write(f">{sid}\n{sample_seq[sid]}\n")
    else:
        for sid in sample_seq:
            with open(f'{outd}/{sid}.fasta', 'w') as fout:
                fout.write(f">{sid}\n{sample_seq[sid]}\n")


if __name__ == "__main__":

    # def run(ref, rid, tabd, outd, consensus):
    # run("../../examples/codoncounter/NC_045512.2.fasta", "test",
    # "../../examples/nucsubs/changes", "../../resulst",
    # "../../examples/nucsubs/consensus", True)
    run()
