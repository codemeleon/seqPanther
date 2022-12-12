#!/usr/bin/env python3

import click
from Bio import SeqIO
import pandas as pd
from glob import glob
from os import path, makedirs


@click.command(context_settings={'help_option_names': ["-h", "--help"]},
               no_args_is_help=True)
@click.option('-r',
              '--ref',
              "ref",
              help="Fasta file or folder",
              default=None,
              show_default=True,
              type=str,
              required=True)
@click.option(
    '-i',
    '--sid',
    "sid",
    help="Reference sequence id. If not give, first will be selected",
    default=None,
    show_default=True,
    type=str)
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
def run(ref, sid, tabd, outd):
    """
    Integrate changes in nucleotide sequences.\n
    fasta: file for folder containing fasta files name ending with .fasta\n
    tab: file containing changes in format of\n
        SampleA\tC:21301:A,C:23063:C,T:25312:G,G:22188:T\n
    outfile: Output fasta file. All the sequences will be in one file

    """

    # important_ids = parse_and_sort(tab)
    makedirs(outd, exist_ok=True)

    seq = None

    if path.isfile(ref):
        for rec in SeqIO.parse(ref, 'fasta'):
            if not sid:
                seq = rec.seq
                break
            if rec.id == sid:
                seq = rec.seq
                break
    if not seq:
        exit(f"No sequence found in {ref} file")

    for fl in glob(f"{tabd}/*.tsv"):
        df = pd.read_table(fl)

        tdf = df[df["type"] == "sub"]
        for _, row in tdf.iterrows():
            seq = seq[:row['coor']] + \
                row['sub'].split(":")[1]+seq[row['coor']+1:]

        tdf = df[df["type"] == "del"]
        for _, row in tdf.iterrows():
            frag = row['sub'].split(':')[1]
            seq = seq[:row["coor"]] + frag + seq[row['coor'] + len(frag):]

        tdf = df[df["type"] == "ins"].sort_values("coor", ascending=False)
        for _, row in tdf.iterrows():
            frag = row['sub'].split(':')[1]
            seq = seq[:row["coor"]] + frag + seq[row['coor']:]
        flb = path.split(fl)[1].split(".tsv")[0]
        with open(f"{outd}/{flb}.fasta", "w") as fout:
            fout.write(f">{flb}\n{seq}\n")


if __name__ == "__main__":
    run()
    # res = parse_and_sort("changes.txt")
    # print(res)
