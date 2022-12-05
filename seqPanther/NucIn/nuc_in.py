#!/usr/bin/env python3

import click
from Bio import SeqIO
import pandas as pd
from glob import glob
from os import path


def parse_and_sort(coor_and_changes):
    """
    Parse the changes for each sequence id and sort them by position in decreasing order.
    example input: SampleA	C:21301:A,C:23063:C,T:25312:G,G:22188:T
    """
    to_return = {}
    with open(coor_and_changes) as f:
        for line in f:
            dt1 = line[:-1].split('\t')
            changes = dt1[1].split(",")
            changes_with_coor = {"coor": [], "from": [], "to": []}
            for change in changes:
                dt2 = change.split(":")
                changes_with_coor["coor"].append(int(dt2[1]))
                changes_with_coor["from"].append(dt2[0])
                changes_with_coor["to"].append(dt2[2])
            changes_with_coor = pd.DataFrame(changes_with_coor).sort_values(
                by="coor", ascending=False)

            to_return[dt1[0]] = changes_with_coor
    return to_return


@click.command()
@click.option('--fasta',
              '-f',
              help="Fasta file or folder",
              type=click.File('r'),
              required=True)
@click.option('--tab',
              '-t',
              help="Nucletide substitution table",
              type=click.File('r'),
              required=True)
@click.option('--outfile',
              '-o',
              help="Output fasta file",
              type=click.File('w'),
              default="output.fasta",
              required=True)
def run(fasta, tab, outfile):
    """
    Integrate changes in nucleotide sequences.\n
    fasta: file for folder containing fasta files name ending with .fasta\n
    tab: file containing changes in format of\n
        SampleA\tC:21301:A,C:23063:C,T:25312:G,G:22188:T\n
    outfile: Output fasta file. All the sequences will be in one file

    """

    important_ids = parse_and_sort(tab)

    sequences = {}

    if path.isfile(fasta):
        for rec in SeqIO.parse(fasta, 'fasta'):
            if rec.id in important_ids:
                seq = rec.seq
                for _, row in important_ids[rec.id].iterrows():
                    seq = seq[:row["coor"]] + row["to"] + seq[
                        row["coor"] + 1:]  # TODO: check if this is correct
                sequences[rec.id] = seq
            else:

                sequences[rec.id] = list(rec.seq)
    else:
        for fl in glob("*.fasta"):
            for rec in SeqIO.parse(fl, 'fasta'):
                if rec.id in important_ids:
                    seq = rec.seq
                    for _, row in important_ids[rec.id].iterrows():
                        seq = seq[:row["coor"]] + row["to"] + seq[
                            row["coor"] + 1:]  # TODO: check if this is correct
                else:
                    sequences[rec.id] = list(rec.seq)
    for sid in sequences:
        outfile.write(f">{sid}\n{sequences[sid]}\n")
    # TODO: Merge the details


if __name__ == "__main__":
    run()
    # res = parse_and_sort("changes.txt")
    # print(res)
