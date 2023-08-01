#!/usr/bin/env python
from typing import List
import tempfile

from os import path, makedirs
from glob import glob
from functools import partial
import click
import numpy as np
import pandas as pd

import pyfaidx
import matplotlib.backends.backend_pdf as bpdf
from pylab import figure, fill_between, scatter, title, xlabel, ylabel, legend, yscale, xlim

from . import auto_cpu, bammer, coors_with_changes, gff_reader

__author__ = "Anmol Kiran, San Emmanuel James"
__organisation__ = (
    "Malawi-Liverpool-Wellcome Trust, Malawi; University of Liverpool, UK")
__github__ = "codemeleon"
__email__ = "akiran@mlw.mw, sanemmanueljames@gmail.com"
__version__ = "0.0.2"


def str2coors(coorstr: str) -> List[List[int]]:
    """Converts comma separated values to coordinates and coordinate ranges."""
    coorslist = [x.strip() for x in coorstr.split(',')]
    try:
        coorrange = []
        for coor in coorslist:
            if '-' in coor:
                start, end = coor.split("-")
                coorrange.append([int(start), int(end)])
                pass
            else:
                coorrange.append([int(coor), int(coor) + 1])
        return coorrange
    except ValueError as err:
        print(err)
        exit("Coordinate accept only , and - as alpha numeric values."
             "check your coordinate input")


@click.command(context_settings={'help_option_names': ["-h", "--help"]},
               no_args_is_help=True)
@click.option(
    "-bam",
    help="Bam files",
    default="test_data/InDel/Insert/SRR17051909.sorted.bam",  # Inserttion
    type=click.Path(exists=True),
    required=True,
    show_default=True,
)
@click.option(
    "-rid",
    help="Reference ID",
    type=str,
    default="NC_045512.2",
    required=True,
    show_default=True,
)
@click.option(
    "-ref",
    help="Reference fasta files",
    type=str,  # click.File("r"),
    default="./test_data/NC_045512.2.fasta",
    required=True,
    show_default=True,
)
@click.option(
    "-coor_range",
    help="Coordinates in the reference, zero index based, end exclusive",
    type=str,
    default=None,  # Forward insert
    show_default=True,
)
@click.option(
    "-gff",
    help="Gff Annotation File",
    type=click.File("r"),
    default="/home/devil/Documents/Tools/BitterBits/src/test_data/genemap.gff",
    required=True,
    show_default=True,
)
@click.option(
    "--ignore_orphans",
    help="Ignore orphaned (Unpaired) reads",
    type=bool,
    default=False,
    show_default=True,
)
@click.option(
    "--min_mapping_quality",
    help="Mapping quality of reads",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--min_base_quality",
    help="Minimum base quality for correct base call",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--ignore_overlaps",
    help="Ignore paired overlapping reads",
    type=bool,
    default=False,
    show_default=True,
)
@click.option("-a",
              "--alt_codon_frac",
              "alt_codon_frac",
              help="Alternative amino acid fraction",
              type=click.FloatRange(0.03, 0.97),
              default=0.03,
              show_default=True)
@click.option(
    "--min_seq_depth",
    help="Minimum sequencing depth at position to be considred",
    type=int,
    default=20,
    show_default=True,
)
@click.option(
    "--max_seq_depth",
    help="Maximum sequencing depth at position to be considred",
    type=int,
    default=1000000,
    show_default=True,
)
@click.option(
    "--alt_nuc_count",
    help="Minimum alternate nucleotide/indel read count fraction",
    type=click.FloatRange(0.003, 0.97),
    default=0.03,
    show_default=True,
)
@click.option(
    "-n",
    "--cpu",
    "cpu",
    help="Number of CPUs to use",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-c",
    "--codoncountfile",
    "codoncountfile",
    help="Ouput codon counting CSV File",
    type=click.File("w"),
    default="codon_output.csv",
    show_default=True,
)
@click.option("-e",
              "--endlen",
              "endlen",
              help="Ingnore mismached around the end of reads",
              type=int,
              default=5,
              show_default=True)
@click.option(
    "-s",
    "--subcountfile",
    "subcountfile",
    help="Ouput  subsubstitution counting CSV File",
    type=click.File("w"),
    default="sub_output.csv",
    show_default=True,
)
@click.option(
    "-i",
    "--indelcountfile",
    "indelcountfile",
    help="Ouput  subsubstitution counting CSV File",
    type=click.File("w"),
    default="indel_output.csv",
    show_default=True,
)
def run(bam: str, rid: str, ref: str, coor_range: str, gff: str,
        ignore_orphans: bool, alt_codon_frac: float, min_mapping_quality: int,
        min_base_quality: int, ignore_overlaps: bool, min_seq_depth: int,
        alt_nuc_count: float, cpu: int, endlen: int, codoncountfile: str,
        subcountfile: str, indelcountfile: str, max_seq_depth: int) -> None:
    """Expected to that bam file is sorted based on coordinate and indexed."""
    for fl in [codoncountfile, subcountfile, indelcountfile]:
        try:
            tp = path.split(fl.name)[0]
            makedirs(tp, exist_ok=True) if tp else None
        except OSError as err:
            exit(str(err))

    gff_data = gff_reader.gff2tab(gff)  # gff to pandas dataframe
    if rid not in gff_data["seq_id"].unique():  # Checking presence of given id
        exit("Reference sequence is not in gff file.\n"
             f"References in the gff file {gff_data['seq_id'].unique()}\n"
             "Exiting")
    else:
        gff_data = gff_data[gff_data["seq_id"] ==
                            rid]  # Filtering gff dataframe

    # reference sequence
    try:
        ref_seq = pyfaidx.Fasta(ref)
    except IOError as err:
        exit(str(err))

    try:
        ref_seq = ref_seq[rid]
    except Exception as err:
        exit(f"Reference id {err} not found in fasta file.")

    # Listing bam files
    t_bam_files = None

    if path.isdir(bam):
        t_bam_files = glob(f"{bam}/*.bam")
        if not t_bam_files:
            exit("No bam files found in the given directory.\n"
                 f"Directory: {bam}\n"
                 "Exiting")
    elif path.isfile(bam) and bam.endswith(".bam"):
        t_bam_files = [bam]
    else:
        exit("Bam file is not in the correct format.\n"
             "Exiting")

    # Sorting and indexing bam files
    tmp_dir = tempfile.mkdtemp()
    bam_files = []
    for bam in t_bam_files:
        sorted_bam = bammer.check_sort_and_index_bam(bam, rid, tmp_dir=tmp_dir)
        if sorted_bam:
            bam_files.append(sorted_bam)
    del t_bam_files

    # NOTE: genomic range
    if not coor_range:
        coor_range = f"1-{len(ref_seq)}"
    coor_ranges = str2coors(coor_range)
    del coor_range

    pool = auto_cpu.cpus(cpu)  # CPU Selection

    # Parameter to select reads
    for start, end in coor_ranges:
        # changes = []
        # TODO: Shift the bottom part and merge in single table

        params = {
            'ref': ref,
            "rid": rid,
            'tmp_dir': tmp_dir,
            "start": start,
            "end": end,
            "gff_data": gff_data,
            "endlen": endlen,
            "ignore_orphans": ignore_orphans,
            "min_mapping_quality": min_mapping_quality,
            "min_seq_depth": min_seq_depth,
            'max_seq_depth': max_seq_depth,
            "min_base_quality": min_base_quality,
            "ignore_overlaps": ignore_overlaps,
            "alt_nuc_count": alt_nuc_count,
            'alt_codon_frac': alt_codon_frac
        }

        changes = partial(coors_with_changes.coor_with_changes_run, params)
        changes = pool.map(changes, bam_files)

        pdf = bpdf.PdfPages("output.pdf")
        all_sub = []
        all_indel = []
        all_codon = []

        for sample, merged_table, depth, subs, indel in changes:

            if depth.empty:
                continue
            all_sub.append(subs)
            all_indel.append(indel)
            all_codon.append(merged_table)
            change_types = pd.DataFrame(
                merged_table.apply(lambda x: x["Nucleotide Change"].split(':'),
                                   axis=1).values.tolist(),
                columns=["coor", "changes"])
            change_types[['from', 'to'
                          ]] = change_types['changes'].str.split('>',
                                                                 expand=True)
            change_types.drop(columns=['changes'], inplace=True)
            # change column type
            change_types['coor'] = change_types['coor'].astype(int)
            change_types = change_types.merge(depth, on="coor", how="inner")
            change_types["type"] = 's'
            change_types.loc[change_types["from"].apply(len) >
                             change_types["to"].apply(len), "type"] = 'd'
            change_types.loc[change_types["from"].apply(len) <
                             change_types["to"].apply(len), "type"] = 'i'

            fig = figure(figsize=(8, 6))
            depth.index = depth.coor
            depth = depth.reindex(
                np.arange(depth.coor.min(),
                          depth.coor.max() + 1)).fillna(0)
            fill_between(depth.index,
                         y1=depth.depth,
                         y2=0,
                         alpha=0.5,
                         color='gray',
                         linewidth=0)
            colors = {"s": "green", "d": "red", "i": "blue"}
            for change_type in colors:
                change = change_types[change_types["type"] == change_type]
                if not change.empty:
                    scatter(change["coor"],
                            change["depth"],
                            color=colors[change_type],
                            label=change_type,
                            alpha=0.4,
                            s=10)

            title(sample)
            xlabel("Position in the reference")
            ylabel("Read coverage")
            legend()

            yscale('log')
            pdf.savefig(fig)

        pdf.close()

        pd.concat(all_codon).to_csv(
            codoncountfile,
            index=False,
            sep="\t" if codoncountfile.name.endswith(".tsv") else ",")
        pd.concat(all_sub).to_csv(
            subcountfile,
            index=False,
            sep="\t" if subcountfile.name.endswith(".tsv") else ",")
        pd.concat(all_indel).to_csv(
            indelcountfile,
            index=False,
            sep="\t" if indelcountfile.name.endswith(".tsv") else ",")


if __name__ == "__main__":
    run()
