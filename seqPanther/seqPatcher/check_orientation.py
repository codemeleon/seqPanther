from Bio import SeqIO
import click
from os import system
import pandas as pd


@click.command()
@click.option(
    "-r",
    "--ref",
    "ref",
    help="Refernce fasta file",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-a",
    "--ab",
    "ab",
    help="ab1 sanger file.",
    type=str,
    default=None,
    show_default=True,
)
def run(ref, ab):
    """TODO: Docstring for run.

    :ref: TODO
    :ab: TODO
    :returns: TODO

    """
    record = SeqIO.read(ab, "abi")
    seq = "".join(chr(x) for x in record.annotations["abif_raw"]["PBAS1"])
    with open("test.fasta", "w") as fout:
        fout.write(f">query\n{seq}\n")
    system(f"blat -noHead {ref} test.fasta test.psl")
    data = (
        pd.read_table("test.psl", header=None).sort_values(
            0, ascending=False).head(1)
    )
    for _, row in data.iterrows():
        if row[8] == "-":
            strand = "reverse"
        else:
            strand = "forward"
        break

    system("rm test.fasta test.psl")

    print(f"{ab} is on {strand}.")


if __name__ == "__main__":
    run()
