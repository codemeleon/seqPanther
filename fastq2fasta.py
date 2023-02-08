from Bio import SeqIO
import click
from os import path, makedirs


@click.command()
@click.option('-i',
              '--infile',
              'infile',
              type=click.File('r'),
              help='input fastq file',
              required=True,
              default=None,
              show_default=True)
@click.option('-o',
              '--outfile',
              'outfile',
              type=click.File('w'),
              help='output fastq file',
              required=True,
              default=None,
              show_default=True)
def run(infile, outfile):
    """Converts fastq consensus file to fasta file, using file name as sequences, 
    consider fastq file contains only one sequence."""
    flbase = path.splitext(path.basename(infile.name))[0]
    outfold = path.split(outfile.name)[0]
    if outfold:
        makedirs(outfold, exist_ok=True)

    for record in SeqIO.parse(infile, 'fastq'):
        outfile.write(f'>{flbase}\n{record.seq}\n')


if __name__ == "__main__":
    run()
