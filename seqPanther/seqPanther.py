#!/usr/bin/env python3
import click
from seqPanther.NucIn import organise, nuc_in
from seqPanther.seqPatcher import seqpatcher
from seqPanther.CodonCounter import CodonCounter

__author__ = "Anmol Kiran"
__organisation__ = ("University College Cork, Ireland")
__github__ = "codemeleon"
__email__ = "akiran@mlw.mw"
__version__ = "0.0.1"


@click.group(context_settings={'help_option_names': ["-h", "--help"]},
             no_args_is_help=True)
def run():
    """A toolset for sequence exploration and manipulation."""
    pass


run.add_command(seqpatcher.run, "seqpatcher")
run.add_command(CodonCounter.run, "codoncounter")
run.add_command(nuc_in.run, "nucsubs")
run.add_command(organise.run, "cc2ns")

if __name__ == "__main__":
    run()
