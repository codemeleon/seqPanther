#!/usr/bin/env python3
import click
from seqPanther.NucIn import nuc_in, organise
from seqPanther.seqPatcher import seqpatcher
from seqPanther.CodonCounter import CodonCounter

__author__ = "Anmol Kiran"
__organisation__ = ("University College Cork, Ireland")
__github__ = "codemeleon"
__email__ = "akiran@mlw.mw"
__version__ = "0.0.1"


@click.group()
@click.option("-a", help="first value")
def run(a):
    """Wrapper for commands below"""
    pass


run.add_command(seqpatcher.run, "seqpatcher")
run.add_command(CodonCounter.run, "codoncounter")
run.add_command(nuc_in.run, "nucsubs")
run.add_command(organise.run, "cc2ns")

if __name__ == "__main__":
    run()
