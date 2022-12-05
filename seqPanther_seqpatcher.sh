#!/bin/sh

# This script runs seqptcher command on test data.

seqpanther seqptcher \
	-s examples/seqpatcher/ab1 \
	-a examples/seqpatcher/assemblies \
	-o examples/seqpatcher/results \
	-t examples/seqpatcher/results/mmf.csv \
	-R # If gene is known
-O anmol \
	-c True \
	-g 10 \
	-3 True \
	-x del
