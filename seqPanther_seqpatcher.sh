#!/bin/sh

# This script runs seqptcher command on test data.

seqpanther seqpatcher \
	-s examples/seqpatcher/ab1 \
	-a examples/seqpatcher/assemblies \
	-o examples/seqpatcher/results \
	-t examples/seqpatcher/results/mmf.csv \
	-O anmol \
	-c True \
	-g 10 \
	-3 True \
	-x del
