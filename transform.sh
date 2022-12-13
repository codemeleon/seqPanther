#!/bin/sh
# For transforming data from codoncounter to nucleotide integrator

seqpanther cc2ns -s sub_output.csv \
	-i indel_output.csv \
	-o newout.csv
