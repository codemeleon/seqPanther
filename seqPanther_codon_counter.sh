#!/bin/env bash

seqpanther codoncounter \
	-bam examples/codoncounter/N47215.bam \
	-rid sars-cov-2 \
	-ref examples/codoncounter/NC_045512.2.fasta \
	-gff examples/codoncounter/genemap.gff \
	--output_type both \
	-coor_range 21000-25000
