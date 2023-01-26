#!/bin/env bash
pip uninstall -y seqPanther
python setup.py clean
python setup.py build
python setup.py install

#python /seqPanther/CodonCounter/CodonCounter.py \
seqpanther codoncounter \
	-bam ./examples/codoncounter/ \
	-rid NC_045512.2 \
	-ref examples/codoncounter/NC_045512.2.fasta \
	-gff examples/codoncounter/genemap.gff \
	-coor_range 21000-25000 \
	-n 5 \
	-e 0
