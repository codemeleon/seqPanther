#!/bin/sh
# For transforming data from codoncounter to nucleotide integrator
pip uninstall -y seqPanther
python setup.py clean
python setup.py build
python setup.py install

seqpanther cc2ns \
	-s sub_output.csv \
	-i indel_output.csv \
	-o ./examples/nucsubs/changes
#python seqPanther/NucIn/organise.py \
