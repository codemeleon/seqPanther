# python ./seqPanther/NucIn/nuc_in.py \
pip uninstall -y seqPanther
python setup.py clean
python setup.py build
python setup.py install

seqpanther nucsubs \
	-i NC_045512.2 \
	-r examples/codoncounter/NC_045512.2.fasta \
	-c examples/nucsubs/consensus \
	-t examples/nucsubs/changes \
	-o examples/nucsubs/results
