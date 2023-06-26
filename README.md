# SeqPanther

`SeqPanther` is a Python application that provides the user with a suite of tools to further interrogate the circumstance under which certain mutations occur or are missed and to enables the user to modify the consensus as needed. The tool is applicable to non-segmented bacterial and viral genomes where reads are mapped to a reference sequence. `SeqPanther` generates detailed reports of mutations identified within a genomic segment or positions of interest, including visualization of the genome coverage and depth. The tool is particularly useful in the examination of multiple next-generation sequencing (NGS) short-read samples. Additionally, we have integrated `Seqpatcher` [Singh et. al.](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000774) which supports the merging of Sanger sequences, or consensus thereof into their respective NGS consensus.

![Seqpanther workflow, genome coverage plot and codon count report](https://github.com/codemeleon/seqPanther/blob/master/manuscript/banner.jpg)

SeqPanther consists of the following set of commands:

- **codoncounter**: performs variant calling and generates nucleotide stats at variant sites and reports the impacts of nucleotide changes on amino acids in the translated proteins.

- **cc2ns** (codon counter 2 nucleotide substitution): performs transformation of codon counter output in a format where the user can select a variant to integrate into the reference or assemblies using the **cc2ns** command. The user is encourage to inspect and edit the output from **codoncounter** to only keep the desired changes before running the **cc2ns** command.

- **nucsubs**: accepts as input the changes file generated by **cc2ns** and modifies the consensus sequences accordingly i.e. adds or changes the mutations as specified in the file. <!--Providing options to users to select changes of their interests-->

- **seqpatcher**: integrates sanger sequencing of missing regions of an incomplete assembly to the assembly. This command is a modification of the [SeqPatcher tool](https://github.com/krisp-kwazulu-natal/seqPatcher).

## Operating system compatibility

Unix and OS X Commandline Application. 

## Dependencies

The tool relies on multiple external open source programs and python modules as listed below:

### External Tools

- python>=3.7

1. Samtools

   - For sorting bam files and indexing them.
   - `conda install -c bioconda samtools`

2. Bcftools

   - For parsing variant calling and generating consensus sequences.
   - `conda install -c bioconda bcftools` to install.

3. Muscle (v. 3.8.31)

   - To perform multiple sequence alignment.
   - `conda install -c bioconda muscle=3.8.31` to install.

4. BLAT

   - To query sequence location in the genome.
   - `conda install -c bioconda ucsc-blat` to install.

5. MAFFT

   - To align consensus sequence against the reference.
   - `conda install -c bioconda mafft` to install.

# Installation

## Option 1: Clone repo and install locally

1. `git clone https://github.com/codemeleon/seqPanther`
2. `cd seqPanther`
3. `pip install .`

## Option 2: Install directly from Git

To install directly from the Github repo, run the command: 

`pip install git+https://github.com/codemeleon/seqPanther.git`

# Usages

seqPanther contains four commands. The commands are listed in the help menu. To view the help menu, type: `seqpanther --help` or just `seqpanther`.

## codoncounter

This command help is accessible using `seqpanther codoncounter` or `seqpanther codoncounter --help`.

## cc2ns

This command help is accessible using `seqpanther cc2ns` or `seqpanther cc2ns --help`.

## nucsubs

This command help is accessible using `seqpanther nucsubs` or `seqpanther nucsubs --help`.

## SeqPatcher

This command help is accessible at `seqpanther seqpatcher` or `seqpanther seqpatcher --help`.

# Example

1. A sample dataset of five BAM files from South African SARS-CoV-2 samples can be downloaded from [SeqPanther zenodo](https://zenodo.org/record/7601513) page. The BAM files were generated by mapping reads mapped against reference (NC_045512.2) using [BWA](https://github.com/lh3/bwa). BWA was also used to generate the SAM files which were converted to BAM and sorted using samtools.

2. Reference sequences in Fasta format and genome annotation GFF files were downloaded from [NCBI Genome Fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz) and [ NCBI Genome GFF](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz) respectively. The downloaded files need to un-compressed.

3. Consensus sequences can be calculated from the BAM files using the command `bcftools mpileup -f GCF_009858895.2_ASM985889v3_genomic.fna bam/K032258-consensus_alignment_sorted.bam| bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > K032258-consensus_alignment_sorted.fastq` and `python fastq2fasta.py -i K032258-consensus_alignment_sorted.fastq -o consensus/K032258-consensus_alignment_sorted.fasta `. `fastq2fasta.py` can be downloaded from the project repo.

4. Store the downloaded BAM files in a folder, e.g. called `bam`.

5. To generate a summary of the codons, nucleotides and indel changes in the Spike gene of these samples, use the command `seqpanther codoncounter -bam bam -rid NC_045512.2 -ref GCF_009858895.2_ASM985889v3_genomic.fna -gff GCF_009858895.2_ASM985889v3_genomic.gff -coor_range 21563-25384`. This will generate the summaries for all the files in bam folder.

6. If you only want to generate the results for a single BAM file, run the command as `seqpanther codoncounter -bam ./bam/K032282-consensus_alignment_sorted.bam -rid NC_045512.2 -ref GCF_009858895.2_ASM985889v3_genomic.fna -gff GCF_009858895.2_ASM985889v3_genomic.gff -coor_range 21563-25384` replacing the BAM file name with your specific bam file name in the command. 

The command will generate four outputs in the current folder including: `sub_output.csv` containing details of the nucleotide substitutions, `indel_output.csv` containing details of the indel events, `codon_output.csv` containing details of the codon changes and `output.pdf` which is a plot of genome depth and breadth of coverage annotated with the positions with mutations and indels.

7. Outputs can be explored using a text file reader (for the text files) and pdf reader (e.g Adobe Reader) for the PDFs. An example command to view the text files would be: `cat sub_output.csv | sed 's/,/ ,/g' | column -t -s, | less -S`. The user needs to explore those files and remove the changes they would like not to be integrated. A text editor of your choice e.g. bbedit or notepad++ can be used to edit the files.

8. In case you decide that there are certain mutations that you need to change, you will have to convert the outputs from `codoncounter` to the format required by the `nucsubs` command and run the command `seqpanther cc2ns -s sub_output.csv -i sub_output.csv -o changes`. It generates a CSV file for each sample in the `./change` folder. 

9. Then execute seqpanther as follows: `seqpanther nucsubs -i NC_045512.2 -r NC_045512.2.fasta -c consensus -t changes -o results` to integrate relevant changes to the consensus sequences. The output will be generated in a folder named `results`.

10. Example data for `seqpanther seqpatcher` is provided in example folder of this project.

11. To run `seqpanther seqpatcher` on the example data, use `seqpanther seqpatcher -s examples/seqpatcher/ab1 -a examples/seqpatcher/assemblies -o Results -t mmf.csv -O sanger.fasta -g 10 -3 True -x del`. It will generate consensus sequences from the Sanger trace files and integrate them into the ngs consensus. The mmf.csv file tells the application whether the sanger data were paired or single ab1, or in single fasta, and sanger.fasta containing all sanger sequences.

# Warning

Recursive use of newly generated consensus sequences might result in incorrect final sequence due to the integration of indel events.

# Bug reporting

To report a bug, request support or propose a new feature, please open an issue.

# Licence

GNU GPL v3.0

# Citation

If you use this software please cite:

SeqPanther: Sequence manipulation and mutation statistics toolset. James Emmanuel San, Stephanie Van Wyk, Houriiyah Tegally, Simeon Eche, Eduan Wilkinson, Aquillah M. Kanzi, Tulio de Oliveira, Anmol M. Kiran. bioRxiv 2023. doi: [10.1101/2023.01.26.525629](https://doi.org/10.1101/2023.01.26.525629)
