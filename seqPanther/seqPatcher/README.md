# SARS-CoV-2 Sequencing: Merge Sanger

The script seqpatcher.py integrates the nucleotides sequenced with the Sanger platform into the corresponding HTS-based (**SARS-CoV-2**) genome assembly.

<!--The script seqpatcher.py integrates the Sanger sequenced segment into the corresponding HTS based assembled incomplete genomes.-->

# Dependencies

## External tools (non-Python)

- [MUSCLE](https://www.drive5.com/muscle/downloads.htm)
  - To perform multiple sequence alignment
- [BLAT](https://hgdownload.soe.ucsc.edu/admin/exe/)
  - To query sequence location in the genome

Both tools can either be downloaded from the links provided, or installed using Anaconda or Ubuntu Linux repositories. Tools must be in the system path and accessible via the terminal.

### Anaconda command

```bash
conda install blat muscle
```

### Ubuntu command

Note: Only system administrators are allowed.

```bash
sudo apt install blat muscle
```

## Python (v3.6+) modules

- pandas==1.2.4
- numpy==1.19.1
- biopython==1.78
- click==7.1.2

### Installation

- Install dependencies with `pip install -r requirements.txt`, download `seqpatcher.py` to your working directory and execute it as instructed in the [Execution](#execution) section.
- Install with `pip install git+https://github.com/krisp-kwazulu-natal/seqPatcher` or `git clone https://github.com/krisp-kwazulu-natal/seqPatcher && cd seqPatcher && python setup.py install`. `seqpatcher.py` should be installed in the system path and accessible from the command-line in the terminal.

## File and sequence naming convention.

1. File name or sample name must not have special characters.
2. Assembly files should be written as <sample_name>.fasta. Sequence's name should be identical to <sample_name>.
3. If assembly files contain multiple sequences, assembly file is generated for each sample ID as explained in point 2.
4. If all the assemblies are in one file, name can be anything with fasta extension. However, the sequence names must match <sample_name>.
5. ab1 file should be named as <sample_name>.\<xyz>.ab1. Where xyz can be any string.
6. If user already has preprocessed Sanger trace files to fasta, it should be named as <sample_name>.\<xyz>.fasta. Where xyz can be any character string.

<a name='cmdoptions'>

# Commandline options

```
Usage: seqpatcher.py [OPTIONS]

  Reports nucleotide sequence from Sanger chromatogram data based on user
  provided parameters and integrate that in assembly generated using NGS
  data

Options:
  -s, --sanger-ab1 TEXT           Folder containing Sanger sequencing trace
                                  files or Fasta files generated from trace
                                  files.  [default: sanger_ab1]

  -a, -assemblies-foder TEXT      Folder containing HTS generate incomplete
                                  assembly Fasta files  [default: assemblies]

  -o, --out-dir TEXT              Result output Folder  [default: Results]
  -t, --tab TEXT                  CSV file for overlapping assemblies and
                                  Sanger IDs. If not specified, stdout.

  -O, --output-fasta TEXT         Ouput file name for Fasta from sanger ab1
  -R, --ref-gene-fasta-file TEXT  Refence gene in fasta format
  -c, --clean-intermediate BOOLEAN
                                  Remove intermediate files  [default: True]
  -g, --gap-allowed INTEGER       Minimum gap length between aligned fragment
                                  to consider the alignment continuous
                                  [default: 10]

  -3, --only-3-nuc BOOLEAN        Allow multiple 3 nucleotide InDels else
                                  replace with reference nucleotides or Ns
                                  [default: True]

  -x, --indel-selection [del|ins|both]
                                  Replace Insertion, Deletion or Both
                                  [default: del]

  -m, --allow-ambi-nuc BOOLEAN    Allow ambigious nucleotide integration, if
                                  present in both forward and reverse
                                  sequences. Else nucleotide will be
                                  calculated.  [default: False]

  --version                       Show the version and exit.
  --help                          Show this message and exit.

```

<!-- -M, --ambigious-base-table BOOLEAN
                                  Generate table of ambigious nucletide in
                                  reads and their replacement in concestion
                                  alongth with the position in consesus
                                  [default: False]--> <!--Include this in next iteration-->

## Analysis requiremenrs

1. HTS assembled incomplete genome.
2. Reference for missing segments that overlap> 50nt on both sides of the missing area in the HTS assembly.
3. Sanger ab1 file (forward, reverse or both), or Fasta sequence generated using Sanger ab1 trace files.
   <!--2. Reference for missing segment overlapping >50nt on either side of missing area in the HTS assembly. Explain in it is it is CDS and reading frame. If information in not included, should be considered as non-coding.--><!--For Future-->
   <a name="execution" />

# Execution

- Execute `python seqpatcher.py --help` for help
- Execute `seqpatcher.py --help` for help, if installed in system path.
- To execute in current folder you must have the following folders
  - ab1: Containing ab1 files single/paired or pre-compiled fasta.
  - assemblies: The folder contains assembly files generated using HTS data. These can be merged in one file.
  - Sanger sequences integrated in assembled genome will be saved in `Results` Folder.
  - Additional files might be generated to user provided path.
    - Multi assembly Fasta file.
    - Multi Fasta for ab1 to fasta converted sequences.

# Base selection

- **Note**: Below cases are valid when InDel are smaller than 10 and not multiple of three
- If average peak value for a given ab1 sequence less than 50, the script displays warnings related to that file.
<!-- If it is not reported in the reference file that it is a coding region, codons will not be considered. Zero based index. Example [[Example_CoVid_S_Gene.fasta](include github link)-->

## Paired ab1 - Preferred

If ambiguious bases are not the same, intersection of nucleotides representing ambigious bases will be considered. In case of no overlap, ambious based represented by union will be considered.

InDels muliple of 3 nucleotides can be permitted else below rule will be applied.

| Ref    | Forward    | Reverse    | Final Outcome                                                                                                                            |
| ------ | ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| -      | Any base/- | Any base/- | -                                                                                                                                        |
| Base A | Base B     | Base B     | Base B                                                                                                                                   |
| Base A | Base A     | Base B     | Based on peak height or ambiguous for both nucleotides                                                                                   |
| Base A | Base B     | Base A     | Based on peak height or ambiguous for both nucleotides                                                                                   |
| Base A | Base B     | -          | Based on neighbouring nucleotide and their peak heights                                                                                  |
| Base A | Base A     | -          | Based on neighbouring nucleotide and their peak heights                                                                                  |
| Base A | -          | -          | Based on neighbouring nucleotide and their peak heights                                                                                  |
| Base A | Base B     | Base Ambi  | Overlapping nucleotide between forward and reverse if any, else ambiguous nucleotide representing all at the position will be considered |
| Base A | Base Ambi  | Base Ambi  | Overlapping nucleotide between forward and reverse if any, else ambiguous nucleotide representing all at the position will be considered |
| Base A | Base B     | Base C     | Overlapping nucleotide between forward and reverse if any, else ambiguous nucleotide representing all at the position will be considered |

## Single ab1 - Not Preferred

| Reference | Forward/Reverse | Selected                                                          |
| --------- | --------------- | ----------------------------------------------------------------- |
| -         | Base A          | Based on neighbouring nucleotide and their peak heights           |
| Base A    | -               | Based on neighbouring nucleotide and their peak heights           |
| Base A    | Base B          | Base B                                                            |
| Base A    | Ambi            | highest peak, or ambiguous nucleotide based on use provided input |

## Single Fasta - Preferred

| Reference | FastaNucleotides | Selected                                                                                  |
| --------- | ---------------- | ----------------------------------------------------------------------------------------- |
| -         | Base A           | Base A                                                                                    |
| Base A    | -                | Base A                                                                                    |
| Base A    | Base B           | Base B                                                                                    |
| Base A    | Ambi             | Ambiguous or common base between the sequence and reference based on user provided option |

## Additional Features

1. Generate table of overlapping IDs of assemblies and sanger sequence data.
2. Converts sanger ab1 to fasta based on user provided criteria.

# License

GPLv3
