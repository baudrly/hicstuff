# hicstuff

[![PyPI version](https://badge.fury.io/py/hicstuff.svg)](https://badge.fury.io/py/hicstuff)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/hicstuff.svg)
[![Read the docs](https://readthedocs.org/projects/hicstuff/badge)](https://hicstuff.readthedocs.io)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://opensource.org/licenses/GPL-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

A lightweight library that generates and handles Hi-C contact maps in either CSV or [instaGRAAL](https://github.com/koszullab/instaGRAAL) format. It is essentially a merge of the [yahcp](https://github.com/baudrly/yahcp) pipeline, the [hicstuff](https://github.com/baudrly/hicstuff) library and extra features illustrated in the [3C tutorial](https://github.com/axelcournac/3C_tutorial) and the [DADE pipeline](https://github.com/scovit/dade), all packaged together for extra convenience.

## Installation

```sh
   pip3 install hicstuff
```

or, for the latest version:

```sh
    pip3 install -e git+https://github.com/koszullab/hicstuff.git@master#egg=hicstuff
```

## Usage

### Full pipeline

All components of the pipelines can be run at once using the `hicstuff pipeline` commands. This allows to generate a contact matrix from reads in a single command. By default, the output sparse matrix is in GRAAL format, but it can be a 2D bedgraph file if required.

    Usage:
        hicstuff pipeline -1 reads_forward.fastq -2 reads_reverse.fastq -f genome.fa [-s size] [-o output_directory] [-e enzyme] [-q quality_min] [--duplicates] [--clean-up]

    Options:
        -1 or --forward: Forward FASTQ reads
        -2 or --reverse: Reverse FASTQ reads
        -f or --fasta: Reference genome to map against in FASTA format
        -o or --output: Output directory. Defaults to the current directory.
        -e or --enzyme: Restriction enzyme if a string, or chunk size (i.e. resolution) if a ber. Defaults to 5000 bp chunks.
        -q or --quality-min: Minimum mapping quality for selecting contacts. Defaults to 30.
        -d or --duplicates: If enabled, trims (10bp) adapters and PCR duplicates prior to ping. Not enabled by default.
        -s or --size: Minimum size threshold to consider contigs. Defaults to 0 (keep all tigs).
        -n or --no-clean-up: If enabled, intermediary BED files will be kept after erating the contact map. Disabled by defaut.
        -p or --pos-matrix: If enabled, generates a sparse matrix with positions (chr,pos) tead of GRAAL-compatible format.
        -t or --threads: Number of threads to use for the aligner and samtools. Defaults to         -T or --tmp: Directory for storing intermediary BED files and temporary sort files. aults to the output directory.
        -m or --minimap: Use the minimap2 aligner instead of bowtie2. Not enabled by default.
        -i or --iterative: Map reads iteratively, by truncating reads to 20bp and then eatedly extending and aligning them.
        -F or --filter: Filter out spurious 3C events (loops and uncuts). Requires -e to be estriction enzyme, not a chunk size.
        -h or --help: Display this help message.

    Output:
         -abs_fragments_contacts_weighted.txt: the sparse contact map
         -fragments_list.txt: information about restriction fragments (or chunks)
         -info_contigs.txt: information about contigs or chromosomes


### Library

All components of the hicstuff program can be used as python modules. See the documentation on [reathedocs](https://hicstuff.readthedocs.io). The expected contact map format for the library is a simple CSV file, and the objects handled by the library are simple ```numpy``` arrays.


### Individual pipeline components

For more advanced usage, different scripts can be used independently on the command line to perform individual parts of the pipeline.

#### Iterative alignment
Truncate reads from a fastq file to 20 basepairs and iteratively extend and re-align the unmapped reads to optimize the proportion of uniquely aligned reads in a 3C library.

    usage:
        hicstuff iteralign [--minimap2] [--threads=1] [--min_len=20] --out_sam=FILE --fasta=FILE <reads.fq>

    arguments:
        reads.fq                Fastq file containing the reads to be aligned

    options:
        -f FILE, --fasta=FILE   Fasta file on which to map the reads.
        -t INT, --threads=INT  Number of parallel threads allocated for the alignment [default: 1].
        -T DIR, --tempdir=DIR  Temporary directory. Defaults to current directory.
        -m, --minimap2     If set, use minimap2 instead of bowtie2 for the alignment.
        -l INT, --min_len=INT  Length to which the reads should be truncated [default: 20].
        -o FILE, --out_sam=FILE Path where the alignment will be written in SAM format.


#### Digestion of the genome

Digests a fasta file into fragments based on a restriction enzyme or a
fixed chunk size. Generates two output files into the target directory
named "info_contigs.txt" and "fragments_list.txt"

    usage:
        digest --enzyme=ENZYME [--size=INT] [--outdir DIR] [--circular] <fasta>

    arguments:
        fasta                     Fasta file to be digested

    options:
        -c, --circular           Whether the genome is circular.
        -e ENZ, --enzyme=ENZ     A restriction enzyme or an integer representing chunk sizes (in bp)
        -s INT, --size=INT       Minimum size threshold to keep fragments [default: 0]
        -o DIR, --outdir=DIR     Directory where the fragments and contigs files will be written

    output:
        -fragments_list.txt: information about restriction fragments (or chunks)
        -info_contigs.txt: information about contigs or chromosomes


#### Filtering of 3C events

Filters spurious 3C events such as loops and uncuts from the library based
on a minimum distance threshold automatically estimated from the library by default. Can also plot 3C library statistics.

    usage:
        hicstuff filter [--interactive | --thresholds INT,INT] [--plot_summary] <input> <output>

    arguments:
        input       2D BED file containing coordinates of Hi-C interacting pairs,
                    the index of their restriction fragment and their strands.
        output      Path to the filtered file, in the same format as the input.

    options:
        -F FILE, --frags=FILE             BED file containing digested genome fragments
        -i, --interactive                 Interactively shows plots and asks for thresholds.
        -t INT-INT, --thresholds=INT-INT  Manually defines integer values for the thresholds in the order [uncut, loop].
        -p, --plot_summary                If set, a plot with library composition informations will be displayed.


#### Viewing the contact map

Visualize a Hi-C matrix file as a heatmap of contact frequencies. Allows to tune visualisation by binning and normalizing the matrix, and to save the output image to disk. If no output is specified, the output is displayed.

    usage:
        hicstuff view [--binning=1] [--normalize] [--max=99] [--output=IMG] <contact_map>

    arguments:
        contact_map             Sparse contact matrix in GRAAL format

    options:
        -b INT, --binning=INT   Subsampling factor to use for binning [default: 1].
        -n, --normalize         Should SCN normalization be performed before rendering the matrix ?
        -m INT, --max=INT       Saturation threshold. Maximum pixel value is set to this percentile [default: 99].
        -o IMG, --output=IMG    Path where the matrix will be stored in PNG format.

### Connecting the modules

All the steps described here are handled automatically when running the `hicstuff pipeline`. But if you want to connect the different modules manually, the intermediary input and output files must be processed using light bash scripting.

#### Extracting contacts from the alignment
The output from iteralign is a SAM file. In order to retrieve Hi-C pairs, you need to run iteralign separately on the two fastq files and process the resulting alignment files processed as follows using bedtools and some bash commands.

1. Convert the SAM files into BED format

```bash
samtools view -bS -F 260 -@ $t -q 30 "for.sam" |
  bedtools bamtobed -i - |
  awk 'OFS="\t" { print $1,$2,$3,$4,$6 }' \
    > contacts_for.bed

samtools view -bS -F 260 -@ $t -q 30 "rev.sam" |
  bedtools bamtobed -i - |
  awk 'OFS="\t" { print $1,$2,$3,$4,$6 }' \
    > contacts_rev.bed
```

2. Put all forward and reverse reads into a single sorted BED file

```bash
sort -k1,1d -k2,2n contacts_for.bed contacts_rev.bed \
  > total_contacts.bed

```

#### Attributing each read to a restriction fragment
To build a a contact matrix, we need to attribute each read to a fragment in the genome. This is done by intersecting all the reads with the digested genome.

1. Extract restriction fragments from the genome using `hicstuff digest`.

```bash
# Generate fragments_list.txt
hicstuff digest --fasta genome.fa \
  --enzyme DpnII \
  --output-dir .

# Make a BED from it
awk 'NR>1 { print $2"\t"$3"\t"$4"\t"(NR-2) }' fragments_list.txt
  >fragments_list.bed

```

2. Intersect the BED files of reads and restriction fragments.

```bash
bedtools intersect -a total_contacts.bed \
  -b fragments_list.bed -wa -wb |
  sort -k4d  \
  > contact_intersect_sorted.bed

```
4. Get reads into a paired BED file (1 pair per line)

```bash
# Note: F is fragment and R is read
# This awk snippet allows to convert a sorted BED file with fields:
#   Rchr Rstart Rend Rname Rstrand Fchr Fstart Fend Fidx
# to a "2D BED" file with:
#   F1chr F1start F1end F1idx R1strand F2chr F2start F2end F2idx R2strand
bed2pairs='
BEGIN{dir="for"; OFS="\t"}
{
  if(dir=="for") {
    fw["name"]=$4; fw["coord"]=$6"\t"$7"\t"$8"\t"$9"\t"$5; dir="rev" }
  else {
    if(fw["name"] == $4) {
        print fw["coord"],$6,$7,$8,$9,$5; dir="for"}
    else {
        dir="rev"; fw["coord"]=$6"\t"$7"\t"$8"\t"$9"\t"$5; fw["name"]=$4}
    }
}
'
awk "$bed2pairs" > contact_intersect_sorted.bed2D
```
The resulting 2D BED file can then be filtered by the `hicstuff filter` module if needed, otherwise, the matrix can be built directly from it. To generate a GRAAL-compatible sparse matrix from the 2D bed file:

```bash
# Remove strand information, sort by fragment combination,
# Count occurrences of each fragment combination and format into csv.
echo -e "id_fragment_a\tid_fragment_b\tn_contact" > matrix.tsv
cut -f4,9 "$tmp_dir/contact_intersect_sorted.bed" |
  sort -V |
  uniq -c |
  sed 's/^\s*//' |
  tr ' ' '\t' |
  awk '{print $0,$1}' |
  cut -f1 --complement >> matrix.tsv
```
