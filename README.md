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

### Fastq to contact map pipeline

    Usage:
        yahcp -1 reads_forward.fastq -2 reads_reverse.fastq -f genome.fa [-s size] [-o output_directory] [-e enzyme] [-q quality_min] [--duplicates] [--clean-up]

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


### Contact map visualization

    Usage:
        vizmap <contact_map> [--binning=1] [--normalize] [--output=<output.png>]
                                [--max=99]

    Options:
        -h, --help              Display this help message.
        --version               Display the program's current version.
        -b 1, --binning 1       Subsampling factor for binning. [default: 1]
        -N, --normalize         Normalize the contact map prior to binning.
                                [default: False]
        -o file, --output file  Save image to output instead of plotting it.
                                [default: output.png]
          -m 99, --max 99         Saturation percentile threshold in the contact map.
                                  [default: 99]

  ### Iterative mapping

  ### Filtering 3C events

    usage: filter3C.py [-h] [-i | -t THRESHOLDS THRESHOLDS] -F FRAGS [-p]
                       input_file [output_file]

    positional arguments:
      input_file            The 2D BED file containing the coordinates of Hi-C
                            interacting pairs and, the indices of their
                            restriction fragments and their strands.
      output_file           Path to the output file (filtered dat.indices).
                            Defaults to stdout.

    optional arguments:
      -h, --help            show this help message and exit
      -i, --interactive     Interactively shows plots and asks the user for
                            thresholds. Overrides predefined thresholds. Disabled
                            by default.
      -t THRESHOLDS THRESHOLDS, --thresholds THRESHOLDS THRESHOLDS
                            The minimum number of restriction fragments between
                          reads to consider loops and uncut events. Estimated
                          automatically by default. Must be two integers
                          separated by a space (-t <loop> <uncut>).
    -F FRAGS, --frags FRAGS
                          Path to bed file containing the sorted list of
                          fragments, with fields: chr start end.
    -p, --plot_summary    Output a piechart summarizing library composition.


### Library

See the documentation on [reathedocs](https://hicstuff.readthedocs.io). The expected contact map format for the library is a simple CSV file, and the objects handled by the library are simple ```numpy``` arrays.


### Individual components

Different scripts can be used independently to perform individual parts of the pipeline.

#### Event filtering
