# hicstuff

[![PyPI version](https://badge.fury.io/py/hicstuff.svg)](https://badge.fury.io/py/hicstuff)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/hicstuff.svg)
[![Read the docs](https://readthedocs.org/projects/hicstuff/badge)](https://hicstuff.readthedocs.io)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://opensource.org/licenses/GPL-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

A lightweight library that generates and handles Hi-C contact maps in either CSV or [instaGRAAL](https://github.com/koszullab/instaGRAAL) format.

## Installation

```sh
   pip3 install hicstuff
```

## Usage

### Contact map generation

```sh
    yahcp -1 reads_forward.fastq -2 reads_reverse.fastq -f genome.fa [-s size] [-o output_directory] [-e enzyme] [-q quality_min] [--duplicates] [--clean-up]
```

See [yahcp](https://github.com/baudrly/yahcp) for more information.

### Contact map visualization

```sh
    vizmap -h
```

### Library

See the documentation on [reathedocs](https://hicstuff.readthedocs.io).