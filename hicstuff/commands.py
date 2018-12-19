# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# cmdoret, 20181412
from docopt import docopt
from hicstuff.hicstuff import bin_sparse, normalize_sparse, bin_kb_sparse
from hicstuff.iteralign import *
from hicstuff.fraglist import write_frag_info, write_sparse_matrix
from hicstuff.filter import get_thresholds, filter_events, process_read_pair
from hicstuff.vizmap import load_raw_matrix, raw_cols_to_sparse, sparse_to_dense, plot_matrix
from re import findall
import sys, os
import shutil as st
import pdb


class AbstractCommand:
    """
    Base class for the commands
    """

    def __init__(self, command_args, global_args):
        """Initialize the commands"""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError


class Iteralign(AbstractCommand):
    """
    Truncate reads from a fastq file to 20 basepairs and iteratively extend and
    re-align the unmapped reads to optimize the proportion of uniquely aligned
    reads in a 3C library.

    usage:
        iteralign [--minimap2] [--threads=1] [--min_len=20] --out_sam=FILE --fasta=FILE <reads.fq>

    arguments:
        reads.fq                Fastq file containing the reads to be aligned

    options:
        -f FILE, --fasta=FILE   Fasta file on which to map the reads.
        -t INT, --threads=INT  Number of parallel threads allocated for the alignment [default: 1].
        -T DIR, --tempdir=DIR  Temporary directory. Defaults to current directory.
        -m, --minimap2     If set, use minimap2 instead of bowtie2 for the alignment.
        -l INT, --min_len=INT  Length to which the reads should be truncated [default: 20].
        -o FILE, --out_sam=FILE Path where the alignment will be written in SAM format.
    """

    def execute(self):
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "."
        if not self.args["--minimap2"]:
            self.args["--minimap2"] = False
        temp_directory = generate_temp_dir(self.args["--tempdir"])
        iterative_align(
            self.args["<reads.fq>"],
            self.args["--tempdir"],
            self.args["--fasta"],
            self.args["--threads"],
            self.args["--out_sam"],
            self.args["--minimap2"],
        )
        # Deletes the temporary folder
        st.rmtree(temp_directory)


class Digest(AbstractCommand):
    """
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

    """

    def execute(self):
        if not self.args["--size"]:
            self.args["--size"] = 0
        if not self.args["--outdir"]:
            self.args["--outdir"] = os.getcwd()
        write_frag_info(
            self.args["<fasta>"],
            self.args["--enzyme"],
            self.args["--size"],
            output_dir=self.args["--outdir"],
            circular=self.args["--circular"],
        )


class Filter(AbstractCommand):
    """
    Filters spurious 3C events such as loops and uncuts from the library based
    on a minimum distance threshold automatically estimated from the library by
    default. Can also plot 3C library statistics.

    usage:
        filter [--interactive | --thresholds INT,INT] [--plot_summary] <input> <output>

    arguments:
        input       2D BED file containing coordinates of Hi-C interacting pairs,
                    the index of their restriction fragment and their strands.
        output      Path to the filtered file, in the same format as the input.

    options:
        -F FILE, --frags=FILE             BED file containing digested genome fragments
        -i, --interactive                 Interactively shows plots and asks for thresholds.
        -t INT-INT, --thresholds=INT-INT  Manually defines integer values for the thresholds in the order [uncut, loop].
        -p, --plot_summary                If set, a plot with library composition informations will be displayed.
    """

    def execute(self):
        output_handle = open(self.args["<output>"], "w")
        if self.args["--thresholds"]:
            # Thresholds supplied by user beforehand
            uncut_thr, loop_thr = self.args["--thresholds"].split("-")
        else:
            # Threshold defined at runtime
            with open(self.args["<input_file>"]) as handle_in:
                uncut_thr, loop_thr = get_thresholds(
                    handle_in, interactive=self.args["--interactive"]
                )
        # Filter library and write to output file
        with open(args.input_file) as handle_in:
            filter_events(
                handle_in,
                output_handle,
                uncut_thr,
                loop_thr,
                plot_events=self.args["--plot_summary"],
            )


class View(AbstractCommand):
    """
    Visualize a Hi-C matrix file as a heatmap of contact frequencies. Allows to
    tune visualisation by binning and normalizing the matrix, and to save the
    output image to disk. If no output is specified, the output is displayed.

    usage:
        view [--binning=1] [--normalize] [--max=99] [--output=IMG] <contact_map>

    arguments:
        contact_map             Sparse contact matrix in GRAAL format

    options:
        -b INT, --binning=INT   Subsampling factor to use for binning [default: 1].
        -n, --normalize         Should SCN normalization be performed before rendering the matrix ?
        -m INT, --max=INT       Saturation threshold. Maximum pixel value is set to this percentile [default: 99].
        -o IMG, --output=IMG    Path where the matrix will be stored in PNG format.
    """

    def execute(self):

        input_map = self.args["<contact_map>"]
        binsuffix = {"bp": 1, "kb": 1000, "Mb": 1_000_000}
        try:
            binning = int(self.args["--binning"])
        except ValueError:
            print("Please provide an integer for binning.", file=sys.stderr)
            raise

        vmax = float(self.args["--max"])

        output_file = self.args["--output"]

        raw_map = load_raw_matrix(input_map)

        sparse_map = raw_cols_to_sparse(raw_map)

        if self.args["--normalize"]:
            sparse_map = normalize_sparse(sparse_map, norm="SCN")

        if binning > 1:
            binned_map = bin_sparse(M=sparse_map, subsampling_factor=binning)
        else:
            binned_map = sparse_map

        try:
            dense_map = sparse_to_dense(binned_map)
            plot_matrix(dense_map, filename=output_file, vmax=vmax)
        except MemoryError:
            print("Contact map is too large to load, try binning more")


class Pipeline(AbstractCommand):
    """
    Entire Pipeline to process fastq files into a Hi-C matrix. Uses all the individual components of hicstuff.

    usage:
        pipeline --fasta FILE --enzyme ENZYME [--size SIZE] [--outdir DIR] <fq1> <fq2>

    arguments:
        fq1:             Forward fastq file
        fq2:             Reverse fastq file

    options:
        -f FILE, --fasta=FILE      Fasta file to be digested
        -e ENZYME, --enzyme=ENZYME A restriction enzyme or an integer representing chunk sizes (in bp)
        -s INT, --size=INT         Minimum size threshold to keep fragments
        -o DIR, --outdir=DIR       Directory where the fragments and contigs files will be written
    """

    # WIP: Parse arguments in docopt instead of Bash ?
    def execute(self):
        pass
