# -*- coding: utf-8 -*-
"""
Created on Jan 26 2018
@author: Remi Montagne & cmdoret
Aligns iteratively reads from a 3C fastq file
"""

import argparse
import os
import subprocess as sp
import pysam as ps
import shutil as st
from random import getrandbits
import compressed_utils as ct
import contextlib

##############################
#        ARGUMENTS
##############################


def parse_args():
    """ Gets the arguments from the command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fq", help="The fastq file containing Hi-C reads.")
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Path to the reference genome, in FASTA format.",
    )
    parser.add_argument(
        "-p",
        "--nb_processors",
        default=1,
        type=int,
        help="number of CPUs used for alignment.",
    )
    parser.add_argument(
        "-o",
        "--out_sam",
        help="Path to the output SAM file for the alignment of in_fq.",
    )
    parser.add_argument(
        "-T",
        "--tempdir",
        default=".",
        help="Directory to write temporary files. Defaults to current directory.",
    )
    parser.add_argument(
        "-m",
        "--minimap2",
        default=False,
        action="store_true",
        help="Use minimap2 instead of bowtie for the alignment.",
    )
    parser.add_argument(
        "-l",
        "--min_len",
        type=int,
        default=20,
        help="Minimum length to which reads should be truncated.",
    )
    return parser.parse_args()


##############################
#        FUNCTIONS
##############################
def generate_temp_dir(path):
    """
    Generates a temporary file with a random name at the input path.
    Parameters
    ----------
    path : str
        The path at which the temporary directory will be created.
    Returns
    -------
    str
        The path of the newly created temporary directory.
    """
    exist = True
    while exist:
        # Keep trying random directory names if they already exist
        directory = str(hex(getrandbits(32)))[2:]
        full_path = os.path.join(path, directory)
        if not os.path.exists(full_path):
            exist = False
    try:
        os.makedirs(full_path)
    except PermissionError:
        raise PermissionError(
            "The temporary directory cannot be created in {}. "
            "Make sure you have write permission.".format(path)
        )
    return full_path


def iterative_align(
    fq_in, tmp_dir, ref, n_cpu, sam_out, minimap2=False, min_len=20
):
    """
    Aligns reads iteratively reads of fq_in with bowtie2 or minimap2. Reads are
    truncated to the 20 first nucleotides and unmapped reads are extended by 20
    nucleotides and realigned on each iteration.
    Parameters
    ----------
    fq_in : str
        Path to input fastq file to align iteratively.
    tmp_dir : str
        Path where temporary files should be written.
    ref : str
        Path to the reference genome.
    n_cpu : int
        The number of CPUs to use for the iterative alignment.
    sam_out : str
        Path where the final alignment should be written in SAM format.
    minimap2 : bool
        If True, use minimap2 instead of bowtie2 for the alignment.
    """
    # initial length of the fragments to align
    n = min_len
    # set with the name of the unaligned reads :
    remaining_reads = set()
    total_reads = 0
    # Store path of SAM containing aligned reads at each iteration.
    iter_out = []

    # If there is already a file with the same name as the output file,
    # remove it. Otherwise, ignore.
    with contextlib.suppress(FileNotFoundError):
        os.remove(sam_out)

    # Bowtie only accepts uncompressed fastq: uncompress it into a temp file
    if not minimap2 and ct.is_compressed(fq_in):
        uncomp_path = os.path.join(tmp_dir, os.path.basename(fq_in) + ".tmp")
        with ct.read_compressed(fq_in) as inf:
            with open(uncomp_path, "w") as uncomp:
                st.copyfileobj(inf, uncomp)
    else:
        uncomp_path = fq_in

    # Index genome if using bowtie2 and index does not exist
    index = os.path.splitext(ref)[0]
    if not minimap2 and not os.path.isfile(index):
        cmd = "bowtie2-build {0} {1}".format(ref, index)
        sp.call(cmd, shell=True)

    # Counting reads
    with ct.read_compressed(uncomp_path) as inf:
        for line in inf:
            total_reads += 1
    total_reads /= 4

    # Use first read to guess read length.
    with ct.read_compressed(uncomp_path) as inf:
        size = inf.readline()
        # Stripping newline.
        size = len(inf.readline().rstrip())

    print("{0} reads to parse".format(total_reads))

    # iterative alignment per se
    while n <= size:
        print("\n" + "-" * 10 + "\nn = {0}".format(n))
        iter_out += [os.path.join(tmp_dir, "trunc_{0}.sam".format(str(n)))]
        # Generate a temporary input fastq file with the n first nucleotids
        # of the reads.
        print("Generating truncated reads")
        truncated_reads = truncate_reads(
            tmp_dir, uncomp_path, remaining_reads, n, min_len
        )

        # Align the truncated reads on reference genome
        print("Aligning reads")
        temp_alignment = "{0}/temp_alignment.sam".format(tmp_dir)
        map_args = {
            "fa": ref,
            "threads": n_cpu,
            "sam": temp_alignment,
            "fq": truncated_reads,
            "idx": index,
        }
        if minimap2:
            cmd = "minimap2 -x sr -a -t {threads} {fa} {fq} > {sam}".format(
                **map_args
            )
        else:
            cmd = "bowtie2 -x {idx} -p {threads} --rdg 500,3 --rfg 500,3 --quiet --very-sensitive -S {sam} {fq}".format(
                **map_args
            )
        sp.call(cmd, shell=True)

        # filter the reads: the reads whose truncated end was aligned are written
        # to the output file.
        # The reads whose truncated end was not aligned are kept for the next round.
        print("Reporting aligned reads")
        remaining_reads = filter_samfile(temp_alignment, iter_out[-1])

        n += 20

    # one last round without trimming
    print("\n" + "-" * 10 + "\nn = {0}".format(size))
    print("Generating truncated reads")
    truncated_reads = truncate_reads(
        tmp_dir, uncomp_path, remaining_reads, size, min_len
    )
    print("Aligning reads")
    if minimap2:
        cmd = "minimap2 -x sr -a -t {1} {0} {3} > {2}".format(
            ref, n_cpu, temp_alignment, truncated_reads
        )
    else:
        cmd = "bowtie2 -x {0} -p {1} --rdg 500,3 --rfg 500,3 --quiet --very-sensitive -S {2} {3}".format(
            index, n_cpu, temp_alignment, truncated_reads
        )
    sp.call(cmd, shell=True)
    print("Reporting aligned reads")
    iter_out += [os.path.join(tmp_dir, "trunc_{0}.sam".format(str(n)))]
    remaining_reads = filter_samfile(temp_alignment, iter_out[-1])
    n_remaining = len(remaining_reads)

    # Report unaligned reads as well
    iter_out += [os.path.join(tmp_dir, "unaligned.sam")]
    temp_sam = ps.AlignmentFile(temp_alignment, "r")
    unmapped = ps.AlignmentFile(iter_out[-1], "w", template=temp_sam)
    for r in temp_sam:
        # Do not write supplementary alignments (keeping 1 alignment/read)
        if r.query_name in remaining_reads and not r.is_supplementary:
            unmapped.write(r)
    unmapped.close()
    temp_sam.close()

    # Merge all aligned reads and unmapped reads into a single sam
    ps.merge("-O", "SAM", "-@", str(n_cpu), sam_out, *iter_out)
    print(
        "{0} reads aligned / {1} total reads.".format(
            total_reads - len(remaining_reads), total_reads
        )
    )

    return 0


def truncate_reads(tmp_dir, infile, unaligned_set, n, min_len):
    """
    Writes the n first nucleotids of each sequence in infile to an auxialiary.
    file in the temporary folder.
    Parameters
    ----------
    tmp_dir : str
        Path to the temporary folder.
    infile : str
        Path to the fastq file to truncate.
    unaligned_set : set
        Contains the names of all reads that did not map unambiguously in
        previous rounds.
    n : int
        The number of basepairs to keep in each truncated sequence.
    str
        Path to the output fastq file containing truncated reads.
    """

    outfile = "{0}/truncated.fastq".format(tmp_dir)
    with ps.FastxFile(infile, "r") as inf, open(outfile, "w") as outf:
        for entry in inf:
            if entry.name in unaligned_set or n == min_len:
                entry.sequence = entry.sequence[:n]
                entry.quality = entry.quality[:n]
                outf.write(str(entry) + "\n")
    return outfile


def filter_samfile(temp_alignment, filtered_out):
    """
    Reads all the reads in the input SAM alignment file.
    Write reads to the output file if they are aligned with a good
    quality, otherwise add their name in a set to stage them for the next round
    of alignment.
    Parameters
    ----------
    temp_alignment : str
        Path to the input temporary alignment.
    outfile : str
        Path to the output filtered temporary alignment.
    Returns
    set:
        Contains the names reads that did not align.
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the final output file.
    # Keep those that do not map unambiguously for the next round.

    unaligned = set()
    temp_sam = ps.AlignmentFile(temp_alignment, "r")
    outf = ps.AlignmentFile(filtered_out, "w", template=temp_sam)
    for r in temp_sam:
        if r.flag in [0, 16] and r.mapping_quality >= 30:
            outf.write(r)
        else:
            unaligned.add(r.query_name)

    print("{0} reads left to map.".format(len(unaligned)))
    temp_sam.close()
    outf.close()

    return unaligned


##############################
#            MAIN
##############################


if __name__ == "__main__":
    # Get arguments
    args = parse_args()

    # Generates temporary folders
    temp_directory = generate_temp_dir(args.tempdir)

    # Aligns iteritatively the fastq file
    print(
        "\nIterative alignment of: {}\n".format(os.path.basename(args.in_fq))
    )
    iterative_align(
        args.in_fq,
        temp_directory,
        args.reference,
        args.nb_processors,
        args.out_sam,
        args.minimap2,
        args.min_len,
    )

    # Deletes the temporary folder
    st.rmtree(temp_directory)
