# -*- coding: utf-8 -*-
"""
Script to analyse the contents of a 3C library in terms of loops, uncuts,
weirds events, and filter those events.
@author: cmdoret (reimplementation of Axel KournaK's code)
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import contextlib
import argparse
import pysam as ps
import os


def parse_args():
    """ Gets the arguments from the command line."""
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "input_file",
        type=str,
        help="The file containing the coordinates of Hi-C interacting pairs "
        "and, the indices of their restriction fragments (.dat.indices format).",
    )
    parser.add_argument(
        "output_file",
        type=str,
        nargs="?",
        help="Path to the output file (filtered dat.indices). Defaults to stdout.",
    )
    group.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        required=False,
        help="Interactively shows plots and asks the user for thresholds. "
        "Overrides predefined thresholds. Disabled by default.",
    )
    group.add_argument(
        "-t",
        "--thresholds",
        nargs=2,
        type=int,
        required=False,
        help="The minimum number of restriction fragments between reads to "
        "consider loops and uncut events. Estimated automatically by default. "
        "Must be two integers separated by a space (-t <loop> <uncut>).",
    )
    parser.add_argument(
        "-p",
        "--plot_summary",
        required=False,
        action="store_true",
        default=False,
        help="Output a piechart summarizing library composition.",
    )
    return parser.parse_args()


def process_read_pair(line, sens2type):
    """
    Takes a read pair (line) from a dat.indices file as input, reorders the pair
    so that read 1 in intrachromosomal pairs always has the smallest genomic
    coordinate.
    Parameters
    ----------
    line : str
        Read pair from a dat.indices file.
    sens2type : dict
        Dictionary with entries mapping a tuple of int (sens1, sens2) to the
        name of the event type. e.g ('-', '-'): '-- (weirds)'.
    Returns
    -------
    dict
        Dictionary with reordered pair where each column from the line as an
        entry. The number of restriction sites separating reads in the pair and
        the type of event are also stored in the dictionary, under keys
        "nsites" and "type".

    Example
    -------
        >>> process_read_pair("a 1 - 3 b 2 - 4", {('-','-'):'--'})
        {'chr1': 'a', 'pos1': 1, 'sens1': '-', 'indice1': 3, 'chr2': 'b', 'pos2': 2, 'sens2': '-', 'indice2': 4, 'nsites': 1, 'type': 'inter'}
    """
    # Split line by whitespace
    p = line.split()
    if len(p) == 6:
        raise ValueError(
            "Your input file has only 6 columns instead of 8. "
            "Did you use a dat file instead of dat.indices ?"
        )
    # Saving each column as a dictionary entry with meaningful key
    cols = [
        "chr1",
        "pos1",
        "sens1",
        "indice1",
        "chr2",
        "pos2",
        "sens2",
        "indice2",
    ]
    p = {cols[i]: p[i] for i in range(len(cols))}
    # Transforming numeric columns to int
    for col in ["pos1", "indice1", "pos2", "indice2"]:
        p[col] = int(p[col])

    # invert records for intrachromosomal pairs where rec2 comes before
    # rec1 in genomic coordinates
    if p["chr1"] == p["chr2"] and p["pos2"] < p["pos1"]:
        p["sens1"], p["sens2"] = p["sens2"], p["sens1"]
        p["pos1"], p["pos2"] = p["pos2"], p["pos1"]
        p["indice1"], p["indice2"] = p["indice2"], p["indice1"]

    # Number of restriction sites separating reads in the pair
    p["nsites"] = p["indice2"] - p["indice1"]
    # Get event type
    if p["chr1"] == p["chr2"]:
        p["type"] = sens2type[(p["sens1"], p["sens2"])]
    else:
        p["type"] = "inter"

    return p


def get_thresholds(in_dat, interactive=False):
    """
    Analyses the events in the first million of Hi-C pairs in the library, plots
    the occurrences of each event type according to number of restriction
    fragments, and asks user interactively for the minimum threshold for uncuts
    and loops.

    Parameters
    ----------
    in_dat: file object
        File handle in read mode to the input .dat file containing Hi-C pairs.
    interactive: bool
        If True, plots are diplayed and thresholds are require
    Returns
    -------
    dictionary
        dictionary with keys "uncuts" and "loops" where the values are the
        corresponding thresholds entered by the user.
    """
    thr_uncut = None
    thr_loop = None
    max_sites = 500
    n_events = {
        event: np.zeros(max_sites)
        for event in [
            "++ (weirds)",
            "-- (weirds)",
            "+- (uncuts)",
            "-+ (loops)",
        ]
    }
    i = 0
    # Map of sense -> name of event for intrachromosomal pairs.
    sens2type = {
        ("+", "+"): "++ (weirds)",
        ("-", "-"): "-- (weirds)",
        ("+", "-"): "+- (uncuts)",
        ("-", "+"): "-+ (loops)",
    }
    # open the file for reading (just the first 1 000 000 lines)
    for line in in_dat:
        i += 1
        if i == 1000000:
            break
        # Process Hi-C pair into a dictionary
        p = process_read_pair(line, sens2type)
        # Type of event and number of restriction site between reads
        etype = p["type"]
        nsites = p["nsites"]
        # Count number of events for intrachrom pairs
        if etype != "inter" and nsites < max_sites:
            n_events[etype][nsites] += 1

    def plot_event(n_events, name):
        plt.xlim([-0.5, 15])
        plt.plot(
            range(n_events[name].shape[0]),
            n_events[name],
            "o-",
            label=name,
            linewidth=2.0,
        )

    if interactive:
        # PLot:
        plot_event(n_events, "++ (weirds)")
        plot_event(n_events, "-- (weirds)")
        plot_event(n_events, "+- (uncuts)")
        plot_event(n_events, "-+ (loops)")
        plt.grid()
        plt.xlabel("Number of restriction fragment(s)")
        plt.ylabel("Number of events")
        plt.yscale("log")
        plt.legend()
        plt.show(block=False)

        # Asks the user for appropriate thresholds
        print(
            "Please enter the number of restriction fragments separating "
            "reads in a Hi-C pair below which \033[91mloops\033[0m and "
            "\033[92muncuts\033[0m events will be excluded\n",
            file=sys.stderr,
        )
        thr_uncut = int(
            input("Enter threshold for the \033[92muncuts\033[0m events (+-):")
        )
        thr_loop = int(
            input("Enter threshold for the \033[91mloops\033[0m events (-+):")
        )
        plt.clf()
    else:
        # Estimate thresholds from data
        all_events = np.log(np.array(list(n_events.values())))
        # Compute median occurences at each restriction sites
        event_med = np.median(all_events, axis=0)
        # Compute MAD, to have a robust estimator of the expected deviation
        # from median at long distances
        mad = np.median(abs(all_events - event_med))
        exp_stdev = mad / 0.67449
        # Iterate over sites, from furthest to closest
        for site in range(max_sites)[::-1]:
            # For uncuts and loops, keep the last (closest) site where the
            # deviation from other events <= expected_stdev
            if (
                abs(np.log(n_events["+- (uncuts)"][site]) - event_med[site])
                <= exp_stdev
            ):
                thr_uncut = site
            if (
                abs(np.log(n_events["-+ (loops)"][site]) - event_med[site])
                <= exp_stdev
            ):
                thr_loop = site
        if thr_uncut is None or thr_loop is None:
            raise ValueError(
                "The threshold for loops or uncut could not be estimated. "
                "Please try running with -i to investigate the problem."
            )
        print(
            "Inferred thresholds: uncuts={0} loops={1}".format(
                thr_uncut, thr_loop
            ),
            file=sys.stderr,
        )
    return thr_uncut, thr_loop


def filter_events(
    in_dat, out_filtered, thr_uncut, thr_loop, plot_events=False
):
    """
    Filter out spurious intrachromosomal Hi-C pairs from input file. +- pairs
    that do not exceed the uncut threshold and -+ pairs that do not exceed the
    loop thresholds are excluded from the ouput file. All others are written.
    Parameters
    ----------
    in_dat : file object
        File handle in read mode to the input "dat" file containing Hi-C pairs.
    out_filtered : file object
        File handle in write mode the output filtered "dat" file.
    thr_uncut : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal +- pair.
    thr_loop : int
        Minimum number of restriction sites between reads to keep an
        intrachromosomal -+ pair.
    plot_events : bool
        If True, a plot summarising the proportion of each type of event will be
        shown after filtering.
    """
    n_uncuts = 0
    n_loops = 0
    n_weirds = 0
    n_int = 0
    lrange_intra = 0
    lrange_inter = 0

    # Map of sense -> name of event for intrachromosomal pairs.
    sens2type = {
        ("+", "+"): "++ (weirds)",
        ("-", "-"): "-- (weirds)",
        ("+", "-"): "+- (uncuts)",
        ("-", "+"): "-+ (loops)",
    }
    i = 0
    # open the files for reading and writing
    for line in in_dat:  # iterate over each line
        p = process_read_pair(line, sens2type)
        if p["chr1"] == p["chr2"]:
            if p["indice1"] == p["indice2"] and p["sens1"] == p["sens2"]:
                n_weirds += 1
            elif p["nsites"] < thr_loop and p["type"] == "-+ (loops)":
                n_loops += 1
            elif p["nsites"] < thr_uncut and p["type"] == "+- (uncuts)":
                n_uncuts += 1
            else:
                lrange_intra += 1
                out_filtered.write(
                    str(p["chr1"])
                    + "\t"
                    + str(p["pos1"])
                    + "\t"
                    + str(p["sens1"])
                    + "\t"
                    + str(p["indice1"])
                    + "\t"
                    + str(p["chr2"])
                    + "\t"
                    + str(p["pos2"])
                    + "\t"
                    + str(p["sens2"])
                    + "\t"
                    + str(p["indice2"])
                    + "\n"
                )
        if p["chr1"] != p["chr2"]:
            lrange_inter += 1
            out_filtered.write(
                str(p["chr1"])
                + "\t"
                + str(p["pos1"])
                + "\t"
                + str(p["sens1"])
                + "\t"
                + str(p["indice1"])
                + "\t"
                + str(p["chr2"])
                + "\t"
                + str(p["pos2"])
                + "\t"
                + str(p["sens2"])
                + "\t"
                + str(p["indice2"])
                + "\n"
            )

    if lrange_inter > 0:
        ratio_inter = round(
            100 * lrange_inter / float(lrange_intra + lrange_inter), 2
        )
    else:
        ratio_inter = 0

    # Dump quick summary of operation results into stderr
    kept = lrange_intra + lrange_inter
    discarded = n_loops + n_uncuts + n_weirds
    print(
        "Proportion of inter contacts: {}%".format(ratio_inter),
        file=sys.stderr,
    )
    print(
        "{0} pairs discarded: Loops: {1}, Uncuts: {2}, Weirds: {3}".format(
            discarded, n_loops, n_uncuts, n_weirds
        ),
        file=sys.stderr,
    )
    print(
        "{0} pairs kept ({1}%)".format(
            kept, round(100 * kept / (kept + discarded), 2)
        ),
        file=sys.stderr,
    )

    # Visualize summary interactively if requested by user
    if plot_events:

        # Plot: make a square figure and axes to plot a pieChart:
        plt.figure(1, figsize=(6, 6))
        ax = plt.axes([0.1, 0.1, 0.8, 0.8])
        # The slices will be ordered and plotted counter-clockwise.
        labels = "Uncuts", "Loops", "Weirds", "3D intra", "3D inter"
        fracs = [n_uncuts, n_loops, n_weirds, lrange_intra, lrange_inter]
        colors = ["salmon", "lightskyblue", "lightcoral", "palegreen", "plum"]
        plt.pie(
            fracs,
            labels=labels,
            colors=colors,
            autopct="%1.1f%%",
            shadow=True,
            startangle=90,
        )
        plt.title(
            "Distribution of library events",
            bbox={"facecolor": "1.0", "pad": 5},
        )
        plt.text(
            0.3,
            1.15,
            "Threshold Uncuts =" + str(thr_uncut),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            0.3,
            1.05,
            "Threshold Loops =" + str(thr_loop),
            fontdict=None,
            withdash=False,
        )

        plt.text(
            -1.5,
            -1.2,
            "Total number of reads =" + str(i),
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.5,
            -1.3,
            "Ratio inter/(intra+inter) =" + str(ratio_inter) + "%",
            fontdict=None,
            withdash=False,
        )
        plt.text(
            -1.5,
            -1.4,
            "selected reads = {0}%".format(
                float(lrange_inter + lrange_intra)
                / (n_loops + n_uncuts + n_weirds + lrange_inter + lrange_intra)
            ),
            fontdict=None,
            withdash=False,
        )
        plt.show()


if __name__ == "__main__":
    args = parse_args()
    # Open connection for writing
    output_handle = (
        sys.stdout if not args.output_file else open(args.output_file, "w")
    )
    if args.thresholds:
        # Thresholds supplied by user beforehand
        uncut_thr, loop_thr = args.thresholds
    else:
        # Threshold defined at runtime
        with open(args.input_file) as handle_in:
            uncut_thr, loop_thr = get_thresholds(
                handle_in, interactive=args.interactive
            )
    # Filter library and write to output file
    with open(args.input_file) as handle_in:
        filter_events(
            handle_in,
            output_handle,
            uncut_thr,
            loop_thr,
            plot_events=args.plot_summary,
        )

