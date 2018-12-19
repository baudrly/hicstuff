#!/usr/bin/env python3
# coding: utf-8

"""Quickly visualize contact maps.

Usage:
    vizmap.py <contact_map> [--binning=1] [--normalize] [--output=<output.png>]
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

"""

import numpy as np
import functools
from matplotlib import pyplot as plt
from scipy import sparse
import docopt

SEABORN = False

try:
    import seaborn as sns

    SEABORN = True
except ImportError:
    pass

try:
    import hicstuff.hicstuff as hcs
except ImportError:
    print("Warning, hicstuff was not found - normalizations won't work")

DEFAULT_DPI = 500
DEFAULT_SATURATION_THRESHOLD = 99

load_raw_matrix = functools.partial(np.genfromtxt, skip_header=True, dtype=np.float64)


def raw_cols_to_sparse(M):
    n = int(np.amax(M[:, :-1]) + 1)

    row = M[:, 0]
    col = M[:, 1]
    data = M[:, 2]
    S = sparse.coo_matrix((data, (row, col)), shape=(n, n))
    return S


def sparse_to_dense(M):

    D = M.todense()
    E = D + np.transpose(D) - 2 * np.diag(np.diag(D))
    return E


def plot_matrix(array, filename, vmax=None, dpi=DEFAULT_DPI):
    """A function that performs all the tedious matplotlib
    magic to draw a 2D array with as few parameters and
    as little whitespace as possible.

    Adjusted from https://github.com/koszullab/metaTOR
    """

    if vmax is None:
        vmax = np.percentile(array, DEFAULT_SATURATION_THRESHOLD)
    # plt.gca().set_axis_off()
    # plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    # plt.margins(0, 0)
    # plt.gca().xaxis.set_major_locator(plt.NullLocator())
    # plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.figure()
    if SEABORN:
        sns.heatmap(array, vmax=vmax, cmap="Reds")
    else:
        plt.imshow(array, vmax=vmax, cmap="Reds", interpolation="none")
        plt.colorbar()
    plt.axis("off")
    if filename:
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.0, dpi=dpi)
        del filename
    else:
        plt.show()


def normalize(M, norm="SCN"):
    """Attempt to normalize if hicstuff is found, does nothing otherwise.
    """
    try:
        return hcs.normalize_sparse(M, norm=norm)
    except NameError:
        return M
