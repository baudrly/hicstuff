#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# cmdoret, 20181214
"""
Simple Hi-C pipeline for generating and manipulating contact matrices.

usage:
    hicstuff [-hv] <command> [<args>...]
options:
    -h, --help                  shows the help
    -v, --version               shows the version
The subcommands are:
    filter      filters Hi-C pairs to exclude spurious events
    digest      digest genome into a list of fragments
    iteralign   iteratively aligns reads to a reference genome
    view        visualize a Hi-C matrix
    pipeline    Hi-C pipeline to generate contact matrix from fastq files
"""

from docopt import docopt
from docopt import DocoptExit

import commands

if __name__ == "__main__":
    args = docopt(__doc__, version="1.0.0", options_first=True)
    # Retrieve the command to execute.
    command_name = args.pop("<command>").capitalize()

    # Retrieve the command arguments.
    command_args = args.pop("<args>")
    if command_args is None:
        command_args = {}

    # After 'popping' '<command>' and '<args>', what is left in the
    # args dictionary are the global arguments.

    # Retrieve the class from the 'commands' module.
    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command.")
        raise DocoptExit()

    # Create an instance of the command.
    command = command_lass(command_args, args)

    # Execute the command.
    command.execute()
