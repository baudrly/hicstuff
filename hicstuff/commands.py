# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# cmdoret, 20181412
from docopt import docopt
import hicstuff
import sys


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


class Digest(AbstractCommand):
    """
    Digests a fasta file into fragments.

    usage:
        digest --fasta FILE --ezyme ENZYME [--size SIZE]

    options:
        --fasta      Fasta file to be digested
        --enzyme     A restriction enzyme or an integer representing chunk sizes (in bp)
        --size       Minimum size threshold to keep fragments
        --output
    """

    def execute(self):
        if not self.args["--size"]:
            self.args["--size"] = 0
        if not self.args["--output"]:
            self.args["--output"] = sys.stderr
        hicstuff.fraglist(
            self.args["--fasta"],
            self.args["--enzyme"],
            self.args["--size"],
            output_frags=self.args["--output"],
        )
