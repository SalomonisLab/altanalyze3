"""Enable ``python3 -m altanalyze3 <subcommand> ...``.

Mirrors the console entry point in ``bin/altanalyze3`` so the package can be run directly from a
source checkout (e.g. on a cluster with ``PYTHONPATH`` set to the source dir, no install needed).
"""
import sys

from altanalyze3.utilities.parser import ArgsParser
from altanalyze3.utilities.helpers import TimeIt


def main(args=None):
    with TimeIt():
        parsed = ArgsParser(sys.argv[1:] if args is None else args)
        parsed.func(parsed)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
