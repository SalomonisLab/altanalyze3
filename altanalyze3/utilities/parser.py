import logging
import argparse
from altanalyze3.utilities.helpers import (
    get_absolute_path,
    get_version
)
from altanalyze3.components.junction_count.main import count_junctions
from altanalyze3.components.intron_count.main import count_introns
from altanalyze3.utilities.io import get_all_bam_chr
from altanalyze3.utilities.constants import IntRetCat


class ArgsParser():

    def __init__(self, args):
        args = args + [""] if len(args) == 0 else args
        self.args, _ = self.get_parser().parse_known_args(args)
        self.normalize_args(["func", "span", "threads", "cpus", "chr", "strandness", "loglevel", "savereads"])
        self.assert_args()
        self.set_args_as_attributes()
        logging.debug("Loaded parameters")
        logging.debug(self.args)

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def get_parser(self):
        """
        Defines arguments for parser. Inlcudes subparsers for
        each component of the tool.
        """
        parent_parser = argparse.ArgumentParser(add_help=False)
        general_parser = argparse.ArgumentParser(description="AltAnalyze3")
        subparsers = general_parser.add_subparsers()
        subparsers.required = True
        # Global parameters for all components of the tool
        general_parser.add_argument(                       
            "--version",
            action="version",
            version=get_version(),
            help="Print current version and exit"
        )
        # Junction count parameters
        junction_parser = subparsers.add_parser(
            "juncount",
            parents=[parent_parser],
            help="Count reads for junctions"
        )
        junction_parser.set_defaults(func=count_junctions)
        junction_parser.add_argument(
            "--bam",
            help="Path to the coordinate-sorted indexed BAM file",
            type=str,
            required=True
        )
        junction_parser.add_argument(
            "--ref",
            help="Path to the coordinate-sorted indexed gene model reference BED file",
            type=str,
            required=True
        )
        junction_parser.add_argument(
            "--span",
            help="5' and 3' overlap that read should have over a splice-site to be counted",
            type=int,
            default=10
        )
        junction_parser.add_argument(
            "--strandness",
            help=" ".join(
                [
                    "Strand specificty of the RNA library."
                    "'unstranded' - reads from the left-most end of the fragment",
                    "(in transcript coordinates) map to the transcript strand, and",
                    "the right-most end maps to the opposite strand.",
                    "'forward' - same as 'unstranded' except we enforce the rule that",
                    "the left-most end of the fragment (in transcript coordinates) is",
                    "the first sequenced (or only sequenced for single-end reads).",
                    "Equivalently, it is assumed that only the strand generated",
                    "during second strand synthesis is sequenced. Used for Ligation and",
                    "Standard SOLiD.",
                    "'reverse' - same as 'unstranded' except we enforce the rule that",
                    "the right-most end of the fragment (in transcript coordinates) is",
                    "the first sequenced (or only sequenced for single-end reads).",
                    "Equivalently, it is assumed that only the strand generated during",
                    "first strand synthesis is sequenced. Used for dUTP, NSR, and NNSR.",
                    "Default: first 'auto' (try to detect strand from the XS tag",
                    "of the read), then downgrade to 'unstranded'"
                ]
            ),
            type=str,
            default="auto",
            choices=["auto", "forward", "reverse", "unstranded"]
        )
        junction_parser.add_argument(
            "--chr",
            help="Select chromosomes to process. Default: all available",
            type=str,
            nargs="*",
            default=[]
        )
        junction_parser.add_argument(
            "--savereads",
            help="Export processed reads into the BAM file. Default: False",
            action="store_true"
        )
        junction_parser.add_argument(
            "--loglevel",
            help="Logging level. Default: info",
            type=str,
            default="info",
            choices=["fatal", "error", "warning", "info", "debug"]
        )
        junction_parser.add_argument(
            "--threads",
            help="Number of threads to run in parallel where applicable",
            type=int,
            default=1
        )
        junction_parser.add_argument(
            "--cpus",
            help="Number of processes to run in parallel where applicable",
            type=int,
            default=1
        )
        junction_parser.add_argument(
            "--output",
            help="Output prefix",
            type=str,
            default="results"
        )

        # Intron count parameters
        intron_parser = subparsers.add_parser(
            "intcount",
            parents=[parent_parser],
            help="Count reads for introns"
        )
        intron_parser.set_defaults(func=count_introns)
        intron_parser.add_argument(
            "--bam",
            help="Path to the coordinate-sorted indexed BAM file",
            type=str,
            required=True
        )
        intron_parser.add_argument(
            "--chr",
            help="Select chromosomes to process. Default: all available",
            type=str,
            nargs="*",
            default=[]
        )
        intron_parser.add_argument(
            "--savereads",
            help="Export processed reads into the BAM file. Default: False",
            action="store_true"
        )
        intron_parser.add_argument(
            "--loglevel",
            help="Logging level. Default: info",
            type=str,
            default="info",
            choices=["fatal", "error", "warning", "info", "debug"]
        )
        intron_parser.add_argument(
            "--threads",
            help="Number of threads to run in parallel where applicable",
            type=int,
            default=1
        )
        intron_parser.add_argument(
            "--cpus",
            help="Number of processes to run in parallel where applicable",
            type=int,
            default=1
        )
        intron_parser.add_argument(
            "--output",
            help="Output prefix",
            type=str,
            default="results"
        )
        return general_parser

    def normalize_args(self, skip=None):
        """
        Converts all relative path arguments to absolute ones relatively
        to the current working directory. Skipped arguments and None will
        be returned unchanged.
        """
        skip = [] if skip is None else skip
        normalized_args = {}
        for key,value in self.args.__dict__.items():
            if key not in skip and value is not None:
                if isinstance(value, list):
                    for v in value:
                        normalized_args.setdefault(key, []).append(
                            get_absolute_path(v)
                        )
                else:
                    normalized_args[key] = get_absolute_path(value)
            else:
                normalized_args[key] = value
        self.args = argparse.Namespace(**normalized_args)

    def assert_args(self):
        """
        Should be used to assert and fix parameters.
        Also can be used to set default values for not
        set parameters in case the later ones depend on other
        parameters that should be first parsed by argparser
        """
        self.assert_common_args()
        if self.args.func == count_junctions:
            self.assert_args_for_count_junctions()
        elif self.args.func == count_introns:
            self.assert_args_for_count_introns()
        else:
            pass

    def assert_args_for_count_junctions(self):
        self.args.chr = get_all_bam_chr(self.args.bam, self.args.threads) if len(self.args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in self.args.chr]
        self.args.strandness = {
            "auto": IntRetCat.AUTO,
            "forward": IntRetCat.FORWARD,
            "reverse": IntRetCat.REVERSE,
            "unstranded": IntRetCat.UNSTRANDED
        }[self.args.strandness]

    def assert_args_for_count_introns(self):
        self.args.chr = get_all_bam_chr(self.args.bam, self.args.threads) if len(self.args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in self.args.chr]

    def assert_common_args(self):
        self.args.loglevel = {
            "fatal": logging.FATAL,
            "error": logging.ERROR,
            "warning": logging.WARNING,
            "info": logging.INFO,
            "debug": logging.DEBUG
        }[self.args.loglevel]