import logging
import pathlib
import argparse
from altanalyze3.utilities.helpers import get_version
from altanalyze3.components.junction_count.main import count_junctions
from altanalyze3.components.intron_count.main import count_introns
from altanalyze3.utilities.io import get_all_bam_chr
from altanalyze3.utilities.constants import IntRetCat


class ArgsParser():

    def __init__(self, args):
        args = args + [""] if len(args) == 0 else args
        self.args, _ = self.get_parser().parse_known_args(args)
        self.resolve_path(["bam", "ref", "output"])
        self.assert_args()
        self.set_args_as_attributes()
        logging.debug("Loaded parameters")
        logging.debug(self.args)

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def add_common_arguments(self, parser):
        self.common_arguments = [
            ("--loglevel", "Logging level. Default: info", str, "info", ["fatal", "error", "warning", "info", "debug"]),
            ("--threads", "Number of threads to run in parallel where applicable", int, 1, None),
            ("--cpus", "Number of processes to run in parallel where applicable", int, 1, None),
            ("--output", "Output prefix", str, "results", None)
        ]
        for param in self.common_arguments:
            parser.add_argument(
                param[0],
                help=param[1],
                type=param[2],
                default=param[3],
                choices=param[4]
            )

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
        self.add_common_arguments(junction_parser)

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
        self.add_common_arguments(intron_parser)
        return general_parser

    def resolve_path(self, selected=None):
        """
        Resolves path of the "selected" parameters.
        The other parameters remain unchanged
        """
        selected = [] if selected is None else selected
        normalized_args = {}
        for key, value in self.args.__dict__.items():
            if key in selected:
                if isinstance(value, list):
                    for v in value:
                        normalized_args.setdefault(key, []).append(
                            pathlib.Path(v).resolve()
                        )
                else:
                    normalized_args[key] = pathlib.Path(value).resolve()
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
        self.args.strandness = IntRetCat[self.args.strandness.upper()]

    def assert_args_for_count_introns(self):
        pass                                                 # nothing to check yet

    def assert_common_args(self):
        self.args.chr = get_all_bam_chr(self.args.bam, self.args.threads) \
            if len(self.args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in self.args.chr]
        self.args.loglevel = getattr(logging, self.args.loglevel.upper())