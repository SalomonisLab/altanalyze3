import os
import sys
import pysam
import logging
import pathlib
import argparse
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.helpers import get_version
from altanalyze3.components.intron_count.main import count_introns
from altanalyze3.components.junction_count.main import count_junctions
from altanalyze3.components.aggregate.main import aggregate
from altanalyze3.utilities.io import (
    get_indexed_references,
    is_bam_indexed
)
from altanalyze3.utilities.constants import (
    IntRetCat,
    MAIN_CRH,
    ChrConverter
)


class ArgsParser():

    def __init__(self, args):
        args = args + [""] if len(args) == 0 else args
        self.args, _ = self.get_parser().parse_known_args(args)
        self.resolve_path(["bam", "ref", "tmp", "output", "juncounts", "intcounts"])
        self.assert_args()
        self.set_args_as_attributes()

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def add_common_arguments(self, parser):
        self.common_arguments = [
            ("--loglevel", "Logging level. Default: info", str, "info", ["fatal", "error", "warning", "info", "debug"]),
            ("--threads", "Number of threads to run in parallel where applicable. Default: 1", int, 1, None),
            ("--cpus", "Number of processes to run in parallel where applicable. Default: 1", int, 1, None),
            ("--tmp", "Temporary files location. Default: tmp", str, "tmp", None),
            ("--output", "Output prefix. Default: results", str, "results", None)
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

        # Intron count parameters
        intron_parser = subparsers.add_parser(
            "intcount",
            parents=[parent_parser],
            help="Count introns reads"
        )
        intron_parser.set_defaults(func=count_introns)
        intron_parser.add_argument(
            "--bam",
            help="Path to the coordinate-sorted indexed BAM file",
            type=str,
            required=True
        )
        intron_parser.add_argument(
            "--ref",
            help="Path to the gene model reference file. Coordinates are treated as 1-based.",
            type=str,
            required=True
        )
        intron_parser.add_argument(
            "--span",
            help="5' and 3' overlap that read should have over a splice-site to be counted",
            type=int,
            default=10
        )
        intron_parser.add_argument(
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
        intron_parser.add_argument(
            "--chr",
            help="Select chromosomes to process. Default: only main chromosomes",
            type=str,
            nargs="*",
            default=MAIN_CRH
        )
        intron_parser.add_argument(
            "--savereads",
            help="Export processed reads into the BAM file. Default: False",
            action="store_true"
        )
        self.add_common_arguments(intron_parser)

        # Junction count parameters
        junction_parser = subparsers.add_parser(
            "juncount",
            parents=[parent_parser],
            help="Count junctions reads"
        )
        junction_parser.set_defaults(func=count_junctions)
        junction_parser.add_argument(
            "--bam",
            help="Path to the coordinate-sorted indexed BAM file",
            type=str,
            required=True
        )
        junction_parser.add_argument(
            "--chr",
            help="Select chromosomes to process. Default: only main chromosomes",
            type=str,
            nargs="*",
            default=MAIN_CRH
        )
        junction_parser.add_argument(
            "--savereads",
            help="Export processed reads into the BAM file. Default: False",
            action="store_true"
        )
        self.add_common_arguments(junction_parser)

        # Aggregate parameters
        aggregate_parser = subparsers.add_parser(
            "aggregate",
            parents=[parent_parser],
            help="Aggregate read counts produced by intcount and juncount"
        )
        aggregate_parser.set_defaults(func=aggregate)
        aggregate_parser.add_argument(
            "--juncounts",
            help=" ".join(
                [
                    "Path the junction counts files. Coordinates are treated as 0-based.",
                    "If provided with --intcounts, the number and the order should",
                    "correspond to --intcounts."
                ]
            ),
            type=str,
            nargs="*"
        )
        aggregate_parser.add_argument(
            "--intcounts",
            help=" ".join(
                [
                    "Path the intron counts files. Coordinates are treated as 0-based.",
                    "If provided with --juncounts, the number and the order should",
                    "correspond to --juncounts."
                ]
            ),
            type=str,
            nargs="*"
        )
        aggregate_parser.add_argument(
            "--aliases",
            help=" ".join(
                [
                    "Column names to be used for the loaded counts. The number of provided",
                    "aliases should be equal to the number of --intcounts and/or --juncounts",
                    "files. Default: rootname of --intcounts and/or --intcounts files."
                ]
            ),
            type=str,
            nargs="*"
        )
        aggregate_parser.add_argument(
            "--ref",
            help=" ".join(
                [
                    "Path to the gene model reference file. Coordinates are treated as 1-based.",
                    "Required if --juncounts parameter was provided."
                ]
            ),
            type=str,
            required=False
        )
        aggregate_parser.add_argument(
            "--chr",
            help="Select chromosomes to process. Default: only main chromosomes",
            type=str,
            nargs="*",
            default=MAIN_CRH
        )
        aggregate_parser.add_argument(
            "--bed",
            help="Export annotated 0-based coordinates as BED file. Default: False",
            action="store_true"
        )
        self.add_common_arguments(aggregate_parser)

        return general_parser

    def resolve_path(self, selected=None):
        """
        Resolves path of the "selected" parameters.
        The other parameters remain unchanged
        """
        selected = [] if selected is None else selected
        normalized_args = {}
        for key, value in self.args.__dict__.items():
            if key in selected and value is not None:
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
        elif self.args.func == aggregate:
            self.assert_args_for_aggregate()

    def assert_args_for_count_junctions(self):
        pass

    def assert_args_for_count_introns(self):
        self.args.strandness = IntRetCat[self.args.strandness.upper()]
        self.args.ref = get_indexed_references(
            location=self.args.ref,
            tmp_location=self.args.tmp,
            selected_chr=self.args.chr,
            only_introns=True
        )

    def assert_args_for_aggregate(self):
        if self.args.juncounts is None and self.args.intcounts is None:
            logging.error("At least one of the --juncounts or --intcounts inputs should be provided")
            sys.exit(1)
        if self.args.juncounts is not None and self.args.intcounts is not None and len(self.args.juncounts) != len(self.args.intcounts):
            logging.error("Number of the provided files for --juncounts and --intcounts should be equal")
            sys.exit(1)
        if self.args.aliases is not None:
            if self.args.juncounts is not None and len(self.args.aliases) != len(self.args.juncounts):
                logging.error("Number of the provided --aliases and --juncounts files should be equal")
                sys.exit(1)
            if self.args.intcounts is not None and len(self.args.aliases) != len(self.args.intcounts):
                logging.error("Number of the provided --aliases and --intcounts files should be equal")
                sys.exit(1)
        if self.args.aliases is None:
            zipped = zip(self.args.juncounts) if self.args.juncounts is not None else zip(self.args.intcounts)
            zipped = zip(*zip(*zipped), self.args.intcounts if self.args.intcounts is not None else self.args.juncounts)
            self.args.aliases = [os.path.commonprefix([a.stem, b.stem]) for a, b in zipped]
        if self.args.juncounts is not None:
            if self.args.ref is None:
                logging.error("--ref parameter is required when using with --intcounts")
                sys.exit(1)
            self.args.ref = get_indexed_references(
                location=self.args.ref,
                tmp_location=self.args.tmp,
                selected_chr=self.args.chr,
                only_introns=False
            )

    def assert_common_args(self):
        self.args.loglevel = getattr(logging, self.args.loglevel.upper())
        setup_logger(logging.root, self.args.loglevel)
        self.args.tmp.mkdir(parents=True, exist_ok=True)                          # safety measure, shouldn't fail
        self.args.output.parent.mkdir(parents=True, exist_ok=True)                # safety measure, shouldn't fail
        self.args.chr = list(map(ChrConverter, self.args.chr))
        if hasattr(self.args, "bam"):
            if not is_bam_indexed(self.args.bam):
                try:
                    pysam.index(str(self.args.bam))                               # attempts to create bai index (will raise if bam is not sorted)
                except pysam.SamtoolsError:
                    logging.error(f"""Failed to index {self.args.bam}. Exiting""")
                    sys.exit(1)