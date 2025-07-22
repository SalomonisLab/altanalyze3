import os
import sys
import pysam
import logging
import pathlib
from pathlib import Path
import argparse
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.helpers import get_version
from altanalyze3.components.intron_count.main import count_introns
from altanalyze3.components.junction_count.main import count_junctions
from altanalyze3.components.aggregate.main import aggregate
from altanalyze3.components.gene_model.gene_model_index import build_index
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
        self.resolve_path(["bam", "ref", "tmp", "output", "juncounts", "intcounts", "index"])
        self.assert_args()
        self.set_args_as_attributes()

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def add_common_arguments(self, parser, exclude=None):
        """
        Add common arguments to the parser.
        Args:
            parser: The parser to add arguments to
            exclude: List of argument names to exclude
        """
        exclude = exclude or []
        self.common_arguments = [
            ("--loglevel", "Logging level. Default: info", str, "info", ["fatal", "error", "warning", "info", "debug"]),
            ("--threads", "Number of threads to run in parallel where applicable. Default: 1", int, 1, None),
            ("--cpus", "Number of processes to run in parallel where applicable. Default: 1", int, 1, None),
            ("--tmp", "Temporary files location. Default: tmp", str, "tmp", None),
            ("--output", "Output prefix. Default: results", str, "results", None)
        ]
        for param in self.common_arguments:
            if param[0].lstrip('--') not in exclude:
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

        # Index parameters
        index_parser = subparsers.add_parser(
            "index",
            parents=[parent_parser],
            help="Create an index from a GTF/GFF file"
        )
        index_parser.set_defaults(func=build_index)
        index_parser.add_argument(
            "--gtf",
            help="Path to the GTF/GFF file to index",
            type=str,
            required=True
        )
        index_parser.add_argument(
            "--output",
            help="Output directory for the index files",
            type=str,
            required=True
        )
        # Add common arguments except output
        self.add_common_arguments(index_parser, exclude=['output'])

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
            default=0
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
            help="Aggregate junction and intron BED files into a single .h5ad matrix"
        )

        aggregate_parser.add_argument(
            "--juncounts", type=str, required=False,
            help="Junction count BED file or directory containing BED files"
        )
        aggregate_parser.add_argument(
            "--intcounts", type=str, required=False,
            help="Intron count BED file or directory containing BED files"
        )
        aggregate_parser.add_argument(
            "--output", type=str, required=True,
            help="Path prefix for output h5ad file (suffix .h5ad will be added)"
        )
        aggregate_parser.add_argument(
            "--chr", nargs="+", default=[],
            help="Optional list of chromosomes to retain (e.g. chr1 chr2 chrX)"
        )

        aggregate_parser.add_argument(
            "--ref", type=str, required=True,
            help="Reference exon BED file used for annotation"
        )

        aggregate_parser.set_defaults(func=aggregate)

        aggregate_parser.add_argument(
            "--loglevel", type=str, default="INFO",
            help="Logging level: DEBUG, INFO, WARNING, ERROR"
        )

        aggregate_parser.add_argument(
            "--tmp", type=Path, default=Path("/tmp/altanalyze_tmp"),
            help="Temporary directory for intermediate files"
        )
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

            from altanalyze3.components.aggregate import annotate  # assumes annotate.py is in that directory
            import anndata

            aggregated_path = Path(self.args.output).with_suffix(".h5ad")
            adata = anndata.read_h5ad(aggregated_path)
            # Use the full exon_file not the compressed bed version
            exon_file = self.args.ref.with_name(self.args.ref.name.replace(".bed.gz", "_all.tsv"))
            annotate.annotate_junctions(adata, exon_file)
            annotated_path = aggregated_path.with_name(aggregated_path.stem + "_annotated.h5ad")
            adata.write(annotated_path)

            # Export dense TSV
            annotate.export_dense_matrix(adata, annotated_path.with_suffix(".tsv"))

        elif self.args.func == build_index:
            self.assert_args_for_index()

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
        # Skip length check if input is Path (new usage: file or dir)
        if isinstance(self.args.juncounts, list) and isinstance(self.args.intcounts, list):
            if len(self.args.juncounts) != len(self.args.intcounts):
                raise ValueError("The number of junction and intron count files must match.")
        def assert_args_for_aggregate(self):
            # Only check list-lengths if both are lists (for backward compatibility)
            if isinstance(self.args.juncounts, list) and isinstance(self.args.intcounts, list):
                if len(self.args.juncounts) != len(self.args.intcounts):
                    raise ValueError("The number of junction and intron count files must match.")

        if self.args.juncounts is not None:
            if self.args.ref is None:
                logging.error("--ref parameter is required when using with --intcounts")
                sys.exit(1)

    def assert_args_for_index(self):
        """
        Assert arguments for the index command
        """
        if not os.path.exists(self.args.gtf):
            logging.error(f"GTF/GFF file not found: {self.args.gtf}")
            sys.exit(1)
        
        # Create output directory if it doesn't exist
        os.makedirs(self.args.output, exist_ok=True)

    def assert_common_args(self):
        self.args.loglevel = getattr(logging, self.args.loglevel.upper())
        setup_logger(logging.root, self.args.loglevel)
        self.args.tmp.mkdir(parents=True, exist_ok=True)                          # safety measure, shouldn't fail
        self.args.output.parent.mkdir(parents=True, exist_ok=True)                # safety measure, shouldn't fail
        
        # Only process chr attribute if it exists
        if hasattr(self.args, 'chr'):
            self.args.chr = list(map(ChrConverter, self.args.chr))
            
        if hasattr(self.args, "bam"):
            if not is_bam_indexed(self.args.bam):
                try:
                    pysam.index(str(self.args.bam))                               # attempts to create bai index (will raise if bam is not sorted)
                except pysam.SamtoolsError:
                    logging.error(f"""Failed to index {self.args.bam}. Exiting""")
                    sys.exit(1)