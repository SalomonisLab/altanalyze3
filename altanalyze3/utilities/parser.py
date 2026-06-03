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
from altanalyze3.components.fastCNV.main import bundled_gene_coordinates
from altanalyze3.components.fastCNV.main import run_from_altanalyze_args as run_fastcnv
from altanalyze3.components.long_read.cli import (
    run_sclr,
    run_sclr_junctions,
    run_sclr_isoforms,
    run_sclr_isoquant,
    run_sclr_diff,
)
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
        self.resolve_path(["bam", "ref", "tmp", "output", "juncounts", "intcounts", "index", "h5ad", "gene_coordinates",
                           "metadata", "ref_gff", "exon_annot", "gene_symbol", "cell_annot", "genome_fasta"])
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

        # fastCNV parameters
        fastcnv_parser = subparsers.add_parser(
            "fastcnv",
            parents=[parent_parser],
            help="Run first-pass state-aware CNV calls from a cellHarmony h5ad"
        )
        fastcnv_parser.set_defaults(func=run_fastcnv)
        fastcnv_parser.add_argument("--h5ad", required=True, type=str, help="Input cellHarmony h5ad file")
        fastcnv_parser.add_argument("--gene-coordinates", default=None, type=str, help="Gene coordinate TSV/CSV")
        fastcnv_parser.add_argument("--species", choices=["human", "mouse"], default=None, help="Use a bundled coordinate resource")
        fastcnv_parser.add_argument("--state-key", required=True, help="obs column with cellHarmony cell-state labels")
        fastcnv_parser.add_argument("--sample-key", default=None, help="Optional obs column with sample labels")
        fastcnv_parser.add_argument("--layer", default="auto", help="'auto', 'X', or a layer name. Auto prefers layers['counts']")
        fastcnv_parser.add_argument("--input-normalized", action="store_true", help="Use selected matrix values as already log-normalized")
        fastcnv_parser.add_argument("--window-genes", type=int, default=41, help="Genes per adaptive genome window")
        fastcnv_parser.add_argument("--stride-genes", type=int, default=7, help="Stride between adaptive genome windows")
        fastcnv_parser.add_argument("--min-chr-genes", type=int, default=25, help="Minimum genes required to use a chromosome")
        fastcnv_parser.add_argument("--min-state-cells", type=int, default=30, help="Minimum cells required to call a cell state")
        fastcnv_parser.add_argument("--anchor-fraction", type=float, default=0.25, help="Lowest-burden state fraction used as WT anchors")
        fastcnv_parser.add_argument("--min-anchor-cells", type=int, default=20, help="Minimum WT-anchor cells per state")
        fastcnv_parser.add_argument("--high-threshold", type=float, default=2.6, help="Window score required to seed an interval")
        fastcnv_parser.add_argument("--low-threshold", type=float, default=1.6, help="Window score used to extend an interval")
        fastcnv_parser.add_argument("--min-run-windows", type=int, default=3, help="Minimum contiguous windows per CNV interval")
        fastcnv_parser.add_argument("--min-interval-genes", type=int, default=60, help="Minimum genes per CNV interval")
        fastcnv_parser.add_argument("--min-mean-score", type=float, default=1.8, help="Minimum mean signed score for interval calls")
        fastcnv_parser.add_argument("--burden-quantile", type=float, default=0.95, help="Upper quantile used for per-cell CNV burden")
        fastcnv_parser.add_argument("--cnv-burden-threshold", type=float, default=1.8, help="Minimum cell burden for CNV-positive status")
        fastcnv_parser.add_argument("--max-clones-per-state", type=int, default=10, help="Maximum NMF clones per cell state")
        fastcnv_parser.add_argument("--max-global-clones", type=int, default=10, help="Maximum clones after cross-state merging")
        fastcnv_parser.add_argument("--min-clone-cells", type=int, default=5, help="Minimum cells required for an NMF clone")
        fastcnv_parser.add_argument("--clone-similarity-threshold", type=float, default=0.88, help="Cosine similarity threshold for cross-state clone merging")
        fastcnv_parser.add_argument("--clone-consensus-fraction", type=float, default=0.45, help="Minimum clone-cell fraction supporting a consensus interval")
        fastcnv_parser.add_argument("--nmf-max-iter", type=int, default=300, help="Maximum sparse NMF iterations")
        fastcnv_parser.add_argument("--skip-clones", action="store_true", help="Only write first-pass cell-level CNV calls")
        fastcnv_parser.add_argument("--skip-pdf", action="store_true", help="Skip clone-level genome PDF export")
        fastcnv_parser.add_argument("--write-h5ad", action="store_true", help="Write an h5ad copy with fastCNV obs columns")
        fastcnv_parser.add_argument("--max-cells-per-state", type=int, default=None, help="Optional cap for fast exploratory runs")
        fastcnv_parser.add_argument("--random-state", type=int, default=0)
        self.add_common_arguments(fastcnv_parser)

        # ---------------------------------------------------------------------
        # Long-read single-cell workflow (sclr): an ALTERNATIVE, parallelizable
        # way to run the same analyses as the driver scripts. One per-sample
        # command (sclr) plus three integration commands. See
        # components/long_read/PARALLEL_CLUSTER_DESIGN.md.
        # ---------------------------------------------------------------------

        # Phase 1 (per sample): BAM -> junction h5ad (+ optional cellHarmony labels)
        sclr_parser = subparsers.add_parser(
            "sclr",
            parents=[parent_parser],
            help="Long-read single-cell per-sample processing (BAM -> junction h5ad); parallelizable"
        )
        sclr_parser.set_defaults(func=run_sclr)
        sclr_parser.add_argument("--metadata", required=True, type=str, help="Tab-delimited metadata file (uid, bam, library, reverse, groups)")
        sclr_parser.add_argument("--sample", default=None, type=str, help="Process ONE uid. Omit to loop serially over all uids in the metadata.")
        sclr_parser.add_argument("--ref_gff", required=True, type=str, help="GENCODE/Ensembl reference GFF/GTF (processed via gff_process)")
        sclr_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Selects bundled annotation defaults and the cellHarmony reference species. Default: human")
        sclr_parser.add_argument("--exon_annot", default=None, type=str, help="Exon annotation file. Default: bundled gzipped <species> Ensembl exon file")
        sclr_parser.add_argument("--gene_symbol", default=None, type=str, help="Ensembl-id -> symbol table. Default: bundled gzipped <species> annotations")
        sclr_parser.add_argument("--cellHarmony_ref", default=None, type=str, help="Align to this reference: a cellHarmony registry id (e.g. hs_bm_reference) or a centroid .txt path. Mutually exclusive with --cell_annot.")
        sclr_parser.add_argument("--cell_annot", default=None, type=str, help="Use existing barcode->cluster annotations (cellHarmony format) instead of aligning. Mutually exclusive with --cellHarmony_ref.")
        self.add_common_arguments(sclr_parser)

        # Phase 2 (integration, 1 job): combine per-sample junction pseudobulks + PSI + splice diff
        sclr_junc_parser = subparsers.add_parser(
            "sclr-junctions",
            parents=[parent_parser],
            help="Long-read P2 (1 job): combine junction pseudobulks + PSI + splicing differentials"
        )
        sclr_junc_parser.set_defaults(func=run_sclr_junctions)
        sclr_junc_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_junc_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_junc_parser.add_argument("--exon_annot", default=None, type=str, help="Exon annotation file. Default: bundled gzipped <species> Ensembl exon file")
        sclr_junc_parser.add_argument("--cell_annot", default=None, type=str, help="Optional explicit barcode->cluster file/dir (else discovered from the sclr cellHarmony outputs)")
        sclr_junc_parser.add_argument("--conditions", default=None, type=str, help="Optional group pairs for SPLICE differentials, e.g. 'young,aged' (semicolon-separate multiple). Omit to skip splice diff.")
        sclr_junc_parser.add_argument("--method", default="mwu", choices=["mwu", "limma"], help="Splice differential test: mwu (default) or limma")
        sclr_junc_parser.add_argument("--gene_symbol", default=None, type=str, help="Ensembl-id -> symbol table. Default: bundled gzipped <species> annotations")
        self.add_common_arguments(sclr_junc_parser)

        # Phase 3 (integration): isoform collapse + per-sample isoform h5ads + protein
        sclr_iso_parser = subparsers.add_parser(
            "sclr-isoforms",
            parents=[parent_parser],
            help="Long-read integration: two-tier isoform collapse + isoform h5ads + protein"
        )
        sclr_iso_parser.set_defaults(func=run_sclr_isoforms)
        sclr_iso_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_iso_parser.add_argument("--ref_gff", required=True, type=str, help="GENCODE/Ensembl reference GFF/GTF")
        sclr_iso_parser.add_argument("--genome_fasta", required=True, type=str, help="Genome FASTA for protein sequence prediction")
        sclr_iso_parser.add_argument("--collapse_method", default="wta", choices=["wta", "em"], help="Isoform read collapse: wta (winner-takes-all, default) or em (fractional EM allocation)")
        sclr_iso_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_iso_parser.add_argument("--exon_annot", default=None, type=str, help="Exon annotation file. Default: bundled gzipped <species> Ensembl exon file")
        sclr_iso_parser.add_argument("--cell_annot", default=None, type=str, help="Optional explicit barcode->cluster file/dir (else discovered from the sclr cellHarmony outputs)")
        self.add_common_arguments(sclr_iso_parser)

        # Phase 4 (per-sample, parallel): isoform re-key + per-sample isoform pseudobulk vs P3 catalog
        sclr_isoquant_parser = subparsers.add_parser(
            "sclr-isoquant",
            parents=[parent_parser],
            help="Long-read P4 (per sample): isoform re-key + pseudobulk against the collapse catalog; parallelizable"
        )
        sclr_isoquant_parser.set_defaults(func=run_sclr_isoquant)
        sclr_isoquant_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_isoquant_parser.add_argument("--sample", default=None, type=str, help="Process ONE uid. Omit to loop over all uids.")
        sclr_isoquant_parser.add_argument("--collapse_method", default="wta", choices=["wta", "em"], help="Must match the method used for sclr-isoforms (P3). wta (default) or em.")
        sclr_isoquant_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_isoquant_parser.add_argument("--cell_annot", default=None, type=str, help="Optional explicit barcode->cluster file/dir (else discovered from the sclr cellHarmony outputs)")
        self.add_common_arguments(sclr_isoquant_parser)

        # Phase 4 combine (1 job): combine per-sample isoform pseudobulks + isoform differentials
        sclr_diff_parser = subparsers.add_parser(
            "sclr-diff",
            parents=[parent_parser],
            help="Long-read integration: differential isoform/junction analysis between groups"
        )
        sclr_diff_parser.set_defaults(func=run_sclr_diff)
        sclr_diff_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_diff_parser.add_argument("--conditions", required=True, type=str, help="Group pairs from the metadata 'groups' column, e.g. 'young,AML-NPM1' (semicolon-separate multiple)")
        sclr_diff_parser.add_argument("--analyses", default="junction,isoform", type=str, help="Comma-separated subset of junction,isoform,isoform-ratio. Default: junction,isoform")
        sclr_diff_parser.add_argument("--method", default="mwu", choices=["mwu", "limma"], help="Differential test: mwu (Mann-Whitney, default) or limma (eBayes moderated t)")
        sclr_diff_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_diff_parser.add_argument("--gene_symbol", default=None, type=str, help="Ensembl-id -> symbol table. Default: bundled gzipped <species> annotations")
        sclr_diff_parser.add_argument("--cell_annot", default=None, type=str, help="Optional explicit barcode->cluster file/dir (else discovered from the sclr cellHarmony outputs)")
        self.add_common_arguments(sclr_diff_parser)

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
        elif self.args.func == run_fastcnv:
            self.assert_args_for_fastcnv()

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

    def assert_args_for_fastcnv(self):
        if not self.args.h5ad.exists():
            logging.error(f"h5ad file not found: {self.args.h5ad}")
            sys.exit(1)
        if self.args.gene_coordinates is None:
            if self.args.species is None:
                logging.error("Either --gene-coordinates or --species is required for fastcnv")
                sys.exit(1)
            self.args.gene_coordinates = bundled_gene_coordinates(self.args.species)
        if not self.args.gene_coordinates.exists():
            logging.error(f"Gene coordinate file not found: {self.args.gene_coordinates}")
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
