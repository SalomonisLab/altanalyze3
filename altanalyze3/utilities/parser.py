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
from altanalyze3.components.psi.psi_single import run_psi
from altanalyze3.components.psi.differential import run_bulk_differential_cli
from altanalyze3.components.gene_model.gene_model_index import build_index
from altanalyze3.components.fastCNV.main import bundled_gene_coordinates
from altanalyze3.components.fastCNV.main import run_from_altanalyze_args as run_fastcnv
from altanalyze3.components.long_read.cli import (
    run_sclr,
    run_sclr_junctions,
    run_sclr_isoforms,
    run_sclr_isoquant,
    run_sclr_gene_aggregate,
    run_sclr_diff,
    run_sclr_iso2func_network,
)
from altanalyze3.components.bam.variant_impact import run_variant_impact
# Optional subcommand handlers. iso2function and snaf may be absent from a given checkout; import
# them lazily so a missing OPTIONAL component never breaks the core CLI (sclr, sclr-junctions, etc.).
# This does NOT swallow the error: the stub re-raises a clear ImportError only if that subcommand is
# actually invoked, so the failure stays loud at the point of use.
def _missing_subcommand(_name, _err):
    def _stub(_args, _n=_name, _e=_err):
        raise ImportError(
            "The '%s' subcommand is unavailable in this altanalyze3 checkout: %s" % (_n, _e))
    return _stub
try:
    from altanalyze3.components.iso2function.cli import run_iso2function
except Exception as _e:
    run_iso2function = _missing_subcommand("sclr-iso2func", _e)
try:
    from altanalyze3.components.snaf.cli import (run_snaf, run_snaf_ts, run_snaf_b,
                                                 run_snaf_precompute_control,
                                                 run_snaf_build_surface_db)
except Exception as _e:
    run_snaf = run_snaf_ts = run_snaf_b = run_snaf_precompute_control = \
        run_snaf_build_surface_db = _missing_subcommand("snaf", _e)
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

        aggregate_parser.add_argument(
            "--novel-gene-mode", dest="novel_gene_mode", type=str,
            choices=["corrected", "legacy"], default="corrected",
            help=("Gene assignment for junctions whose two splice sites both miss the reference. "
                  "'corrected' (default) tests the chromosome, makes minus-strand genes reachable, "
                  "and annotates each splice site with its own coordinate. 'legacy' reproduces the "
                  "pre-2026-08 scan exactly, which ignores the chromosome and can never return a "
                  "minus-strand gene; use it only to reproduce an older result.")
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

        # PSI parameters
        psi_parser = subparsers.add_parser(
            "psi",
            help="Compute PSI per splicing event from an annotated junction matrix or h5ad"
        )
        psi_parser.set_defaults(func=run_psi)
        psi_parser.add_argument(
            "--junctions", dest="junctions", type=str, required=True,
            help="Annotated junction matrix (the *_annotated.tsv from aggregate) or a junction h5ad"
        )
        psi_parser.add_argument(
            "--output", type=str, required=True,
            help="Output PSI file"
        )
        psi_parser.add_argument(
            "--gene", dest="query_gene", type=str, default=None,
            help="Restrict the run to one gene. Default: all genes"
        )
        psi_parser.add_argument(
            "--min-reads", dest="min_reads", type=int, default=20,
            help=("Per-sample read floor an event must clear on BOTH its inclusion and its "
                  "exclusion side. Default: 20, the value that was previously hard-coded. Lower it "
                  "for shallow libraries, where 20 can reject every event.")
        )
        psi_parser.add_argument(
            "--min-denominator", dest="min_denominator", type=int, default=5,
            help="Below this total read count, that sample's PSI is null. Default: 5"
        )
        psi_parser.add_argument(
            "--min-dpsi-range", dest="min_dpsi_range", type=float, default=0.1,
            help="An event must vary by at least this much PSI across samples. Default: 0.1"
        )
        psi_parser.add_argument(
            "--min-junction-reads", dest="min_read", type=int, default=None,
            help="Drop a junction unless some sample exceeds this. h5ad input only. Default: off"
        )
        psi_parser.add_argument(
            "--tmp", type=Path, default=Path("/tmp/altanalyze_tmp"),
            help="Temporary directory. assert_common_args creates it for every subcommand"
        )
        psi_parser.add_argument(
            "--loglevel", type=str, default="INFO",
            help="Logging level: DEBUG, INFO, WARNING, ERROR"
        )

        # Differential splicing parameters
        diff_parser = subparsers.add_parser(
            "diff-splice",
            help="Differential splicing between two groups of bulk samples, using the long-read engine"
        )
        diff_parser.set_defaults(func=run_bulk_differential_cli)
        diff_parser.add_argument("--psi", type=str, required=True, help="PSI file from 'altanalyze3 psi'")
        diff_parser.add_argument(
            "--groups", type=str, required=True,
            help="Two-column file: sample<TAB>group. A header row is optional"
        )
        diff_parser.add_argument("--condition1", type=str, required=True, help="First group name")
        diff_parser.add_argument("--condition2", type=str, required=True, help="Second group name")
        diff_parser.add_argument(
            "--method", type=str, default="limma", choices=["limma", "mwu"],
            help=("'limma' (default) runs the empirical-Bayes moderated t-test; 'mwu' runs the "
                  "Mann-Whitney rank test. Both write the same columns.")
        )
        diff_parser.add_argument(
            "--output", type=str, required=True,
            help="Directory for <cond1>-<cond2>-bulk_stats.txt"
        )
        diff_parser.add_argument("--tmp", type=Path, default=Path("/tmp/altanalyze_tmp"), help="Temporary directory")
        diff_parser.add_argument("--loglevel", type=str, default="INFO", help="Logging level")

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
        sclr_parser.add_argument("--force", action="store_true", help="Reprocess and OVERWRITE existing per-sample outputs in place (re-extract / rebuild the junction h5ad with the cellHarmony barcode match / re-aggregate the gene h5ad) instead of reloading or skipping them. Nothing is deleted out-of-band.")
        sclr_parser.add_argument("--skip-bam-extract", dest="skip_bam_extract", action="store_true", help="RESUME phase 1 without re-extracting the BAMs: reuse the already-written <library>.gff.gz / <library>.h5ad and pick up at cellHarmony + junction export. Combine with --force to rebuild the junction h5ad / redo the barcode match (e.g. after fixing the annotation file). Skips the multi-hour extraction.")
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
        sclr_isoquant_parser.add_argument("--force", action="store_true", help="Reprocess + overwrite existing per-sample outputs instead of skipping.")
        self.add_common_arguments(sclr_isoquant_parser)

        # Gene aggregation (per sample, parallelizable): molecule h5ad -> <library>-gene.h5ad. Needed
        # before gene-level differentials on runs whose -gene.h5ad were not produced during extraction.
        sclr_gene_agg_parser = subparsers.add_parser(
            "sclr-gene-aggregate",
            parents=[parent_parser],
            help="Long-read (per sample): aggregate the molecule h5ad to gene level (<library>-gene.h5ad); parallelizable. Idempotent."
        )
        sclr_gene_agg_parser.set_defaults(func=run_sclr_gene_aggregate)
        sclr_gene_agg_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_gene_agg_parser.add_argument("--sample", default=None, type=str, help="Process ONE uid (parallel fan-out). Omit to loop over all uids.")
        sclr_gene_agg_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_gene_agg_parser.add_argument("--cell_annot", default=None, type=str, help="Barcode->cluster file/dir (cellHarmony format), same source the isoform/junction diffs use; joined in-memory to build the per-sample gene pseudobulk")
        sclr_gene_agg_parser.add_argument("--force", action="store_true", help="Re-aggregate + overwrite the existing <library>-gene.h5ad instead of skipping.")
        self.add_common_arguments(sclr_gene_agg_parser)

        # Phase 4 combine (1 job): combine per-sample isoform pseudobulks + isoform differentials
        sclr_diff_parser = subparsers.add_parser(
            "sclr-diff",
            parents=[parent_parser],
            help="Long-read integration: differential isoform/junction analysis between groups"
        )
        sclr_diff_parser.set_defaults(func=run_sclr_diff)
        sclr_diff_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr)")
        sclr_diff_parser.add_argument("--conditions", required=True, type=str, help="Group pairs from the metadata 'groups' column, e.g. 'young,AML-NPM1' (semicolon-separate multiple)")
        sclr_diff_parser.add_argument("--analyses", default="junction,isoform", type=str, help="Comma-separated subset of junction,isoform,isoform-ratio,gene. 'gene' runs gene-level differential expression (log2 CP10k pseudobulks, same metaDataAnalysis stats, stats-only -> diff-cluster-gene/diff-covariate-gene). Default: junction,isoform")
        sclr_diff_parser.add_argument("--method", default="mwu", choices=["mwu", "limma"], help="Differential test: mwu (Mann-Whitney, default) or limma (eBayes moderated t)")
        sclr_diff_parser.add_argument("--species", default="human", choices=["human", "mouse"], help="Default: human")
        sclr_diff_parser.add_argument("--gene_symbol", default=None, type=str, help="Ensembl-id -> symbol table. Default: bundled gzipped <species> annotations")
        sclr_diff_parser.add_argument("--cell_annot", default=None, type=str, help="Optional explicit barcode->cluster file/dir (else discovered from the sclr cellHarmony outputs)")
        sclr_diff_parser.add_argument("--network-viz", action="store_true", help="After differentials, generate iso2function condition interaction-network figures (cross-cell-state + group-contrast rewiring) for the top differential TF isoforms")
        sclr_diff_parser.add_argument("--network-viz-group", default=None, type=str, help="Reference group for the cross-cell-state network analysis (default: control of the first --conditions pair)")
        sclr_diff_parser.add_argument("--network-viz-edges", default="all,PDI", type=str, help="Interaction types for the network viz: any of all,PDI,PPI (default all,PDI)")
        self.add_common_arguments(sclr_diff_parser)

        # Standalone (remote-submittable) network visualization -- same step as sclr-diff --network-viz.
        sclr_netviz_parser = subparsers.add_parser(
            "sclr-iso2func-network",
            parents=[parent_parser],
            help="iso2function condition interaction-network figures for a finished long-read analysis dir (cross-cell-state + group-contrast rewiring)"
        )
        sclr_netviz_parser.set_defaults(func=run_sclr_iso2func_network)
        sclr_netviz_parser.add_argument("--metadata", required=True, type=str, help="Metadata file (same one used for sclr); its directory is the default dataset dir")
        sclr_netviz_parser.add_argument("--dataset", default=None, type=str, help="Directory with the combined pseudobulk h5ads (default: the metadata file's directory)")
        sclr_netviz_parser.add_argument("--out_dir", default=None, type=str, help="Output dir (default: <dataset>/iso2function_network)")
        sclr_netviz_parser.add_argument("--reference-group", default="young", type=str, help="Reference group for the cross-cell-state analysis (default: young)")
        sclr_netviz_parser.add_argument("--conditions", default=None, type=str, help="Optional group-contrast pairs, e.g. 'AML-SRSF2,aged;MDS-post-SRSF2,MDS-pre-SRSF2'")
        sclr_netviz_parser.add_argument("--edge-types", default="all,PDI", type=str, help="Interaction types: any of all,PDI,PPI (default all,PDI)")
        sclr_netviz_parser.add_argument("--top-tfs", default=10, type=int, help="Number of top rewiring TF isoforms to plot per group (default 10)")
        self.add_common_arguments(sclr_netviz_parser)

        variant_impact_parser = subparsers.add_parser(
            "variant-impact",
            parents=[parent_parser],
            help="Cell-state-resolved expression impact of called variants: per variant, run MarkerFinder (MUT vs WT) within each cell state with >= --min-cells of each, and optionally impute UNK-genotype cells"
        )
        variant_impact_parser.set_defaults(func=run_variant_impact)
        variant_impact_parser.add_argument("--metadata", required=True, type=str, help="sclr metadata file (uid, bam, library, ...); samples, BAMs and per-sample h5ad paths are taken from it")
        variant_impact_parser.add_argument("--variants", required=True, type=str, help="Variant positions, genotyped straight from the BAM's =/X CIGAR (no reference FASTA): 'chr:pos:label' semicolon-separated (e.g. 'chr21:34799432:RUNX1_W279*;chr17:76736877:SRSF2_P95R'), OR a mutations file (chrom start end label type)")
        variant_impact_parser.add_argument("--level", default="gene", choices=["gene", "isoform", "both"], help="Expression matrix to test: gene (<library>-gene.h5ad), isoform (<library>-isoform.h5ad), or both. Default: gene")
        variant_impact_parser.add_argument("--cell_annot", dest="cell_annot", default=None, type=str, help="cellHarmony barcode->cluster file for cell states. If omitted, obs['cluster'] in the h5ad is used")
        variant_impact_parser.add_argument("--mut-min", dest="mut_min", default=1, type=int, help="Min mutant (X CIGAR) reads in a cell to call MUT for the COMBINED confident set. Default: 1")
        variant_impact_parser.add_argument("--wt-min", dest="wt_min", default=2, type=int, help="Min wild-type (= CIGAR) reads in a cell to call WT for the COMBINED confident set. Default: 2")
        variant_impact_parser.add_argument("--min-cells", dest="min_cells", default=20, type=int, help="Minimum called MUT and WT cells in a cell state to run MarkerFinder. Default: 20")
        variant_impact_parser.add_argument("--top-n", dest="top_n", default=50, type=int, help="Top markers per group (MUT/WT). Default: 50")
        variant_impact_parser.add_argument("--min-mapq", dest="min_mapq", default=None, type=int, help="Minimum read mapping quality for genotyping. Default: 20 for discovery genotyping, 1 for --supervised/--validate-supervised (matching the callers, which accept low-MAPQ reads in repeats/homopolymers). Pass a value to override either.")
        variant_impact_parser.add_argument("--indel-window", dest="indel_window", default=5, type=int, help="bp window around an indel locus to absorb left/right-alignment when picking the major indel allele. Default: 5")
        variant_impact_parser.add_argument("--supervised", action="store_true", help="Non-default: FULLY SUPERVISED genotyping -- confirm the EXACT pre-defined REF>ALT allele from the caller (no allele discovery). Requires REF and ALT per --variants entry ('chr:pos:label:type:REF:ALT' or file columns 6,7). The default discovery genotyper is unchanged.")
        variant_impact_parser.add_argument("--validate-supervised", dest="validate_supervised", action="store_true", help="Non-default: validate the supervised genotyper on ALL four-caller variants of the requested sample(s) (discovery outputs beside the BAM) and write the specific variants it fails to detect; skips the impact/MarkerFinder pipeline.")
        variant_impact_parser.add_argument("--passenger-association", dest="passenger_association", action="store_true", help="Non-default: passenger->driver association. Per sample, build cell x (driver + four-caller-union nuclear+mt passenger) matrix with the supervised exact-allele genotyper, run the mt_variants.py statistics (sample + patient), and write the hits table. No caller re-runs.")
        variant_impact_parser.add_argument("--driver-variants", dest="driver_variants", default=None, type=str, help="Driver (pathogenic) panel TSV: name<TAB>chrom<TAB>pos<TAB>ref<TAB>alt. Required for --passenger-association.")
        variant_impact_parser.add_argument("--results-extraction", dest="results_extraction", default=None, type=str, help="results_variant_extraction directory (has {uid}_variant_readcounts.tsv); enables the variant_extract recovery assessment (step 7).")
        variant_impact_parser.add_argument("--patient-map", dest="patient_map", default=None, type=str, help="sample<TAB>patient TSV for the patient-level association.")
        variant_impact_parser.add_argument("--max-q", dest="max_q", default=0.05, type=float, help="Max BH q-value for a significant passenger association. Default: 0.05")
        variant_impact_parser.add_argument("--max-background", dest="max_background", default=0.15, type=float, help="Max background (fraction of wild-type cells carrying the marker) for the specificity gate. Default: 0.15")
        variant_impact_parser.add_argument("--build-only", dest="build_only", action="store_true", help="--passenger-association: only build the per-sample matrices (for a compute-node array); skip nomination.")
        variant_impact_parser.add_argument("--nominate-only", dest="nominate_only", action="store_true", help="--passenger-association: skip matrix building; run nomination + hits table on matrices already built.")
        variant_impact_parser.add_argument("--impute", action="store_true", help="Impute genotype for undetected cells from the per-cell-state MarkerFinder signature, then expand and re-render the heatmap")
        variant_impact_parser.add_argument("--gene_symbol", default=None, type=str, help="Optional Ensembl-id -> symbol table for gene-level heatmap labels")
        variant_impact_parser.add_argument("--samples", nargs="*", default=None, help="Optional subset of uids/libraries to process (default: all in metadata)")
        variant_impact_parser.add_argument("--output", default=None, type=str, help="Output directory (default: ./variant_impact)")
        self.add_common_arguments(variant_impact_parser, exclude=["output"])

        # iso2function: isoform-resolved PPI/PDI/function annotation + interaction networks (TFIso atlas)
        iso2func_parser = subparsers.add_parser(
            "sclr-iso2func",
            parents=[parent_parser],
            help="Annotate isoform PPI/PDI/function (TFIso atlas) and build structure-keyed interaction networks"
        )
        iso2func_parser.set_defaults(func=run_iso2function)
        iso2func_parser.add_argument("action", choices=["ingest", "crosswalk", "associate", "network", "enrich", "all"],
                                     help="pipeline stage to run")
        iso2func_parser.add_argument("--out-dir", default=None, type=str, help="Output dir for parsed/crosswalk/product tables (default: component data/ dir)")
        iso2func_parser.add_argument("--artifacts-dir", default=None, type=str, help="Output dir for network/figure artifacts (default: component artifacts/)")
        iso2func_parser.add_argument("--clone", action="append", default=None, help="Restrict the network to these atlas clone_ids (repeatable)")
        iso2func_parser.add_argument("--detected-only", action="store_true", help="Drop assayed-negative edges from the exported network")
        self.add_common_arguments(iso2func_parser)

        # SNAF: tumor-specificity-only (DB-free, offline, cross-platform)
        snaf_ts_parser = subparsers.add_parser(
            "snaf-ts",
            parents=[parent_parser],
            help="SNAF tumor-specificity scoring only (mean/mle/BayesTS). Needs a junction-count matrix + normal control h5ad; no Alt91_db, no MHC binder, no internet."
        )
        snaf_ts_parser.set_defaults(func=run_snaf_ts)
        snaf_ts_parser.add_argument("--juncounts", required=True, type=str, help="Junction-count matrix (TSV; rows=junction UIDs, cols=tumor samples)")
        snaf_ts_parser.add_argument("--control_h5ad", required=True, type=str, help="Normal-tissue control h5ad (e.g. GTEx_junction_counts.h5ad), obs=junctions/var=samples with var['tissue'],var['total_count'],obs['mean']")
        snaf_ts_parser.add_argument("--add_control", default=None, type=str, help="Extra controls 'name=path.txt,name2=path.h5ad' (TSV DataFrame or h5ad)")
        snaf_ts_parser.add_argument("--methods", nargs="*", default=["mle"], choices=["mle", "bayesian"], help="Tumor-specificity methods to add beyond the always-on 'mean'. Default: mle")
        snaf_ts_parser.add_argument("--filter_mode", default="maxmin", choices=["maxmin", "prevalance"], help="Neojunction sifting mode. Default: maxmin")
        snaf_ts_parser.add_argument("--min_samples", default=1, type=int, help="Keep only neojunctions expressed in >= this many tumor samples (recurrence filter). Default: 1 (no filtering)")
        snaf_ts_parser.add_argument("--t_min", default=20, type=int, help="maxmin: min (tumor - normal-mean) count. Default: 20")
        snaf_ts_parser.add_argument("--n_max", default=3, type=int, help="maxmin: max normal-mean count. Default: 3")
        snaf_ts_parser.add_argument("--normal_cutoff", default=5, type=int, help="prevalance: normal count cutoff. Default: 5")
        snaf_ts_parser.add_argument("--tumor_cutoff", default=20, type=int, help="prevalance: tumor count cutoff. Default: 20")
        snaf_ts_parser.add_argument("--normal_prevalance_cutoff", default=0.01, type=float, help="prevalance: max normal fraction. Default: 0.01")
        snaf_ts_parser.add_argument("--tumor_prevalance_cutoff", default=0.1, type=float, help="prevalance: min tumor fraction. Default: 0.1")
        snaf_ts_parser.add_argument("--bayes_mode", default="XY", choices=["XY", "Y"], help="BayesTS model (XY=tissue-dist+expression). Default: XY")
        snaf_ts_parser.add_argument("--bayes_epoch", default=2000, type=int, help="BayesTS SVI steps. Default: 2000")
        snaf_ts_parser.add_argument("--control_stats", default=None, type=str, help="Precomputed control-stats table from `snaf-precompute-control` (default: auto-detect <control_h5ad>.snaf_stats.tsv.gz). When present, the full control matrix is NOT loaded.")
        snaf_ts_parser.add_argument("--max_bayests_percentile", default=0.9, type=float, help="Drop sifted neojunctions whose precomputed BayesTS percentile exceeds this (0-1; lower=more tumor-specific). DEFAULT 0.9 (BayesTS filtering ON); auto-skips when the control has no BayesTS. Pass 1.0 (or a value >=1) to disable.")
        self.add_common_arguments(snaf_ts_parser)

        # SNAF: full MHC-bound T-antigen pipeline
        snaf_parser = subparsers.add_parser(
            "snaf",
            parents=[parent_parser],
            help="SNAF full T-antigen pipeline: sifting -> in-silico translation -> MHC binding (MHCflurry/netMHCpan) -> immunogenicity -> burden/frequency. Needs --db_dir (Alt91_db + controls) and --hla."
        )
        snaf_parser.set_defaults(func=run_snaf)
        snaf_parser.add_argument("--juncounts", required=True, type=str, help="Junction-count matrix (TSV; rows=junction UIDs, cols=tumor samples)")
        snaf_parser.add_argument("--db_dir", required=True, type=str, help="SNAF reference dir containing Alt91_db/ and controls/")
        snaf_parser.add_argument("--hla", required=True, type=str, help="Per-sample HLA table (sample<TAB>hla1,hla2,... or sample + HLA columns)")
        snaf_parser.add_argument("--gtex_db", default=None, type=str, help="Override control h5ad path (default: <db_dir>/controls/GTEx_junction_counts.h5ad)")
        snaf_parser.add_argument("--gtex_mode", default="count", choices=["count", "psi"], help="Control mode. Default: count")
        snaf_parser.add_argument("--binding_method", default="MHCflurry", choices=["MHCflurry", "netMHCpan"], help="MHC binder. Default: MHCflurry (cross-platform)")
        snaf_parser.add_argument("--software_path", default=None, type=str, help="netMHCpan executable path (only for --binding_method netMHCpan)")
        snaf_parser.add_argument("--genome_fasta", default=None, type=str, help="Local indexed genome FASTA for offline UTR/novel-exon sequence retrieval (else UCSC DAS)")
        snaf_parser.add_argument("--download_ref", action="store_true", help="If the reference bundle is missing under --db_dir, download it (~2.7 GB) from altanalyze.org and extract it there instead of erroring")
        snaf_parser.add_argument("--add_control", default=None, type=str, help="Extra controls 'name=path.txt,name2=path.h5ad'")
        snaf_parser.add_argument("--not_in_db", action="store_true", help="Drop junctions documented in any Ensembl transcript")
        snaf_parser.add_argument("--strict", action="store_true", help="Require canonical start-codon support during translation")
        snaf_parser.add_argument("--filter_mode", default="maxmin", choices=["maxmin", "prevalance"], help="Neojunction sifting mode. Default: maxmin")
        snaf_parser.add_argument("--min_samples", default=1, type=int, help="Keep only neojunctions expressed in >= this many tumor samples (recurrence filter). Default: 1 (no filtering)")
        snaf_parser.add_argument("--t_min", default=20, type=int, help="maxmin: min (tumor - normal-mean) count. Default: 20")
        snaf_parser.add_argument("--n_max", default=3, type=int, help="maxmin: max normal-mean count. Default: 3")
        snaf_parser.add_argument("--normal_cutoff", default=5, type=int, help="prevalance: normal count cutoff. Default: 5")
        snaf_parser.add_argument("--tumor_cutoff", default=20, type=int, help="prevalance: tumor count cutoff. Default: 20")
        snaf_parser.add_argument("--normal_prevalance_cutoff", default=0.01, type=float, help="prevalance: max normal fraction. Default: 0.01")
        snaf_parser.add_argument("--tumor_prevalance_cutoff", default=0.1, type=float, help="prevalance: min tumor fraction. Default: 0.1")
        snaf_parser.add_argument("--control_stats", default=None, type=str, help="Precomputed control-stats table from `snaf-precompute-control` (default: auto-detect <gtex_db>.snaf_stats.tsv.gz). When present, the full control matrix is NOT loaded and BayesTS is not re-run.")
        snaf_parser.add_argument("--max_bayests_percentile", default=0.9, type=float, help="Drop sifted neojunctions whose precomputed BayesTS percentile exceeds this (0-1; lower=more tumor-specific). DEFAULT 0.9 (BayesTS filtering ON); auto-skips when the control has no BayesTS. Pass 1.0 (or a value >=1) to disable.")
        self.add_common_arguments(snaf_parser)

        # SNAF-B: surface / B-antigen pipeline (pure-python: tmhmm.py + Biopython; no REST)
        snaf_b_parser = subparsers.add_parser(
            "snaf-b",
            parents=[parent_parser],
            help="SNAF surface/B-antigen pipeline: surface-gene neojunctions -> ORF recovery -> transmembrane topology (pure-python tmhmm.py) -> novelty vs UniProt -> stringency-gated candidates. Needs --db_dir; runs standalone (no SNAF-T, no HLA) when --freq_path is omitted."
        )
        snaf_b_parser.set_defaults(func=run_snaf_b)
        snaf_b_parser.add_argument("--juncounts", required=True, type=str, help="Junction-count matrix (TSV; rows=junction UIDs, cols=tumor samples)")
        snaf_b_parser.add_argument("--db_dir", required=True, type=str, help="SNAF reference dir containing Alt91_db/ and controls/")
        snaf_b_parser.add_argument("--freq_path", default=None, type=str, help="T-antigen frequency table from a prior `snaf` run (frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt). OPTIONAL: omit it and SNAF-B derives the equivalent table over its membrane neojunctions, so no SNAF-T run and no HLA types are needed.")
        snaf_b_parser.add_argument("--surface_db", default=None, type=str, help="Custom cell-surface gene database REPLACING the built-in Alt91_db surfaceome: a directory from `snaf-build-surface-db`, or a bare gene table with an Ensembl-gene-ID column. Genes with no reference protein are excluded and counted.")
        snaf_b_parser.add_argument("--mode", default="short_read", choices=["short_read", "long_read", "find_full_length"], help="Surface prediction mode. Default: short_read")
        snaf_b_parser.add_argument("--validation_gtf", default=None, type=str, help="Long-read/EST GTF or GFF, plain or gzipped (e.g. SQANTI, or a long-read combined.gff.gz) enabling stringency-4/5 support gates; omit for stringency-3-only (fully offline)")
        snaf_b_parser.add_argument("--no_tmhmm", action="store_true", help="Disable the transmembrane-topology gate (otherwise pure-python tmhmm.py is used)")
        snaf_b_parser.add_argument("--tmhmm_path", default=None, type=str, help="Path to a legacy TMHMM 2.0c binary (Linux); if unset, the pure-python tmhmm.py is used")
        snaf_b_parser.add_argument("--n_stride", default=2, type=int, help="ORF-check stride. Default: 2")
        snaf_b_parser.add_argument("--gtex_db", default=None, type=str, help="Override control h5ad path (default: <db_dir>/controls/GTEx_junction_counts.h5ad)")
        snaf_b_parser.add_argument("--control_stats", default=None, type=str, help="Precomputed control-stats table from `snaf-precompute-control` (default: auto-detect <gtex_db>.snaf_stats.tsv.gz). When present, the full control matrix is NOT loaded.")
        snaf_b_parser.add_argument("--gtex_mode", default="count", choices=["count", "psi"], help="Control mode. Default: count")
        snaf_b_parser.add_argument("--genome_fasta", default=None, type=str, help="Local indexed genome FASTA for offline UTR-event sequence retrieval (else UCSC DAS)")
        snaf_b_parser.add_argument("--add_control", default=None, type=str, help="Extra controls 'name=path.txt,name2=path.h5ad'")
        snaf_b_parser.add_argument("--not_in_db", action="store_true", help="Drop junctions documented in any Ensembl transcript")
        snaf_b_parser.add_argument("--download_ref", action="store_true", help="If the reference bundle is missing under --db_dir, download it (~2.7 GB) instead of erroring")
        snaf_b_parser.add_argument("--filter_mode", default="maxmin", choices=["maxmin", "prevalance"], help="Neojunction sifting mode. Default: maxmin")
        snaf_b_parser.add_argument("--min_samples", default=1, type=int, help="Keep only neojunctions expressed in >= this many tumor samples (recurrence filter). Default: 1 (no filtering)")
        snaf_b_parser.add_argument("--t_min", default=20, type=int, help="maxmin: min (tumor - normal-mean) count. Default: 20")
        snaf_b_parser.add_argument("--n_max", default=3, type=int, help="maxmin: max normal-mean count. Default: 3")
        snaf_b_parser.add_argument("--normal_cutoff", default=5, type=int, help="prevalance: normal count cutoff. Default: 5")
        snaf_b_parser.add_argument("--tumor_cutoff", default=20, type=int, help="prevalance: tumor count cutoff. Default: 20")
        snaf_b_parser.add_argument("--normal_prevalance_cutoff", default=0.01, type=float, help="prevalance: max normal fraction. Default: 0.01")
        snaf_b_parser.add_argument("--tumor_prevalance_cutoff", default=0.1, type=float, help="prevalance: min tumor fraction. Default: 0.1")
        snaf_b_parser.add_argument("--max_bayests_percentile", default=0.9, type=float, help="Drop sifted neojunctions whose precomputed BayesTS percentile exceeds this (0-1; lower=more tumor-specific). DEFAULT 0.9 (BayesTS filtering ON); auto-skips when the control has no BayesTS. Pass 1.0 (or a value >=1) to disable.")
        self.add_common_arguments(snaf_b_parser)

        # SNAF-B surface database builder: format a user surfaceome list (e.g. SURFY) into the
        # whitelist + reference-protein pair that SNAF-B needs.
        snaf_sdb_parser = subparsers.add_parser(
            "snaf-build-surface-db",
            parents=[parent_parser],
            help="Format a user cell-surface gene list (e.g. SURFY) into a SNAF-B surface database: the gene whitelist plus the reference protein sequences SNAF-B compares novel ORFs against. Feed the result to `snaf-b --surface_db`."
        )
        snaf_sdb_parser.set_defaults(func=run_snaf_build_surface_db)
        snaf_sdb_parser.add_argument("--gene_table", required=True, type=str, help="User surfaceome list: any TSV/CSV with a column of Ensembl gene IDs (a symbol column is used for labels). CRLF and a header row are handled.")
        snaf_sdb_parser.add_argument("--db_dir", required=True, type=str, help="SNAF reference dir containing Alt91_db/ (supplies the built-in reference sequences and the Ensembl-91 gene models)")
        snaf_sdb_parser.add_argument("--output", required=True, type=str, help="Directory to write the surface database into (surface_genes.txt, surface_reference.fasta, surface_db_params.json)")
        snaf_sdb_parser.add_argument("--uniprot_dir", default=None, type=str, help="AltAnalyze UniProt dir (e.g. .../EnsMart100/uniprot/Hs) with Hs_Ensembl-UniProt.txt + uniprot_sequence.txt. Without it, listed genes that have no built-in reference protein are excluded (and counted) instead of being built.")
        snaf_sdb_parser.add_argument("--mode", default="replace", choices=["replace", "union"], help="replace: the database is exactly the user's list (default). union: the user's list plus the built-in Alt91_db surfaceome.")
        snaf_sdb_parser.add_argument("--name", default=None, type=str, help="Label recorded in surface_db_params.json")
        snaf_sdb_parser.add_argument("--download_ref", action="store_true", help="If the reference bundle is missing under --db_dir, download it (~2.7 GB) instead of erroring")
        self.add_common_arguments(snaf_sdb_parser, exclude=["output"])

        # SNAF precompute-control: build the small per-junction control-stats table ONCE
        # (mean/std/mle/normal_prevalence + BayesTS sigma/percentile) so subsequent runs never
        # load the multi-GB control matrix or re-run BayesTS.
        snaf_pc_parser = subparsers.add_parser(
            "snaf-precompute-control",
            parents=[parent_parser],
            help="SNAF one-time control preprocessing: read a control h5ad and write a small per-junction stats table (mean/std/mle/normal_prevalence + BayesTS sigma/percentile). Runs consume this instead of the full matrix -- no multi-GB load, no per-run BayesTS."
        )
        snaf_pc_parser.set_defaults(func=run_snaf_precompute_control)
        snaf_pc_parser.add_argument("--control_h5ad", required=True, type=str, help="Control h5ad to summarize (e.g. GTEx_junction_counts.h5ad)")
        snaf_pc_parser.add_argument("--stats_output", default=None, type=str, help="Output stats-table path (default: <control_h5ad>.snaf_stats.tsv.gz, auto-detected at runtime)")
        snaf_pc_parser.add_argument("--normal_cutoff", default=5, type=int, help="Count cutoff for normal_prevalence. Default: 5 (MUST match the runtime --normal_cutoff)")
        snaf_pc_parser.add_argument("--bayes_mode", default="XY", choices=["XY", "Y"], help="BayesTS model. Default: XY")
        snaf_pc_parser.add_argument("--bayes_epoch", default=2000, type=int, help="BayesTS SVI steps. Default: 2000")
        snaf_pc_parser.add_argument("--bayes_batch", default=50000, type=int, help="BayesTS junctions per batched joint run (bounds memory). Default: 50000")
        snaf_pc_parser.add_argument("--bayes_cores", default=None, type=int, help="Parallel workers for the (independent, identical-result) BayesTS batches. Default: cpu_count-2. Set 1 to force serial. Bit-identical to serial; ~Nx faster.")
        snaf_pc_parser.add_argument("--no_bayes", action="store_true", help="Skip BayesTS (write only mean/std/mle/normal_prevalence)")
        snaf_pc_parser.add_argument("--bayes_juncounts", default=None, type=str, help="Restrict BayesTS to the junctions in this cohort count matrix (intersected with the control). ~50x faster than all-GTEx; percentile ranks tumor-specificity among the tested cohort. Basic stats (mean/mle/prevalence) are still written for ALL control junctions.")
        self.add_common_arguments(snaf_pc_parser)

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
