"""Command-line entry points for the long-read single-cell workflow.

These are an ALTERNATIVE way to run exactly the same analyses as the driver scripts
(``run_full_workflow.py`` / ``run_workflow_simple.py``). They call the identical
``isoform_automate`` / ``cellHarmony_automate`` / ``comparisons`` functions, but slice the work so
the embarrassingly-parallel per-sample step (``sclr``) can be fanned out as one cluster job per
sample, with the cross-sample integration steps (``sclr-junctions`` / ``sclr-isoforms`` /
``sclr-diff``) run once afterwards. Outputs are identical to the single-process driver.

Phase mapping (see PARALLEL_CLUSTER_DESIGN.md):
  sclr            -> per-sample BAM extract -> junction h5ad (+ optional cellHarmony cluster labels)
  sclr-junctions  -> return_cluster_order + pre_process_samples   (junction agg + PSI + diff splicing)
  sclr-isoforms   -> combine_processed_samples                    (two-tier collapse + isoform h5ads)
  sclr-diff       -> compute_differentials                        (differential isoform/junction)

Per-sample BAM directory and reverse-complement are read from the metadata row, so nothing about
either appears on the command line. ``--sample`` selects ONE uid; omitting it loops serially over
all uids in the metadata (the small-N convenience path -- parallelism comes from the submit script
launching one ``sclr --sample <uid>`` per job).
"""

import os
import glob
import gzip
import json
import shutil
import logging
from pathlib import Path

# Bundled, gzipped default annotation files (packaged in the repo). Resolved per --species.
_HERE = os.path.dirname(os.path.abspath(__file__))
_RESOURCE_DIR = os.path.join(_HERE, "resources")
_DEFAULT_ANNOT = {
    "human": {
        "exon_annot": os.path.join(_RESOURCE_DIR, "Hs_Ensembl_exon.txt.gz"),
        "gene_symbol": os.path.join(_RESOURCE_DIR, "Hs_Ensembl-annotations.txt.gz"),
    },
    "mouse": {
        "exon_annot": os.path.join(_RESOURCE_DIR, "Mm_Ensembl_exon.txt.gz"),
        "gene_symbol": os.path.join(_RESOURCE_DIR, "Mm_Ensembl-annotations.txt.gz"),
    },
}

# cellHarmony reference registry shipped with cellHarmony-web (reused, not duplicated).
_CH_REGISTRY = os.path.join(_HERE, "..", "cellHarmony", "flask", "reference_config.json")


# --------------------------------------------------------------------------- helpers


def _maybe_gunzip(path, work_dir):
    """The long-read readers (``parse_exon_file`` etc.) open these as plain text. If the resolved
    default is gzipped, materialize a plain-text copy under ``work_dir`` once and return that.
    A user-supplied uncompressed path is returned unchanged."""
    path = str(path)
    if not path.endswith(".gz"):
        return path
    os.makedirs(work_dir, exist_ok=True)
    out = os.path.join(work_dir, os.path.basename(path)[:-3])
    if not os.path.exists(out):
        with gzip.open(path, "rt") as src, open(out, "w") as dst:
            shutil.copyfileobj(src, dst)
    return out


def _resolve_default(value, species, key, work_dir):
    """Return an explicit ``--exon_annot``/``--gene_symbol`` path, or the bundled default for the
    species (gunzipped on demand). Errors clearly if the bundled default is missing."""
    if value:
        return str(value)
    default = _DEFAULT_ANNOT.get(species, {}).get(key)
    if not default or not os.path.exists(default):
        raise FileNotFoundError(
            f"No --{key} given and no bundled default for species '{species}' "
            f"(expected {default}). Pass --{key} explicitly or build one with `altanalyze3 index`."
        )
    return _maybe_gunzip(default, work_dir)


def _resolve_cellharmony_ref(value, species):
    """Resolve ``--cellHarmony_ref`` to a centroid/states TSV path.

    Accepts EITHER a reference ``id`` from the cellHarmony registry (reference_config.json) for the
    selected species (e.g. ``hs_bm_reference``, ``hs_lung_hlca_reference``, ``mm_bm_reference``), OR
    a path to a user-provided centroid ``.txt``. Reuses the web's resolution so every stored
    reference is available; raises with the valid ids if the value matches neither.
    """
    value = str(value)
    if os.path.exists(value):
        return value  # user-provided centroid file
    registry_path = Path(_CH_REGISTRY)
    if not registry_path.exists():
        raise FileNotFoundError(
            f"--cellHarmony_ref '{value}' is not a file and the reference registry is missing "
            f"({registry_path}). Pass a centroid .txt path instead."
        )
    with registry_path.open("r", encoding="utf-8") as handle:
        registry = json.load(handle)
    valid = []
    for entry in registry.get("species", []):
        if entry.get("id") != species:
            continue
        for ref in entry.get("references", []):
            valid.append(ref.get("id"))
            if ref.get("id") == value:
                states = ref.get("states_tsv")
                p = Path(states)
                if not p.is_absolute():
                    p = (registry_path.parent / p).resolve()
                return str(p)
    raise ValueError(
        f"--cellHarmony_ref '{value}' is neither an existing file nor a known reference id for "
        f"species '{species}'. Valid ids: {', '.join(valid) if valid else '(none)'}."
    )


def _cellharmony_outdir(metadata_file):
    return os.path.join(os.path.dirname(os.path.abspath(str(metadata_file))) or os.getcwd(),
                        "cellHarmony")


def _discover_barcode_cluster_dirs(metadata_file, cell_annot=None):
    """For the integration steps: gather the per-sample barcode->cluster TSVs produced by the
    phase-1 ``sclr`` jobs (under ``<metadata dir>/cellHarmony/<library>/<library>_barcode_clusters.txt``),
    or use an explicit ``--cell_annot`` file/dir if the user supplied their own annotations."""
    if cell_annot:
        cell_annot = str(cell_annot)
        if os.path.isdir(cell_annot):
            return sorted(glob.glob(os.path.join(cell_annot, "*barcode_clusters.txt")) or
                          glob.glob(os.path.join(cell_annot, "*.txt")))
        return [cell_annot]
    ch = _cellharmony_outdir(metadata_file)
    found = sorted(glob.glob(os.path.join(ch, "*", "*_barcode_clusters.txt")))
    if not found:
        raise FileNotFoundError(
            f"No barcode->cluster TSVs found under {ch}. Run `altanalyze3 sclr` for each sample "
            f"first (or pass --cell_annot)."
        )
    return found


def _assert_phase1_complete(sample_dict):
    """Validation gate: every sample must have its per-sample junction h5ad before integration."""
    import altanalyze3.components.long_read.isoform_automate as isoa  # noqa: F401 (path setup)
    missing = []
    for uid, samples in sample_dict.items():
        for s in samples:
            gff = str(s["gff"])
            h5 = Path(f"{gff.split('.g')[0]}-junction.h5ad")
            if not h5.exists():
                missing.append(str(h5))
    if missing:
        raise FileNotFoundError(
            "Phase-1 outputs missing for {} sample(s); run `altanalyze3 sclr` for them first:\n  {}"
            .format(len(missing), "\n  ".join(missing))
        )


def _subset_metadata_to_sample(sample_dict, sample):
    """Restrict the import_metadata sample_dict to a single uid (the --sample value)."""
    if sample not in sample_dict:
        raise ValueError(
            f"--sample '{sample}' not found in metadata. Available uids: "
            f"{', '.join(sample_dict.keys())}"
        )
    return {sample: sample_dict[sample]}


def _write_single_sample_metadata(metadata_file, sample, work_dir):
    """Write a one-uid metadata file (same header + only that uid's row(s)) and return its path.

    ``import_metadata(extract_from_bams=True)`` performs BAM extraction inline while it parses EVERY
    row, so filtering the returned sample_dict afterwards would still extract all samples. To
    process a single ``--sample`` we must hand import_metadata a metadata file containing only that
    uid. This also matches the per-sample-metadata cluster pattern (one job per file).
    """
    import pandas as pd
    df = pd.read_csv(str(metadata_file), sep="\t", dtype=str)
    sub = df[df["uid"].astype(str) == str(sample)]
    if sub.empty:
        raise ValueError(
            f"--sample '{sample}' not found in {metadata_file}. Available uids: "
            f"{', '.join(sorted(df['uid'].astype(str).unique()))}"
        )
    os.makedirs(work_dir, exist_ok=True)
    out = os.path.join(work_dir, f"_metadata_{sample}.txt")
    sub.to_csv(out, sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- commands


def run_sclr(args):
    """Phase 1 (per sample): BAM -> molecule h5ad + per-sample junction h5ad (+ cluster labels).

    Idempotent and parallel-safe: processes one --sample, or loops serially over all uids when
    --sample is omitted. cellHarmony alignment (--cellHarmony_ref) or pre-existing annotations
    (--cell_annot) are optional; outputs land beside each sample's BAM and under
    <metadata dir>/cellHarmony/, ready for the integration steps.
    """
    import altanalyze3.components.long_read.isoform_automate as isoa
    import altanalyze3.components.cellHarmony.cellHarmony_automate as cha

    species = args.species
    work_dir = _cellharmony_outdir(args.metadata)
    exon_annot = _resolve_default(getattr(args, "exon_annot", None), species, "exon_annot", work_dir)
    gene_symbol = _resolve_default(getattr(args, "gene_symbol", None), species, "gene_symbol", work_dir)

    if getattr(args, "cellHarmony_ref", None) and getattr(args, "cell_annot", None):
        raise ValueError("--cellHarmony_ref and --cell_annot are mutually exclusive.")

    # Step 1: BAM -> molecule-level h5ad + gff (per-sample junction h5ad written by extraction path).
    # import_metadata extracts inline while parsing EVERY row, so to process a single --sample we
    # hand it a one-uid metadata file (otherwise all samples would be extracted).
    metadata_for_extract = str(args.metadata)
    if getattr(args, "sample", None):
        metadata_for_extract = _write_single_sample_metadata(args.metadata, args.sample, work_dir)
    sample_dict = isoa.import_metadata(
        metadata_for_extract, include_hashed_samples=True,
        extract_from_bams=True, reference_model=exon_annot,
    )

    # Step 2 (optional): gene aggregation + cellHarmony per sample -> barcode->cluster TSV.
    if getattr(args, "cellHarmony_ref", None):
        ref = _resolve_cellharmony_ref(args.cellHarmony_ref, species)
        cha.annotate_samples(
            sample_dict,
            cellharmony_ref=ref,
            gene_translation_file=gene_symbol,
            output_dir=work_dir,
        )
        logging.info("sclr: wrote cellHarmony barcode->cluster labels under %s", work_dir)
    elif getattr(args, "cell_annot", None):
        logging.info("sclr: using user-provided annotations %s (no alignment run)", args.cell_annot)
    else:
        logging.info("sclr: no cluster labeling requested (supply at integration time).")

    # Step 3 (P1): per-sample junction quantification + junction pseudobulk + Tier-1 collapse inputs.
    # Everything here is per-sample (no other sample needed), so it belongs in this parallel job --
    # the combine (P2) then only concatenates. Needs this sample's barcode->cluster labels.
    barcode_cluster_dirs = _discover_barcode_cluster_dirs(args.metadata, getattr(args, "cell_annot", None))
    import altanalyze3.components.long_read.isoform_matrix as iso
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)
    for uid, libs in sample_dict.items():
        isoa.export_sample_junctions(libs, exon_annot, barcode_sample_dict, uid=uid)
    logging.info("sclr: per-sample junction h5ad + pseudobulk written.")

    logging.info("sclr: phase-1 complete for %s",
                 args.sample if getattr(args, "sample", None) else "all samples")


def run_sclr_junctions(args):
    """Phase 2 (one job, cross-sample): combine the per-sample junction pseudobulks built in P1 ->
    filter -> PSI. No per-sample work (that was P1).

    NOTE: splice DIFFERENTIALS are NOT run here. Their annotation step
    (junction_isoform.annotate) loads gff-output/protein_summary.txt + transcript_associations.txt,
    which are produced by P3 (sclr-isoforms, the translation/collapse job). P2 only depends on P1, so
    protein_summary.txt does not exist yet when P2 runs. Splice + isoform differentials therefore both
    run in the post-P3 combine job (sclr-diff)."""
    import altanalyze3.components.long_read.isoform_automate as isoa
    import altanalyze3.components.long_read.isoform_matrix as iso

    species = args.species
    work_dir = _cellharmony_outdir(args.metadata)
    exon_annot = _resolve_default(getattr(args, "exon_annot", None), species, "exon_annot", work_dir)

    sample_dict = isoa.import_metadata(str(args.metadata), include_hashed_samples=True)
    _assert_phase1_complete(sample_dict)

    barcode_cluster_dirs = _discover_barcode_cluster_dirs(args.metadata, getattr(args, "cell_annot", None))
    iso.return_cluster_order(barcode_cluster_dirs)  # validates/orders clusters across samples
    isoa.combine_junctions(str(args.metadata), barcode_cluster_dirs, exon_annot)
    logging.info("sclr-junctions: junction combine + PSI complete (differentials run in sclr-diff, post-P3).")


def run_sclr_isoforms(args):
    """Phase 3 (one job, cross-sample): isoform COLLAPSE catalog + translation/FASTA only. The heavy
    per-sample re-key is P4 (sclr-isoquant). Persists the catalog + structure->exemplar (+EM soft)
    maps so P4 can re-key from disk."""
    import altanalyze3.components.long_read.isoform_automate as isoa

    species = args.species
    work_dir = _cellharmony_outdir(args.metadata)
    exon_annot = _resolve_default(getattr(args, "exon_annot", None), species, "exon_annot", work_dir)

    sample_dict = isoa.import_metadata(str(args.metadata), include_hashed_samples=True)
    _assert_phase1_complete(sample_dict)

    collapse_method = getattr(args, "collapse_method", "wta")
    isoa.build_isoform_catalog(
        str(args.metadata), exon_annot, str(args.ref_gff), str(args.genome_fasta),
        collapse_method=collapse_method, force_recollapse=True,
    )
    logging.info("sclr-isoforms: collapse catalog (%s) + translation/FASTA complete (re-key is P4).",
                 collapse_method)


def run_sclr_isoquant(args):
    """Phase 4 (per sample, parallel): re-key each sample's molecules onto the P3 catalog + build the
    per-sample isoform pseudobulk. Loads the catalog maps from disk (no dependency on P3 memory).
    With --sample, processes one uid; without, loops all. The final isoform-differential combine is a
    separate one-job call (sclr-diff)."""
    import altanalyze3.components.long_read.isoform_automate as isoa
    import altanalyze3.components.long_read.isoform_matrix as iso

    species = args.species
    work_dir = _cellharmony_outdir(args.metadata)
    gff_output_dir = os.path.join(os.path.dirname(os.path.abspath(str(args.metadata))) or os.getcwd(),
                                  "gff-output")
    if not os.path.exists(os.path.join(gff_output_dir, "FINAL_structure_to_exemplar.tsv")):
        raise FileNotFoundError(
            f"Collapse catalog not found in {gff_output_dir}. Run `altanalyze3 sclr-isoforms` (P3) first."
        )

    sample_dict = isoa.import_metadata(str(args.metadata), include_hashed_samples=True)
    if getattr(args, "sample", None):
        sample_dict = _subset_metadata_to_sample(sample_dict, args.sample)

    barcode_cluster_dirs = _discover_barcode_cluster_dirs(args.metadata, getattr(args, "cell_annot", None))
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)
    collapse_method = getattr(args, "collapse_method", "wta")
    for uid, libs in sample_dict.items():
        isoa.export_sample_isoform(libs, gff_output_dir, barcode_sample_dict, uid=uid,
                                   collapse_method=collapse_method, compute_tpm=True)
    logging.info("sclr-isoquant: per-sample isoform re-key + pseudobulk (%s) complete for %s.",
                 collapse_method, args.sample if getattr(args, "sample", None) else "all samples")


def _parse_conditions(value):
    """Parse '--conditions groupA,groupB' (semicolon-separate multiple) -> [(a,b), ...]. None/'' -> []."""
    conditions = []
    if not value:
        return conditions
    for pair in str(value).split(";"):
        a, _, b = pair.partition(",")
        if a and b:
            conditions.append((a.strip(), b.strip()))
    return conditions


def run_sclr_diff(args):
    """Phase 4 combine (post-P3): run differential analysis. Pass --analyses junction, isoform, or
    both (default). For isoform analyses it first concatenates the per-sample isoform pseudobulks
    (built by sclr-isoquant) into the combined isoform matrix; the junction (splice) analysis reuses
    the PSI matrix produced in P2.

    All differentials run HERE (P4, post-P3) rather than in P2 because their annotation
    (junction_isoform.annotate) requires gff-output/protein_summary.txt + transcript_associations.txt,
    which are P3 (sclr-isoforms) products. P2 is gated only on P1, so those files don't exist yet when
    P2 runs."""
    import altanalyze3.components.long_read.isoform_automate as isoa
    import altanalyze3.components.long_read.isoform_matrix as iso
    import altanalyze3.components.long_read.comparisons as comp

    species = args.species
    work_dir = _cellharmony_outdir(args.metadata)
    gene_symbol = _resolve_default(getattr(args, "gene_symbol", None), species, "gene_symbol", work_dir)

    sample_dict = isoa.import_metadata(str(args.metadata), include_hashed_samples=True)
    barcode_cluster_dirs = _discover_barcode_cluster_dirs(args.metadata, getattr(args, "cell_annot", None))
    cluster_order = iso.return_cluster_order(barcode_cluster_dirs)

    conditions = _parse_conditions(getattr(args, "conditions", None))
    if not conditions:
        raise ValueError("--conditions must be like 'young,AML-NPM1' (semicolon-separate multiple).")
    analyses = [a.strip() for a in str(getattr(args, "analyses", "junction,isoform")).split(",") if a.strip()]
    if not analyses:
        raise ValueError("--analyses must include at least one of: junction, isoform, isoform-ratio.")
    wants_isoform = any(a.startswith('isoform') for a in analyses)

    # Guard: differential annotation needs P3 outputs. Fail loudly with a clear message rather than
    # the cryptic FileNotFoundError on 'gff-output/protein_summary.txt' deep in annotate().
    protein_summary = os.path.join('gff-output', 'protein_summary.txt')
    if not os.path.exists(protein_summary):
        raise FileNotFoundError(
            f"{protein_summary} not found -- run `altanalyze3 sclr-isoforms` (P3, translation) before "
            f"sclr-diff. Differential annotation depends on the P3 protein summary."
        )

    # Combine the per-sample isoform pseudobulks (built per-sample in P4 sclr-isoquant) -- only when an
    # isoform analysis is requested. Junction-only diffs reuse the P2 PSI matrix and need no combine.
    if wants_isoform:
        isoa.combine_isoforms(str(args.metadata), barcode_cluster_dirs)

    comp.compute_differentials(
        sample_dict, conditions, cluster_order, gene_symbol,
        analyses=analyses, method=getattr(args, "method", "mwu"),
    )
    logging.info("sclr-diff: %s differentials complete for %s", analyses, conditions)
