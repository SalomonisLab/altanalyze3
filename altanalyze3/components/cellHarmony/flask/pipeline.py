from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import re
import zipfile
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.cellHarmony import cellHarmony_differential, cellHarmony_lite
from altanalyze3.components.visualization import NetPerspective, approximate_umap as approx_mod

from .job_manager import JobStore


def _flatten_matrix(x) -> np.ndarray:
    if sp.issparse(x):
        return np.asarray(x.todense()).ravel()
    arr = np.asarray(x)
    return arr.ravel()


def _split_uploads(files: List[Dict[str, str]], uploads_dir: Path) -> Tuple[List[Tuple[str, str]], Optional[str]]:
    h5_files: List[Tuple[str, str]] = []
    h5ad_entries: List[str] = []
    for record in files:
        path = uploads_dir / record["filename"]
        sample_name = str(record.get("sample_name", "")).strip() or path.stem
        suffix = path.suffix.lower()
        if suffix == ".h5ad":
            h5ad_entries.append(str(path))
        else:
            h5_files.append((str(path), sample_name))
    h5ad_file = None
    if h5ad_entries:
        if len(h5ad_entries) > 1:
            raise ValueError("Multiple .h5ad uploads detected; please upload a single combined .h5ad or 10x files.")
        if h5_files:
            raise ValueError("Mixing .h5ad with other file types is not supported.")
        h5ad_file = h5ad_entries[0]
    if not h5_files and not h5ad_file:
        raise ValueError("No compatible input files detected.")
    return h5_files, h5ad_file


def _lookup_reference(species: str, reference_id: str, registry_path: Path) -> Dict:
    with registry_path.open("r", encoding="utf-8") as handle:
        registry = json.load(handle)
    for entry in registry.get("species", []):
        if entry.get("id") != species:
            continue
        for ref in entry.get("references", []):
            if ref.get("id") == reference_id:
                resolved = dict(ref)
                base_dir = registry_path.parent
                for key in (
                    "states_tsv",
                    "reference_clusters_tsv",
                    "reference_coords_tsv",
                    "reference_h5ad",
                ):
                    value = resolved.get(key)
                    if not value:
                        continue
                    path = Path(value)
                    if not path.is_absolute():
                        resolved[key] = str((base_dir / path).resolve())
                return resolved
    raise ValueError(f"Reference '{reference_id}' for species '{species}' not found in registry.")


def _ensure_reference_fields(reference_entry: Dict) -> None:
    required = ["states_tsv", "reference_clusters_tsv", "reference_coords_tsv"]
    missing = [key for key in required if key not in reference_entry]
    if missing:
        raise ValueError(f"Reference metadata missing required keys: {', '.join(missing)}")


def _first_reference_gene(states_tsv: str) -> Optional[str]:
    try:
        reference_df = pd.read_csv(states_tsv, sep="\t", index_col=0, nrows=1)
    except Exception:
        return None
    if reference_df.empty:
        return None
    return str(reference_df.index[0])


def _job_sample_names(meta: Dict) -> List[str]:
    sample_names: List[str] = []
    for record in meta.get("files", []):
        sample_name = str(record.get("sample_name", "")).strip()
        if sample_name:
            sample_names.append(sample_name)
    return sample_names


def _sample_name_aliases(meta: Dict) -> Dict[str, List[str]]:
    aliases: Dict[str, List[str]] = {}
    for record in meta.get("files", []):
        sample_name = str(record.get("sample_name", "")).strip()
        filename = str(record.get("filename", "")).strip()
        if not sample_name or not filename:
            continue
        path = Path(filename)
        candidates = [sample_name, path.stem]
        if path.suffix.lower() == ".h5":
            candidates.append(path.name.replace(".h5", ""))
        elif path.suffix.lower() == ".h5ad":
            candidates.append(path.name.replace(".h5ad", ""))
        match = re.match(r"(.+)_\d+$", path.stem)
        if match:
            candidates.append(match.group(1))
        deduped = []
        seen = set()
        for candidate in candidates:
            candidate = str(candidate).strip()
            if candidate and candidate not in seen:
                deduped.append(candidate)
                seen.add(candidate)
        aliases[sample_name] = deduped
    return aliases


def _resolve_samples_for_adata(
    adata: ad.AnnData,
    meta: Dict,
    selected_samples: List[str],
) -> Tuple[str, Dict[str, str]]:
    candidate_columns = [column for column in ("Library", "group", "sample") if column in adata.obs.columns]
    if not candidate_columns:
        raise ValueError("Aligned AnnData is missing the sample metadata required for group comparison.")

    aliases = _sample_name_aliases(meta)
    selected_order = [str(sample).strip() for sample in selected_samples if str(sample).strip()]
    best_column = None
    best_resolved: Dict[str, str] = {}
    best_score = -1

    for column in candidate_columns:
        available = set(adata.obs[column].astype(str))
        resolved: Dict[str, str] = {}
        score = 0
        for sample in selected_order:
            candidates = aliases.get(sample, [])
            if sample not in candidates:
                candidates = [sample] + candidates
            match = next((candidate for candidate in candidates if candidate in available), None)
            if match is not None:
                resolved[sample] = match
                score += 1
        if score > best_score:
            best_column = column
            best_resolved = resolved
            best_score = score
        if score == len(selected_order):
            return column, resolved

    missing = [sample for sample in selected_order if sample not in best_resolved]
    details = []
    if best_column is not None:
        values = sorted(set(adata.obs[best_column].astype(str)))
        preview = ", ".join(values[:10])
        suffix = ", …" if len(values) > 10 else ""
        details.append(f"best matching column '{best_column}' contains: {preview}{suffix}")
    if candidate_columns:
        details.append(f"available sample columns: {', '.join(candidate_columns)}")
    extra = f" ({'; '.join(details)})" if details else ""
    raise ValueError(f"Selected samples were not found in the aligned AnnData: {', '.join(missing)}{extra}")


def _candidate_population_columns(h5ad_path: Path, preferred: Optional[List[str]] = None) -> List[Dict[str, object]]:
    if not h5ad_path.exists():
        return []

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs = adata.obs.copy()
        n_obs = max(int(adata.n_obs), 1)
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()

    max_categories = max(8, min(200, n_obs // 2 if n_obs > 1 else 2))
    excluded_columns = {
        "sample",
        "group",
        "Library",
        "pct_counts_mt",
        "UMAP-X",
        "UMAP-Y",
        "UMAP_leiden_X",
        "UMAP_leiden_Y",
    }

    options: List[Dict[str, object]] = []
    for column in obs.columns:
        if column in excluded_columns:
            continue
        series = obs[column]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            continue
        non_null = series[series.notna()].astype(str).str.strip()
        non_null = non_null[non_null != ""]
        unique_values = pd.Index(non_null).unique().tolist()
        n_unique = len(unique_values)
        if n_unique < 2 or n_unique > max_categories:
            continue
        options.append(
            {
                "value": str(column),
                "label": str(column),
                "n_categories": n_unique,
            }
        )

    preferred_order = {value: idx for idx, value in enumerate(preferred or [])}
    options.sort(key=lambda item: (0 if item["value"] in preferred_order else 1, preferred_order.get(item["value"], 999), item["label"]))
    return options


def _default_differential_state(enabled: bool, default_population_col: Optional[str]) -> Dict[str, object]:
    if enabled:
        message = "Select a cell-state field and numerator/denominator sample groups to run cellHarmony-differential."
    else:
        message = "Differential analysis is available when the job includes two or more samples."
    return {
        "status": "idle",
        "progress": 0,
        "message": message,
        "config": {},
        "artifacts": {},
        "networks": [],
        "go_terms_included": False,
        "default_population_col": default_population_col,
    }


def _update_differential_state(store: JobStore, job_id: str, **changes) -> Dict[str, object]:
    meta = store.get_job(job_id)
    current = dict(meta.get("differential") or {})
    current.update(changes)
    store.update_job(job_id, differential=current)
    return current


def _group_display_label(samples: List[str]) -> str:
    ordered = [str(sample).strip() for sample in samples if str(sample).strip()]
    if len(ordered) <= 3:
        return " + ".join(ordered)
    head = " + ".join(ordered[:3])
    return f"{head} + {len(ordered) - 3} more"


def _comparison_tag(group1_samples: List[str], group2_samples: List[str]) -> str:
    lhs = "__".join(NetPerspective.safe_component(sample, fallback="sample") for sample in group1_samples)
    rhs = "__".join(NetPerspective.safe_component(sample, fallback="sample") for sample in group2_samples)
    return f"{lhs}_vs_{rhs}"


def _write_zip_archive(source_dir: Path, archive_path: Path) -> Path:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for file_path in sorted(source_dir.rglob("*")):
            if not file_path.is_file():
                continue
            if file_path.suffix.lower() == ".h5ad":
                continue
            if "_goelite_downloads" in file_path.parts:
                continue
            handle.write(file_path, arcname=str(file_path.relative_to(source_dir)))
    return archive_path


def run_cellharmony_pipeline(job_id: str, store: JobStore, registry_path: Path) -> Dict[str, Path]:
    """
    Execute the production cellHarmony-lite + approximate UMAP workflow for a job.

    Returns
    -------
    Dict[str, Path]
        Mapping of artifact identifiers to their on-disk locations.
    """

    meta = store.get_job(job_id)
    reference_entry = _lookup_reference(meta["species"], meta["reference"], registry_path)
    _ensure_reference_fields(reference_entry)
    store.append_log(job_id, "Reference metadata loaded.")
    default_gene = _first_reference_gene(reference_entry["states_tsv"])

    uploads_dir = store.uploads_dir(job_id)
    outputs_dir = store.outputs_dir(job_id)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    h5_files, h5ad_file = _split_uploads(meta["files"], uploads_dir)
    qc = meta.get("qc", {})

    store.append_log(job_id, "Running cellHarmony_lite pipeline.")
    cellHarmony_lite.combine_and_align_h5(
        h5_files=h5_files,
        h5ad_file=h5ad_file,
        cellharmony_ref=reference_entry["states_tsv"],
        output_dir=str(outputs_dir),
        export_cptt=False,
        export_h5ad=False,
        min_genes=int(qc.get("min_genes", 200)),
        min_cells=int(qc.get("min_cells", 3)),
        min_counts=int(qc.get("min_counts", 500)),
        mit_percent=int(qc.get("mit_percent", 10)),
        generate_umap=False,
        save_adata=True,
        unsupervised_cluster=False,
        alignment_mode="classic",
        min_alignment_score=float(qc.get("align_cutoff", 0.1)),
        gene_translation_file=None,
        metacell_align=False,
        ambient_correct_cutoff=None,
    )

    assignments_path = outputs_dir / "cellHarmony_lite_assignments.txt"
    combined_h5ad_path = outputs_dir / "combined_with_umap_and_markers.h5ad"
    if not assignments_path.exists():
        raise FileNotFoundError("cellHarmony-lite assignments file missing.")
    if not combined_h5ad_path.exists():
        raise FileNotFoundError("Combined AnnData output missing; ensure save_adata=True.")

    query_cluster_key = Path(reference_entry["states_tsv"]).stem
    reference_cluster_key = reference_entry.get("cluster_key", query_cluster_key)

    store.append_log(job_id, "Running approximate UMAP placement.")
    reference_adata = approx_mod._load_reference_from_tsv(
        reference_entry["reference_coords_tsv"],
        reference_entry["reference_clusters_tsv"],
        umap_key="X_umap",
        cluster_key=reference_cluster_key,
    )
    query_adata = approx_mod._load_query_from_tsv(
        assignments_path,
        cluster_key=query_cluster_key,
    )
    approx_result = approx_mod.approximate_umap(
        query=query_adata,
        reference=reference_adata,
        query_cluster_key=query_cluster_key,
        reference_cluster_key=reference_cluster_key,
        umap_key="X_umap",
        jitter=0.05,
        num_reference_cells=1,
    )

    prefix = outputs_dir / "approximate"
    prefix.mkdir(parents=True, exist_ok=True)
    text_outputs = approx_result.write_text_outputs(prefix)
    annotated_pdf, plain_pdf = approx_result.write_comparison_pdf(
        reference_adata=reference_adata,
        umap_key="X_umap",
        reference_cluster_key=reference_cluster_key,
        query_cluster_key=query_cluster_key,
        output_path=prefix.with_name(prefix.name + "-comparison.pdf"),
    )

    artifacts = {
        "assignments": assignments_path,
        "combined_h5ad": combined_h5ad_path,
        "umap_coordinates": Path(text_outputs["coordinates"]),
        "umap_augmented": Path(text_outputs["augmented"]),
        "umap_placeholder_expression": Path(text_outputs["placeholder_expression"]),
        "umap_pdf": Path(annotated_pdf),
        "umap_pdf_plain": Path(plain_pdf),
    }
    for key, path in artifacts.items():
        store.add_artifact(job_id, key, path)

    sample_names = _job_sample_names(meta)
    population_columns = _candidate_population_columns(
        combined_h5ad_path,
        preferred=[query_cluster_key, reference_cluster_key],
    )
    default_population_col = next(
        (
            candidate["value"]
            for candidate in population_columns
            if candidate["value"] in {query_cluster_key, reference_cluster_key}
        ),
        population_columns[0]["value"] if population_columns else None,
    )
    differential_options = {
        "enabled": len(sample_names) >= 2,
        "sample_names": sample_names,
        "population_columns": population_columns,
        "default_population_col": default_population_col,
    }

    store.update_job(
        job_id,
        cluster_key=query_cluster_key,
        reference_cluster_key=reference_cluster_key,
        reference_coords_tsv=reference_entry["reference_coords_tsv"],
        reference_clusters_tsv=reference_entry["reference_clusters_tsv"],
        default_gene=default_gene,
        differential_options=differential_options,
        differential=_default_differential_state(len(sample_names) >= 2, default_population_col),
        message="Approximate UMAP completed.",
    )
    store.append_log(job_id, "cellHarmony-lite pipeline finished.")
    return artifacts


def run_cellharmony_differential(job_id: str, store: JobStore) -> Dict[str, object]:
    meta = store.get_job(job_id)
    differential_meta = dict(meta.get("differential") or {})
    config = dict(differential_meta.get("config") or {})

    group1_samples = [str(sample).strip() for sample in config.get("group1_samples", []) if str(sample).strip()]
    group2_samples = [str(sample).strip() for sample in config.get("group2_samples", []) if str(sample).strip()]
    population_col = str(config.get("population_col", "")).strip()

    if not population_col:
        raise ValueError("Differential analysis requires a cell-state field.")
    if not group1_samples or not group2_samples:
        raise ValueError("Both differential groups must include at least one sample.")
    overlap = sorted(set(group1_samples) & set(group2_samples))
    if overlap:
        raise ValueError(f"Samples cannot appear in both groups: {', '.join(overlap)}")

    combined_h5ad_value = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
    if not combined_h5ad_value:
        raise FileNotFoundError("Combined AnnData output unavailable for differential analysis.")
    combined_h5ad = Path(combined_h5ad_value)
    if not combined_h5ad.exists():
        raise FileNotFoundError("Combined AnnData output unavailable for differential analysis.")

    comparison_tag = _comparison_tag(group1_samples, group2_samples)
    run_id = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    run_root = store.outputs_dir(job_id) / "cellHarmony_differential" / f"{run_id}_{NetPerspective.safe_component(comparison_tag)}"
    heatmap_dir = run_root / "heatmaps"
    deg_dir = run_root / "DEGs"
    cellfreq_dir = run_root / "cell-frequency"
    goelite_dir = run_root / "GeneSetEnrichment"
    network_dir = run_root / "interaction-plots"
    for path in (heatmap_dir, deg_dir, cellfreq_dir, goelite_dir, network_dir):
        path.mkdir(parents=True, exist_ok=True)

    case_label = _group_display_label(group1_samples)
    control_label = _group_display_label(group2_samples)

    _update_differential_state(
        store,
        job_id,
        run_id=run_id,
        comparison_tag=comparison_tag,
        case_label=case_label,
        control_label=control_label,
        message="Loading aligned AnnData for differential analysis.",
        progress=15,
    )
    print(f"[INFO] Running cellHarmony-differential for {case_label} vs {control_label}")

    adata = ad.read_h5ad(combined_h5ad)
    if population_col not in adata.obs.columns:
        raise ValueError(f"'{population_col}' is not present in the aligned AnnData observations.")

    sample_col, resolved_samples = _resolve_samples_for_adata(
        adata,
        meta,
        group1_samples + group2_samples,
    )
    resolved_group1 = [resolved_samples[sample] for sample in group1_samples]
    resolved_group2 = [resolved_samples[sample] for sample in group2_samples]
    sample_values = adata.obs[sample_col].astype(str)
    print(
        f"[INFO] Differential sample matching uses obs['{sample_col}'] with "
        f"group1={resolved_group1} and group2={resolved_group2}"
    )

    subset_mask = sample_values.isin(resolved_group1 + resolved_group2)
    if int(np.asarray(subset_mask).sum()) == 0:
        raise ValueError("No cells were found for the selected sample groups.")
    adata = adata[subset_mask].copy()
    web_group_col = "__cellharmony_web_group__"
    sample_values = adata.obs[sample_col].astype(str)
    adata.obs[web_group_col] = np.where(sample_values.isin(resolved_group1), case_label, control_label)

    _update_differential_state(
        store,
        job_id,
        message="Computing differential expression across aligned cell states.",
        progress=35,
    )
    de_store = cellHarmony_differential.run_de_for_comparisons(
        adata=adata,
        population_col=population_col,
        covariate_col=web_group_col,
        case_label=case_label,
        control_label=control_label,
        method="wilcoxon",
        alpha=0.05,
        fc_thresh=1.2,
        min_cells_per_group=20,
        use_rawp=False,
        progress_callback=None,
    )
    pop_order = cellHarmony_differential._extract_lineage_order(adata, population_col)
    if not pop_order:
        detailed = de_store.get("detailed_deg")
        if detailed is not None and "population_or_pattern" in detailed.columns:
            pop_order = sorted(detailed["population_or_pattern"].astype(str).unique().tolist())
        else:
            pop_order = []

    cellfreq_outputs = cellHarmony_differential._write_cell_frequency_plots(
        adata=adata,
        population_col=population_col,
        covariate_col=web_group_col,
        conditions=[control_label, case_label],
        outdir=str(cellfreq_dir),
        comparison_label=f"{case_label} vs {control_label}",
        comparison_prefix=comparison_tag,
        pop_order=pop_order,
    ) or {}

    _update_differential_state(
        store,
        job_id,
        message="Building GO-enriched cellHarmony heatmap.",
        progress=60,
    )
    background_genes = adata.var_names.astype(str)
    if "gene_symbols" in adata.var.columns:
        background_genes = adata.var["gene_symbols"].astype(str)
    elif "features" in adata.var.columns:
        background_genes = adata.var["features"].astype(str)

    goelite_payload = cellHarmony_differential.run_goelite_for_clusters(
        de_store=de_store,
        background_genes=background_genes,
        outdir=str(goelite_dir),
        comparison_tag=comparison_tag,
        species=meta.get("species"),
    )
    if goelite_payload and goelite_payload.get("runner") and goelite_payload.get("prepared"):
        cellHarmony_differential.run_goelite_for_clusters_directional(
            de_store=de_store,
            background_genes=background_genes,
            outdir=str(goelite_dir / "up-regulated"),
            comparison_tag=f"{comparison_tag}_up",
            species=meta.get("species"),
            direction="up",
            runner=goelite_payload["runner"],
            prepared=goelite_payload["prepared"],
        )
        cellHarmony_differential.run_goelite_for_clusters_directional(
            de_store=de_store,
            background_genes=background_genes,
            outdir=str(goelite_dir / "down-regulated"),
            comparison_tag=f"{comparison_tag}_down",
            species=meta.get("species"),
            direction="down",
            runner=goelite_payload["runner"],
            prepared=goelite_payload["prepared"],
        )
    show_go_terms = bool(goelite_payload and goelite_payload.get("results") is not None and not goelite_payload["results"].empty)

    heatmap_pdf_name = f"heatmap_{comparison_tag}_by_{population_col}.pdf"
    heatmap_tsv_name = f"heatmap_{comparison_tag}_by_{population_col}.tsv"
    heatmap_pdf, heatmap_tsv = cellHarmony_differential.build_fixed_order_heatmap(
        de_store,
        str(heatmap_dir),
        population_col,
        1.2,
        heatmap_png=heatmap_pdf_name,
        heatmap_tsv=heatmap_tsv_name,
        goelite_payload=goelite_payload,
        show_go_terms=show_go_terms,
    )
    heatmap_pdf = Path(heatmap_pdf)
    heatmap_png = heatmap_pdf.with_suffix(".png")
    heatmap_tsv = Path(heatmap_tsv)
    heatmap_svg = heatmap_pdf.with_suffix(".svg")

    detailed_deg = de_store.get("detailed_deg", pd.DataFrame())
    summary_deg = de_store.get("summary_per_population", pd.DataFrame())
    assigned_deg = de_store.get("assigned_groups", pd.DataFrame())
    pooled_deg = de_store.get("pooled_overall")
    coreg_deg = de_store.get("coreg_pooled")

    deg_artifacts: Dict[str, Path] = {}
    named_tables = {
        f"DEG_detailed_{comparison_tag}.tsv": detailed_deg,
        f"DEG_summary_{comparison_tag}.tsv": summary_deg,
        f"DEG_assigned_groups_{comparison_tag}.tsv": assigned_deg,
    }
    for filename, frame in named_tables.items():
        if frame is None or frame.empty:
            continue
        path = deg_dir / filename
        frame.to_csv(path, sep="\t", index=False)
        deg_artifacts[path.stem] = path
        print(f"[INFO] Wrote {path}")
    if pooled_deg is not None and not pooled_deg.empty:
        path = deg_dir / f"DEG_pooled_overall_{comparison_tag}.tsv"
        pooled_deg.to_csv(path, sep="\t")
        deg_artifacts[path.stem] = path
        print(f"[INFO] Wrote {path}")
    if coreg_deg is not None and not coreg_deg.empty:
        path = deg_dir / f"DEG_coreg_pooled_{comparison_tag}.tsv"
        coreg_deg.to_csv(path, sep="\t", index=False)
        deg_artifacts[path.stem] = path
        print(f"[INFO] Wrote {path}")
    differential_h5ad_path = deg_dir / f"differentials_only_{comparison_tag}.h5ad"
    cellHarmony_differential.write_differentials_only_h5ad(adata, de_store, str(differential_h5ad_path))

    _update_differential_state(
        store,
        job_id,
        message="Rendering interaction networks for aligned cell states.",
        progress=80,
    )
    networks: List[Dict[str, str]] = []
    try:
        interactions_df = NetPerspective.load_interaction_data()
    except Exception as exc:
        interactions_df = None
        print(f"[WARN] Interaction network export skipped: {exc}")

    if interactions_df is not None and detailed_deg is not None and not detailed_deg.empty:
        used_ids: set[str] = set()
        for index, (population, pop_df) in enumerate(detailed_deg.groupby("population"), start=1):
            stats_df = pop_df.copy().reset_index(drop=True)
            stats_df["gene"] = stats_df["gene"].astype(str)
            stats_df["log2fc"] = pd.to_numeric(stats_df.get("log2fc"), errors="coerce")
            keep_columns = [col for col in ("gene", "log2fc", "fdr", "pval") if col in stats_df.columns]
            selected = (
                stats_df.loc[:, keep_columns]
                .dropna(subset=["gene", "log2fc"])
                .drop_duplicates(subset=["gene"])
            )
            if selected.empty or selected["gene"].nunique() < 2:
                continue
            network_id = NetPerspective.safe_component(str(population), fallback=f"network_{index}")
            if network_id in used_ids:
                network_id = f"{network_id}_{index}"
            used_ids.add(network_id)
            output_prefix = network_dir / network_id
            try:
                pdf_path, png_path, tsv_path = NetPerspective.generate_network_for_genes(
                    selected,
                    interactions_df,
                    str(output_prefix),
                    gene_column="gene",
                    fold_change_column="log2fc",
                    pval_column="fdr" if "fdr" in selected.columns else ("pval" if "pval" in selected.columns else None),
                    max_genes=None,
                )
            except NetPerspective.NetworkGenerationError as exc:
                print(f"[WARN] No interaction edges found for {population}; skipping network plot. {exc}")
                continue
            except ImportError as exc:
                print(f"[WARN] Interaction network export disabled: {exc}")
                break
            networks.append(
                {
                    "id": network_id,
                    "population": str(population),
                    "pdf": pdf_path,
                    "png": png_path,
                    "tsv": tsv_path,
                }
            )

    manifest = {
        "run_id": run_id,
        "population_col": population_col,
        "case_label": case_label,
        "control_label": control_label,
        "group1_samples": group1_samples,
        "group2_samples": group2_samples,
        "go_terms_included": show_go_terms,
        "network_count": len(networks),
    }
    manifest_path = run_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    archive_path = run_root.parent / f"{run_root.name}.zip"
    _write_zip_archive(run_root, archive_path)
    print(f"[INFO] Wrote differential results archive: {archive_path}")

    differential_artifacts: Dict[str, Path] = {
        "archive": archive_path,
        "heatmap_pdf": heatmap_pdf,
        "heatmap_png": heatmap_png,
        "heatmap_svg": heatmap_svg,
        "heatmap_tsv": heatmap_tsv,
        "manifest": manifest_path,
        "differentials_only_h5ad": differential_h5ad_path,
    }
    goelite_tsv = goelite_dir / f"GOElite_{NetPerspective.safe_component(comparison_tag)}.tsv"
    if goelite_tsv.exists():
        differential_artifacts["goelite_tsv"] = goelite_tsv
    for key, path_str in cellfreq_outputs.items():
        differential_artifacts[f"cell_frequency_{key}"] = Path(path_str)
    differential_artifacts.update(deg_artifacts)

    _update_differential_state(
        store,
        job_id,
        run_id=run_id,
        comparison_tag=comparison_tag,
        case_label=case_label,
        control_label=control_label,
        config={
            "population_col": population_col,
            "group1_samples": group1_samples,
            "group2_samples": group2_samples,
        },
        artifacts={key: str(path) for key, path in differential_artifacts.items()},
        networks=networks,
        go_terms_included=show_go_terms,
        message="cellHarmony-differential finished.",
        progress=95,
    )
    return {
        "artifacts": differential_artifacts,
        "networks": networks,
        "run_id": run_id,
        "go_terms_included": show_go_terms,
    }
