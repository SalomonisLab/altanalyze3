#!/usr/bin/env python3
"""GO-Elite enrichment for UDON clusters. Takes the per-cluster marker genes
(uns['udon_marker_genes_top_n']: columns 'marker','top_cluster') and runs the
altanalyze3 GO-Elite engine (same API cellHarmony uses), writing one combined
TSV (and a scatter PDF if available) to the UDON results folder."""
import os, sys, pandas as pd

def _import_goelite():
    try:
        from altanalyze3.components.goelite import GOEliteRunner, EnrichmentSettings, prepare_species_resources
        return GOEliteRunner, EnrichmentSettings, prepare_species_resources
    except Exception:
        here = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.abspath(os.path.join(here, "..", "..", ".."))  # dir containing the altanalyze3 package
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
        from altanalyze3.components.goelite import GOEliteRunner, EnrichmentSettings, prepare_species_resources
        return GOEliteRunner, EnrichmentSettings, prepare_species_resources

_SPECIES = {"Hs": "human", "Mm": "mouse", "human": "human", "mouse": "mouse"}


def run_goelite_on_udon(adata, outdir, species="Hs", background_genes=None,
                        marker_key="udon_marker_genes_top_n", min_term_size=5,
                        max_term_size=2000, cache_dir=None, logger=print):
    """Run GO-Elite per UDON cluster. Returns the combined results DataFrame (or None)."""
    try:
        GOEliteRunner, EnrichmentSettings, prepare_species_resources = _import_goelite()
    except Exception as e:
        logger(f"[goelite] SKIPPED — could not import GO-Elite engine: {e}")
        return None

    markers = adata.uns.get(marker_key)
    if markers is None or len(markers) == 0:
        logger("[goelite] SKIPPED — no marker genes found")
        return None
    markers = pd.DataFrame(markers)
    if background_genes is None:
        background_genes = list(adata.varm["pseudobulk_folds"].index) if "pseudobulk_folds" in adata.varm \
            else list(adata.var_names)

    species_key = _SPECIES.get(species, "human")
    logger(f"[goelite] preparing {species_key} resources ...")
    parsed = prepare_species_resources(species_key, cache_dir=cache_dir)
    settings = EnrichmentSettings(min_term_size=int(min_term_size), max_term_size=int(max_term_size))
    runner = GOEliteRunner(parsed, settings=settings)
    prepared = runner.prepare_background([str(g) for g in pd.Index(background_genes).dropna().unique()])

    os.makedirs(outdir, exist_ok=True)
    genes_dir = os.path.join(outdir, "cluster_gene_lists")
    os.makedirs(genes_dir, exist_ok=True)

    rows = []
    for cluster, grp in markers.groupby("top_cluster"):
        genes = [str(g) for g in grp["marker"].tolist()]
        with open(os.path.join(genes_dir, f"cluster_{cluster}_markers.txt"), "w") as fh:
            fh.write("\n".join(genes) + "\n")
        try:
            results = runner.run_prepared(genes, prepared, apply_prioritization=True)
        except Exception as e:
            logger(f"[goelite] cluster {cluster}: enrichment error {e}")
            continue
        n_sel = 0
        for res in results:
            node = runner.go_tree.get(res.term_id) if hasattr(runner, "go_tree") else None
            sel = bool(getattr(res, "selected", False))
            n_sel += int(sel)
            rows.append({
                "cluster": cluster,
                "term_id": res.term_id,
                "term_name": getattr(node, "name", "") if node else "",
                "namespace": getattr(node, "namespace", "") if node else "",
                "p_value": getattr(res, "p_value", None),
                "fdr": getattr(res, "fdr", None),
                "z_score": getattr(res, "z_score", None),
                "overlap": getattr(res, "overlap", None),
                "term_size": getattr(res, "total_genes", getattr(res, "term_genes", None)),
                "query_size": len(genes),
                "selected": sel,
                "overlap_genes": ",".join(getattr(res, "overlap_genes", []) or []),
            })
        logger(f"[goelite] cluster {cluster}: {len(genes)} genes -> {n_sel} prioritized terms")

    if not rows:
        logger("[goelite] no enrichment results")
        return None
    df = pd.DataFrame(rows)
    out_tsv = os.path.join(outdir, "GOElite_UDON.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    sel_tsv = os.path.join(outdir, "GOElite_UDON_selected.tsv")
    df[df["selected"]].to_csv(sel_tsv, sep="\t", index=False)
    logger(f"[goelite] wrote {out_tsv} ({len(df)} rows, {int(df['selected'].sum())} selected)")

    # optional scatter PDF (same util cellHarmony uses), best-effort
    try:
        from altanalyze3.components.goelite.plotting import write_goelite_scatter_pdf
        pdf_df = df.rename(columns={"cluster": "population"})
        write_goelite_scatter_pdf(pdf_df, os.path.join(outdir, "GOElite_UDON.pdf"), title_prefix="GO-Elite (UDON)")
        logger("[goelite] wrote GOElite_UDON.pdf")
    except Exception as e:
        logger(f"[goelite] PDF skipped: {e}")
    return df
