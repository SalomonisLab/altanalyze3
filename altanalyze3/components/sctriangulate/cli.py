#!/usr/bin/env python3
"""
scTriangulate CLI for AltAnalyze3 -- reconcile conflicting cluster annotations.

Runs the real entry point, ``ScTriangulate(...).lazy_run()``. Nothing here
reimplements the workflow; the flags map onto lazy_run's own parameters.

Example:
  python3 -m altanalyze3.components.sctriangulate.cli \\
      --h5ad   /path/to/input.h5ad \\
      --query  sctri_rna_leiden_1,sctri_rna_leiden_2,sctri_rna_leiden_3 \\
      --outdir /path/to/output

Outputs land in --outdir:
  sctriangulate.h5ad              annotated object (obs carries 'pruned', 'raw',
                                  'final_annotation', 'confidence', 'prefixed')
  sctri_barcode2cellmetadata.txt  per-cell metadata table
  raw_cluster_goodness.txt        per raw cluster: size, win fraction
  stability_<annotation>.txt      per-cluster stability metrics
  sctri_gene_to_df_*.txt          marker and exclusive gene tables
  umap_sctriangulate_*.pdf        UMAPs of the final and pruned labels
  after_metrics.p / after_shapley.p / after_rank_pruning.p   resumable checkpoints
  run_parameters.json             every parameter and version used by this run

FLAGS THAT CHANGE THE RESULT
----------------------------
Everything else only changes how the same answer is computed. These four change
what the answer is, so each one records what it did in run_parameters.json and
each one belongs in a methods section:

  --downsample N          cap cells per cluster. Cell counts are the denominator
                          of every stability metric.
  --prefilter-markers N   cap genes. Selects genes with the labels it then
                          scores, so it is circular; measured to move the demo's
                          raw labels to ARI 0.79 at N=200.
  --scale-sccaf           mean-centre before the SCCAF fit. Also forces the
                          legacy dense SCCAF path.
  --criterion 1..6        which gene classes count as artifacts.

--keep-layers changes the OUTPUT FILE but not the result: without it, layers
that --layer does not name are dropped from memory and from sctriangulate.h5ad.

WHAT run_parameters.json RECORDS
--------------------------------
Every flag above, plus: n_cells_input / n_genes_input against
n_cells_analysed / n_genes_analysed, which differ when a reduction ran;
downsample_info and prefilter_info, the full per-annotation record of each
reduction; dropped_layers; the resolved reference annotation; cluster counts per
annotation; wall_seconds; n_pruned_clusters; the exact command line; and the
version of python, numpy, scipy, pandas, scanpy, anndata, sklearn and gseapy.

DEFAULTS THAT DIFFER FROM UPSTREAM scTriangulate
------------------------------------------------
  enrichment              off here, always on upstream. It feeds no metric and
                          no Shapley decision. Turn it on with --enrichment.
  sccaf-mode              'optimized' here. Same accuracies as upstream's dense
                          procedure, and seeded, which upstream was not.
  sparse compute          on. Upstream densified X for every metric.
  layers                  dropped unless --layer names them. Use --keep-layers.
See README.md for the measurement behind each one.
"""

import os
import sys
import json
import time
import argparse
import platform


def build_parser():
    p = argparse.ArgumentParser(
        description='scTriangulate: mix and match conflicting single cell cluster annotations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--h5ad', required=True,
                   help='Input AnnData .h5ad. adata.obs must carry one column per annotation.')
    p.add_argument('--query', required=True,
                   help='Comma-separated obs column names to reconcile.')
    p.add_argument('--outdir', required=True, help='Output directory (created if absent).')
    p.add_argument('--reference', default=None,
                   help='Annotation used as the reference. Default: the first --query entry.')
    p.add_argument('--species', default='human', choices=['human', 'mouse'])
    p.add_argument('--criterion', type=int, default=2, choices=[1, 2, 3, 4, 5, 6],
                   help='Which artifact gene classes count as artifacts. See metrics.read_artifact_genes.')
    p.add_argument('--layer', default=None,
                   help='adata layer holding raw counts, for tf-idf when X is skewed.')
    p.add_argument('--cores', type=int, default=None,
                   help='Worker processes for the per-annotation metric step. Default: min(n_annotations, n_cpus).')
    p.add_argument('--sequential', action='store_true',
                   help='Compute metrics one annotation at a time (lower peak memory).')
    p.add_argument('--shapley-mode', default=None,
                   choices=['shapley_all_or_none', 'shapley', 'rank', 'rank_all_or_none'],
                   help='Default: shapley_all_or_none for <=15 annotations, else rank.')
    p.add_argument('--shapley-bonus', type=float, default=0.01)
    p.add_argument('--win-fraction-cutoff', type=float, default=0.25)
    p.add_argument('--reassign-abs-thresh', type=int, default=10)
    p.add_argument('--keep-layers', action='store_true',
                   help='Keep every adata layer in memory and in the output h5ad. By default '
                        'only the layer named by --layer is kept, because nothing else reads '
                        'them and they sit in memory through the whole run.')
    p.add_argument('--downsample', type=int, default=None, metavar='N',
                   help='Keep at most N cells per cluster, sharing one whitelist across all '
                        'annotations: the first --query annotation draws its full N per cluster, '
                        'and every later annotation draws only the shortfall its clusters still '
                        'have. Off by default. Cell counts are the denominator of every stability '
                        'metric, so this changes the result; use it to make a run fit, not to '
                        'make it faster.')
    p.add_argument('--prefilter-markers', type=int, default=None, metavar='N',
                   help='Before any metric, keep only the union of the top N MarkerFinder genes '
                        'per cluster of each annotation. Off by default, and not needed for '
                        'memory since the sparse path replaced the dense one. It selects genes '
                        'using the same labels it then scores, and on the demo N=200 moved the '
                        'raw labels to ARI 0.79; see benchmarks/prefilter_validate.py.')
    p.add_argument('--prefilter-cells', type=int, default=1000, metavar='N',
                   help='Cells per cluster used to rank markers for --prefilter-markers. '
                        'Ranking needs enough cells per state, not every cell.')
    p.add_argument('--subsample-seed', type=int, default=0,
                   help='Seed for --downsample and --prefilter-markers. Same seed, same cells.')
    p.add_argument('--sccaf-mode', default='optimized', choices=['optimized', 'legacy'],
                   help="'optimized' keeps the matrix sparse and fits the per-class "
                        "problems on threads (about 3x faster, about 20x less memory). "
                        "'legacy' is upstream's dense single-threaded procedure. Both are "
                        "seeded; upstream was not.")
    p.add_argument('--scale-sccaf', action='store_true',
                   help='Scale the matrix before the SCCAF logistic regression.')
    p.add_argument('--predict-doublet', action='store_true',
                   help='Run scrublet. Off by default; doublet score is then a constant 0.5, as upstream.')
    p.add_argument('--assess-raw', action='store_true')
    p.add_argument('--assess-pruned', action='store_true')
    p.add_argument('--viewer-cluster', action='store_true')
    p.add_argument('--viewer-heterogeneity', action='store_true')
    a = p.add_mutually_exclusive_group()
    a.add_argument('--annotate', dest='annotate', action='store_true', default=True,
                   help='After pruning, run MarkerFinder on the pruned clusters, name every '
                        'cluster by GO-Elite BioMarkers cell-state enrichment, order the '
                        'centroids with HOPACH, and redraw MarkerFinder in that order. On by '
                        'default. Writes obs["cluster_name"], obs["hopach_cluster"] and '
                        'uns["lineage_order"], and an antibody-only heatmap when the object '
                        'carries AB_ features.')
    a.add_argument('--no-annotate', dest='annotate', action='store_false',
                   help='Stop after pruning, as the pipeline did before this step became a '
                        'default. obs["pruned"] is then the last output.')
    p.add_argument('--annotate-lead', default=None, metavar='OBS_COLUMN',
                   help='obs column holding a trusted annotation, e.g. a published cell type '
                        'that also competes in --query. Its dominant label then names the '
                        'cluster and the enriched term becomes the fallback.')
    p.add_argument('--annotate-top-n', type=int, default=60,
                   help='Markers per cluster for both MarkerFinder passes.')
    p.add_argument('--annotate-cells-per-cluster', type=int, default=100,
                   help='Cells drawn per cluster in the heatmaps.')
    p.add_argument('--annotate-hopach-distance', default='cor',
                   choices=['cor', 'abscor', 'cosangle', 'abscosangle', 'euclid', 'manhattan'],
                   help='HOPACH distance over the cluster centroids.')
    p.add_argument('--annotate-max-fdr', type=float, default=1e-5,
                   help='Reject a cell-state term at a worse FDR than this and name the cluster '
                        'from its strongest marker instead. The enrichment always returns its '
                        'best term, however weak, which invents cell types. Measured over 79 '
                        'clusters: every term 0.05 accepted and 1e-5 rejected overlapped by only '
                        '2-8 of about 36-60 markers.')
    p.add_argument('--annotate-min-overlap', type=int, default=5,
                   help='Reject a cell-state term overlapping the cluster markers by fewer genes '
                        'than this, however small its FDR.')
    p.add_argument('--annotate-biomarker-file', default=None,
                   help='Override the GO-Elite BioMarkers table. Defaults to the EnsMart72 table '
                        'for --species.')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--enrichment', dest='enrichment', action='store_true', default=False,
                   help='Add the artifact enrichr/GSEA annotation of marker genes. Off by '
                        'default: it feeds no stability metric and no Shapley decision, it '
                        'is read only by plot_cluster_feature, plot_multi_modal_feature_rank, '
                        'penalize_artifact(mode=cellcycle) and the viewer, and the GSEA '
                        'permutations cost about 3.7 s of a 21 s demo run. Upstream always '
                        'computed it. Required by --viewer-cluster.')
    g.add_argument('--no-enrichment', dest='enrichment', action='store_false',
                   help='Explicitly skip it. Already the default.')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    import matplotlib
    matplotlib.use('Agg')
    import scanpy as sc

    from . import ScTriangulate
    from . import metrics as M

    query = [q.strip() for q in args.query.split(',') if q.strip()]
    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.h5ad)

    # sc.read materialises every layer. Only --layer is ever read, by
    # tf_idf10_for_cluster; the rest sit in memory through the whole run, and the
    # peak comes later in the metric step, not here. Dropping them does not lower
    # the read peak, it lowers everything after it. It also drops them from the
    # written h5ad, which is why --keep-layers exists.
    dropped_layers = []
    if not args.keep_layers:
        dropped_layers = [k for k in list(adata.layers.keys()) if k != args.layer]
        for k in dropped_layers:
            del adata.layers[k]
        if dropped_layers:
            print('freed {} unused layer(s), also absent from the output h5ad: {}'.format(
                len(dropped_layers), ', '.join(dropped_layers)))

    missing = [q for q in query if q not in adata.obs.columns]
    if missing:
        raise SystemExit('these --query columns are not in adata.obs: {}\navailable: {}'.format(
            missing, list(adata.obs.columns)))

    enrichment = args.enrichment
    if args.viewer_cluster and not enrichment:
        # the viewer's per-cluster page plots the enrichment panel
        print('--viewer-cluster needs the enrichment columns; turning --enrichment on')
        enrichment = True
    M.RUN_ENRICHMENT = enrichment

    params = {
        'input_h5ad': os.path.abspath(args.h5ad),
        'outdir': os.path.abspath(args.outdir),
        'query': query,
        'reference': args.reference or query[0],
        'species': args.species,
        'criterion': args.criterion,
        'layer': args.layer,
        'cores': args.cores,
        'compute_metrics_parallel': not args.sequential,
        'shapley_mode': args.shapley_mode,
        'shapley_bonus': args.shapley_bonus,
        'win_fraction_cutoff': args.win_fraction_cutoff,
        'reassign_abs_thresh': args.reassign_abs_thresh,
        'scale_sccaf': args.scale_sccaf,
        'predict_doublet': args.predict_doublet,
        'run_enrichment': enrichment,
        'sccaf_mode': args.sccaf_mode,
        'downsample': args.downsample,
        'prefilter_markers': args.prefilter_markers,
        'prefilter_cells': args.prefilter_cells,
        'subsample_seed': args.subsample_seed,
        'keep_layers': args.keep_layers,
        'dropped_layers': dropped_layers,
        'n_cells_input': int(adata.shape[0]),
        'n_genes_input': int(adata.shape[1]),
        'n_clusters_per_annotation': {q: int(adata.obs[q].nunique()) for q in query},
        'python': sys.version.split()[0],
        'platform': platform.platform(),
        'command': ' '.join(sys.argv),
    }
    for mod in ['numpy', 'scipy', 'pandas', 'scanpy', 'anndata', 'sklearn', 'gseapy']:
        try:
            params['version_' + mod] = __import__(mod).__version__
        except Exception:
            params['version_' + mod] = 'not installed'

    t0 = time.perf_counter()
    sctri = ScTriangulate(dir=args.outdir, adata=adata, query=query,
                          reference=args.reference, species=args.species,
                          criterion=args.criterion, predict_doublet=args.predict_doublet)
    sctri.run_enrichment = enrichment
    sctri.sccaf_mode = args.sccaf_mode

    # both reductions run before lazy_run, and both change what the metrics see
    if args.downsample is not None:
        info = sctri.downsample_cells(max_per_cluster=args.downsample,
                                      seed=args.subsample_seed)
        params['downsample_info'] = info
        print('downsample: {} of {} cells kept ({:.1f}x fewer)'.format(
            info['cells_after'], info['cells_before'], info['reduction']))
    if args.prefilter_markers is not None:
        info = sctri.prefilter_marker_genes(n_per_cluster=args.prefilter_markers,
                                            max_cells_per_cluster=args.prefilter_cells,
                                            seed=args.subsample_seed)
        params['prefilter_info'] = info
        print('prefilter: {} of {} genes kept ({:.1f}x fewer)'.format(
            info['genes_after'], info['genes_before'], info['reduction']))
    params['n_cells_analysed'] = int(sctri.adata.shape[0])
    params['n_genes_analysed'] = int(sctri.adata.shape[1])

    sctri.lazy_run(compute_metrics_parallel=not args.sequential,
                   compute_shapley_parallel=not args.sequential,
                   scale_sccaf=args.scale_sccaf,
                   layer=args.layer,
                   cores=args.cores,
                   shapley_mode=args.shapley_mode,
                   shapley_bonus=args.shapley_bonus,
                   win_fraction_cutoff=args.win_fraction_cutoff,
                   reassign_abs_thresh=args.reassign_abs_thresh,
                   assess_raw=args.assess_raw,
                   assess_pruned=args.assess_pruned,
                   viewer_cluster=args.viewer_cluster,
                   viewer_heterogeneity=args.viewer_heterogeneity)
    params['n_pruned_clusters'] = int(sctri.adata.obs['pruned'].nunique())

    # Steps 4-6, on by default. lazy_run stops at obs['pruned'], which names a source cluster and
    # no biology. annotate_pruned_clusters runs MarkerFinder, names each cluster by GO-Elite
    # BioMarkers cell-state enrichment, orders the centroids with HOPACH and redraws MarkerFinder
    # in that order. A failure here must not destroy the pruning result the run already earned, so
    # the error is recorded and re-raised only after sctriangulate.h5ad and run_parameters.json
    # are on disk.
    params['annotate'] = bool(args.annotate)
    annotate_error = None
    if args.annotate:
        from .annotate import annotate_pruned_clusters
        try:
            summary = annotate_pruned_clusters(
                sctri.adata,
                outdir=args.outdir,
                cluster_key='pruned',
                species=args.species,
                top_n=args.annotate_top_n,
                cells_per_cluster=args.annotate_cells_per_cluster,
                layer=args.layer,
                lead_annotation=args.annotate_lead,
                biomarker_file=args.annotate_biomarker_file,
                enrichment_max_fdr=args.annotate_max_fdr,
                enrichment_min_overlap=args.annotate_min_overlap,
                hopach_distance=args.annotate_hopach_distance,
                covariate_columns=query,
            )
            params['annotate_summary'] = {
                'n_named_clusters': len(summary['names']),
                'hopach_k': summary['hopach_k'],
                'hopach_mss': summary['hopach_mss'],
                'hopach_distance': args.annotate_hopach_distance,
                'lead_annotation': args.annotate_lead,
                'cluster_table': summary['cluster_table'],
                'name_source_counts': {
                    src: sum(1 for v in summary['name_source'].values() if v == src)
                    for src in sorted(set(summary['name_source'].values()))
                },
            }
            sctri.adata.write(os.path.join(args.outdir, 'sctriangulate.h5ad'))
        except Exception as exc:  # recorded, re-raised below; never silently swallowed
            annotate_error = exc
            params['annotate_error'] = '{}: {}'.format(type(exc).__name__, exc)
            print('[annotate] FAILED: {}: {}'.format(type(exc).__name__, exc))
            print('[annotate] the pruning result is intact; rerun with --no-annotate to skip '
                  'this stage, or fix the cause and rerun.')

    params['wall_seconds'] = time.perf_counter() - t0

    with open(os.path.join(args.outdir, 'run_parameters.json'), 'w') as f:
        json.dump(params, f, indent=1)

    print('\nfinished in {:.1f} s'.format(params['wall_seconds']))
    print('pruned clusters: {}'.format(params['n_pruned_clusters']))
    if args.annotate and annotate_error is None:
        print('named clusters: {} in {} HOPACH groups'.format(
            params['annotate_summary']['n_named_clusters'], params['annotate_summary']['hopach_k']))
    print('outputs in: {}'.format(os.path.abspath(args.outdir)))
    print('parameters: {}'.format(os.path.join(os.path.abspath(args.outdir), 'run_parameters.json')))
    if annotate_error is not None:
        raise annotate_error


if __name__ == '__main__':
    main()
