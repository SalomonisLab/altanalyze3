"""Command-line entry points for the SNAF neoantigen workflow.

Two subcommands, both cross-platform (macOS / Windows / Linux) and offline:

  altanalyze3 snaf-ts   Tumor-specificity scoring only. Needs just a junction-count
                        matrix and a normal-tissue control h5ad -- NO Alt91_db, NO
                        MHC binder, NO internet. Scores every neojunction with
                        'mean', 'mle' (closed-form) and/or 'bayesian' (BayesTS).

  altanalyze3 snaf      Full MHC-bound T-antigen pipeline: sifting -> in-silico
                        translation -> MHC binding (MHCflurry by default; netMHCpan
                        optional) -> immunogenicity (deepimmuno) -> burden/frequency
                        tables. Needs the reference DB (--db_dir with Alt91_db +
                        controls) and per-sample HLA types.

SNAF is imported lazily inside each function so `--help` and CLI startup stay fast
and don't require the optional heavy deps (tensorflow, torch/pyro, dash).
"""
import os
import sys
import json
import logging


def _ensure_fork_safety():
    """macOS: guarantee OBJC_DISABLE_INITIALIZE_FORK_SAFETY is in the environment by
    re-exec'ing ONCE before any heavy (TensorFlow/Accelerate) import, so the fork-based
    parallel binding/immunogenicity workers don't abort or deadlock. No-op elsewhere or if
    already set. Kept dependency-free (os/sys only) so it can run before importing snaf."""
    if sys.platform == 'darwin' and os.environ.get('OBJC_DISABLE_INITIALIZE_FORK_SAFETY') != 'YES':
        os.environ['OBJC_DISABLE_INITIALIZE_FORK_SAFETY'] = 'YES'
        os.execv(sys.executable, [sys.executable] + sys.argv)

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- helpers
def _read_count_matrix(path, label='junction matrix'):
    """Read an AltAnalyze junction-count matrix, compressed or not, and report what it lost.

    * Compression is inferred from the extension, so `.txt`, `.txt.gz`, `.zip`, `.bz2`,
      `.xz` and `.zst` all work -- a cohort matrix can stay in the archive it shipped in.
    * A row whose field count exceeds the header's is skipped rather than aborting the
      whole read, and a row with too few fields is padded with NaN by the parser; both
      kinds are counted and reported here, never dropped in silence.
    """
    import pandas as pd
    import numpy as np
    kw = {'sep': '\t', 'index_col': 0}
    if str(path).endswith('.zip'):
        # pandas only reads a single-member zip; say so plainly instead of failing deep
        import zipfile
        with zipfile.ZipFile(path) as zf:
            members = [n for n in zf.namelist() if not n.endswith('/')]
        if len(members) != 1:
            raise ValueError('{}: a zipped count matrix must hold exactly one file, {} has '
                             '{}: {}'.format(label, path, len(members), members[:5]))
    try:
        df = pd.read_csv(path, on_bad_lines='warn', **kw)
    except TypeError:                      # pandas < 1.3
        df = pd.read_csv(path, error_bad_lines=False, warn_bad_lines=True, **kw)
    n_read = df.shape[0]
    bad = df.isna().any(axis=1)
    n_bad = int(bad.sum())
    if n_bad:
        logger.warning('%s %s: %d / %d rows carry a missing value (a short or corrupt line) '
                       'and are dropped; the first are %s', label, path, n_bad, n_read,
                       list(df.index[bad][:5]))
        df = df.loc[~bad, :]
    # AltAnalyze junction UIDs sometimes carry a trailing '=...'; strip and dedup
    df.index = [str(item).split('=')[0] for item in df.index]
    dup = df.index.duplicated()
    n_dup = int(np.count_nonzero(dup))
    if n_dup:
        logger.warning('%s %s: %d / %d rows are duplicate UIDs after stripping the coordinate '
                       'suffix; keeping the first of each', label, path, n_dup, df.shape[0])
    df = df.loc[np.logical_not(dup), :]
    logger.info('%s %s: %d rows x %d samples', label, path, df.shape[0], df.shape[1])
    return df


def _read_junction_matrix(path):
    return _read_count_matrix(path, label='junction matrix')


def _read_hla(path, samples):
    """Parse per-sample HLA types into the list-of-lists SNAF expects (aligned to
    the junction-matrix columns). Accepts a 2-column table `sample<TAB>hla1,hla2,...`
    or a table whose first column is the sample and remaining columns are HLAs."""
    import pandas as pd
    raw = pd.read_csv(path, sep='\t', index_col=0, header=0)
    hlas = []
    for s in samples:
        if s not in raw.index:
            raise ValueError("sample '{}' from the junction matrix is missing in --hla".format(s))
        row = raw.loc[s]
        if raw.shape[1] == 1:
            vals = str(row.iloc[0]).replace(';', ',').split(',')
        else:
            vals = [str(v) for v in row.values if str(v) not in ('nan', '')]
        hlas.append([v.strip() for v in vals if v.strip()])
    return hlas


def _require_file(path, flag):
    """Fail loudly and early if a required input path is missing."""
    if path is None:
        raise ValueError('{} is required'.format(flag))
    p = str(path)
    if not os.path.exists(p):
        raise FileNotFoundError('{} file not found: {}'.format(flag, p))
    return p


def _load_add_control(spec):
    """--add_control 'name=path.txt,name2=path.h5ad' -> {name: DataFrame|AnnData}.

    Text controls go through _read_count_matrix, so they may be gzipped/zipped and a
    corrupt row is reported and dropped instead of aborting the run."""
    if not spec:
        return None
    import anndata as ad
    out = {}
    for item in spec.split(','):
        name, path = item.split('=', 1)
        _require_file(path, "--add_control '{}'".format(name))
        if path.endswith('.h5ad'):
            out[name] = ad.read_h5ad(path)
        else:
            out[name] = _read_count_matrix(path, label="control '{}'".format(name))
    return out


# --------------------------------------------------------------------------- snaf-ts
def run_snaf_ts(args):
    """Tumor-specificity-only entry (DB-free, offline)."""
    _ensure_fork_safety()
    import numpy as np
    import pandas as pd
    from altanalyze3.components import snaf
    from altanalyze3.components.snaf import gtex as snaf_gtex

    outdir = str(args.output)
    os.makedirs(outdir, exist_ok=True)
    _require_file(args.juncounts, '--juncounts')
    control_h5ad = str(args.control_h5ad)
    if not os.path.exists(control_h5ad):
        from altanalyze3.components.snaf.reference import control_filename, download_command
        raise FileNotFoundError(
            "--control_h5ad not found: {}\nThe GTEx junction control ({}) ships in the SNAF "
            "reference bundle. Fetch it with:\n\n    {}\n".format(
                control_h5ad, control_filename('count'),
                download_command(os.path.dirname(control_h5ad) or 'snaf_reference')))
    df = _read_junction_matrix(args.juncounts)
    logger.info("junction matrix: %s", df.shape)

    add_control = _load_add_control(getattr(args, 'add_control', None))

    # configure controls (sets the gtex globals used by sifting + tumor_specificity). When a
    # precomputed stats table exists (auto-detected or --control_stats), the full control
    # matrix is NOT loaded and BayesTS is not re-run.
    snaf.gtex_configuration(
        df, control_h5ad, args.t_min, args.n_max, args.normal_cutoff, args.tumor_cutoff,
        args.normal_prevalance_cutoff, args.tumor_prevalance_cutoff, add_control=add_control,
        control_stats_path=(str(args.control_stats) if getattr(args, 'control_stats', None) else None))

    jcmq = snaf.JunctionCountMatrixQuery(
        junction_count_matrix=df, cores=args.cpus, add_control=add_control,
        not_in_db=False, outdir=outdir, filter_mode=args.filter_mode,
        min_samples=getattr(args, 'min_samples', 1),
        max_bayests_percentile=getattr(args, 'max_bayests_percentile', None))
    valid = list(jcmq.valid)
    logger.info("neojunctions passing sifting: %d", len(valid))

    # Vectorized tumor-specificity over all neojunctions (identical results to the
    # per-uid loop, but sparse/batched -- removes the 405k-iteration Python bottleneck).
    want = ['mean'] + (['mle'] if 'mle' in args.methods else [])
    ts = snaf.tumor_specificity_batch(valid, methods=tuple(want))
    out = pd.DataFrame(index=valid)
    out['mean_gtex_count'] = ts['mean'].values
    if 'mle' in args.methods:
        out['tumor_specificity_mle'] = ts['mle'].values
    if 'bayesian' in args.methods:
        sig = snaf.compute_bayests_sigma(snaf_gtex.adata, uids=valid, mode=args.bayes_mode,
                                         epoch=args.bayes_epoch)
        out['tumor_specificity_bayesian'] = out.index.map(sig['mean_sigma'].to_dict())
    out_path = os.path.join(outdir, 'neojunction_tumor_specificity.csv')
    out.to_csv(out_path)
    logger.info("wrote %s (%d neojunctions)", out_path, out.shape[0])
    print('SNAF tumor-specificity written to {}'.format(out_path))
    return out_path


# --------------------------------------------------------------------------- snaf (full T pipeline)
def run_snaf(args):
    """Full MHC-bound T-antigen pipeline entry."""
    _ensure_fork_safety()
    from altanalyze3.components import snaf

    outdir = str(args.output)
    os.makedirs(outdir, exist_ok=True)
    _require_file(args.juncounts, '--juncounts')
    _require_file(args.hla, '--hla')
    df = _read_junction_matrix(args.juncounts)
    samples = list(df.columns)
    hlas = _read_hla(args.hla, samples)
    add_control = _load_add_control(getattr(args, 'add_control', None))

    genome_fasta = str(args.genome_fasta) if getattr(args, 'genome_fasta', None) else None
    # Honor --gtex_db (initialize() now accepts it; it was previously discarded).
    gtex_db = str(args.gtex_db) if getattr(args, 'gtex_db', None) else None
    # Loudly flag the one step that can silently reach for the network: novel-exon /
    # UTR sequence retrieval falls back to the UCSC DAS web service when no local
    # genome FASTA is given and SNAF_OFFLINE != 1.
    if genome_fasta is None and os.environ.get('SNAF_OFFLINE', '0') != '1':
        logger.warning(
            'No --genome_fasta and SNAF_OFFLINE != 1: novel-exon/UTR sequence retrieval '
            'will query the UCSC DAS web service over the network (slow, non-reproducible). '
            'Pass --genome_fasta <hg38.fa> for a fully offline run, or set SNAF_OFFLINE=1 '
            'to forbid the network.')

    # initialize() resolves/validates (and optionally downloads) the reference and
    # configures everything.
    snaf.initialize(
        df=df, db_dir=str(args.db_dir), gtex_mode=args.gtex_mode,
        software_path=(str(args.software_path) if getattr(args, 'software_path', None) else None),
        binding_method=args.binding_method, t_min=args.t_min, n_max=args.n_max,
        normal_cutoff=args.normal_cutoff, tumor_cutoff=args.tumor_cutoff,
        normal_prevalance_cutoff=args.normal_prevalance_cutoff,
        tumor_prevalance_cutoff=args.tumor_prevalance_cutoff, add_control=add_control,
        genome_fasta=genome_fasta, gtex_db=gtex_db,
        download_ref=getattr(args, 'download_ref', False),
        control_stats_path=(str(args.control_stats) if getattr(args, 'control_stats', None) else None))

    jcmq = snaf.JunctionCountMatrixQuery(
        junction_count_matrix=df, cores=args.cpus, add_control=add_control,
        not_in_db=args.not_in_db, outdir=outdir, filter_mode=args.filter_mode,
        min_samples=getattr(args, 'min_samples', 1),
        max_bayests_percentile=getattr(args, 'max_bayests_percentile', None))
    import time as _t
    _tj = _t.time()
    jcmq.run(hlas=hlas, strict=args.strict, outdir=outdir, name='after_prediction.p')
    print('[STAGE-TIMING] jcmq.run total (translate+bind+immuno): {:.1f}s'.format(_t.time()-_tj)); _tg = _t.time()
    snaf.JunctionCountMatrixQuery.generate_results(
        path=os.path.join(outdir, 'after_prediction.p'), outdir=outdir)
    print('[STAGE-TIMING] generate_results (burden+freq+symbols): {:.1f}s'.format(_t.time()-_tg))
    print('SNAF T-antigen results written to {}'.format(outdir))
    return outdir


# --------------------------------------------------------------------------- snaf-precompute-control
def run_snaf_precompute_control(args):
    """Precompute the per-junction control-stats table (mean/std/mle/normal_prevalence/
    BayesTS sigma+percentile) from a control h5ad. One-time per control; the result is
    reused at runtime as a small lookup instead of loading the full control + BayesTS."""
    from altanalyze3.components.snaf import control_stats
    _require_file(args.control_h5ad, '--control_h5ad')
    out = str(args.stats_output) if getattr(args, 'stats_output', None) else None
    bayes_uids = None
    if getattr(args, 'bayes_juncounts', None):
        _require_file(args.bayes_juncounts, '--bayes_juncounts')
        import pandas as _pd
        # SNAF cohort UIDs carry a '=chr...' coordinate suffix; the control h5ad keys are the
        # AltAnalyze part before '=' (same transform gtex_configuration matches on).
        _idx = _pd.read_csv(str(args.bayes_juncounts), sep='\t', index_col=0, usecols=[0]).index.astype(str)
        bayes_uids = sorted({u.split('=')[0] for u in _idx})
        print('BayesTS restricted to {} cohort junctions (from {} rows in {})'.format(len(bayes_uids), len(_idx), args.bayes_juncounts))
    path = control_stats.precompute_control_stats(
        str(args.control_h5ad), out_path=out, cutoff=args.normal_cutoff,
        bayes_mode=args.bayes_mode, bayes_epoch=args.bayes_epoch,
        bayes_batch=args.bayes_batch, compute_bayes=(not args.no_bayes), bayes_uids=bayes_uids,
        bayes_cores=getattr(args, 'bayes_cores', None))
    print('SNAF control-stats table written to {}'.format(path))
    return path


def _derive_surface_frequency_table(jcmq, uids, outdir):
    """Write the frequency table SNAF-B needs WITHOUT running the T-antigen pipeline.

    SNAF-B's `report_candidates` reads exactly three things out of the SNAF-T table: the
    uid (from the 'uid,uid' index), `tumor_specificity_mean` and `tumor_specificity_mle`.
    None of them come from MHC binding -- they are the sifting condition matrix and the
    control-side tumor specificity, both already computed by the time SNAF-B has its
    membrane tuples. So this builds the same table over the membrane neojunctions using
    the SAME functions the T pipeline uses (`stage0_compatible_results`'s per-junction
    expressed-sample logic, then `enhance_frequency_table`), rather than approximating it.

    SCOPE: the rows cover the membrane neojunctions only (SNAF-T's stage-0 table covers
    every neojunction). That is complete for SNAF-B, whose candidates are a subset of the
    membrane set, and is stated in the file name.
    """
    import pandas as pd
    from altanalyze3.components.snaf.snaf import enhance_frequency_table

    cond = jcmq.cond_df
    keep = [u for u in uids if u in cond.index]
    missing = len(uids) - len(keep)
    if missing:
        logger.warning('%d / %d membrane neojunctions are absent from cond_df and are not in '
                       'the derived frequency table', missing, len(uids))
    sub = cond.loc[keep]
    records = []
    for index, series in sub.iterrows():
        samples = series.loc[series].index.tolist()
        records.append((index, samples, len(samples)))
    df = pd.DataFrame.from_records(records, columns=['junction', 'samples', 'n_sample'])
    df = df.set_index(keys='junction').sort_values(by='n_sample', ascending=False)
    raw = os.path.join(outdir, 'frequency_surface_stage0.txt')
    df.to_csv(raw, sep='\t')
    df = pd.read_csv(raw, sep='\t', index_col=0)
    df.index = [','.join([item, item]) for item in df.index]
    uid_path = os.path.join(outdir, 'frequency_surface_stage0_verbosity1_uid.txt')
    df.to_csv(uid_path, sep='\t')
    df = pd.read_csv(uid_path, sep='\t', index_col=0)
    name = 'frequency_surface_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt'
    enhance_frequency_table(df, True, True, outdir, name)
    out = os.path.join(outdir, name)
    logger.info('derived SNAF-B frequency table (%d membrane neojunctions) -> %s', len(keep), out)
    return out


# --------------------------------------------------------------------------- snaf-build-surface-db
def run_snaf_build_surface_db(args):
    """Format a user surfaceome gene list into a SNAF-B surface database."""
    from altanalyze3.components.snaf.surface.surface_db import build_surface_db
    import json
    _require_file(args.gene_table, '--gene_table')
    from altanalyze3.components.snaf.reference import ensure_reference
    root = ensure_reference(str(args.db_dir), gtex_mode='count',
                            download=getattr(args, 'download_ref', False))
    uniprot_dir = str(args.uniprot_dir) if getattr(args, 'uniprot_dir', None) else None
    if uniprot_dir and not os.path.isdir(uniprot_dir):
        raise FileNotFoundError('--uniprot_dir not found: {}'.format(uniprot_dir))
    outdir, params = build_surface_db(
        gene_table=str(args.gene_table), db_dir=root, outdir=str(args.output),
        uniprot_dir=uniprot_dir, mode=args.mode, name=getattr(args, 'name', None))
    print(json.dumps(params['counts'], indent=2))
    print('SNAF-B surface database written to {}'.format(os.path.abspath(outdir)))
    return outdir


# --------------------------------------------------------------------------- snaf-b (surface / B-antigen)
def run_snaf_b(args):
    """SNAF surface / B-antigen pipeline entry.

    Pure-Python and offline-first: transmembrane topology via tmhmm.py (no TMHMM binary),
    pairwise alignment via Biopython (no EMBOSS-needle REST), UTR retrieval via the local
    genome FASTA when given. Any missing optional dependency (tmhmm.py, mygene) is noted
    and the pipeline proceeds with the remaining steps.

    Needs --db_dir (Alt91_db + controls). --freq_path (the frequency table from a prior
    `altanalyze3 snaf` run) is optional: when omitted, the equivalent table is derived here
    over the membrane neojunctions, so SNAF-B runs WITHOUT SNAF-T and without HLA types.
    --validation_gtf (e.g. a long-read SQANTI GTF, plain or gzipped) enables the
    stringency-4/5 EST/long-read support gates; without it, stringency-3 candidates are
    still produced fully offline. --surface_db replaces the built-in Alt91_db surfaceome
    whitelist with a user database (e.g. SURFY).
    """
    _ensure_fork_safety()
    from altanalyze3.components import snaf
    from altanalyze3.components.snaf import surface
    from altanalyze3.components.snaf.reference import ensure_reference

    outdir = str(args.output)
    os.makedirs(outdir, exist_ok=True)
    _require_file(args.juncounts, '--juncounts')
    freq_path = _require_file(args.freq_path, '--freq_path') \
        if getattr(args, 'freq_path', None) else None
    surface_db = str(args.surface_db) if getattr(args, 'surface_db', None) else None
    if surface_db is not None and not os.path.exists(surface_db):
        raise FileNotFoundError('--surface_db not found: {}'.format(surface_db))
    validation_gtf = str(args.validation_gtf) if getattr(args, 'validation_gtf', None) else None
    if validation_gtf:
        _require_file(validation_gtf, '--validation_gtf')
    add_control = _load_add_control(getattr(args, 'add_control', None))
    genome_fasta = str(args.genome_fasta) if getattr(args, 'genome_fasta', None) else None
    gtex_db = str(args.gtex_db) if getattr(args, 'gtex_db', None) else None
    mode = args.mode

    if genome_fasta is None and os.environ.get('SNAF_OFFLINE', '0') != '1':
        logger.warning(
            'No --genome_fasta and SNAF_OFFLINE != 1: UTR-event sequence retrieval may query '
            'the UCSC DAS web service. Pass --genome_fasta <hg38.fa> or set SNAF_OFFLINE=1 for '
            'a fully offline run.')

    df = _read_junction_matrix(args.juncounts)
    # Resolve/validate the reference (tolerates the tarball data/ nesting).
    root = ensure_reference(str(args.db_dir), gtex_mode=args.gtex_mode,
                            download=getattr(args, 'download_ref', False))
    # Configure gtex + Alt91_db (needed by get_membrane_tuples) and the surface globals.
    snaf.initialize(
        df=df, db_dir=str(args.db_dir), gtex_mode=args.gtex_mode, binding_method=None,
        t_min=args.t_min, n_max=args.n_max, normal_cutoff=args.normal_cutoff,
        tumor_cutoff=args.tumor_cutoff, normal_prevalance_cutoff=args.normal_prevalance_cutoff,
        tumor_prevalance_cutoff=args.tumor_prevalance_cutoff, add_control=add_control,
        genome_fasta=genome_fasta, gtex_db=gtex_db, download_ref=getattr(args, 'download_ref', False),
        control_stats_path=(str(args.control_stats) if getattr(args, 'control_stats', None) else None))
    surface.initialize(db_dir=root, surface_db=surface_db)

    membrane_tuples, jcmq = snaf.JunctionCountMatrixQuery.get_membrane_tuples(
        df, return_jcmq=True, add_control=add_control, not_in_db=args.not_in_db, outdir=outdir,
        filter_mode=args.filter_mode, cores=args.cpus,
        min_samples=getattr(args, 'min_samples', 1),
        max_bayests_percentile=getattr(args, 'max_bayests_percentile', None))
    logger.info("neojunctions passing sifting: %d | on a surface gene: %d",
                len(jcmq.valid), len(membrane_tuples))
    if not membrane_tuples:
        raise RuntimeError(
            'no neojunction fell on a surface gene: {} junctions passed sifting but none is in '
            'the surface database ({}). Check --surface_db and the junction UID gene ids.'.format(
                len(jcmq.valid), surface_db or 'built-in Alt91_db surfaceome'))

    if freq_path is None:
        logger.warning('no --freq_path given: deriving the frequency table from this run '
                       '(SNAF-T is not required for SNAF-B)')
        freq_path = _derive_surface_frequency_table(
            jcmq, [t[0] for t in membrane_tuples], outdir)

    # short_read/find_full_length: validation_gtf feeds the str4/5 support gate in
    # generate_full_results. long_read: the gtf is the primary transcript source in run().
    run_gtf = validation_gtf if mode == 'long_read' else None
    surface.run(uids=membrane_tuples, outdir=outdir, prediction_mode=mode, n_stride=args.n_stride,
                gtf=run_gtf, tmhmm=(not args.no_tmhmm),
                software_path=(str(args.tmhmm_path) if getattr(args, 'tmhmm_path', None) else None))
    surface.generate_full_results(outdir=outdir, freq_path=freq_path, mode=mode,
                                  validation_gtf=validation_gtf)
    print('SNAF-B surface results written to {}'.format(os.path.join(outdir, 'B_candidates')))
    return outdir
