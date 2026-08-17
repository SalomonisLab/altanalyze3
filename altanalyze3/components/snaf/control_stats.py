"""Precomputed control-reference statistics for SNAF.

BayesTS and the GTEx control summaries depend ONLY on the control reference, not on the
tumor cohort. So they should be computed ONCE per control (shipped for the default
Ensembl91 GTEx control, or computed on first use of a new control and cached) and then
reused at runtime as a small lookup table -- instead of loading the multi-GB control
h5ad and re-running BayesTS on every analysis.

The table has one row per control junction with columns:
    mean                 : mean count across control samples  (== adata.obs['mean'])
    std                  : std of count across control samples
    mle                  : closed-form half-normal tumor-specificity MLE (control-side)
    normal_prevalence    : fraction of control samples with count > cutoff
    bayests_sigma        : BayesTS per-junction sigma in [0,1]; LOWER == more tumor-specific
    bayests_percentile   : global rank of bayests_sigma in (0,1]; the SNAF "p-value"

At runtime SNAF filters junctions on a default `bayests_percentile` (or sigma) threshold
(user-tunable); junctions ABSENT from the control are maximally tumor-specific and always
kept. This eliminates the large fraction of junctions expressed in normal tissue up front,
before the expensive translation/binding steps.
"""
import os
import logging
import numpy as np
import pandas as pd
from scipy.sparse import issparse

logger = logging.getLogger(__name__)

STATS_COLUMNS = ['mean', 'std', 'mle', 'normal_prevalence', 'bayests_sigma', 'bayests_percentile']


def default_stats_path(control_h5ad):
    """Canonical location of the precomputed table for a given control h5ad."""
    base = str(control_h5ad)
    for ext in ('.h5ad', '.h5'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    return base + '.snaf_stats.tsv.gz'


# --------------------------------------------------------------------------- bundled BayesTS
# BayesTS is expensive (hours of pyro SVI) and its answer depends only on the REFERENCE
# cohort, never on the tumour cohort being analysed. So the finished per-junction scores for
# the two standard normal references ship with SNAF and are loaded directly -- SNAF must
# never recompute them. Built by BayesTS/build_bundled_bayests.py; see BayesTS/README.md.
BAYESTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BayesTS')

BUNDLED_BAYESTS = {
    'gtex':          'GTEx_BayesTS.tsv.gz',            # 2,476,734 junctions, 2,629 samples / 51 tissues
    'tabulasapiens': 'TabulaSapiens_BayesTS.tsv.gz',   # 2,852,713 junctions, 762 pseudobulks / 112 cell types
}

# Both references are used by default. A junction must look tumour-specific in EVERY
# reference to survive, so the scores are combined by taking the MAXIMUM percentile: broad
# expression in either normal cohort disqualifies it. That is the conservative direction --
# it can only reject candidates, never invent them.
DEFAULT_BAYESTS_REFERENCES = ('gtex', 'tabulasapiens')


def bundled_bayests_path(name):
    """Absolute path of a bundled BayesTS table, or None when it is not installed."""
    fn = BUNDLED_BAYESTS.get(str(name).lower())
    if fn is None:
        raise KeyError('unknown BayesTS reference {!r}; known: {}'.format(
            name, sorted(BUNDLED_BAYESTS)))
    p = os.path.join(BAYESTS_DIR, fn)
    return p if os.path.exists(p) else None


def load_bundled_bayests_per_reference(references=DEFAULT_BAYESTS_REFERENCES, uids=None):
    """{reference_name: DataFrame(bayests_sigma, bayests_percentile)} -- kept SEPARATE.

    Use this, not load_bundled_bayests, whenever the answer will be interpreted. GTEx and
    Tabula Sapiens measure different normal populations: a junction silent in GTEx bulk
    tissue but expressed across Tabula Sapiens cell types is a different candidate from the
    reverse, and collapsing them to one number destroys exactly that distinction.
    """
    out = OrderedDict()
    want = set(uids) if uids is not None else None
    for name in references:
        try:
            path = bundled_bayests_path(name)
        except KeyError:
            logger.warning('unknown BayesTS reference %r; skipped', name)
            continue
        if path is None:
            logger.warning('bundled BayesTS reference %r is not installed under %s; skipped',
                           name, BAYESTS_DIR)
            continue
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df.loc[~df.index.duplicated()]
        if want is not None:
            df = df.loc[df.index.intersection(want)]
        logger.warning('BayesTS reference %s: %d junctions loaded from %s (no recomputation)',
                       name, df.shape[0], os.path.basename(path))
        out[name] = df
    return out


def load_bundled_bayests(references=DEFAULT_BAYESTS_REFERENCES, uids=None):
    """Per-junction BayesTS sigma/percentile from the bundled references.

    :param references: which bundled references to combine.
    :param uids: restrict to these junction UIDs (the cohort's tested set) so only the
        intersecting rows are held -- the tables are ~2.5-2.9M rows each.
    NOTE: this COMBINES references into one number (max percentile) and is only a
    convenience for a single go/no-go filter. For anything that will be interpreted, use
    load_bundled_bayests_per_reference so each normal cohort stays visible.

    :return: DataFrame indexed by UID with 'bayests_sigma' and 'bayests_percentile',
        or None when no bundled table is installed. A junction present in only some
        references is scored on those; a junction absent from all is simply not in the
        result, and callers must treat "absent" as maximally tumour-specific (never as 0).
    """
    frames = []
    used = []
    want = set(uids) if uids is not None else None
    for name in references:
        try:
            path = bundled_bayests_path(name)
        except KeyError:
            logger.warning('unknown BayesTS reference %r; skipped', name)
            continue
        if path is None:
            logger.warning('bundled BayesTS reference %r is not installed under %s; skipped',
                           name, BAYESTS_DIR)
            continue
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df.loc[~df.index.duplicated()]
        if want is not None:
            df = df.loc[df.index.intersection(want)]
        logger.warning('BayesTS reference %s: %d junctions loaded from %s (no recomputation)',
                       name, df.shape[0], os.path.basename(path))
        frames.append(df)
        used.append(name)
    if not frames:
        return None
    if len(frames) == 1:
        out = frames[0]
    else:
        # maximum percentile across references = least tumour-specific verdict wins
        out = pd.concat(frames, axis=1, keys=used)
        pct = out.xs('bayests_percentile', axis=1, level=1).max(axis=1, skipna=True)
        sig = out.xs('bayests_sigma', axis=1, level=1).max(axis=1, skipna=True)
        out = pd.DataFrame({'bayests_sigma': sig, 'bayests_percentile': pct})
    logger.warning('BayesTS combined over %s: %d junctions (percentile = max across '
                   'references, so broad expression in ANY normal cohort disqualifies)',
                   '+'.join(used), out.shape[0])
    return out


def load_control_stats(path):
    """Load a stats table (indexed by junction UID); returns None if absent."""
    if path is None or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep='\t', index_col=0)
    return df


def _basic_stats(adata, cutoff):
    """mean, std, mle, normal_prevalence per junction from an in-memory control AnnData."""
    X = adata.X
    n = adata.shape[1]
    # Prefer the stored obs['mean']/['std'] so the summary is BIT-IDENTICAL to the full-matrix
    # runtime (gtex_configuration/sifting read the stored obs['mean']); recomputing in float64
    # differs from the float32-stored value by ~1e-5 and can flip boundary junctions.
    if 'mean' in adata.obs.columns:
        mean = np.asarray(adata.obs['mean'].values, dtype=np.float64)
    else:
        mean = np.asarray(X.mean(axis=1)).ravel()
    if 'std' in adata.obs.columns:
        std = np.asarray(adata.obs['std'].values, dtype=np.float64)
    else:
        if issparse(X):
            mean_sq = np.asarray(X.multiply(X).mean(axis=1)).ravel()
        else:
            mean_sq = (np.asarray(X) ** 2).mean(axis=1)
        std = np.sqrt(np.maximum(mean_sq - mean ** 2, 0.0))
    # per-sample library size for CPM: use the stored var['total_count'] when present so the
    # MLE matches tumor_specificity_batch exactly (it reads the same stored value); only fall
    # back to recomputing the column sum if the control lacks it.
    if 'total_count' in adata.var.columns:
        total_count = np.asarray(adata.var['total_count'].values, dtype=np.float64)
    else:
        total_count = np.asarray(X.sum(axis=0)).ravel() / 1e6
    inv = np.divide(1.0, total_count, out=np.zeros_like(total_count), where=total_count > 0)
    if issparse(X):
        Xs = X.multiply(inv[None, :])
        mle = np.sqrt(np.asarray(Xs.multiply(Xs).sum(axis=1)).ravel() / n)
        prev = np.asarray((X > cutoff).sum(axis=1)).ravel() / n
    else:
        cpm = np.asarray(X) * inv[None, :]
        mle = np.sqrt((cpm ** 2).mean(axis=1))
        prev = (np.asarray(X) > cutoff).sum(axis=1) / n
    mle = np.minimum(mle, 1.0)
    return pd.DataFrame({'mean': mean, 'std': std, 'mle': mle, 'normal_prevalence': prev},
                        index=adata.obs_names)


_BAYES_ADATA = None   # set before the fork; workers inherit it copy-on-write (never pickled)


def _bayes_batch_worker(payload):
    '''Fork-worker: BayesTS sigma for ONE batch of junctions. Independent + seeded, so the
    result is identical to the serial path. torch is imported HERE (not in the parent) so there
    is no torch-after-fork hazard; threads are capped to avoid N_workers x cores oversubscription.'''
    import os as _os
    batch, mode, epoch, threads = payload
    try:
        import torch; torch.set_num_threads(max(1, threads))
    except Exception:
        pass
    for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        _os.environ[_v] = str(max(1, threads))
    from . import bayests
    sig = bayests.compute_bayests_sigma(_BAYES_ADATA[batch, :].to_memory(), uids=batch,
                                        mode=mode, epoch=epoch)
    return sig['mean_sigma']


def precompute_control_stats(control_h5ad, out_path=None, cutoff=5, bayes_mode='XY',
                             bayes_epoch=2000, bayes_batch=50000, compute_bayes=True,
                             bayes_uids=None, bayes_cores=None):
    """Read a control h5ad once and write the per-junction stats table. One-time job.

    bayes_uids: if given, compute BayesTS ONLY for these junctions (the cohort's tested set,
    intersected with the control) instead of all ~2.5M -- ~50x less work. The percentile then
    ranks tumor-specificity among the tested cohort rather than globally; junctions outside the
    set get NaN BayesTS (they are never queried by that cohort's run)."""
    import anndata as ad
    from . import gtex as _gtex
    out_path = out_path or default_stats_path(control_h5ad)
    logger.warning('Precomputing SNAF control stats from %s -> %s (one-time)', control_h5ad, out_path)
    adata = ad.read_h5ad(control_h5ad)
    # ensure the fields BayesTS/basic stats rely on (mirrors gtex_configuration)
    if 'mean' not in adata.obs:
        adata.obs['mean'] = np.asarray(adata.X.mean(axis=1)).ravel()
    if 'total_count' not in adata.var:
        adata.var['total_count'] = np.asarray(adata.X.sum(axis=0)).ravel() / 1e6
    if 'tissue' not in adata.var:
        adata.var['tissue'] = 'unknown'

    stats = _basic_stats(adata, cutoff)

    if compute_bayes:
        from . import bayests
        uids = adata.obs_names.tolist()
        if bayes_uids is not None:
            wanted = set(bayes_uids)
            uids = [u for u in uids if u in wanted]   # preserve h5ad order; intersect with control
            logger.warning('BayesTS restricted to cohort: %d of %d control junctions', len(uids), adata.n_obs)
        sigmas = pd.Series(index=uids, dtype=float)
        # Each 50k-junction batch is an independent, seeded pyro SVI -> the batches parallelize
        # across cores with a result identical to running them serially. This is the difference
        # between ~5.5 h (serial) and ~1 h on a 10-core box for the full ~2.5M-junction GTEx.
        batches = [uids[s0:s0 + bayes_batch] for s0 in range(0, len(uids), bayes_batch)]
        import multiprocessing as _mp, os as _os
        from .snaf import _mp_context
        n_cpu = _os.cpu_count() or 2
        n_workers = bayes_cores if bayes_cores else max(1, n_cpu - 2)
        n_workers = max(1, min(n_workers, len(batches)))
        ctx = _mp_context() if n_workers > 1 else None   # None on Windows / macOS-no-forkflag -> serial
        if ctx is None or n_workers == 1:
            for i, batch in enumerate(batches):
                sig = bayests.compute_bayests_sigma(adata[batch, :].to_memory(), uids=batch,
                                                    mode=bayes_mode, epoch=bayes_epoch)
                sigmas.loc[sig.index] = sig['mean_sigma'].values
                logger.warning('  BayesTS %d/%d junctions (serial)', min((i + 1) * bayes_batch, len(uids)), len(uids))
        else:
            threads = max(1, n_cpu // n_workers)   # total threads ~= cores, no oversubscription
            global _BAYES_ADATA
            _BAYES_ADATA = adata                    # fork children inherit this (copy-on-write)
            payloads = [(b, bayes_mode, bayes_epoch, threads) for b in batches]
            logger.warning('BayesTS: %d batches across %d workers (%d threads each)', len(batches), n_workers, threads)
            try:
                pool = ctx.Pool(processes=n_workers)
                try:
                    for i, part in enumerate(pool.imap(_bayes_batch_worker, payloads)):
                        sigmas.loc[part.index] = part.values
                        logger.warning('  BayesTS batch %d/%d done', i + 1, len(batches))
                    pool.close(); pool.join()
                finally:
                    pool.terminate()
            finally:
                _BAYES_ADATA = None
        stats['bayests_sigma'] = sigmas.reindex(stats.index).values
        # global percentile (rank) of sigma across ALL junctions -- the SNAF "p-value"
        stats['bayests_percentile'] = stats['bayests_sigma'].rank(method='average', pct=True).values
    else:
        stats['bayests_sigma'] = np.nan
        stats['bayests_percentile'] = np.nan

    stats.to_csv(out_path, sep='\t')
    logger.warning('Wrote %d-junction control stats table to %s', stats.shape[0], out_path)
    return out_path
