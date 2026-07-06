"""Self-contained regression tests for the modernized SNAF workflow.

Everything here runs on SYNTHETIC data (generated in a tmp dir), so the suite is
reproducible on any machine / CI with no reference DB and no network. The same
checks were additionally validated against a real 6000-junction GTEx slice during
development.

Covers:
  * vectorized/sparse sifting == legacy dense sifting  (maxmin & prevalance, +controls, +exonlist)
  * closed-form half-normal MLE == scipy minimize_scalar
  * sparse add_control combine == legacy dense pandas join
  * BayesTS tumor-specificity: runs on CPU, sigma in [0,1], monotone with normal breadth  (skipped if pyro missing)
  * netMHCpan stdout parser (pure-python, portable)
  * _parallel_apply fork/serial correctness
  * offline genome-FASTA sequence retrieval
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy.sparse import csr_matrix

from altanalyze3.components.snaf import gtex as G
from altanalyze3.components.snaf import binding as BD
from altanalyze3.components.snaf import snaf as SN

try:
    import pyro  # noqa
    HAS_PYRO = True
except Exception:
    HAS_PYRO = False


# --------------------------------------------------------------------------- fixtures
def _make_synthetic_gtex(path, n_junctions=1200, n_samples=400, n_tissues=12, seed=0):
    rng = np.random.default_rng(seed)
    uids = ['ENSG{:08d}:E{}.1-E{}.1'.format(i, i % 25 + 1, i % 25 + 2) for i in range(n_junctions)]
    # sparse-ish normal expression, some junctions broad, some rare
    dense = (rng.random((n_junctions, n_samples)) < rng.random((n_junctions, 1)) * 0.3) * \
            rng.integers(1, 8, size=(n_junctions, n_samples))
    X = csr_matrix(dense.astype(np.float32))
    tissues = rng.choice(['tis{}'.format(k) for k in range(n_tissues)], n_samples)
    var = pd.DataFrame({'tissue': tissues, 'total_count': rng.uniform(5, 20, n_samples)},
                       index=['g{}'.format(i) for i in range(n_samples)])
    A = ad.AnnData(X=X, obs=pd.DataFrame(index=uids), var=var)
    A.obs['mean'] = np.asarray(X.mean(axis=1)).ravel()
    A.write(path)
    return uids


@pytest.fixture(scope='module')
def synthetic(tmp_path_factory):
    d = tmp_path_factory.mktemp('snaf_synth')
    gtex_path = str(d / 'gtex.h5ad')
    uids = _make_synthetic_gtex(gtex_path)
    rng = np.random.default_rng(1)
    chosen = list(rng.choice(uids, 800, replace=False))
    novel = ['ENSGnovel{:04d}:E1.1-E2.1'.format(i) for i in range(150)]
    tumor_uids = chosen + novel
    base = rng.poisson(2.0, size=(len(tumor_uids), 25)).astype(float)
    spike = (rng.random((len(tumor_uids), 25)) < 0.15) * rng.integers(25, 300, size=(len(tumor_uids), 25))
    jcm = pd.DataFrame(base + spike, index=tumor_uids, columns=['t{}'.format(i) for i in range(25)])
    return {'dir': str(d), 'gtex': gtex_path, 'uids': uids, 'jcm': jcm, 'tumor_uids': tumor_uids}


def _configure(synthetic):
    G.gtex_configuration(synthetic['jcm'], synthetic['gtex'], 20, 3, 5, 20, 0.01, 0.1, add_control=None)


def _complete_exonlist(tumor_uids):
    # Add an entry for every junction that reaches the dict lookup in the sifting
    # (i.e. the branch NOT short-circuited by '_'/'U'/'ENSG'/'I' in the exon field),
    # so the legacy code -- which does dict_exonlist[ensg] with no .get -- never
    # KeyErrors. Half are made "documented" (pair present), half not.
    d = {}
    k = 0
    for u in tumor_uids:
        exons = ':'.join(u.split(':')[1:])
        if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
            continue
        ensg = u.split(':')[0]; e1, e2 = exons.split('-')
        d.setdefault(ensg, []).append('E0.1|{}|{}|E9.1'.format(e1, e2) if k % 2 == 0 else 'E0.1|E9.1')
        k += 1
    return d


# --------------------------------------------------------------------------- sifting
@pytest.mark.parametrize('mode', ['maxmin', 'prevalance'])
def test_sifting_equivalence(synthetic, mode):
    jcm = synthetic['jcm']; out = synthetic['dir']
    legacy = getattr(G, 'multiple_crude_sifting_{}'.format(mode))
    fast = getattr(G, 'multiple_crude_sifting_{}_fast'.format(mode))
    exonlist = _complete_exonlist(synthetic['tumor_uids'])
    ctrl_df = pd.DataFrame(np.random.default_rng(2).poisson(1, size=(len(jcm), 10)).astype(float),
                           index=jcm.index, columns=['c{}'.format(i) for i in range(10)])
    for add_control, dexon in [(None, None), (None, exonlist), ({'c': ctrl_df}, None), ({'c': ctrl_df}, exonlist)]:
        _configure(synthetic)
        vL, iL, cL = legacy(jcm, add_control, dexon, out)
        _configure(synthetic)
        vF, iF, cF = fast(jcm, add_control, dexon, out)
        assert set(vL) == set(vF)
        assert set(iL) == set(iF)
        assert np.array_equal(cL.values.astype(bool), cF.loc[cL.index, cL.columns].values.astype(bool))


def test_add_control_sparse_combine(synthetic):
    jcm = synthetic['jcm']
    rng = np.random.default_rng(4)
    ctrl = pd.DataFrame(rng.poisson(1, size=(500, 8)).astype(float),
                        index=list(rng.choice(jcm.index, 500, replace=False)),
                        columns=['x{}'.format(i) for i in range(8)])
    adata = G.gtex_configuration(jcm, synthetic['gtex'], 20, 3, 5, 20, 0.01, 0.1, add_control={'x': ctrl})
    gt = ad.read_h5ad(synthetic['gtex'])
    gt_sub = gt[[o for o in gt.obs_names if o in set(jcm.index)], :]
    ref = gt_sub.to_df().join(ctrl.loc[list(set(ctrl.index) & set(jcm.index))], how='outer').fillna(0)
    new = adata.to_df().loc[ref.index, ref.columns]
    assert np.allclose(new.values, ref.values, atol=1e-4)
    assert np.allclose(adata.obs['mean'].reindex(ref.index).values, ref.mean(axis=1).values, atol=1e-4)
    assert np.allclose(adata.var['total_count'].reindex(ref.columns).values, (ref.sum(axis=0) / 1e6).values, atol=1e-4)


def test_mle_closed_form_matches_optimizer():
    from scipy.optimize import minimize_scalar
    rng = np.random.default_rng(7)
    for _ in range(200):
        n = rng.integers(5, 200)
        y = np.where(rng.random(n) < 0.2, rng.gamma(1.0, 0.3, n), 0.0)
        legacy = minimize_scalar(G.mle_func, bounds=(0, 1), args=(y,), method='bounded').x
        closed = min(np.sqrt(np.mean(y ** 2)), 1.0)
        assert abs(legacy - closed) < 1e-3


# --------------------------------------------------------------------------- BayesTS
@pytest.mark.skipif(not HAS_PYRO, reason='pyro not installed')
def test_bayests_sigma(synthetic):
    from altanalyze3.components.snaf import bayests as B
    A = ad.read_h5ad(synthetic['gtex'])
    breadth = np.asarray((A.X > 0).sum(axis=1)).ravel()
    order = np.argsort(breadth)
    idx = np.concatenate([order[:60], order[-60:]])
    uids = [A.obs_names[i] for i in idx]
    res = B.compute_bayests_sigma(A, uids=uids, mode='XY', epoch=400, seed=0)
    sig = res['mean_sigma']
    assert (sig >= 0).all() and (sig <= 1).all()
    lo = res.loc[[A.obs_names[i] for i in order[:60] if A.obs_names[i] in res.index], 'mean_sigma'].mean()
    hi = res.loc[[A.obs_names[i] for i in order[-60:] if A.obs_names[i] in res.index], 'mean_sigma'].mean()
    assert hi > lo   # broadly-expressed -> higher sigma (less tumor-specific)


# --------------------------------------------------------------------------- portability
def test_netmhcpan_parser():
    out = (" Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL BindLevel\n"
           "   1  HLA-A*02:01  AAAWYLWEV AAAWYLWEV 0 0 0 0 0 AAAWYLWEV PEP1 0.850 0.120 <= SB\n"
           "   1  HLA-A*01:01  AAGLQDCTM AAGLQDCTM 0 0 0 0 0 AAGLQDCTM PEP2 0.004 8.500\n")
    df = BD._parse_netMHCpan_stdout(out, {'AAAWYLWEV', 'AAGLQDCTM'})
    assert set(df['peptide']) == {'AAAWYLWEV', 'AAGLQDCTM'}       # binder + non-binder both kept
    r = df.set_index('peptide')
    assert abs(r.loc['AAAWYLWEV', 'score'] - 0.120) < 1e-9 and r.loc['AAAWYLWEV', 'identity'] == 'SB'
    assert r.loc['AAGLQDCTM', 'identity'] == '' and r.loc['AAAWYLWEV', 'hla'] == 'HLA-A*02:01'


def _double(x):   # module-level so it is picklable for the fork pool
    return x * 2


def test_parallel_apply_fork_and_serial():
    arg_tuples = [(i,) for i in range(6)]
    assert SN._parallel_apply(_double, arg_tuples, cores=1) == [0, 2, 4, 6, 8, 10]     # serial
    assert SN._parallel_apply(_double, arg_tuples, cores=2) == [0, 2, 4, 6, 8, 10]     # fork or serial


def test_offline_genome_fasta(tmp_path):
    import pysam
    fa = str(tmp_path / 'mini.fa')
    with open(fa, 'w') as f:
        f.write('>chr1\n' + 'ACGTACGTAC' * 20 + '\n>2\n' + 'TTTTGGGGCC' * 20 + '\n')
    pysam.faidx(fa)
    SN.set_genome_fasta(fa)
    assert SN.retrieveSeqFromUCSCapi('chr1', 1, 10) == 'ACGTACGTAC'
    assert SN.retrieveSeqFromUCSCapi('chr2', 1, 10) == 'TTTTGGGGCC'   # name toggle chr2 -> 2
    SN.set_genome_fasta(None)
    os.environ['SNAF_OFFLINE'] = '1'
    try:
        assert SN.retrieveSeqFromUCSCapi('chrX', 1, 10) == '#' * 10   # no fasta + offline -> placeholder
    finally:
        os.environ.pop('SNAF_OFFLINE', None)


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-v']))
