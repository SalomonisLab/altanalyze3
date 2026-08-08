"""Can the SCCAF logistic regression be parallelised without changing its answer?

SCCAF is the one stage I could not speed up. It fits
`LogisticRegression(penalty='l1', solver='liblinear', max_iter=100000)` and costs
1.3 s to 3.0 s per annotation on the demo, 4.2 s to 8.6 s on the 8k x 20k
simulation. liblinear is single-threaded C.

sklearn hands the whole multiclass problem to liblinear, which runs one-versus-rest
INSIDE the C code, one class after another. If the same coefficients came out of C
independent binary fits, those fits could run on C cores at once.

The catch is the random stream. sklearn draws ONE seed and liblinear's coordinate
descent consumes `rand()` across all sub-problems in sequence, so sub-problem k
sees a stream that depends on sub-problems 0..k-1. Fitting them separately restarts
that stream. This script measures whether the difference matters.

What has to hold for the parallel version to ship:
  1. the per-cluster accuracies (the metric the Shapley game reads) are identical;
  2. the confusion matrix is identical;
  3. and it is actually faster.

Usage:
  python3.11 sccaf_parallel_probe.py [--h5ad data/demo_pbmc3k.h5ad]
"""

import os
import sys
import time
import argparse

import numpy as np
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from sklearn.preprocessing import LabelEncoder                      # noqa: E402
from sklearn.model_selection import StratifiedShuffleSplit          # noqa: E402
from sklearn.linear_model import LogisticRegression                 # noqa: E402
from sklearn.metrics import confusion_matrix                        # noqa: E402
from joblib import Parallel, delayed                                # noqa: E402

from altanalyze3.components.sctriangulate import ScTriangulate      # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M       # noqa: E402

SEED = M.SCCAF_RANDOM_STATE


def prepare(adata, key):
    artifact = M._artifact_gene_set(adata.uns.get('_species', 'human'), 2)
    keep = ~adata.var_names.isin(artifact)
    X = M._as_dense(adata[:, keep].X)
    Y = LabelEncoder().fit_transform(adata.obs[key].values)
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=0)
    train, test = next(iter(sss.split(X, Y)))
    return X[train], Y[train], X[test], Y[test]


def fit_multiclass(X_train, Y_train):
    return LogisticRegression(penalty='l1', solver='liblinear',
                              max_iter=100000, random_state=SEED).fit(X_train, Y_train)


def fit_per_class(X_train, Y_train, n_jobs):
    """One binary liblinear fit per class, in parallel, assembled like sklearn does."""
    classes = np.unique(Y_train)

    def one(c):
        m = LogisticRegression(penalty='l1', solver='liblinear',
                               max_iter=100000, random_state=SEED)
        m.fit(X_train, (Y_train == c).astype(int))
        return m.coef_[0], m.intercept_[0]

    out = Parallel(n_jobs=n_jobs, prefer='processes')(delayed(one)(c) for c in classes)
    coef = np.vstack([o[0] for o in out])
    intercept = np.array([o[1] for o in out])
    return classes, coef, intercept


def accuracies(Y_true, Y_pred, n_classes):
    m = confusion_matrix(Y_true, Y_pred, labels=np.arange(n_classes))
    with np.errstate(invalid='ignore', divide='ignore'):
        acc = m.diagonal() / m.sum(axis=1)
    return m, acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', default=os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad'))
    ap.add_argument('--keys', default='sctri_rna_leiden_1,sctri_rna_leiden_2,sctri_rna_leiden_3')
    ap.add_argument('--n-jobs', type=int, default=8)
    args = ap.parse_args()

    keys = args.keys.split(',')
    adata = sc.read(args.h5ad)
    sctri = ScTriangulate(dir=os.path.join(HERE, 'out_tests'), adata=adata,
                          query=keys, predict_doublet=False, verbose=0)
    sctri._to_dense()

    print('{:<26s} {:>7s} {:>9s} {:>9s} {:>9s} {:>10s} {:>10s}'.format(
        'annotation', 'classes', 'sklearn', 'parallel', 'speedup',
        'same acc', 'same conf'))
    print('-' * 90)

    all_same = True
    for key in keys:
        adata_c = M.check_filter_single_cluster(sctri.adata, key)
        X_train, Y_train, X_test, Y_test = prepare(adata_c, key)
        n_classes = len(np.unique(Y_train))

        t0 = time.perf_counter()
        model = fit_multiclass(X_train, Y_train)
        t_sklearn = time.perf_counter() - t0
        m_ref, acc_ref = accuracies(Y_test, model.predict(X_test), n_classes)

        t0 = time.perf_counter()
        classes, coef, intercept = fit_per_class(X_train, Y_train, args.n_jobs)
        t_par = time.perf_counter() - t0
        scores = X_test @ coef.T + intercept
        pred = classes[scores.argmax(axis=1)]
        m_par, acc_par = accuracies(Y_test, pred, n_classes)

        same_acc = np.array_equal(np.nan_to_num(acc_ref, nan=-1),
                                  np.nan_to_num(acc_par, nan=-1))
        same_conf = np.array_equal(m_ref, m_par)
        all_same = all_same and same_acc and same_conf

        print('{:<26s} {:>7d} {:>8.2f}s {:>8.2f}s {:>8.2f}x {:>10s} {:>10s}'.format(
            key, n_classes, t_sklearn, t_par, t_sklearn / t_par,
            'yes' if same_acc else 'NO', 'yes' if same_conf else 'NO'))
        if not same_acc:
            d = np.abs(np.nan_to_num(acc_ref) - np.nan_to_num(acc_par))
            print('      max accuracy difference {:.6f}, {} of {} clusters differ'.format(
                d.max(), int((d > 0).sum()), len(d)))
        if not same_conf:
            print('      {} of {} predicted labels differ'.format(
                int((model.predict(X_test) != pred).sum()), len(Y_test)))
        cd = np.abs(model.coef_ - coef)
        print('      max |coef difference| {:.3e}, coefficients bitwise equal: {}'.format(
            cd.max(), np.array_equal(model.coef_, coef)))

    print()
    print('VERDICT: parallel per-class fitting reproduces the SCCAF metric exactly: {}'.format(
        all_same))
    return 0 if all_same else 1


if __name__ == '__main__':
    sys.exit(main())
