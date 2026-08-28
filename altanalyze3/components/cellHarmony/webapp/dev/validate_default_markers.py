"""Prove the vectorized _default_marker_genes reproduces the old statistic.

Old code: for each group, scan candidate genes, gap = mean(inside) - mean(outside),
keep the largest gap not already chosen.  New code: same gap, computed for every
gene at once with an indicator-matrix product.  Restricting the new one to the
old one's candidate columns must give identical genes and identical gaps.
"""
import time
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

H5 = ("/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/"
      "webapp/jobs/8677866414b84a07a0b66a56fb0d1321/outputs/combined_with_umap_and_markers.h5ad")
CLUSTER_KEY = "Mm-MarrowAtlas-L4"
LIMIT = 12

adata = ad.read_h5ad(H5)
values_of = adata.obs[CLUSTER_KEY].astype(str).to_numpy()
groups = [str(c) for c in adata.obs[CLUSTER_KEY].cat.categories][:LIMIT]
var_names = [str(v) for v in adata.var_names]
X = adata.X

# A small candidate set so the OLD nested loop finishes in seconds.
step = 200
candidates = list(range(0, len(var_names), step))
print(f"cells={adata.n_obs} genes={len(var_names)} groups={len(groups)} candidates={len(candidates)}")

def dense_column(row):
    col = X[:, row]
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()

# ---- OLD ---------------------------------------------------------------
t0 = time.time()
masks = {g: (values_of == g) for g in groups}
old_chosen, old_gaps, seen = [], [], set()
for group in groups:
    inside = masks[group]
    if not inside.any():
        continue
    best, best_gap = None, 0.0
    for row in candidates:
        v = dense_column(row)
        gap = float(v[inside].mean() - v[~inside].mean())
        if gap > best_gap and var_names[row] not in seen:
            best, best_gap = var_names[row], gap
    if best:
        seen.add(best); old_chosen.append(best); old_gaps.append(best_gap)
old_time = time.time() - t0

# ---- NEW (same maths, restricted to the same candidate columns) --------
t0 = time.time()
index = {g: i for i, g in enumerate(groups)}
codes = pd.Series(values_of).map(index).fillna(-1).to_numpy(dtype=np.int64)
rows = np.nonzero(codes >= 0)[0]
indicator = sp.csr_matrix((np.ones(rows.size), (codes[rows], rows)), shape=(len(groups), adata.n_obs))
sums = np.asarray((indicator @ X).todense(), dtype=np.float64)
counts = np.bincount(codes[rows], minlength=len(groups)).astype(np.float64)
total = np.asarray(X.sum(axis=0), dtype=np.float64).ravel()
inside_mean = sums / np.maximum(counts[:, None], 1.0)
outside_mean = (total[None, :] - sums) / np.maximum(float(adata.n_obs) - counts[:, None], 1.0)
gap_all = inside_mean - outside_mean
new_time = time.time() - t0

cand = np.array(candidates)
new_chosen, new_gaps, seen = [], [], set()
for i, group in enumerate(groups):
    if counts[i] <= 0:
        continue
    order = cand[np.argsort(-gap_all[i][cand])]
    for row in order:
        if gap_all[i][int(row)] <= 0:
            break
        name = var_names[int(row)]
        if name in seen:
            continue
        seen.add(name); new_chosen.append(name); new_gaps.append(float(gap_all[i][int(row)]))
        break

print(f"\nold nested loop: {old_time:.1f}s   new one-pass: {new_time:.3f}s "
      f"({old_time/max(new_time,1e-9):.0f}x on {len(candidates)} candidates)")
print(f"{'group':<14}{'old gene':<14}{'new gene':<14}{'old gap':>10}{'new gap':>10}  match")
same = 0
for g, og, ng, ogap, ngap in zip(groups, old_chosen, new_chosen, old_gaps, new_gaps):
    ok = (og == ng) and abs(ogap - ngap) < 1e-6
    same += ok
    print(f"{g:<14}{og:<14}{ng:<14}{ogap:10.5f}{ngap:10.5f}  {'yes' if ok else 'NO'}")
print(f"\nidentical: {same}/{len(old_chosen)} groups")

# Also check the gap value itself for every (group, candidate) pair, not only the winner.
maxdiff = 0.0
for i, group in enumerate(groups):
    inside = masks[group]
    for row in candidates[:40]:
        v = dense_column(row)
        ref = float(v[inside].mean() - v[~inside].mean())
        maxdiff = max(maxdiff, abs(ref - float(gap_all[i][row])))
print(f"max |old gap - new gap| over {len(groups)*40} (group, gene) pairs: {maxdiff:.3e}")
