P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

old = '''    adata = cache["adata"]
    _, groups, values_of = _group_axis(cache, group_by)
    var_names = [str(v) for v in cache["var_names"]]
    if not groups or not var_names:
        return []

    # A sample keeps this responsive on a large matrix; the strongest markers
    # are not hiding in the tail of the gene list.
    step = max(1, len(var_names) // 4000)
    candidates = list(range(0, len(var_names), step))
    masks = {g: (values_of == g) for g in groups}
    chosen, seen = [], set()
    for group in groups[:limit]:
        inside = masks[group]
        if not inside.any():
            continue
        best, best_gap = None, 0.0
        for row in candidates:
            values = _dense_column(adata, row)
            gap = float(values[inside].mean() - values[~inside].mean())
            if gap > best_gap and var_names[row] not in seen:
                best, best_gap = var_names[row], gap
        if best:
            seen.add(best)
            chosen.append(best)
    return chosen'''

new = '''    adata = cache["adata"]
    _, groups, values_of = _group_axis(cache, group_by)
    var_names = [str(v) for v in cache["var_names"]]
    if not groups or not var_names:
        return []
    groups = list(groups[:limit])

    # One pass over the matrix, not one column slice per (gene, group). The
    # previous version sliced `adata.X[:, row]` inside a nested loop: on a
    # 2,797-cell CSR matrix that is 7.7 ms per slice, 4,764 sampled genes x 12
    # groups = 57,168 slices, measured at 7.3 minutes for one figure. The
    # endpoint held the event loop for all of it, so every other panel in the
    # app froze. An indicator matrix multiplied into X gives the same sums for
    # every gene at once, and it reads every gene rather than a sample of them.
    X = adata.X
    n_cells, n_genes = X.shape
    index = {g: i for i, g in enumerate(groups)}
    codes = pd.Series(values_of).map(index).fillna(-1).to_numpy(dtype=np.int64)
    inside_any = codes >= 0
    if not inside_any.any():
        return []
    rows = np.nonzero(inside_any)[0]
    indicator = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.float64), (codes[rows], rows)),
        shape=(len(groups), n_cells))
    sums = indicator @ X
    sums = np.asarray(sums.todense() if sp.issparse(sums) else sums, dtype=np.float64)
    counts = np.bincount(codes[rows], minlength=len(groups)).astype(np.float64)
    # The contrast is against every other cell in the dataset, which is what the
    # previous `values[~inside]` measured, so the statistic is unchanged.
    total = np.asarray(X.sum(axis=0), dtype=np.float64).ravel()
    inside_mean = sums / np.maximum(counts[:, None], 1.0)
    outside_mean = (total[None, :] - sums) / np.maximum(float(n_cells) - counts[:, None], 1.0)
    gap = inside_mean - outside_mean

    chosen, seen = [], set()
    for position, group in enumerate(groups):
        if counts[position] <= 0:
            continue
        for row in np.argsort(-gap[position])[:200]:
            if gap[position][int(row)] <= 0:
                break
            name = var_names[int(row)]
            if name in seen:
                continue
            seen.add(name)
            chosen.append(name)
            break
    return chosen'''

assert src.count(old) == 1, src.count(old)
open(P, "w").write(src.replace(old, new))
print("patch2 applied")
