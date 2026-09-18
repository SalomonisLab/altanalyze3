"""Shared data-backed Chat protocols for uploaded and precomputed scALABLE data.

Executors originate in scalable_viewer. Both interfaces call these same methods.
Donor aggregation is supplied by the adapter and cached independently of questions.
"""
from typing import Any, Dict, List, Optional
from pathlib import Path
import numpy as np
from altanalyze3.components.visualization.scalable_viewer import bundle_meta


def _donor_axis(ds):
    return ds.donor_axis()


def _donor_pseudobulk(ds, rows, state='', min_cells=5):
    return ds.donor_pseudobulk(rows,state,min_cells)


def _candidate_rows(ds, limit=3000):
    return list(range(0,len(ds.symbols),max(1,len(ds.symbols)//limit)))[:limit]


def _bundle_default_genes(ds, groups):
    rows=ds.markers()
    return list(dict.fromkeys(r['gene'] for r in rows if r['cluster'] in groups))[:6] or list(ds.symbols[:6])

_REFERENCE_WORDS = {'non','no','not','control','controls','ctrl','healthy','normal','never','none','negative','unaffected','ref'}

def _run_severity_gradient(ds, state: str, covariate: str) -> Dict[str, Any]:
    """Genes whose per-donor level tracks a numeric clinical variable.

    A case-control contrast says whether a gene differs between two labelled
    groups. It cannot say whether the gene follows severity. Spearman against a
    measured variable separates a step at the diagnosis boundary from a decline
    with lung function, and only the second is a graded mechanism.
    """
    from scipy import stats as _st

    try:
        kind, values, _ = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded in this dataset.",
                "status": "not_covered"}
    if kind != "numeric":
        return {"answer": f"{covariate} is categorical, so it has no gradient. "
                          "Ask for a group comparison instead.",
                "status": "not_covered"}

    rows = _candidate_rows(ds)
    matrix, donors, counts = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state} to correlate anything.",
                "status": "not_covered"}

    # One covariate value per donor: the mean over that donor's cells, which is
    # exact here because these covariates are donor-level and constant within a
    # donor.
    donor_code, donor_labels = _donor_axis(ds)
    per_donor = []
    values = np.asarray(values, dtype=np.float64)
    for name in donors:
        mask = donor_code == donor_labels.index(name)
        per_donor.append(float(np.nanmean(values[mask])) if mask.any() else np.nan)
    covariate_values = np.asarray(per_donor, dtype=np.float64)
    good = np.isfinite(covariate_values)
    if int(good.sum()) < 8:
        return {"answer": f"{covariate} is missing for too many donors in {state}.",
                "status": "not_covered"}

    # spearmanr returns a scalar for one gene, a matrix for many. Indexing the
    # scalar case as a matrix raised "IndexError: invalid index to scalar
    # variable" for EC general capillary, where only one gene survived the
    # per-donor filter.
    # spearmanr returns a scalar for a single pair and a correlation matrix for
    # many, so the RESULT decides how to read it, not the input shape. Indexing
    # a scalar as a matrix raised IndexError for EC general capillary, where the
    # per-donor filter left one gene.
    # Spearman is Pearson on ranks, so the whole scan is one matrix product.
    # `spearmanr` on the stacked matrix builds the full gene-by-gene correlation
    # matrix, 3022 x 3022 here, to use one column of it: that cost 5.07 s.
    scanned = matrix[:, good]
    n = int(good.sum())

    def rank_rows(block):
        return _st.rankdata(block,axis=1,method="average")

    gene_ranks = rank_rows(scanned)
    covariate_ranks = rank_rows(covariate_values[good][None, :])[0]
    gene_centred = gene_ranks - gene_ranks.mean(axis=1, keepdims=True)
    covariate_centred = covariate_ranks - covariate_ranks.mean()
    denominator = (np.linalg.norm(gene_centred, axis=1)
                   * np.linalg.norm(covariate_centred))
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where(denominator > 0,
                       gene_centred @ covariate_centred / np.maximum(denominator, 1e-12),
                       0.0)
    rho = np.clip(np.nan_to_num(rho, nan=0.0), -1.0, 1.0)
    # The same t approximation scipy uses for Spearman.
    with np.errstate(invalid="ignore", divide="ignore"):
        tstat = rho * np.sqrt((n - 2) / np.maximum(1e-12, 1 - rho ** 2))
    pval = np.nan_to_num(2 * _st.t.sf(np.abs(tstat), max(1, n - 2)), nan=1.0)

    varying = scanned.std(axis=1) > 0
    rho = np.where(varying[:len(rho)], rho, 0.0)
    order = [int(i) for i in np.argsort(-np.abs(rho)) if varying[int(i)]][:25]
    # A linear fit alongside the rank correlation. Spearman says a relationship
    # is monotonic; the slope says how much the gene moves per unit of the
    # variable, which is what "regulated with age" actually claims. Both are
    # reported, because a high rho with a slope near zero is a real ordering of
    # noise.
    x = covariate_values[good]
    x_centred = x - x.mean()
    denominator = float((x_centred ** 2).sum())
    table_rows = []
    points: Dict[str, Any] = {}
    for i in order:
        y = scanned[int(i)]
        slope = float((x_centred * (y - y.mean())).sum() / denominator) if denominator else 0.0
        intercept = float(y.mean() - slope * x.mean())
        predicted = slope * x + intercept
        ss_res = float(((y - predicted) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        gene = ds.symbols[rows[int(i)]]
        table_rows.append([gene, round(float(rho[int(i)]), 4), float(pval[int(i)]),
                           round(slope, 6), round(r2, 4), int(good.sum()),
                           "up with " + covariate if slope > 0 else "down with " + covariate])
        if len(points) < 6:
            points[gene] = {"x": [round(float(v), 4) for v in x],
                            "y": [round(float(v), 5) for v in y],
                            "slope": round(slope, 6), "intercept": round(intercept, 6),
                            "rho": round(float(rho[int(i)]), 4), "r2": round(r2, 4)}
    if not table_rows:
        # Every scanned gene was flat across these donors, so no correlation is
        # defined. An empty table would read as "nothing correlates", which is a
        # different and unsupported claim.
        return {
            "answer": (f"No gene varies across the {int(good.sum())} donors that "
                       f"contribute cells to {state}, so no correlation with "
                       f"{covariate} is defined. Too few cells per donor in this "
                       "state, not an absence of signal."),
            "status": "not_covered",
            "table": {"columns": ["donors with cells in " + state],
                      "rows": [[d] for d in donors[:40]]},
        }
    return {
        "answer": (f"Genes in {state} whose per-donor level tracks {covariate}, "
                   f"Spearman across {int(good.sum())} donors, "
                   f"{len(rows)} genes scanned."),
        "status": "answered",
        "table": {"columns": ["gene", "rho", "p", "slope", "r2", "n_donors", "direction"],
                  "rows": table_rows},
        # One panel per gene: each donor's level against the variable, with the
        # fitted line. A CombPlot here drew all 39 cell states for a correlation
        # computed inside one of them, which showed the wrong thing entirely.
        "plot": {"kind": "gradient", "covariate": covariate, "state": state,
                 "donors": [d for d, ok in zip(donors, good) if ok], "points": points},
    }


def _run_coexpression(ds, state: str, seed: str) -> Dict[str, Any]:
    """Genes co-varying with a seed gene across donors, inside one cell state."""
    seed_row = ds.resolve_gene(seed)
    if seed_row is None:
        return {"answer": f"{seed} is not in this dataset.", "status": "not_covered"}
    rows = _candidate_rows(ds)
    if seed_row not in rows:
        rows = [seed_row] + rows
    matrix, donors, _ = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    seed_index = rows.index(seed_row)
    seed_values = matrix[seed_index]
    centred = matrix - matrix.mean(axis=1, keepdims=True)
    seed_centred = seed_values - seed_values.mean()
    denominator = (np.linalg.norm(centred, axis=1) * np.linalg.norm(seed_centred))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denominator > 0, centred @ seed_centred / np.maximum(denominator, 1e-12), 0.0)
    if np.linalg.norm(seed_centred)==0:
        return {'status':'not_covered','answer':f'{seed} is constant across the contributing samples; correlation is undefined.'}
    r[seed_index] = 1.0
    order = [int(i) for i in np.argsort(-r) if int(i) != seed_index and denominator[int(i)]>0][:25]
    return {
        "answer": (f"Genes co-varying with {seed} across {len(donors)} donors in {state}. "
                   "Correlation is over per-donor pseudobulk, the same axis disease varies along."),
        "status": "answered",
        "table": {"columns": ["gene", "r with " + seed, "n_donors"],
                  "rows": [[ds.symbols[rows[i]], round(float(r[i]), 4), len(donors)]
                           for i in order]},
        "plot": {"kind": "gradient", "covariate":seed, "state":state,"donors":donors,
                 "points":{ds.symbols[rows[i]]:{"x":seed_values.tolist(),"y":matrix[i].tolist(),
                     "rho":float(r[i]),"r2":float(r[i]**2),
                     "slope":float(np.polyfit(seed_values,matrix[i],1)[0]),
                     "intercept":float(np.polyfit(seed_values,matrix[i],1)[1])} for i in order[:6]}},
    }


def _categorical_twin(ds, covariate: str) -> str:
    """A categorical covariate standing for a numeric one, or "".

    Asking which states are depleted in GOLD IV routes the covariate slot to
    `gold_ordinal`, which is numeric, and composition needs groups. The same
    grading is also stored categorically in `Group`, whose levels are named
    GOLD I, II / GOLD III / GOLD IV. Matching on the words shared between the
    numeric covariate's name and a categorical covariate's levels finds it.
    """
    wanted = {w for w in bundle_meta._label_tokens(covariate) if len(w) > 2}
    if not wanted:
        return ""
    for name, info in (ds.covariate_names() or {}).items():
        if (info or {}).get("kind") != "categorical":
            continue
        try:
            kind, _v, labels = ds.covariate_values(name)
        except KeyError:
            continue
        if kind != "categorical" or labels is None or len(labels) < 2:
            continue
        words = {w for label in labels for w in bundle_meta._label_tokens(str(label))}
        words |= set(bundle_meta._label_tokens(name))
        if wanted & words:
            return name
    return ""


def _named_level(question: str, labels) -> int:
    """Index of the level the question names, or -1.

    `Group` carries three GOLD levels. A question about GOLD IV must compare
    GOLD IV, not whichever two levels happen to be encoded first.
    """
    asked = set(bundle_meta._label_tokens(question or ""))
    best, best_score = -1, 0
    for index, label in enumerate(labels):
        tokens = set(bundle_meta._label_tokens(str(label)))
        if tokens and tokens <= asked and len(tokens) > best_score:
            best, best_score = index, len(tokens)
    return best


def _run_composition(ds, covariate: str, question: str = "") -> Dict[str, Any]:
    """Which cell states change in abundance between the levels of a variable."""
    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded here.", "status": "not_covered"}
    swapped = ""
    if kind != "categorical" or labels is None:
        twin = _categorical_twin(ds, covariate)
        if twin:
            swapped, covariate = covariate, twin
            kind, values, labels = ds.covariate_values(covariate)
    if kind != "categorical" or labels is None:
        return {"answer": f"{covariate} is numeric; composition needs groups.",
                "status": "not_covered"}

    donor_code, donor_labels = _donor_axis(ds)
    if donor_code is None:
        return {"answer": "This bundle records no donor column.", "status": "not_covered"}
    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)

    # Fraction of each donor's cells in each state, then the mean per group.
    per_donor = np.zeros((len(donor_labels), len(ds.states)), dtype=np.float64)
    valid=(donor_code>=0)&(states>=0)
    np.add.at(per_donor, (donor_code[valid], states[valid]), np.asarray(ds.cell_weights)[valid])
    totals = per_donor.sum(axis=1, keepdims=True)
    fractions = np.divide(per_donor, np.maximum(totals, 1))
    donor_group = np.full(len(donor_labels), -1, dtype=np.int64)
    for donor in range(len(donor_labels)):
        mask = donor_code == donor
        if mask.any():
            known=np.unique(values[mask][values[mask]>=0])
            if len(known)==1:donor_group[donor]=int(known[0])

    present = [g for g in range(len(labels)) if (donor_group == g).any()]
    if len(present) < 2:
        return {"answer": f"{covariate} has fewer than two groups with donors.",
                "status": "not_covered"}
    # Reference group first, so it takes the sky-blue bar and the disease the
    # red one. Stored order put COPD first and inverted the colours.
    named = _named_level(question, labels)
    if named in present and len(present) > 2:
        # Three GOLD levels: compare the one asked about against the mildest
        # remaining level rather than an arbitrary pair.
        rest = [g for g in present if g != named]
        first, second = _control_first(labels, [rest[0], named])
    else:
        first, second = _control_first(labels, present)
    a = fractions[donor_group == first].mean(axis=0)
    b = fractions[donor_group == second].mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log2((a + 1e-6) / (b + 1e-6))
    # A Mann-Whitney per state: donor fractions are bounded and skewed, so a
    # rank test is the defensible one, and the two groups have unequal sizes.
    from scipy import stats as _st

    a_donors = fractions[donor_group == first]
    b_donors = fractions[donor_group == second]
    if min(len(a_donors),len(b_donors))<2:
        return {'status':'not_covered','answer':'At least two biological samples per group are required for composition summaries.'}
    pvals = []
    for index in range(len(ds.states)):
        x, y = a_donors[:, index], b_donors[:, index]
        if x.size < 3 or y.size < 3 or (x.std() == 0 and y.std() == 0):
            pvals.append(1.0)
            continue
        try:
            pvals.append(float(_st.mannwhitneyu(x, y, alternative="two-sided").pvalue))
        except ValueError:
            pvals.append(1.0)
    # Rank by the test, not by the log ratio. The ratio carries a pseudocount,
    # so a state present at 0.00012 against 0.00000 topped the list on a
    # difference of one ten-thousandth of a donor's cells, which is noise. The
    # rank test asks whether the donors actually separate.
    order = sorted(range(len(ds.states)),
                   key=lambda i: (pvals[i], -abs(float(a[i] - b[i]))))[:25]
    return {
        "answer": ((f"{swapped} is numeric, so the grouped variable {covariate} "
                    f"was used instead. " if swapped else "")
                   + f"Cell-state abundance, {labels[first]} against {labels[second]}, "
                   f"as each donor's share of their own cells. "
                   f"{int((donor_group == first).sum())} and "
                   f"{int((donor_group == second).sum())} donors. "
                   "Abundance and expression are confounded: a depleted state looks "
                   "changed in any pooled expression comparison."),
        "status": "answered",
        "table": {"columns": ["cell state", f"mean fraction {labels[first]}",
                              f"mean fraction {labels[second]}", "log2 ratio", "p"],
                  "rows": [[ds.states[int(i)], round(float(a[int(i)]), 5),
                            round(float(b[int(i)]), 5), round(float(ratio[int(i)]), 4),
                            float(pvals[int(i)])]
                           for i in order]},
        # A paired bar per cell state, one bar per group, so the two frequencies
        # are read directly. A log ratio alone hides whether a state is common
        # in both groups or rare in both.
        # Every donor's own fraction travels with the summary, so the figure can
        # show the biological replicates behind each bar rather than only the
        # mean. A mean of 141 donors and a mean of 4 look identical otherwise.
        "plot": {"kind": "frequency",
                 "groups": [str(labels[first]), str(labels[second])],
                 "n_donors": [int((donor_group == first).sum()),
                              int((donor_group == second).sum())],
                 "states": [ds.states[int(i)] for i in order],
                 "a": [round(float(a[int(i)]), 6) for i in order],
                 "b": [round(float(b[int(i)]), 6) for i in order],
                 # Standard error of the mean, which is what an error bar on a
                 # group mean should show.
                 "a_sem": [round(float(a_donors[:, int(i)].std(ddof=1)
                                       / max(1.0, np.sqrt(a_donors.shape[0]))), 6)
                           for i in order],
                 "b_sem": [round(float(b_donors[:, int(i)].std(ddof=1)
                                       / max(1.0, np.sqrt(b_donors.shape[0]))), 6)
                           for i in order],
                 "a_points": [[round(float(v), 6) for v in a_donors[:, int(i)]]
                              for i in order],
                 "b_points": [[round(float(v), 6) for v in b_donors[:, int(i)]]
                              for i in order],
                 "p": [float(pvals[int(i)]) for i in order]},
    }


def _run_concordance(ds, state: str = "") -> Dict[str, Any]:
    """Where the two cell annotations agree, and where they do not."""
    primary=ds.sv.get('cluster_key') or getattr(ds,'cluster_key','cell_state')
    second = getattr(ds,'annotation_column','') or next((name for name in
        ("TGEN-IPF","Population","celltype_level3","original_cell_type","annotation","reference_cluster","predicted_label")
        if name != primary and name in ds.covariate_names()), "")
    if not second:
        return {"answer": "This bundle carries only one cell annotation.",
                "status": "not_covered"}
    kind, values, labels = ds.covariate_values(second)
    if kind != "categorical" or labels is None:
        return {"answer": f"{second} is not categorical.", "status": "not_covered"}

    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)
    wanted = range(len(ds.states)) if not state else [ds.states.index(state)]
    table_rows = []
    for index in wanted:
        mask = states == index
        total = int(mask.sum())
        if not total:
            continue
        counts = np.bincount(values[mask & (values>=0)], minlength=len(labels))
        for other in np.argsort(-counts)[:3]:
            if counts[int(other)] == 0:
                continue
            table_rows.append([ds.states[int(index)], str(labels[int(other)]),
                               int(counts[int(other)]),
                               round(float(counts[int(other)]) / total, 4)])
    scope = state or f"all {len(ds.states)} cell states"
    return {
        "answer": (f"How {ds.sv.get('cluster_key') or 'cell_state'} maps onto {second}, "
                   f"for {scope}. A state whose top match holds well under 1.0 is a "
                   "boundary case, and a differential there may be an artefact of "
                   "which labelling was used."),
        "status": "answered",
        "table": {"columns": ["cell state", second, "n cells", "fraction"],
                  "rows": table_rows[:60]},
        "plot": {"kind": "heatmap"},
    }


def _run_dose_response(ds, state: str, covariate: str, genes: List[str]) -> Dict[str, Any]:
    """Whether a gene changes stepwise across the levels of an ordered variable."""
    from scipy import stats as _st

    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return {"answer": f"{covariate} is not recorded here.", "status": "not_covered"}
    if kind != "categorical" or labels is None:
        # "across GOLD stages" resolves to `gold_ordinal`, which is numeric.
        # The same variable exists categorically as `Group`, and an ordered
        # question wants the levels, so swap rather than refuse.
        swapped = ""
        for name, info in (ds.covariate_names() or {}).items():
            if info.get("kind") != "categorical":
                continue
            stem = covariate.lower().split("_")[0][:4]
            if stem and (stem in name.lower() or name.lower()[:4] == stem):
                swapped = name
                break
        if not swapped and "Group" in (ds.covariate_names() or {}):
            swapped = "Group"
        if not swapped:
            return {"answer": f"{covariate} is numeric; ask for a gradient instead.",
                    "status": "not_covered"}
        covariate = swapped
        kind, values, labels = ds.covariate_values(covariate)

    kind,values,labels=ds.ordered_covariate_values(covariate)

    wanted = [g for g in (genes or []) if ds.resolve_gene(g) is not None]
    if not wanted:
        wanted = _bundle_default_genes(ds, [state])[:6]
    rows = [ds.resolve_gene(g) for g in wanted]
    matrix, donors, _ = _donor_pseudobulk(ds, rows, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    donor_code, donor_labels = _donor_axis(ds)
    values = np.asarray(values, dtype=np.int64)
    donor_level = []
    for name in donors:
        mask = donor_code == donor_labels.index(name)
        known=np.unique(values[mask][values[mask]>=0])
        donor_level.append(int(known[0]) if len(known)==1 else -1)
    level = np.asarray(donor_level, dtype=np.float64)
    good = level >= 0

    if int(good.sum())<8 or len(set(level[good]))<2:
        return {'status':'not_covered','answer':'At least eight samples across two ordered levels are required.'}
    table_rows = []; points={}
    for position, gene in enumerate(wanted):
        per_level = []
        for index in sorted(set(level[good].astype(int))):
            selected = good & (level == index)
            per_level.append(round(float(matrix[position][selected].mean()), 4)
                             if selected.any() else None)
        rho, pval = _st.spearmanr(level[good], matrix[position][good])
        clean = [v for v in per_level if v is not None]
        monotonic = (all(x <= y for x, y in zip(clean, clean[1:]))
                     or all(x >= y for x, y in zip(clean, clean[1:])))
        x=level[good];y=matrix[position][good]
        slope,intercept=np.polyfit(x,y,1)
        spread=float(((y-y.mean())**2).sum())
        points[gene]={'x':x.tolist(),'y':y.tolist(),'slope':round(float(slope),6),'intercept':round(float(intercept),6),
                      'rho':round(float(rho),4) if np.isfinite(rho) else None,
                      'r2':round(1-float(((y-(slope*x+intercept))**2).sum())/spread,4) if spread>0 else None}
        table_rows.append([gene, str(per_level), round(float(rho), 4) if np.isfinite(rho) else None,
                           float(pval) if np.isfinite(pval) else None, "yes" if monotonic else "no"])
    order = [str(l) for l in labels]
    return {
        "answer": (f"{', '.join(wanted)} across {covariate} in {state}, "
                   f"{int(good.sum())} donors. Levels in order: {', '.join(order)}. "
                   "A stepwise change behaves like a progression marker; a change only "
                   "at the extreme behaves like an end-stage marker."),
        "status": "answered",
        "table": {"columns": ["gene", "mean per level", "trend rho", "p", "monotonic"],
                  "rows": table_rows},
        "plot": {"kind": "gradient", "state":state, "covariate":covariate, "levels":order,
                 "donors":[d for d,ok in zip(donors,good) if ok],"points":points},
    }


def _run_pathway_program(assets: Dict[str, Any], ds, state: str, contrast: str,
                         direction: str = "both", limit: int = 25) -> Dict[str, Any]:
    """Enriched processes for one cell state, read from the shipped GO-Elite table.

    A gene list is not a mechanism. GO-Elite is precomputed for each comparison,
    so the enriched terms are read directly rather than re-derived, and the
    genes driving each term come with them.

    The table keys its rows as `<cell state>__<direction>`, so a question about
    what is up in AT2 reads `AT2__up`.
    """
    import csv as _csv

    # The assets record the file directly; an earlier version guessed at the
    # directory layout with rglob and found nothing.
    head = (contrast or "").split("::")[0]
    differential = (assets or {}).get("differential") or {}
    entry = differential.get(contrast) or {}
    if not entry:
        entry = {}  # Never substitute another comparison's enrichment.
    path = str((entry or {}).get("goelite_tsv") or "")
    if not path or not Path(path).exists():
        return {"answer": f"No GO-Elite results are shipped for {head}.",
                "status": "not_covered"}
    candidates = [Path(path)]

    wanted = set()
    if direction in ("up", "both"):
        wanted.add(f"{state}__up")
    if direction in ("down", "both"):
        wanted.add(f"{state}__down")

    rows = []
    with candidates[0].open() as handle:
        for row in _csv.DictReader(handle, delimiter="\t"):
            if row.get("population") in wanted:
                try:
                    z = float(row.get("z_score") or 0)
                    fdr = float(row["fdr"]) if row.get("fdr") not in (None, "") else 1.
                except ValueError:
                    continue
                rows.append((z, fdr, row))
    if not rows:
        return {
            "answer": (f"GO-Elite is shipped for {head} but holds no terms for {state}. "
                       "It covers the cell states with enough differential genes to test."),
            "status": "not_covered",
        }

    # A term overlapping one gene can reach Z=27 and says nothing: the Z score
    # rewards a small expected count. Terms are ranked by Z but must overlap at
    # least two genes, and the count that was dropped is reported rather than
    # hidden, because a state whose only terms are one-gene terms has no
    # enrichment worth reading.
    single = [item for item in rows if int(float(item[2].get("overlap") or 0)) < 2]
    rows = [item for item in rows if int(float(item[2].get("overlap") or 0)) >= 2]
    if not rows:
        return {
            "answer": (f"GO-Elite found {len(single)} terms for {state} in {head}, and every "
                       "one overlaps a single gene. A one-gene term reaches a high Z score "
                       "because the expected count is tiny, so there is no enrichment here "
                       "worth reporting."),
            "status": "not_covered",
        }
    rows.sort(key=lambda item: -item[0])
    table_rows, chart = [], []
    for z, fdr, row in rows[:limit]:
        way = "up" if str(row.get("population", "")).endswith("__up") else "down"
        genes = str(row.get("overlap_genes") or "").replace("|", ", ")[:70]
        table_rows.append([row.get("term_name"), way, round(z, 3), fdr,
                           int(float(row.get("overlap") or 0)), genes])
        chart.append([f"{row.get('term_name')} ({way})", z if way == "up" else -z])

    return {
        "answer": (f"Enriched processes in {state} for {head}: {len(rows)} terms "
                   f"overlapping two or more genes, ranked by Z score. "
                   f"{len(single)} one-gene terms were left out, because a term with one "
                   "overlapping gene reaches a high Z on a tiny expected count. The genes "
                   "driving each term are listed beside it."),
        "status": "answered",
        "table": {"columns": ["term", "direction", "z", "fdr", "n overlap", "genes"],
                  "rows": table_rows},
        # Name the column explicitly: the last column here is the gene list, a
        # string, so a bar chart guessing "rightmost column" drew nothing.
        "plot": {"kind": "barchart", "value_column": "z",
                 "label_column": "term", "sign_column": "direction"},
        "chart": chart,
    }


def _control_first(labels, present, reference: str = ""):
    """Order two levels so the reference group comes first.

    The bundle stores levels in whatever order the categories were encoded, and
    for `copd_status` that is COPD before non-COPD. Drawn as stored, the disease
    takes the control colour. `reference`, when given, is the right-hand side of
    a contrast id (`COPD_vs_non-COPD` -> `non-COPD`) and decides it outright;
    otherwise a level whose label carries a negation or health word is the
    reference. Neither test firing leaves the stored order alone.
    """
    first, second = present[0], present[1]
    words = lambda index: {w.lower() for w in
                           str(labels[index]).replace("-", " ").replace("_", " ").split()}
    if reference:
        wanted = {w.lower() for w in
                  reference.replace("-", " ").replace("_", " ").split() if w}
        score_first = len(wanted & words(first))
        score_second = len(wanted & words(second))
        if score_second > score_first:
            return second, first
        if score_first > score_second:
            return first, second
    if (words(second) & _REFERENCE_WORDS) and not (words(first) & _REFERENCE_WORDS):
        return second, first
    return first, second


def _frequency_payload(ds, covariate: str, states_in_order: Optional[List[str]] = None,
                       limit: int = 25, case_label: str = "",
                       control_label: str = "") -> Optional[Dict[str, Any]]:
    """Per-donor cell-state frequency for the two groups of a categorical variable.

    Shared by `composition_shift`, which ranks by the rank test, and by
    `most_affected_state`, which keeps its own transcriptional ranking and uses
    this only for the figure. Returns None when the variable has fewer than two
    groups with donors.
    """
    from scipy import stats as _st

    try:
        kind, values, labels = ds.covariate_values(covariate)
    except KeyError:
        return None
    if kind != "categorical" or labels is None:
        return None
    donor_code, donor_labels = _donor_axis(ds)
    if donor_code is None:
        return None

    states = np.asarray(ds.state_code, dtype=np.int64)
    values = np.asarray(values, dtype=np.int64)
    per_donor = np.zeros((len(donor_labels), len(ds.states)), dtype=np.float64)
    valid=(donor_code>=0)&(states>=0)
    np.add.at(per_donor, (donor_code[valid], states[valid]), np.asarray(ds.cell_weights)[valid])
    totals = per_donor.sum(axis=1, keepdims=True)
    fractions = np.divide(per_donor, np.maximum(totals, 1))

    donor_group = np.full(len(donor_labels), -1, dtype=np.int64)
    for donor in range(len(donor_labels)):
        mask = donor_code == donor
        if mask.any():
            known=np.unique(values[mask][values[mask]>=0])
            if len(known)==1:donor_group[donor]=int(known[0])
    present = [g for g in range(len(labels)) if (donor_group == g).any()]
    if len(present) < 2:
        return None
    # When the caller knows which two levels the contrast compared, use exactly
    # those. `gold_ordinal` carries four levels, and taking the first two
    # present drew GOLD III for a GOLD IV contrast.
    first = second = None
    if case_label and control_label:
        by_key = {bundle_meta._label_tokens(str(labels[g])): g for g in present}
        first = by_key.get(bundle_meta._label_tokens(control_label))
        second = by_key.get(bundle_meta._label_tokens(case_label))
    if first is None or second is None or first == second:
        first, second = _control_first(labels, present)

    a_donors = fractions[donor_group == first]
    b_donors = fractions[donor_group == second]
    if min(len(a_donors),len(b_donors))<2:return None
    a, b = a_donors.mean(axis=0), b_donors.mean(axis=0)
    pvals = []
    for index in range(len(ds.states)):
        x, y = a_donors[:, index], b_donors[:, index]
        if x.size < 3 or y.size < 3 or (x.std() == 0 and y.std() == 0):
            pvals.append(1.0)
            continue
        try:
            pvals.append(float(_st.mannwhitneyu(x, y, alternative="two-sided").pvalue))
        except ValueError:
            pvals.append(1.0)

    if states_in_order:
        order = [ds.states.index(s) for s in states_in_order if s in ds.states][:limit]
    else:
        order = sorted(range(len(ds.states)),
                       key=lambda i: (pvals[i], -abs(float(a[i] - b[i]))))[:limit]

    sem = lambda block, i: float(block[:, i].std(ddof=1) / max(1.0, np.sqrt(block.shape[0])))
    return {
        "kind": "frequency",
        "groups": [str(labels[first]), str(labels[second])],
        "n_donors": [int((donor_group == first).sum()), int((donor_group == second).sum())],
        "states": [ds.states[int(i)] for i in order],
        "a": [round(float(a[int(i)]), 6) for i in order],
        "b": [round(float(b[int(i)]), 6) for i in order],
        "a_sem": [round(sem(a_donors, int(i)), 6) for i in order],
        "b_sem": [round(sem(b_donors, int(i)), 6) for i in order],
        "a_points": [[round(float(v), 6) for v in a_donors[:, int(i)]] for i in order],
        "b_points": [[round(float(v), 6) for v in b_donors[:, int(i)]] for i in order],
        "p": [float(pvals[int(i)]) for i in order],
        "_table": {"a": a, "b": b, "p": pvals, "order": order,
                   "labels": [str(labels[first]), str(labels[second])]},
    }


def _run_most_affected_state(ds, contrast: str, limit: int = 39) -> Dict[str, Any]:
    """Rank cell states by how strongly they respond to a contrast.

    The question is where the disease acts, not which genes move. Ranking the
    states by the number of genes passing FDR, and by the median effect among
    them, says which compartment to look at before any gene list is opened.

    This used to be aliased to the generic differential, so it answered with a
    gene volcano: the right numbers for a different question.
    """
    manifest = ds.deg_manifest() or {}
    comparisons = manifest.get("comparisons", [])
    chosen = next((c for c in comparisons if (c.get("id") or "") == contrast), None)
    if chosen is None:
        return {"answer": "This bundle carries no per-cell-state comparison.",
                "status": "not_covered"}

    table = ds.deg_table(chosen.get("id"), 1000000, 0.05, None) or {}
    rows = table.get("rows") if isinstance(table, dict) else table
    if not rows:
        return {"answer": f"{chosen.get('id')} returned no genes below FDR 0.05.",
                "status": "not_covered"}

    per_state: Dict[str, List[float]] = {}
    top_gene: Dict[str, tuple] = {}
    for row in rows:
        state = str(row.get("population") or row.get("cluster") or "")
        if not state:
            continue
        try:
            fold = float(row.get("log2fc") or 0.0)
        except (TypeError, ValueError):
            continue
        per_state.setdefault(state, []).append(abs(fold))
        best = top_gene.get(state)
        if best is None or abs(fold) > best[1]:
            top_gene[state] = (row.get("gene"), abs(fold), fold)

    ranked = sorted(per_state.items(), key=lambda kv: (-len(kv[1]), -float(np.median(kv[1]))))
    table_rows, chart = [], []
    for state, folds in ranked[:limit]:
        gene, _, signed = top_gene.get(state, ("", 0.0, 0.0))
        table_rows.append([state, len(folds), round(float(np.median(folds)), 4),
                           gene, round(float(signed), 4)])
        chart.append([state, len(folds)])

    covered = len(per_state)
    # The abundance figure must split donors on the same variable the
    # differential compared, so it is derived from the contrast id.
    frequency = None
    try:
        if hasattr(ds,'comparison_covariate'):
            grouping,case_value,control_value=ds.comparison_covariate(contrast)
        else:
            categorical = bundle_meta._categorical_covariates(ds)
            grouping, case_value, control_value = bundle_meta._contrast_group_field(
                ds, chosen, categorical)
    except Exception:                                    # noqa: BLE001
        grouping, case_value, control_value = "", "", ""
    if grouping:
        frequency = _frequency_payload(ds, grouping, [r[0] for r in table_rows],
                                       case_label=case_value,
                                       control_label=control_value)
        if frequency:
            frequency.pop("_table", None)
    return {
        "answer": (f"Cell states ranked by how strongly they respond to "
                   f"{chosen.get('id')}, counting genes below FDR 0.05. "
                   f"{covered} of the {len(ds.states)} states have any, and the "
                   "contrast is only computed for the states with enough cells. "
                   + (f"The table ranks by differential genes. The bars show the "
                      f"same states' abundance in each group, "
                      f"{frequency['groups'][0]} against {frequency['groups'][1]}, "
                      f"as each donor's share of their own cells, with one dot per "
                      f"donor and standard error. A state can change in expression, "
                      f"in abundance, or in both, and the two are confounded: a "
                      f"depleted state looks changed in any pooled comparison."
                      if frequency else
                      "Where the disease acts, before any gene list is opened.")),
        "status": "answered",
        "table": {"columns": ["cell state", "n significant", "median |log2fc|",
                              "top gene", "its log2fc"],
                  "rows": table_rows},
        # The question spans two things: which states change transcriptionally,
        # and which change in abundance. The table ranks by differential genes;
        # the figure shows each group's cell frequency for those same states, in
        # the same order, so both are visible at once.
        "plot": frequency or {"kind": "barchart", "value_column": "n significant",
                              "label_column": "cell state"},
    }


def _run_donor_heterogeneity(ds, state: str, contrast: str,
                             n_genes: int = 40) -> Dict[str, Any]:
    """Whether a signature is carried by most donors or by a few.

    A significant fold change is a statement about group means. With 178 donors
    a gene can clear FDR because a handful carry it strongly while the rest show
    nothing, and that is a subtype rather than a disease mechanism.

    So the signature is scored per donor: the genes up in the contrast minus the
    genes down, each standardised across donors, inside the one cell state. A
    heatmap of those genes with donors ordered by their score shows at a glance
    whether the case donors form a block or scatter through the controls.

    This used to answer with a gene volcano, which shows the group means the
    question is asking to look past.
    """
    manifest = ds.deg_manifest() or {}
    comparisons = manifest.get("comparisons", [])
    chosen = next((c for c in comparisons if (c.get("id") or "") == contrast), None)
    if chosen is None:
        return {"answer": "This bundle carries no per-cell-state comparison.",
                "status": "not_covered"}

    table = ds.deg_table(chosen.get("id"), 100000, 0.05, state) or {}
    rows = (table.get("rows") if isinstance(table, dict) else table) or []
    if not rows:
        return {"answer": f"{chosen.get('id')} is not computed for {state}.",
                "status": "not_covered"}

    ups, downs = [], []
    for row in sorted(rows, key=lambda r: -abs(float(r.get("log2fc") or 0))):
        gene = row.get("gene")
        row_index = ds.resolve_gene(str(gene))
        if row_index is None:
            continue
        if float(row.get("log2fc") or 0) > 0 and len(ups) < n_genes // 2:
            ups.append((gene, row_index))
        elif float(row.get("log2fc") or 0) < 0 and len(downs) < n_genes // 2:
            downs.append((gene, row_index))
        if len(ups) >= n_genes // 2 and len(downs) >= n_genes // 2:
            break
    signature = ups + downs
    if len(signature) < 4:
        return {"answer": f"Too few significant genes in {state} to score a signature.",
                "status": "not_covered"}

    indices = [i for _, i in signature]
    matrix, donors, counts = _donor_pseudobulk(ds, indices, state)
    if matrix is None or len(donors) < 8:
        return {"answer": f"Too few donors contribute cells to {state}.",
                "status": "not_covered"}

    # Standardise each gene across donors so one loud gene cannot set the score.
    centre = matrix.mean(axis=1, keepdims=True)
    spread = matrix.std(axis=1, keepdims=True)
    z = np.divide(matrix - centre, np.where(spread > 0, spread, 1.0))
    up_rows = np.arange(len(ups))
    down_rows = np.arange(len(ups), len(signature))
    score = (z[up_rows].mean(axis=0) if up_rows.size else 0.0) - (z[down_rows].mean(axis=0) if down_rows.size else 0.0)

    donor_code, donor_labels = _donor_axis(ds)
    group_of = {}
    case_name = (chosen.get("id") or "").split("::")[0].split("_vs_")[0].replace("_", " ")
    for name in (getattr(ds,"comparison_group_field",""), "copd_status", "dx_category", "condition", "Group"):
        if name in (ds.covariate_names() or {}):
            kind, values, labels = ds.covariate_values(name)
            if kind == "categorical" and labels is not None:
                values = np.asarray(values, dtype=np.int64)
                for donor in donors:
                    mask = donor_code == donor_labels.index(donor)
                    if mask.any():
                        known=np.unique(values[mask][values[mask]>=0])
                        if len(known)==1:group_of[donor]=str(labels[int(known[0])])
            break

    order = np.argsort(-score)
    ordered = [donors[int(i)] for i in order]
    groups = [group_of.get(d, "") for d in ordered]
    levels = sorted({g for g in groups if g})
    # How far up the ranking the case donors sit: if the signature is universal
    # they are spread through it, if it is a subtype they cluster at the top.
    summary = ""
    if len(levels) == 2:
        case = next((l for l in levels if case_name.lower().split()[0] in l.lower()), levels[0])
        positions = [i for i, g in enumerate(groups) if g == case]
        n_case = len(positions)
        top_half = sum(1 for i in positions if i < len(ordered) / 2)
        summary = (f" {top_half} of the {n_case} {case} donors fall in the top half of the "
                   f"ranking; an even split would put about {n_case // 2} there.")

    return {
        "answer": (f"The {chosen.get('id')} signature in {state}, scored per donor: "
                   f"{len(ups)} up genes minus {len(downs)} down genes, each standardised "
                   f"across the {len(donors)} donors with cells in this state.{summary} "
                   "A fold change is a group mean; this is who carries it."),
        "status": "answered",
        "table": {"columns": ["donor", "signature score", "group", "n cells"],
                  "rows": [[donors[int(i)], round(float(score[int(i)]), 4),
                            group_of.get(donors[int(i)], ""), int(counts[int(i)])]
                           for i in order]},
        "plot": {"kind": "signature",
                 "genes": [g for g, _ in signature],
                 "n_up": len(ups),
                 "donors": ordered,
                 "groups": groups,
                 "score": [round(float(score[int(i)]), 4) for i in order],
                 "z": [[round(float(z[r][int(c)]), 3) for c in order]
                       for r in range(len(signature))]},
    }

