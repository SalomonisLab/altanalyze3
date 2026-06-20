"""Count-aware junction->isoform ATTRIBUTION -- an alternative to isoform_predict.py.

isoform_predict.py answers "which annotated isoform best represents each side of the
event?" by choosing, for each competing junction, the best-supported isoform that
CONTAINS that junction (Ensembl > coding > NMD > length). By construction it can only
ever name an isoform carrying the junction, so when a junction's PSI moves because a
COMPETING isoform changed (the denominator confound), it cannot identify the driver.

This module instead answers "which isoform's regulation DROVE the PSI change?" using a
deterministic one-at-a-time ΔPSI contribution decomposition (method M4 in the
code-dev benchmark): junction counts are modelled as n = A x (A = junction-by-isoform
incidence, x = isoform abundance), so PSI_u = n_u / Σ_{v∈clique(u)} n_v; the
contribution of each isoform is the ΔPSI that isoform alone would produce, holding all
others at baseline. The driver is the isoform with the largest |contribution|, and the
mechanism is 'inclusion_isoform_change' if the driver contains the junction or
'competing_isoform_change' if it does not. An EM deconvolution fallback (M5) is
available when isoform abundances are unavailable (junction counts only).

It writes a NEW '*-attributed.txt' next to the input and never overwrites the original
event table; the current Isoform_1/Isoform_2 columns are preserved for comparison.

Analogous prior work for the decomposition: control-coefficient sensitivity in
metabolic control analysis (Kacser & Burns 1973; Heinrich & Rapoport 1974), additive
index decomposition (Ang 2005), and the one-at-a-time case of Shapley attribution
(Shapley 1953; Lundberg & Lee 2017). The EM fallback is the RSEM/kallisto isoform-
quantification EM (Dempster, Laird & Rubin 1977; Li & Dewey 2011; Bray et al. 2016).
"""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

_NOVEL_RE = re.compile(r"_(\d+)$")


# --------------------------------------------------------------------------- #
# structure / coordinate helpers (self-contained; mirror psi/psi_single cliques)
# --------------------------------------------------------------------------- #
def _base_token(t: str) -> str:
    return _NOVEL_RE.sub("", t)


def _novel_coord(t: str) -> Optional[int]:
    m = _NOVEL_RE.search(t)
    return int(m.group(1)) if m else None


def _structure_junctions(tokens: List[str]) -> List[str]:
    return [f"{tokens[i]}-{tokens[i + 1]}" for i in range(len(tokens) - 1)]


def _strip_gene(j: str) -> str:
    return j.split(":", 1)[1] if ":" in j else j


def _strip_version(x: str) -> str:
    return x.split(".", 1)[0].strip()


def load_exon_coords(exon_path: str, genes: set) -> Tuple[Dict, Dict[str, str]]:
    """(gene, base_token) -> (lo, hi) genomic; gene -> strand. Streams the file."""
    coords: Dict[Tuple[str, str], Tuple[int, int]] = {}
    strand: Dict[str, str] = {}
    with open(exon_path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 or p[0] not in genes:
                continue
            try:
                a, b = int(p[4]), int(p[5])
            except ValueError:
                continue
            coords[(p[0], p[1])] = (min(a, b), max(a, b))
            strand.setdefault(p[0], p[3])
    return coords, strand


def _junction_interval(gene, ta, tb, strand, coords) -> Optional[Tuple[int, int]]:
    a, b = _base_token(ta), _base_token(tb)
    if (gene, a) not in coords or (gene, b) not in coords:
        return None
    a_lo, a_hi = coords[(gene, a)]
    b_lo, b_hi = coords[(gene, b)]
    na, nb = _novel_coord(ta), _novel_coord(tb)
    if strand == "+":
        donor = na if na is not None else a_hi
        acceptor = nb if nb is not None else b_lo
    else:
        donor = na if na is not None else a_lo
        acceptor = nb if nb is not None else b_hi
    lo, hi = min(donor, acceptor), max(donor, acceptor)
    return (lo, hi + 1) if lo == hi else (lo, hi)


def _overlap(a, b) -> bool:
    return max(a[0], b[0]) < min(a[1], b[1])


# --------------------------------------------------------------------------- #
# attribution core (M4 contribution decomposition; M5 EM fallback)
# --------------------------------------------------------------------------- #
def _psi_u(counts, u_idx, clique_idx) -> float:
    d = float(sum(counts[k] for k in clique_idx))
    return counts[u_idx] / d if d > 0 else np.nan


def _solo_contributions(A, x0, x1, u_idx, clique_idx) -> np.ndarray:
    j0 = A @ x0
    p0 = _psi_u(j0, u_idx, clique_idx)
    out = np.zeros(A.shape[1])
    for i in range(A.shape[1]):
        xj = x0.copy()
        xj[i] = x1[i]
        out[i] = _psi_u(A @ xj, u_idx, clique_idx) - p0
    return out


def _poisson_em(A, j, n_iter=200, tol=1e-6) -> np.ndarray:
    x = np.full(A.shape[1], max(j.sum(), 1.0) / A.shape[1])
    col = A.sum(0); col[col == 0] = 1.0
    for _ in range(n_iter):
        lam = A @ x; lam[lam < 1e-9] = 1e-9
        x_new = x * (A.T @ (j / lam)) / col
        if np.max(np.abs(x_new - x)) < tol:
            return x_new
        x = x_new
    return x


class GeneModel:
    """All isoforms of one gene with structures, abundances and incidence."""

    def __init__(self, gene_id, strand, iso_ids, structures, x0, x1, coords):
        self.gene_id = gene_id
        self.strand = strand
        self.iso_ids = iso_ids
        self.structures = structures            # list of token lists
        self.x0 = np.asarray(x0, float)
        self.x1 = np.asarray(x1, float)
        # unique junctions across isoforms
        jset, seen = [], set()
        for toks in structures:
            for jn in _structure_junctions(toks):
                if jn not in seen:
                    seen.add(jn); jset.append(jn)
        self.junctions = jset
        self.jidx = {j: i for i, j in enumerate(jset)}
        A = np.zeros((len(jset), len(iso_ids)))
        for ci, toks in enumerate(structures):
            for jn in _structure_junctions(toks):
                A[self.jidx[jn], ci] = 1.0
        self.A = A
        # genomic intervals for cliques
        self.interval = {}
        for jn in jset:
            ta, tb = jn.rsplit("-", 1)
            iv = _junction_interval(gene_id, ta, tb, strand, coords)
            if iv is not None:
                self.interval[jn] = iv

    def clique_idx(self, u: str, partner: Optional[str]) -> List[int]:
        """Junctions competing with u: genomic-overlap clique (psi_single), else the
        minimal {u, partner} pair when coordinates are unavailable."""
        members = set()
        if u in self.interval:
            ui = self.interval[u]
            for v in self.junctions:
                if v in self.interval and _overlap(ui, self.interval[v]):
                    members.add(v)
        members.add(u)
        if partner and partner in self.jidx:
            members.add(partner)
        return [self.jidx[v] for v in members if v in self.jidx]


def attribute_event(gm: GeneModel, u: str, partner: Optional[str],
                    method: str = "auto", min_iso_reads: float = 5.0) -> Dict:
    if u not in gm.jidx:
        return {"driver": "", "mechanism": "no_junction_in_models", "direction": 0,
                "contribution": 0.0, "method": "", "ranked": ""}
    u_idx = gm.jidx[u]
    clq = gm.clique_idx(u, partner)
    if len(clq) < 2:
        return {"driver": "", "mechanism": "no_competition", "direction": 0,
                "contribution": 0.0, "method": "", "ranked": ""}

    x0, x1, used = gm.x0, gm.x1, "contribution_decomp"
    if method == "em" or (method == "auto" and (x0.sum() < min_iso_reads or x1.sum() < min_iso_reads)):
        # reconstruct abundances from junction counts (n = A x) via EM
        j0 = gm.A @ gm.x0; j1 = gm.A @ gm.x1
        x0, x1, used = _poisson_em(gm.A, j0), _poisson_em(gm.A, j1), "em_fallback"

    contribs = _solo_contributions(gm.A, x0, x1, u_idx, clq)
    order = np.argsort(-np.abs(contribs))
    ranked = [(gm.iso_ids[i], float(contribs[i])) for i in order if abs(contribs[i]) > 1e-4]
    if not ranked:
        return {"driver": "", "mechanism": "no_signal", "direction": 0,
                "contribution": 0.0, "method": used, "ranked": ""}
    top = order[0]
    with_u = gm.A[u_idx, :] == 1.0
    mech = "inclusion_isoform_change" if with_u[top] else "competing_isoform_change"
    if len(order) >= 2:
        c0, c1 = abs(contribs[order[0]]), abs(contribs[order[1]])
        if c1 >= 0.6 * c0 and (with_u[order[0]] != with_u[order[1]]):
            mech = "mixed"
    direction = int(np.sign(x1[top] - x0[top]))
    return {
        "driver": gm.iso_ids[top], "mechanism": mech, "direction": direction,
        "contribution": float(contribs[top]), "method": used,
        "ranked": ";".join(f"{i}:{c:+.3f}" for i, c in ranked[:5]),
    }


# --------------------------------------------------------------------------- #
# loaders for real AltAnalyze3 differential outputs
# --------------------------------------------------------------------------- #
def _detect_mean_cols(cols: List[str]) -> Tuple[str, str]:
    means = [c for c in cols if c.startswith("Mean_")]
    if len(means) < 2:
        raise ValueError(f"need two Mean_<cond> columns, found {means}")
    return means[0], means[1]


def load_isoform_abundances(stats_tsv: str, control_label: Optional[str]):
    """gene -> (iso_ids, structures, x0_linear, x1_linear). Mean_* are log2(CPM+1)
    (comparisons.py); inverted to linear CPM. control_label picks which Mean_* is x0."""
    df = pd.read_csv(stats_tsv, sep="\t", dtype=str).fillna("")
    m1, m2 = _detect_mean_cols(list(df.columns))
    if control_label and control_label in m2:
        ctrl_col, dis_col = m2, m1
    else:
        ctrl_col, dis_col = m1, m2     # default: first Mean_ is control
    genes: Dict[str, dict] = defaultdict(lambda: {"iso": [], "struct": [], "x0": [], "x1": []})
    for _, r in df.iterrows():
        feat = r.get("Feature", "")
        if ":" not in feat:
            continue
        gene, iso = feat.split(":", 1)
        toks = [t for t in r.get("Junctions", "").split("|") if t]
        if len(toks) < 2:
            continue
        try:
            x0 = max(0.0, 2.0 ** float(r[ctrl_col]) - 1.0)
            x1 = max(0.0, 2.0 ** float(r[dis_col]) - 1.0)
        except ValueError:
            continue
        g = genes[gene]
        g["iso"].append(iso); g["struct"].append(toks); g["x0"].append(x0); g["x1"].append(x1)
    return genes, ctrl_col, dis_col


def attribute_event_table(event_tsv: str, isoform_stats_tsv: str, exon_db: str,
                          control_label: Optional[str] = None, output_tsv: Optional[str] = None,
                          method: str = "auto") -> pd.DataFrame:
    print(f"[isoform_attribute] event table : {event_tsv}")
    print(f"[isoform_attribute] isoform stats: {isoform_stats_tsv}")
    events = pd.read_csv(event_tsv, sep="\t", dtype=str).fillna("")
    gene_ab, ctrl_col, dis_col = load_isoform_abundances(isoform_stats_tsv, control_label)
    print(f"[isoform_attribute] control={ctrl_col} disease={dis_col}; "
          f"loaded isoforms for {len(gene_ab)} genes")

    genes_needed = set(gene_ab.keys())
    coords, strand = load_exon_coords(exon_db, genes_needed)

    models: Dict[str, GeneModel] = {}

    def model_for(gene):
        if gene in models:
            return models[gene]
        g = gene_ab.get(gene)
        if not g or len(g["iso"]) < 2:
            models[gene] = None
            return None
        models[gene] = GeneModel(gene, strand.get(gene, "+"), g["iso"], g["struct"],
                                 g["x0"], g["x1"], coords)
        return models[gene]

    out_rows = []
    n_attr = n_confound = n_agree = 0
    for _, row in events.iterrows():
        rec = row.to_dict()
        feat = row.get("Feature", "")
        res = {"attr_driver_isoform": "", "attr_driver_mechanism": "", "attr_driver_direction": "",
               "attr_driver_contribution": "", "attr_method": "", "attr_ranked": "",
               "attr_vs_current": ""}
        if "|" in feat:
            j1, j2 = feat.split("|", 1)
            gene = j1.split(":", 1)[0]
            u = _strip_gene(j1.strip())
            partner = _strip_gene(j2.strip())
            gm = model_for(gene)
            if gm is not None:
                a = attribute_event(gm, u, partner, method=method)
                res["attr_driver_isoform"] = a["driver"]
                res["attr_driver_mechanism"] = a["mechanism"]
                res["attr_driver_direction"] = {1: "up", -1: "down", 0: ""}[a["direction"]]
                res["attr_driver_contribution"] = f"{a['contribution']:+.4f}" if a["driver"] else ""
                res["attr_method"] = a["method"]
                res["attr_ranked"] = a["ranked"]
                if a["driver"]:
                    n_attr += 1
                    cur = {_strip_version(_strip_gene(row.get("Isoform_1|Length", "").split("|")[0])),
                           _strip_version(_strip_gene(row.get("Isoform_2|Length", "").split("|")[0]))}
                    drv = _strip_version(a["driver"])
                    if drv in cur:
                        res["attr_vs_current"] = "agrees"; n_agree += 1
                    else:
                        res["attr_vs_current"] = "new_driver"
                    if a["mechanism"] == "competing_isoform_change":
                        n_confound += 1
        rec.update(res)
        out_rows.append(rec)

    result = pd.DataFrame(out_rows)
    if output_tsv is None:
        base, ext = os.path.splitext(event_tsv)
        output_tsv = f"{base}-attributed{ext or '.txt'}"
    result.to_csv(output_tsv, sep="\t", index=False)
    print(f"[isoform_attribute] attributed {n_attr} events; "
          f"{n_agree} agree with current pair, {n_attr - n_agree} name a different driver; "
          f"{n_confound} are competing-isoform (confound) events the current rule cannot name")
    print(f"[isoform_attribute] wrote: {output_tsv}")
    return result


def _find_exon_db(work_dir: str, exon_db: Optional[str]) -> Optional[str]:
    if exon_db and os.path.exists(exon_db):
        return exon_db
    import glob
    for pat in ("*_Ensembl_exon.txt", "*Ensembl_exon.txt"):
        hits = glob.glob(os.path.join(work_dir, pat))
        if hits:
            return hits[0]
    return None


def attribute_differentials(work_dir: str, exon_db: Optional[str] = None,
                            control_label: str = "young") -> List[str]:
    """Workflow hook: for every dPSI event table produced by compute_differentials,
    write a companion '*-attributed.txt' with the count-aware driver attribution.
    Pairs dPSI-{covariate,cluster}-events/<name> with diff-{covariate,cluster}-isoform/<name>.
    Best-effort and non-destructive (originals untouched). Returns written paths."""
    import glob
    exon = _find_exon_db(work_dir, exon_db)
    if exon is None:
        print("[isoform_attribute] no Ensembl_exon.txt found; skipping attribution step")
        return []
    pairs = [("dPSI-covariate-events", "diff-covariate-isoform"),
             ("dPSI-cluster-events", "diff-cluster-isoform")]
    written = []
    for ev_dir, iso_dir in pairs:
        evp, isop = os.path.join(work_dir, ev_dir), os.path.join(work_dir, iso_dir)
        if not (os.path.isdir(evp) and os.path.isdir(isop)):
            continue
        out_dir = os.path.join(work_dir, "attribution", ev_dir)
        os.makedirs(out_dir, exist_ok=True)
        for ev_file in sorted(glob.glob(os.path.join(evp, "*_stats-annotated.txt"))):
            base = os.path.basename(ev_file)
            iso_file = os.path.join(isop, base)
            if not os.path.exists(iso_file):
                continue
            out = os.path.join(out_dir, base.replace("_stats-annotated.txt", "-attributed.txt"))
            try:
                attribute_event_table(ev_file, iso_file, exon, control_label=control_label, output_tsv=out)
                written.append(out)
            except Exception as e:  # never break the workflow
                print(f"[isoform_attribute] skipped {base}: {e}")
    print(f"[isoform_attribute] attribution complete: {len(written)} event tables -> {work_dir}/attribution/")
    return written


def _build_parser():
    p = argparse.ArgumentParser(
        description="Count-aware junction->isoform attribution (alternative to isoform_predict.py).")
    p.add_argument("event_tsv", help="dPSI event table (with Feature, Isoform_1|Length, Isoform_2|Length).")
    p.add_argument("isoform_stats_tsv", help="matching diff-covariate/cluster isoform stats "
                                             "(Feature=gene:iso, Junctions, Mean_<cond> columns).")
    p.add_argument("exon_db", help="Hs_Ensembl_exon.txt (for genomic-overlap cliques).")
    p.add_argument("--control-label", default="young",
                   help="substring of the control Mean_<cond> column (default: young).")
    p.add_argument("--method", default="auto", choices=["auto", "contribution", "em"],
                   help="auto = contribution decomposition, EM only when isoform abundances are sparse.")
    p.add_argument("-o", "--output", default=None, help="output TSV (default <event>-attributed.txt)")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    # attribute_event: "auto" = contribution unless sparse->EM; "em" = force EM;
    # anything else (e.g. "contribution") = force the contribution decomposition.
    attribute_event_table(args.event_tsv, args.isoform_stats_tsv, args.exon_db,
                          control_label=args.control_label, output_tsv=args.output,
                          method=args.method)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
