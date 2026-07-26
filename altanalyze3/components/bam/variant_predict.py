"""Predict somatic-mutation status from RNA expression, and transfer the model across platforms.

Companion to ``variant_impact.py``. Where ``variant_impact`` asks *"given per-cell genotypes, is the
mutation's expression program recoverable within one sample?"*, this module asks the harder,
generalisable question:

    Can a classifier trained on a large **bulk** RNA-seq cohort (BEAT-AML) predict the mutation
    status of **independent** samples -- and does it transfer down to single-cell data?

Contents
--------
data          load_bulk_expression, load_bulk_labels, representation()
models        MODELS registry -- signature score, elastic net, L1/L2 logistic, linear SVM,
              random forest, gradient boosting, nearest shrunken centroid, k-TSP, rank-space
              elastic net, PCA-logistic, clinical-covariate-only baseline
evaluation    nested_cv (group-aware, leak-free feature filtering), permutation_null, metrics
transfer      fit_full, transfer_scores, quantile_align
per-cell      cardelino_posterior (sample-prior shrinkage), numbat_joint (allele + expression),
              comutation_prior (Monopogen-style co-segregation prior)

Design notes that matter for correctness
----------------------------------------
* **Grouping.** BEAT-AML contains multiple specimens per subject (698 RNA-seq / 632 subjects).
  Splitting at the specimen level leaks a patient across train and test, so every split here is
  ``StratifiedGroupKFold`` on ``dbgap_subject_id``.
* **Labels.** A sample is a negative only if it was actually exome-sequenced (``analysisExomeSeq ==
  'y'``) and the gene is absent from its ``variantSummary``. Samples with no exome are *unknown*,
  not wild-type, and are dropped.
* **Leak-free filtering.** Expression filtering, variance ranking, scaling and gene selection all
  happen inside the training fold (sklearn ``Pipeline``), never on the full matrix.
* **Cross-platform.** Bulk polyA whole-marrow vs 3' single-cell pseudobulk differ in scale,
  composition and gene detection. Three representations are supported: within-dataset gene
  z-score (``z``), within-sample gene rank (``rank``), and raw log2-CPM (``log``). Rank-space and
  k-TSP are invariant to any monotone per-sample transform, which is why they are the candidates
  most likely to survive the platform jump.
"""

from __future__ import annotations

import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------------

BEATAML_ID_COLS = ("stable_id", "display_label", "description", "biotype")


def load_bulk_expression(counts_path, min_cpm=1.0, min_frac=0.05):
    """BEAT-AML gene x sample log2(CPM+1), collapsed to HGNC symbol.

    Duplicate symbols are summed (they are alternative Ensembl IDs for the same locus).
    Genes below ``min_cpm`` in more than ``1 - min_frac`` of samples are dropped -- this is a
    dataset-level (not label-aware) filter, so it cannot leak the outcome.
    """
    df = pd.read_csv(counts_path, sep="\t", low_memory=False)
    sym = df["display_label"].astype(str)
    mat = df.drop(columns=[c for c in BEATAML_ID_COLS if c in df.columns])
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat.index = sym
    mat = mat.groupby(level=0).sum()
    cpm = mat.div(mat.sum(axis=0).replace(0, 1), axis=1) * 1e6
    keep = (cpm >= min_cpm).mean(axis=1) >= min_frac
    return np.log2(cpm.loc[keep] + 1.0)


def parse_variant_summary(value):
    """Gene symbols from a BEAT-AML ``variantSummary`` cell.

    Entries are pipe-separated and take two forms -- ``'SRSF2 (p.P95R; MAF 44%)'`` and the
    lesion shorthand ``'FLT3-ITD'``. The leading symbol token normalises both.
    """
    genes = set()
    if isinstance(value, str) and value.strip():
        for entry in value.split("|"):
            m = re.match(r"\s*([A-Za-z][A-Za-z0-9]*)", entry)
            if m:
                genes.add(m.group(1))
    return genes


#: The allele-fraction field is written three different ways in the same column -- ``MAF 44%``,
#: ``VAF 39%`` and ``MAF <0.5%`` (sometimes HTML-escaped as ``&lt;``). A regex that only matches
#: ``MAF\s*([0-9.]+)%`` silently misses every ``VAF`` entry and every ``<`` entry -- i.e. exactly the
#: lowest-fraction calls a subclonality filter exists to remove. Matching all three forms is what
#: makes ``min_maf`` actually do anything.
_MAF_RE = re.compile(r"\b(?:MAF|VAF)\s*(?:&lt;|&gt;|[<>~=])?\s*([0-9.]+)\s*%")

#: Legacy symbols that must be merged before a gene becomes a phenotype. ``MLL`` and ``KMT2A`` are
#: the same locus; leaving them split creates two labels, two entries in the multiple-testing
#: family, and two different AUROCs for one gene.
HGNC_ALIASES = {"MLL": "KMT2A", "MLL2": "KMT2D", "MLL3": "KMT2C", "MLL4": "KMT2B",
                "PAX5": "PAX5", "C-KIT": "KIT", "FLT-3": "FLT3", "TET-2": "TET2"}

#: Substrings that mark a ``variantSummary`` entry as not a somatic point mutation of the gene body.
#: ``p.G105G`` style synonymous changes and structural shorthands are separated out rather than
#: pooled, because they have different (or no) expression consequences.
_SYNONYMOUS_RE = re.compile(r"p\.([A-Z])(\d+)\1\b")


def canonical_gene(symbol):
    return HGNC_ALIASES.get(symbol.upper(), symbol.upper())


def parse_variant_maf(value):
    """{gene: max allele fraction} from a ``variantSummary`` cell (nan where none is quoted)."""
    out = {}
    if not isinstance(value, str):
        return out
    for entry in value.split("|"):
        m = re.match(r"\s*([A-Za-z][A-Za-z0-9]*)", entry)
        if not m:
            continue
        maf = _MAF_RE.search(entry)
        v = float(maf.group(1)) / 100.0 if maf else np.nan
        g = canonical_gene(m.group(1))
        if g not in out or (not np.isnan(v) and (np.isnan(out[g]) or v > out[g])):
            out[g] = v
    return out


def parse_lesions(value, drop_synonymous=True):
    """Lesion-class labels from a ``variantSummary`` cell.

    A gene symbol is not a phenotype. ``FLT3-ITD`` (an internal tandem duplication driving a STAT5 /
    high-blast program) and ``FLT3-TKD`` (a D835 kinase-domain point mutation) have different, partly
    opposing expression correlates, and pooling them with 35 miscellaneous FLT3 coding variants
    produces a label with no coherent phenotype. The same applies to ``KMT2A-PTD``. This returns
    both the gene-level token and the lesion-level token so either can be modelled.
    """
    out = set()
    if not isinstance(value, str) or not value.strip():
        return out
    for entry in value.split("|"):
        m = re.match(r"\s*([A-Za-z][A-Za-z0-9]*)(-[A-Za-z0-9]+)?", entry)
        if not m:
            continue
        gene = canonical_gene(m.group(1))
        suffix = (m.group(2) or "").upper().lstrip("-")
        if drop_synonymous and _SYNONYMOUS_RE.search(entry):
            continue                                  # e.g. 'IDH1 (p.G105G; MAF 50%)'
        if "PTD" in entry.upper():
            out.add(f"{gene}-PTD")
        elif suffix in ("ITD", "TKD"):
            out.add(f"{gene}-{suffix}")
        else:
            out.add(gene)
        out.add(gene)                                 # gene-level token always emitted too
    return out


#: Genes whose somatic mutations are effectively confined to a handful of codons, so that a
#: hotspot panel with adequate depth and no call is genuine evidence of wild-type. Everything else
#: (RUNX1, TP53, TET2, ASXL1, BCOR, EZH2, STAG2, PHF6, CEBPA, NF1, ...) carries dispersed
#: loss-of-function mutations across the whole coding sequence: a hotspot panel that finds nothing
#: has not shown the gene is wild-type, it has only failed to look. Treating those as negatives is
#: the single easiest way to manufacture a falsely good transfer result, so they are graded weak.
HOTSPOT_COMPLETE_GENES = frozenset({
    "NRAS", "KRAS", "IDH1", "IDH2", "NPM1", "JAK2", "KIT", "SF3B1", "SRSF2", "U2AF1",
    "PTPN11", "CSF3R", "MYD88", "SETBP1", "FLT3", "BRAF", "MPL", "CALR",
})

#: Case-mix covariates: everything that determines an AML expression profile *other than* the
#: mutation being predicted. ``consensusAMLFusions`` matters most -- fusion class (CBFB-MYH11,
#: RUNX1-RUNX1T1, PML-RARA, KMT2A rearrangements) is the strongest and most reproducible expression
#: stratification in AML, and it is tightly linked to the mutations under test (KIT with
#: core-binding-factor AML; NRAS/KRAS with inv(16); FLT3 with PML-RARA and NPM1). ``fabBlastMorphology``
#: carries the differentiation axis and the ontogeny flags carry the secondary-AML axis, which is
#: *defined* by mutations in SRSF2/SF3B1/U2AF1/ASXL1/EZH2/BCOR/STAG2 -- i.e. predicting SRSF2 from
#: bulk expression is, to first order, predicting secondary-type AML ontogeny.
CASEMIX_COLS = ("%.Blasts.in.BM", "%.Blasts.in.PB", "ageAtSpecimenAcquisition",
                "consensus_sex", "specimenType", "cohort", "centerID",
                "diseaseStageAtSpecimenCollection", "consensusAMLFusions",
                "fabBlastMorphology", "isDenovo", "isTransformed", "priorMDS")

#: Columns that must NEVER enter a "case-mix" baseline. ELN2017 is *defined by* NPM1, FLT3-ITD
#: allelic ratio, biallelic CEBPA, RUNX1, ASXL1 and TP53 status; using it as a confounder control
#: for predicting those genes hands the baseline the answer, turning the intended floor into a
#: partial oracle and making "expression beats case mix" untestable in either direction.
GENOTYPE_DERIVED_COLS = ("ELN2017", "CEBPA_Biallelic", "FLT3_ITD", "NPM1", "RUNX1",
                         "variantSummary", "riskGroup", "riskCyto")

COVARIATE_COLS = CASEMIX_COLS          # backwards-compatible alias


def load_bulk_labels(clinical_path, sheet="summary", require_exome=True, min_maf=None,
                     drop_empty_summary=True, lesion_level=True):
    """Return ``(labels, covariates, groups)`` for BEAT-AML RNA-seq specimens.

    ``labels``      sample x gene boolean; a cell is True when the gene appears in the specimen's
                    ``variantSummary``.
    ``covariates``  the clinical confounders used for the covariate-only baseline.
    ``groups``      ``dbgap_subject_id`` per sample -- the CV grouping variable.

    ``require_exome`` drops specimens with ``analysisExomeSeq != 'y'``: without DNA sequencing an
    absent gene is *unknown*, not wild-type, and silently calling it a negative inflates every
    downstream metric. ``min_maf`` additionally demotes calls below a variant-allele fraction to
    unlabelled (dropped per gene) -- subclonal lesions cannot be expected to shift bulk expression.
    """
    d = pd.read_excel(clinical_path, sheet_name=sheet)
    d = d[d["dbgap_rnaseq_sample"].notna()].copy()
    d["dbgap_rnaseq_sample"] = d["dbgap_rnaseq_sample"].astype(str)
    if require_exome:
        d = d[d["analysisExomeSeq"].astype(str).str.lower() == "y"]
    d = d.drop_duplicates(subset="dbgap_rnaseq_sample").set_index("dbgap_rnaseq_sample")

    if drop_empty_summary:
        # 206 of the 658 exome-sequenced specimens have a completely empty variantSummary. An AML
        # exome does not yield zero drivers, so these are curation gaps, not wild-type genomes:
        # keeping them makes ~31% of every negative class "nobody filled in the field".
        nonempty = d["variantSummary"].astype(str).str.strip().ne("") & d["variantSummary"].notna()
        d = d[nonempty]

    gene_sets = d["variantSummary"].map(
        parse_lesions if lesion_level else (lambda v: {canonical_gene(g)
                                                       for g in parse_variant_summary(v)}))
    all_genes = sorted(set().union(*gene_sets)) if len(gene_sets) else []
    labels = pd.DataFrame(False, index=d.index, columns=all_genes)
    for s, gs in gene_sets.items():
        if gs:
            labels.loc[s, sorted(gs)] = True

    if "CEBPA_Biallelic" in d.columns:                 # only biallelic CEBPA has the published
        bi = d["CEBPA_Biallelic"].astype(str).str.lower()   # expression signature; monoallelic
        labels["CEBPA-biallelic"] = (bi == "bi").values     # largely does not
    if min_maf is not None:
        mafs = d["variantSummary"].map(parse_variant_maf)
        for s, mm in mafs.items():
            for g, v in mm.items():
                if g in labels.columns and not np.isnan(v) and v < min_maf:
                    labels.loc[s, g] = np.nan          # unlabelled: subclonal
    cov = d[[c for c in CASEMIX_COLS if c in d.columns]].copy()
    groups = d["dbgap_subject_id"].astype(str)
    return labels, cov, groups


def encode_covariates(cov, max_levels=40, min_level_n=5):
    """Numeric design matrix from the clinical covariate frame.

    Three things here matter, because this matrix is the *baseline* that any claim of "expression
    predicts genotype rather than case mix" has to clear -- and every encoding shortcut weakens the
    baseline, which inflates the apparent incremental value of expression:

    * **Missingness is a level, not a zero.** ``consensusAMLFusions`` is blank for most specimens and
      ``fabBlastMorphology`` for many. Dropping NaN before one-hot encoding makes "no fusion",
      "not assessed" and "blank" the same all-zero row. Each column therefore gets an explicit
      ``=MISSING`` indicator, and numeric columns get a ``_isna`` flag alongside median imputation.
    * **Low-cardinality integer codes are categorical.** ``centerID`` is a nominal site label that
      happens to be written 1-8; treating it as a continuous predictor imposes a meaningless
      ordering and amounts to no site adjustment at all.
    * **No alphabetical truncation.** Keeping the first 12 levels by sort order silently discards
      whichever fusion classes sort late.
    """
    out = {}
    for c in cov.columns:
        s = cov[c]
        num = pd.to_numeric(s, errors="coerce")
        nuniq = num.dropna().nunique()
        is_numeric = num.notna().mean() > 0.5 and not (nuniq <= 12 and
                                                       float(num.dropna().mod(1).abs().max() or 0) == 0)
        if is_numeric:
            out[c] = num.fillna(num.median())
            if num.isna().any():
                out[f"{c}_isna"] = num.isna().astype(float)
        else:
            lab = s.astype(str).where(s.notna() & (s.astype(str).str.strip() != ""), "MISSING")
            counts = lab.value_counts()
            for lev in counts[counts >= min_level_n].index[:max_levels]:
                out[f"{c}={lev}"] = (lab == lev).astype(float)
    return pd.DataFrame(out, index=cov.index)


# --------------------------------------------------------------------------------------------
# representations
# --------------------------------------------------------------------------------------------

def representation(expr, kind="z"):
    """Transform a gene x sample log2-CPM matrix into the space a model expects.

    ``z``     gene-wise z-score **within this dataset** -- removes the platform's per-gene offset,
              which is what makes a bulk-trained coefficient vector applicable to pseudobulk.
    ``rank``  within-sample rank of each gene, scaled to [0, 1] -- invariant to any monotone
              per-sample transform (library size, normalisation, 3'-bias), so it survives the
              bulk -> single-cell jump when z-scoring does not.
    ``log``   unchanged log2-CPM.
    """
    if kind == "log":
        return expr.copy()
    if kind == "z":
        mu = expr.mean(axis=1)
        sd = expr.std(axis=1).replace(0, np.nan)
        z = expr.sub(mu, axis=0).div(sd, axis=0)
        return z.fillna(0.0)
    if kind == "rank":
        r = expr.rank(axis=0, method="average")
        return r.div(float(len(expr.index)))
    raise ValueError(f"unknown representation {kind}")


def quantile_align(target, reference):
    """Map each target sample's expression onto the reference cohort's mean quantile profile.

    Classic cross-platform quantile normalisation: the two datasets are forced to share a marginal
    distribution so that a coefficient vector learnt on one is on-scale for the other. Genes are
    intersected first by the caller.
    """
    ref_profile = np.sort(reference.values, axis=0).mean(axis=1)
    out = np.empty_like(target.values, dtype=float)
    for j in range(target.shape[1]):
        order = np.argsort(np.argsort(target.values[:, j]))
        out[:, j] = np.interp(order / max(len(order) - 1, 1),
                              np.linspace(0, 1, len(ref_profile)), ref_profile)
    return pd.DataFrame(out, index=target.index, columns=target.columns)


# --------------------------------------------------------------------------------------------
# custom estimators
# --------------------------------------------------------------------------------------------

def top_variable_genes(expr, k=4000, min_k=200):
    """The ``k`` most variable genes of a gene x sample **log2-CPM** matrix.

    Gene selection must happen on the untransformed expression, BEFORE any per-gene
    standardisation. Z-scoring sets every gene's variance to 1, so a variance filter applied after
    it has nothing to rank on and silently returns an arbitrary tie-broken subset -- in this dataset
    that dropped XIST, every Y-chromosome gene, HOXA9 and MEIS1 while keeping lincRNAs. The filter
    is unsupervised (labels are never consulted), so applying it at the dataset level does not leak
    the outcome, and applying the identical gene set to both cohorts is what makes a transferred
    coefficient vector meaningful.
    """
    v = expr.var(axis=1)
    return list(v.nlargest(max(min(k, len(v)), min_k)).index)


class VarianceTopK(BaseEstimator):
    """Keep the K most variable features, fitted on the training fold only.

    Refuses to run on a matrix whose feature variances are all equal (the post-z-score case), rather
    than returning an arbitrary subset. Callers must select genes on log2-CPM first, via
    ``top_variable_genes``.
    """

    def __init__(self, k=5000):
        self.k = k

    def fit(self, X, y=None):
        v = np.asarray(X).var(axis=0)
        if X.shape[1] > 1 and np.ptp(v) < 1e-8 * max(np.median(v), 1e-12) and self.k < X.shape[1]:
            raise ValueError(
                "VarianceTopK: all feature variances are identical (input looks z-scored). "
                "Select genes on log2-CPM with top_variable_genes() before standardising.")
        k = min(self.k, X.shape[1])
        self.idx_ = np.argsort(v)[::-1][:k]
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.idx_]


class SignatureScore(BaseEstimator, ClassifierMixin):
    """The naive baseline: differential genes -> mean(up) - mean(down).

    This is the method already applied to these data (MarkerFinder signature scored in bulk). It is
    included so every learned model is judged against what a marker-gene signature alone achieves,
    not against chance.
    """

    def __init__(self, top_n=50):
        self.top_n = top_n

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y).astype(bool)
        t, _ = stats.ttest_ind(X[y], X[~y], axis=0, equal_var=False)
        t = np.nan_to_num(t)
        n = min(self.top_n, X.shape[1] // 2)
        self.up_ = np.argsort(t)[::-1][:n]
        self.dn_ = np.argsort(t)[:n]
        self.mu_ = X.mean(axis=0)
        self.sd_ = X.std(axis=0)
        self.sd_[self.sd_ == 0] = 1.0
        self.classes_ = np.array([False, True])
        return self

    def decision_function(self, X):
        z = (np.asarray(X, dtype=float) - self.mu_) / self.sd_
        return z[:, self.up_].mean(axis=1) - z[:, self.dn_].mean(axis=1)

    def predict(self, X):
        return self.decision_function(X) > 0


class NearestShrunkenCentroid(BaseEstimator, ClassifierMixin):
    """PAM (Tibshirani 2002) with a continuous discriminant, so it can be scored by AUROC.

    sklearn's ``NearestCentroid`` implements the shrinkage but exposes only hard labels; the
    difference of standardised centroid distances is the natural continuous analogue.
    """

    def __init__(self, shrink=2.0):
        self.shrink = shrink

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y).astype(bool)
        self.grand_ = X.mean(axis=0)
        s = X.std(axis=0)
        self.s_ = s + np.median(s[s > 0]) if (s > 0).any() else s + 1.0
        cents = {}
        for lab, mask in ((True, y), (False, ~y)):
            n = max(int(mask.sum()), 1)
            m = 1.0 / np.sqrt(n) - 1.0 / np.sqrt(len(y))
            m = max(m, 1e-6)
            d = (X[mask].mean(axis=0) - self.grand_) / (m * self.s_)
            d = np.sign(d) * np.maximum(np.abs(d) - self.shrink, 0.0)   # soft-threshold
            cents[lab] = self.grand_ + m * self.s_ * d
        self.cents_ = cents
        self.classes_ = np.array([False, True])
        return self

    def decision_function(self, X):
        X = np.asarray(X, dtype=float)
        d1 = (((X - self.cents_[True]) / self.s_) ** 2).sum(axis=1)
        d0 = (((X - self.cents_[False]) / self.s_) ** 2).sum(axis=1)
        return d0 - d1

    def predict(self, X):
        return self.decision_function(X) > 0


class KTSP(BaseEstimator, ClassifierMixin):
    """k Top-Scoring Pairs (Geman 2004; Tan 2005) -- decisions are *within-sample gene comparisons*.

    Every rule has the form ``expr[i] < expr[j]``, so the classifier is invariant to library size,
    normalisation and any monotone per-sample transform. That is precisely the failure mode when
    carrying a bulk-trained model onto single-cell pseudobulk, which is why this method is a
    first-class candidate here rather than a curiosity.
    """

    def __init__(self, k=9, n_pre=200):
        self.k = k
        self.n_pre = n_pre

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y).astype(bool)
        t, _ = stats.ttest_ind(X[y], X[~y], axis=0, equal_var=False)
        pre = np.argsort(np.abs(np.nan_to_num(t)))[::-1][:min(self.n_pre, X.shape[1])]
        A = X[:, pre]
        # P(x_i < x_j) within each class, for all pairs, vectorised over samples
        lt1 = (A[y][:, :, None] < A[y][:, None, :]).mean(axis=0)
        lt0 = (A[~y][:, :, None] < A[~y][:, None, :]).mean(axis=0)
        delta = lt1 - lt0
        np.fill_diagonal(delta, 0.0)
        order = np.argsort(np.abs(delta).ravel())[::-1]
        used, pairs = set(), []
        for flat in order:
            i, j = divmod(int(flat), delta.shape[1])
            if i in used or j in used or i == j:
                continue
            used.update((i, j))
            pairs.append((pre[i], pre[j], 1.0 if delta[i, j] > 0 else -1.0))
            if len(pairs) >= self.k:
                break
        self.pairs_ = pairs
        self.classes_ = np.array([False, True])
        return self

    def decision_function(self, X):
        X = np.asarray(X, dtype=float)
        if not self.pairs_:
            return np.zeros(X.shape[0])
        votes = np.zeros(X.shape[0])
        for i, j, sgn in self.pairs_:
            votes += sgn * ((X[:, i] < X[:, j]).astype(float) * 2 - 1)
        return votes / len(self.pairs_)

    def predict(self, X):
        return self.decision_function(X) > 0


# --------------------------------------------------------------------------------------------
# model registry
# --------------------------------------------------------------------------------------------

def _pipe(*steps):
    return Pipeline(list(steps))


def build_model(name, k_var=5000, seed=0, grid="full"):
    """(estimator, param_grid, representation) for a registry name.

    ``grid='fast'`` shrinks every hyperparameter grid to a single sensible value and loosens the
    saga tolerance. It exists for the 27-gene x 11-model screening pass, where the question is
    *which model family is worth pursuing*, not *what is its best possible AUROC*; the shortlisted
    models are then re-run with the full grid. Screening on a reduced grid biases all families in
    the same direction, so the ranking it produces is still usable.
    """
    fast = grid == "fast"
    tol = 1e-3 if fast else 1e-4
    filt = ("filter", VarianceTopK(k=k_var))
    sc = ("scale", StandardScaler())
    n_iter = 800 if fast else 3000
    spec = {
        "signature": (_pipe(filt, ("clf", SignatureScore())),
                      {"clf__top_n": [25, 50, 100]}, "z"),
        "enet": (_pipe(filt, sc, ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                                            max_iter=n_iter, tol=tol,
                                                            l1_ratio=0.5))),
                 {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.15, 0.5, 0.85]}, "z"),
        "lasso": (_pipe(filt, sc, ("clf", LogisticRegression(penalty="l1", solver="liblinear",
                                                             max_iter=n_iter, tol=tol))),
                  {"clf__C": [0.01, 0.1, 1.0]}, "z"),
        "ridge": (_pipe(filt, sc, ("clf", LogisticRegression(penalty="l2", solver="lbfgs",
                                                             max_iter=n_iter, tol=tol))),
                  {"clf__C": [0.001, 0.01, 0.1, 1.0]}, "z"),
        "svm": (_pipe(filt, sc, ("clf", LinearSVC(max_iter=5000, dual=True, tol=tol))),
                {"clf__C": [0.001, 0.01, 0.1]}, "z"),
        "rf": (_pipe(filt, ("clf", RandomForestClassifier(n_estimators=300 if fast else 500,
                                                          n_jobs=1, random_state=seed,
                                                          class_weight="balanced"))),
               {"clf__max_features": ["sqrt", 0.05]}, "z"),
        "gbm": (_pipe(("filter", VarianceTopK(k=min(k_var, 2000))),
                      ("clf", HistGradientBoostingClassifier(random_state=seed,
                                                             max_iter=100 if fast else 200))),
                {"clf__learning_rate": [0.05, 0.15], "clf__max_leaf_nodes": [7, 31]}, "z"),
        "nsc": (_pipe(filt, ("clf", NearestShrunkenCentroid())),
                {"clf__shrink": [0.5, 1.0, 2.0, 4.0]}, "z"),
        "ktsp": (_pipe(filt, ("clf", KTSP())),
                 {"clf__k": [5, 9, 15], "clf__n_pre": [100, 300]}, "rank"),
        "enet_rank": (_pipe(filt, sc, ("clf", LogisticRegression(penalty="elasticnet",
                                                                 solver="saga", max_iter=n_iter,
                                                                 tol=tol, l1_ratio=0.5))),
                      {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.15, 0.5, 0.85]}, "rank"),
        "pca_logit": (_pipe(filt, sc, ("pca", PCA(n_components=40, random_state=seed)),
                            ("clf", LogisticRegression(max_iter=n_iter, tol=tol))),
                      {"clf__C": [0.01, 0.1, 1.0]}, "z"),
    }
    if name not in spec:
        raise ValueError(f"unknown model {name}")
    est, full_grid, rep = spec[name]
    if not fast:
        return est, full_grid, rep
    fast_grid = {
        "signature": {"clf__top_n": [50]}, "enet": {"clf__C": [0.1], "clf__l1_ratio": [0.5]},
        "lasso": {"clf__C": [0.1]}, "ridge": {"clf__C": [0.01]}, "svm": {"clf__C": [0.01]},
        "rf": {"clf__max_features": ["sqrt"]}, "gbm": {"clf__learning_rate": [0.05]},
        "nsc": {"clf__shrink": [1.0, 2.0]}, "ktsp": {"clf__k": [9]},
        "enet_rank": {"clf__C": [0.1], "clf__l1_ratio": [0.5]}, "pca_logit": {"clf__C": [0.1]},
    }
    return est, fast_grid[name], rep


MODELS = ("signature", "enet", "lasso", "ridge", "svm", "rf", "gbm", "nsc", "ktsp",
          "enet_rank", "pca_logit")

#: representations that are invariant to per-sample monotone transforms, i.e. the ones expected to
#: survive a platform change without recalibration.
PLATFORM_ROBUST = ("ktsp", "enet_rank")


# --------------------------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------------------------

def _scores(est, X):
    if hasattr(est, "decision_function"):
        return np.asarray(est.decision_function(X)).ravel()
    return np.asarray(est.predict_proba(X))[:, 1]


def metrics(y, s):
    y = np.asarray(y).astype(bool)
    if y.sum() == 0 or (~y).sum() == 0:
        return dict(auroc=np.nan, auprc=np.nan, auprc_lift=np.nan, n=len(y), n_pos=int(y.sum()))
    prev = y.mean()
    ap = average_precision_score(y, s)
    return dict(auroc=float(roc_auc_score(y, s)), auprc=float(ap),
                auprc_lift=float(ap / prev), n=int(len(y)), n_pos=int(y.sum()))


def nested_cv(X, y, groups, model_name, n_splits=5, n_repeats=3, inner_splits=3,
              k_var=5000, seed=0, n_jobs=1, fixed_params=None, return_params=False,
              grid="full"):
    """Group-aware repeated nested CV. Returns ``(metrics_dict, oof_scores)``.

    Hyperparameters are chosen inside each outer training fold (inner ``StratifiedGroupKFold``),
    so the reported AUROC contains no model-selection optimism. Out-of-fold scores are averaged
    across repeats and are what downstream calibration/threshold analysis uses.

    ``fixed_params`` skips the inner search and uses the given hyperparameters directly. This
    exists for the permutation null, where re-running model selection inside every one of hundreds
    of permutations is not affordable; the selected-on-real-data values are frozen and applied to
    both the observed and the permuted runs so the comparison stays like-for-like.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(bool)
    groups = np.asarray(groups)
    est, param_grid, _rep = build_model(model_name, k_var=k_var, seed=seed, grid=grid)
    oof = np.zeros((n_repeats, len(y)))
    per_fold, chosen = [], []
    for r in range(n_repeats):
        outer = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed + r)
        for tr, te in outer.split(X, y, groups):
            if y[tr].sum() < 3 or (~y[tr]).sum() < 3:
                oof[r, te] = np.nan
                continue
            if fixed_params is not None:
                best = clone(est).set_params(**fixed_params).fit(X[tr], y[tr])
            else:
                inner = StratifiedGroupKFold(n_splits=inner_splits, shuffle=True,
                                             random_state=seed + r)
                gs = GridSearchCV(clone(est), param_grid, scoring="roc_auc",
                                  cv=inner.split(X[tr], y[tr], groups[tr]), n_jobs=n_jobs,
                                  refit=True, error_score=np.nan)
                gs.fit(X[tr], y[tr])
                best = gs.best_estimator_
                chosen.append(gs.best_params_)
            s = _scores(best, X[te])
            oof[r, te] = s
            per_fold.append(metrics(y[te], s))
    mean_oof = np.nanmean(oof, axis=0)
    ok = ~np.isnan(mean_oof)
    out = metrics(y[ok], mean_oof[ok])
    out["auroc_fold_mean"] = float(np.nanmean([f["auroc"] for f in per_fold])) if per_fold else np.nan
    out["auroc_fold_sd"] = float(np.nanstd([f["auroc"] for f in per_fold])) if per_fold else np.nan
    out["n_folds"] = len(per_fold)
    if return_params:
        return out, mean_oof, chosen
    return out, mean_oof


def make_folds(y, groups, n_splits=5, seed=0):
    """Fold assignment shared by every feature set, so scores can be compared pairwise.

    The incremental-value question (does expression add anything beyond case mix?) requires the
    competing models to be evaluated on identical held-out samples. Generating folds separately per
    model would leave the difference confounded with the split.
    """
    y = np.asarray(y).astype(bool)
    groups = np.asarray(groups)
    skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(skf.split(np.zeros((len(y), 1)), y, groups))


def oof_from_folds(X, y, groups, folds, model_name, k_var=5000, seed=0, inner_splits=3,
                   grid="full", n_jobs=1):
    """Out-of-fold scores on a fixed fold assignment, with hyperparameters chosen inside each fold."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(bool)
    groups = np.asarray(groups)
    est, param_grid, _ = build_model(model_name, k_var=k_var, seed=seed, grid=grid)
    out = np.full(len(y), np.nan)
    for tr, te in folds:
        if y[tr].sum() < 3 or (~y[tr]).sum() < 3:
            continue
        inner = StratifiedGroupKFold(n_splits=inner_splits, shuffle=True, random_state=seed)
        gs = GridSearchCV(clone(est), param_grid, scoring="roc_auc",
                          cv=inner.split(X[tr], y[tr], groups[tr]), n_jobs=n_jobs,
                          refit=True, error_score=np.nan)
        gs.fit(X[tr], y[tr])
        s = _scores(gs.best_estimator_, X[te])
        out[te] = (stats.rankdata(s) / len(s))     # rank-normalise within fold before pooling:
    return out                                     # raw decision values differ in scale by fold


def nested_method_selection(X_by_rep, y, groups, methods, n_splits=5, inner_splits=3,
                            k_var=5000, seed=0, grid="full", n_jobs=1):
    """Choose the model family *inside* each outer fold, then score the held-out fold.

    Taking the best of N methods per gene and comparing it to a single-method null is
    anti-conservative: the reported statistic is ``max over N correlated methods`` while the null is
    the distribution of one method. Nesting the choice makes the outer-fold score an unbiased
    estimate of the whole "pick a method, then fit" procedure -- which is also what actually gets
    transferred to an external cohort. ``X_by_rep`` is {'z': array, 'rank': array}.
    """
    y = np.asarray(y).astype(bool)
    groups = np.asarray(groups)
    folds = make_folds(y, groups, n_splits=n_splits, seed=seed)
    oof = np.full(len(y), np.nan)
    picked = []
    for tr, te in folds:
        if y[tr].sum() < 3 or (~y[tr]).sum() < 3:
            continue
        inner = list(StratifiedGroupKFold(n_splits=inner_splits, shuffle=True,
                                          random_state=seed).split(np.zeros((len(tr), 1)),
                                                                   y[tr], groups[tr]))
        best, best_auc = None, -np.inf
        for meth in methods:
            X = X_by_rep[build_model(meth)[2]]
            s_in = np.full(len(tr), np.nan)
            for itr, ite in inner:
                est, pg, _ = build_model(meth, k_var=k_var, seed=seed, grid=grid)
                try:
                    est.set_params(**{k: v[0] for k, v in pg.items()}).fit(X[tr][itr], y[tr][itr])
                    s_in[ite] = stats.rankdata(_scores(est, X[tr][ite])) / len(ite)
                except Exception:
                    pass
            ok = ~np.isnan(s_in)
            if ok.sum() < 5 or y[tr][ok].sum() < 2 or (~y[tr][ok]).sum() < 2:
                continue
            auc = roc_auc_score(y[tr][ok], s_in[ok])
            if auc > best_auc:
                best, best_auc = meth, auc
        if best is None:
            continue
        picked.append(best)
        X = X_by_rep[build_model(best)[2]]
        est, pg, _ = build_model(best, k_var=k_var, seed=seed, grid=grid)
        gs = GridSearchCV(clone(est), pg, scoring="roc_auc",
                          cv=StratifiedGroupKFold(n_splits=inner_splits, shuffle=True,
                                                  random_state=seed).split(X[tr], y[tr], groups[tr]),
                          n_jobs=n_jobs, refit=True, error_score=np.nan)
        gs.fit(X[tr], y[tr])
        s = _scores(gs.best_estimator_, X[te])
        oof[te] = stats.rankdata(s) / len(s)
    ok = ~np.isnan(oof)
    m = metrics(y[ok], oof[ok])
    m["methods_picked"] = ",".join(picked)
    return m, oof


def permutation_null(X, y, groups, model_name, n_perm=200, seed=0, obs_auroc=None, **kw):
    """Empirical p for the pooled out-of-fold AUROC under permuted labels.

    Labels are permuted **within the group structure** (subject labels shuffled, then broadcast),
    which preserves the multi-specimen-per-patient dependence -- shuffling specimens independently
    would create an easier null and an anti-conservative p.
    """
    rng = np.random.default_rng(seed)
    if obs_auroc is None:
        obs, _ = nested_cv(X, y, groups, model_name, seed=seed, n_repeats=1, **kw)
    else:
        obs = {"auroc": float(obs_auroc)}
    y = np.asarray(y).astype(bool)
    groups = np.asarray(groups)
    uniq = np.unique(groups)
    gy = pd.Series(y).groupby(groups).max().reindex(uniq).values
    null = []
    for i in range(n_perm):
        perm = rng.permutation(gy)
        ymap = dict(zip(uniq, perm))
        yp = np.array([ymap[g] for g in groups])
        if yp.sum() < 5 or (~yp).sum() < 5:
            continue
        m, _ = nested_cv(X, yp, groups, model_name, seed=seed + 1000 + i, n_repeats=1, **kw)
        null.append(m["auroc"])
    null = np.asarray([v for v in null if not np.isnan(v)])
    p = float((1 + (null >= obs["auroc"]).sum()) / (1 + len(null))) if len(null) else np.nan
    return dict(auroc=obs["auroc"], p_perm=p, null_median=float(np.median(null)) if len(null) else np.nan,
                null_p95=float(np.percentile(null, 95)) if len(null) else np.nan, n_perm=len(null))


def bootstrap_ci(y, s, n_boot=2000, seed=0, groups=None):
    """Percentile CI for AUROC; resamples *groups* when given, so patients stay intact."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    vals = []
    if groups is None:
        idx_pool = [np.array([i]) for i in range(len(y))]
    else:
        groups = np.asarray(groups)
        idx_pool = [np.where(groups == g)[0] for g in np.unique(groups)]
    for _ in range(n_boot):
        pick = rng.integers(0, len(idx_pool), len(idx_pool))
        idx = np.concatenate([idx_pool[i] for i in pick])
        if y[idx].sum() == 0 or (~y[idx]).sum() == 0:
            continue
        vals.append(roc_auc_score(y[idx], s[idx]))
    if not vals:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def stratified_auroc(y, s, strata):
    """AUROC computed within each stratum then sample-size weighted.

    Used to ask whether an apparently predictive signature is really tracking a confounder
    (blast fraction, specimen site): if within-stratum AUROC collapses to chance while the pooled
    value is high, the confounder -- not genotype -- is what the model learnt.
    """
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=float)
    strata = np.asarray(strata)
    num = den = 0.0
    parts = {}
    for lev in pd.unique(strata):
        m = strata == lev
        if y[m].sum() < 3 or (~y[m]).sum() < 3:
            continue
        a = roc_auc_score(y[m], s[m])
        parts[str(lev)] = (float(a), int(m.sum()))
        num += a * m.sum()
        den += m.sum()
    return (float(num / den) if den else np.nan), parts


def delta_auroc_ci(y, base, full, n_boot=2000, seed=0, groups=None):
    """Paired bootstrap CI and p for AUROC(full) - AUROC(base) on the SAME samples.

    Printing two AUROCs side by side is not evidence of incremental value. The question "does
    expression add anything beyond case mix and cell composition?" is a statement about the
    *difference*, and the two scores are computed on identical samples, so the bootstrap must be
    paired -- resampling them independently would inflate the interval and hide a real effect (or
    manufacture one).
    """
    y = np.asarray(y).astype(bool)
    base = np.asarray(base, dtype=float)
    full = np.asarray(full, dtype=float)
    obs = roc_auc_score(y, full) - roc_auc_score(y, base)
    rng = np.random.default_rng(seed)
    pool = ([np.array([i]) for i in range(len(y))] if groups is None
            else [np.where(np.asarray(groups) == g)[0] for g in np.unique(groups)])
    d = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(pool), len(pool))
        idx = np.concatenate([pool[i] for i in pick])
        if y[idx].sum() == 0 or (~y[idx]).sum() == 0:
            continue
        d.append(roc_auc_score(y[idx], full[idx]) - roc_auc_score(y[idx], base[idx]))
    d = np.asarray(d)
    if not len(d):
        return dict(delta=obs, lo=np.nan, hi=np.nan, p_two_sided=np.nan)
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return dict(delta=float(obs), lo=float(np.percentile(d, 2.5)),
                hi=float(np.percentile(d, 97.5)), p_two_sided=float(min(p, 1.0)))


def nnls_deconvolve(bulk_log2cpm, ref_log2cpm, n_markers=50, min_genes=200):
    """Non-negative least-squares cell-type fractions for each bulk sample.

    ``ref_log2cpm`` is a gene x cell-state signature matrix (here, cell-state pseudobulk centroids
    from the matched single-cell atlas). Both matrices are returned to linear CPM before fitting,
    because a mixture is linear in expression, not in its logarithm -- deconvolving log values is a
    common and silent error. Markers are the top ``n_markers`` genes per state by the ratio of that
    state's expression to the maximum of the others.

    This exists to answer the question the whole project turns on: how much of a per-gene mutation
    AUROC is really cell-composition prediction?
    """
    from scipy.optimize import nnls
    shared = [g for g in ref_log2cpm.index if g in set(bulk_log2cpm.index)]
    if len(shared) < min_genes:
        raise ValueError(f"only {len(shared)} shared genes for deconvolution")
    R = (2.0 ** ref_log2cpm.loc[shared]) - 1.0
    B = (2.0 ** bulk_log2cpm.loc[shared]) - 1.0
    marks = set()
    for st in R.columns:
        others = R.drop(columns=[st]).max(axis=1)
        score = (R[st] + 1.0) / (others + 1.0)
        marks.update(score.nlargest(n_markers).index)
    marks = sorted(marks)
    Rm, Bm = R.loc[marks].values, B.loc[marks].values
    Rm = Rm / np.maximum(Rm.sum(axis=0, keepdims=True), 1e-9)
    frac = np.zeros((B.shape[1], R.shape[1]))
    for j in range(B.shape[1]):
        b = Bm[:, j] / max(Bm[:, j].sum(), 1e-9)
        w, _ = nnls(Rm, b)
        frac[j] = w / max(w.sum(), 1e-9)
    return pd.DataFrame(frac, index=bulk_log2cpm.columns, columns=R.columns)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.clip(ranked, 0, 1)
    q[ok] = tmp
    return q


# --------------------------------------------------------------------------------------------
# transfer
# --------------------------------------------------------------------------------------------

def fit_full(X, y, groups, model_name, k_var=5000, seed=0, inner_splits=5, n_jobs=1,
             grid="full"):
    """Refit on the entire cohort with hyperparameters chosen by group-aware inner CV."""
    est, param_grid, _ = build_model(model_name, k_var=k_var, seed=seed, grid=grid)
    inner = StratifiedGroupKFold(n_splits=inner_splits, shuffle=True, random_state=seed)
    gs = GridSearchCV(clone(est), param_grid, scoring="roc_auc",
                      cv=inner.split(np.asarray(X), np.asarray(y).astype(bool), np.asarray(groups)),
                      n_jobs=n_jobs, refit=True, error_score=np.nan)
    gs.fit(np.asarray(X, dtype=float), np.asarray(y).astype(bool))
    return gs.best_estimator_, gs.best_params_


def transfer_scores(fitted, X_target):
    return _scores(fitted, np.asarray(X_target, dtype=float))


def align_matrices(source_expr, target_expr, min_genes=1000):
    """Intersect two gene x sample log2-CPM matrices on shared symbols, preserving gene order."""
    shared = [g for g in source_expr.index if g in set(target_expr.index)]
    if len(shared) < min_genes:
        raise ValueError(f"only {len(shared)} shared genes")
    return source_expr.loc[shared], target_expr.loc[shared]


# --------------------------------------------------------------------------------------------
# single-cell layer: adaptations of published clone-assignment methods
# --------------------------------------------------------------------------------------------
#
# None of cardelino / Monopogen / SComatic / Numbat predicts genotype from expression alone -- they
# all read the variant allele directly. What they contribute is *how to combine weak evidence*:
#   cardelino  per-cell posterior = sparse variant reads x a clone-structure prior
#   Numbat     joint allele + expression likelihood
#   Monopogen  co-segregation / linkage across cells to rescue low-coverage sites
# The three functions below port those combination rules to an expression-only score, so a bulk
# classifier can be evaluated against -- and fused with -- the allelic evidence.

def _logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def calibrate_scores(scores, y, scores_target):
    """Map raw decision values onto probabilities with a Platt fit on the source cohort."""
    lr = LogisticRegression(max_iter=1000)
    lr.fit(np.asarray(scores).reshape(-1, 1), np.asarray(y).astype(int))
    return lr.predict_proba(np.asarray(scores_target).reshape(-1, 1))[:, 1]


def cardelino_posterior(p_expr, sample_prior, weight=1.0):
    """cardelino-style shrinkage of a noisy per-cell score toward the sample's clone prior.

    cardelino stabilises sparse per-cell variant evidence with a prior from the bulk-derived clone
    tree. The analogue here: ``sample_prior`` is the mutation's cell fraction in that sample (or the
    bulk classifier's sample-level probability), and each cell's expression-derived probability is
    combined with it on the log-odds scale. ``weight`` scales how much the per-cell evidence is
    allowed to move the prior; weight=0 returns the prior (no per-cell information).
    """
    p_expr = np.clip(np.asarray(p_expr, dtype=float), 1e-6, 1 - 1e-6)
    prior = float(np.clip(sample_prior, 1e-6, 1 - 1e-6))
    lo_prior = np.log(prior / (1 - prior))
    lr_cell = np.log(p_expr / (1 - p_expr)) - np.log(0.5 / 0.5)
    return _logistic(lo_prior + weight * lr_cell)


def numbat_joint(p_expr, alt_reads, ref_reads, err=0.01, het_vaf=0.5, weight=1.0):
    """Numbat-style joint allele + expression posterior for one variant, per cell.

    The allelic likelihood is binomial: P(alt | MUT) uses ``het_vaf`` (or 1.0 for a hemizygous /
    LOH locus), P(alt | WT) uses the sequencing error rate ``err``. Cells with zero coverage get a
    flat allelic likelihood, so their posterior is driven entirely by expression -- which is the
    practical question this addresses: *does expression rescue the cells long-read coverage misses?*
    """
    p_expr = np.clip(np.asarray(p_expr, dtype=float), 1e-6, 1 - 1e-6)
    a = np.asarray(alt_reads, dtype=float)
    r = np.asarray(ref_reads, dtype=float)
    n = a + r
    ll_mut = a * np.log(het_vaf) + r * np.log(1 - het_vaf)
    ll_wt = a * np.log(err) + r * np.log(1 - err)
    ll_mut[n == 0] = 0.0
    ll_wt[n == 0] = 0.0
    lo = (ll_mut - ll_wt) + weight * (np.log(p_expr / (1 - p_expr)))
    return _logistic(lo)


def comutation_prior(labels, gene, min_pairs=10):
    """Monopogen-style co-segregation prior: log-odds of ``gene`` given each co-mutated partner.

    Monopogen exploits linkage between sites to genotype loci that individually lack coverage. The
    expression-space analogue is the cohort's co-mutation structure (SRSF2-ASXL1, NPM1-FLT3, ...):
    knowing one call shifts the prior on another. Returns {partner: log-odds contribution}.
    """
    y = labels[gene].astype(bool)
    out = {}
    for other in labels.columns:
        if other == gene:
            continue
        x = labels[other].astype(bool)
        if x.sum() < min_pairs:
            continue
        a = int((x & y).sum()) + 0.5
        b = int((x & ~y).sum()) + 0.5
        c = int((~x & y).sum()) + 0.5
        d = int((~x & ~y).sum()) + 0.5
        out[other] = float(np.log((a * d) / (b * c)))
    return out


# --------------------------------------------------------------------------------------------
# pseudobulk helpers (single-cell side)
# --------------------------------------------------------------------------------------------

def pseudobulk_from_adata(adata, group_keys, layer="counts", min_cells=10):
    """Sum raw counts within each group -> gene x group log2(CPM+1) plus a cell-count table."""
    import scipy.sparse as sp
    obs = adata.obs
    key = obs[list(group_keys)].astype(str).agg("|".join, axis=1) if len(group_keys) > 1 \
        else obs[group_keys[0]].astype(str)
    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    cols, counts = [], []
    names = []
    for g, idx in pd.Series(np.arange(len(key)), index=key.values).groupby(level=0):
        ii = idx.values
        if len(ii) < min_cells:
            continue
        v = X[ii].sum(axis=0)
        v = np.asarray(v).ravel() if sp.issparse(X) else np.asarray(v).ravel()
        cols.append(v)
        counts.append(len(ii))
        names.append(g)
    if not cols:
        return pd.DataFrame(), pd.Series(dtype=int)
    mat = pd.DataFrame(np.vstack(cols).T, index=adata.var_names.astype(str), columns=names)
    cpm = mat.div(mat.sum(axis=0).replace(0, 1), axis=1) * 1e6
    return np.log2(cpm + 1.0), pd.Series(counts, index=names)
