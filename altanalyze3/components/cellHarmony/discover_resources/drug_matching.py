"""One drug-matching engine: one scored row per source perturbation signature.

WHY THIS REPLACES THE OLD SCORER
Nathan, 2026-09-07: "I want DART or DART like, I don't want to average any drug effects
(each drug transcriptome measurments are important, not consesnsus which is misleading)."

The store this replaces held 21,299 rows that were each a CMap MODZ consensus over many
LINCS signatures -- 312,438 source measurements collapsed into compound averages. A
compound tested in 54 cell lines at 10 doses became one row whose dose read "10 doses".
That is the object the instruction rejects. This engine scores every source signature on
its own and never combines effect vectors across signatures.

  result_unit                  source_signature
  local_effect_aggregation     none
  local_score_aggregation      none
  compound_representative      none

WHAT ELSE CHANGED, AND WHY EACH ONE IS A DEFECT AND NOT A PREFERENCE

  A q of exactly 0 was read as 1.0, because `(r.get("q") or 1.0)` treats 0.0 as absent,
  so the most significant rows were counted as non-significant.

  `opposing` was the first N rows of a list sorted by rho, whatever the sign. When every
  signature mimicked the program, positive rows were still presented as opposing it.

  The reversal percentile was `(values > rho).mean()`, which is 0 for every row when all
  scores tie. The midrank convention gives 50, and 0 reads as "nothing opposes more
  strongly" when in fact nothing differs at all.

  A drug effect of exactly zero was called a direction: with a query gene up and the drug
  effect 0, `(value > 0) == (drug > 0)` is False, so the gene was counted as moving
  oppositely. Zero is neutral and is counted as neutral here.

  rho was rounded to four decimals before sorting, so genuinely different scores tied and
  their order came from list position. Sorting uses full precision and breaks ties on the
  signature key.

  A CTRP sensitivity association ranked beside drug treatments. Its sign says a gene's
  expression tracks sensitivity across untreated cell lines, not that a compound moved it.
  Treatment scoring excludes those partitions.

INFERENCE IS NOT AVAILABLE, AND SAYING SO IS THE POINT
The old scorer reported a Student-t p for the Spearman with genes as observations, then a
Benjamini-Hochberg q over reference rows. Genes are not independent observations, the
LINCS inferred genes are functions of the landmarks, and the query genes were selected on
their own significance. Matching scipy proves the arithmetic, not the null. So this engine
reports no p and no q and sets `inference_mode="not_available"`. Source per-gene p-values,
which are measurements, are still exported with the shared genes.

Directional evidence is reported as counts rather than inferred from the sign of rho. A
rank correlation of -1 does not mean the genes moved in opposite directions: for a query
D = [1..20] and a drug effect P = 100 - D, rho is -1 while every gene is positive in both.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
from functools import lru_cache
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    # Annotation only. `_caches()` imports OrderedDict at first use, so the module
    # costs nothing at import time; without this guard ruff's F82 set reports the
    # two string annotations below as undefined names and CI stops.
    from collections import OrderedDict

import numpy as np
from . import config

logger = logging.getLogger("lungmap.drug_matching")

METHOD_VERSION = "cellref2-dart-like-signature-v1"
DEFAULT_STORE = str(config.cellref2_dir() / "drugstore_v1")
MIN_OVERLAP = 10          # legacy compatibility floor, not a validated power calculation
GENE_CAP = 1000

# WHICH GENES A SCREEN IS ALLOWED TO USE, declared per request and recorded on the result.
#
# Nathan, 2026-09-08: "Our method should run analyses both ways... inferred transcriptome
# and the exact way LungChat does... let's compare both."
#
#   all        every gene the signature reports. For LINCS Level 5 that is 978 measured
#              landmarks plus 11,350 inferred from them, so a 3,563-gene query meets about
#              2,609 of them.
#   measured   only genes the assay measured. For LINCS that is the 978 landmarks, which is
#              what the iLINCS connectivity screen uses: every row of its published output
#              reports num_genes = 284 for the same 3,563-gene query. A partition with no
#              landmark mask is an array measurement throughout and is not filtered.
#
# Both recover saracatinib at rank 1 from the published alveolar fibroblast IPF query.
# `measured` also recovers its exact context, HME1 3.33 uM 3 h, at rho -0.4043 against the
# published -0.4505; `all` returns MCF10A 10 uM 6 h at -0.3239. Neither is the right
# answer on its own, which is why the policy is a choice and travels with the result.
GENE_SPACES = ("all", "measured")


def _store_dir() -> str:
    return os.environ.get("CELLREF2_DRUGSTORE", DEFAULT_STORE)


# --------------------------------------------------------------------- the store


class Partition:
    """One source release on its own gene axis, in its own encoding.

    Partitions are never merged into a single matrix. LINCS Level 5 lives on 12,328 genes
    at int8 scale 10; the iLINCS libraries live on a symbol union at their own encoding.
    Putting them on one axis would waste half the array and would force one scale on two
    sources whose values mean different things.
    """

    def __init__(self, base: str, spec: dict[str, Any]) -> None:
        self.id = str(spec["partition_id"])
        self.source = str(spec.get("source", ""))
        self.release = str(spec.get("release", ""))
        self.layout = str(spec.get("layout", "genes_x_signatures"))
        self.scale = float(spec.get("effect_scale", 1.0))
        self.sentinel = spec.get("effect_missing_sentinel")
        self.effect_type = str(spec.get("effect_type", ""))
        self.signature_type = str(spec.get("signature_type", ""))
        self.gene_axis_id = str(spec.get("gene_axis_id", ""))
        self.gene_axis_sha256 = str(spec.get("gene_axis_sha256", ""))
        self.pvalue_available = bool(spec.get("pvalue_available", False))
        self.n_signatures = int(spec["n_signatures"])
        folder = f"{base}/partitions/{self.id}"
        self.folder = folder
        with open(f"{folder}/genes.json") as handle:
            self.genes: list[str] = json.load(handle)
        self.index = {str(g).upper(): i for i, g in enumerate(self.genes)}
        self.values = np.load(f"{folder}/values.npy", mmap_mode="r")
        p_path = f"{folder}/neglog10p.npy"
        self.neglog10p = (np.load(p_path, mmap_mode="r")
                          if self.pvalue_available and os.path.exists(p_path) else None)
        # An L1000 partition marks which genes the assay measured; the rest are inferred
        # from them by the source and must not be exported as independent measurements.
        mask_path = f"{folder}/is_landmark.npy"
        self.is_landmark = (np.load(mask_path, mmap_mode="r")
                            if os.path.exists(mask_path) else None)
        expected = ((len(self.genes), self.n_signatures)
                    if self.layout == "genes_x_signatures"
                    else (self.n_signatures, len(self.genes)))
        if tuple(self.values.shape) != expected:
            raise ValueError(f"{self.id}: values {self.values.shape} is not {expected}")
        # METADATA IS SEEKED, NOT LOADED. Holding every record as a Python dictionary
        # cost 1,489 MiB for 312,438 signatures before a single one was scored. Only the
        # key list stays resident; a record is read from its own line when a row is built.
        with open(f"{folder}/keys.json") as handle:
            self.keys: list[str] = json.load(handle)
        if len(self.keys) != self.n_signatures:
            raise ValueError(f"{self.id}: {len(self.keys)} keys for "
                             f"{self.n_signatures} signatures")
        self.key_row = {key: i for i, key in enumerate(self.keys)}
        self.offsets = np.load(f"{folder}/record_offsets.npy")
        self._records_path = f"{folder}/records.jsonl"
        self._handle = None

    def record(self, row: int) -> dict[str, Any]:
        """One signature's metadata, read from its own line."""
        if self._handle is None:
            self._handle = open(self._records_path, "rb")
        start = int(self.offsets[row])
        length = int(self.offsets[row + 1]) - start
        self._handle.seek(start)
        return json.loads(self._handle.read(length).decode())

    def rows(self, gene_rows: np.ndarray, lo: int, hi: int) -> np.ndarray:
        """Decoded effects for these genes over signatures [lo, hi), as float64.

        The sentinel becomes NaN here and nowhere else, so no caller has to know the
        encoding, and an int8 partition is never decoded with another partition's scale.
        """
        if self.layout == "genes_x_signatures":
            raw = self.values[gene_rows, lo:hi]
        else:
            raw = self.values[lo:hi, :][:, gene_rows].T
        out = np.asarray(raw, dtype=np.float64)
        if self.sentinel is not None:
            out[np.asarray(raw) == self.sentinel] = np.nan
        if self.scale != 1.0:
            out = out / self.scale
        return out

    def column(self, gene_rows: np.ndarray, column: int) -> np.ndarray:
        """One signature's effects at these genes. Reads one signature, not the library."""
        return self.rows(gene_rows, column, column + 1)[:, 0]

    def pvalues(self, gene_rows: np.ndarray, column: int) -> np.ndarray | None:
        if self.neglog10p is None:
            return None
        if self.layout == "genes_x_signatures":
            nl = np.asarray(self.neglog10p[gene_rows, column], dtype=np.float64)
        else:
            nl = np.asarray(self.neglog10p[column, gene_rows], dtype=np.float64)
        with np.errstate(over="ignore"):
            out = np.where(np.isfinite(nl), 10.0 ** (-nl), np.nan)
        return out


class Store:
    def __init__(self, base: str) -> None:
        self.base = base
        with open(f"{base}/manifest.json") as handle:
            self.manifest = json.load(handle)
        # THE LOADER REFUSES AN AGGREGATED STORE. A store built by the old consensus
        # chain declares local_effect_aggregation="modz"; serving it would silently put
        # compound averages back on the page.
        aggregation = str(self.manifest.get("local_effect_aggregation", "none"))
        if aggregation != "none":
            raise ValueError(
                f"{base}: local_effect_aggregation={aggregation!r}; this engine serves "
                "source signatures only")
        self.partitions = {p["partition_id"]: Partition(base, p)
                           for p in self.manifest["partitions"]}
        self.n_signatures = sum(p.n_signatures for p in self.partitions.values())
        seen: set[str] = set()
        for partition in self.partitions.values():
            clash = seen & set(partition.keys)
            if clash:
                raise ValueError(f"{base}: {len(clash)} signature keys appear in more "
                                 f"than one partition, e.g. {sorted(clash)[:3]}")
            seen |= set(partition.keys)
        self.version = str(self.manifest.get("store_version", ""))

    def resolve(self, signature_key: str) -> tuple[str, int] | None:
        """Which partition and row a namespaced key names."""
        for pid, partition in self.partitions.items():
            row = partition.key_row.get(signature_key)
            if row is not None:
                return pid, row
        return None

    def record(self, partition_id: str, row: int) -> dict[str, Any]:
        return self.partitions[partition_id].record(row)


@lru_cache(maxsize=4)
def _store_cached(base: str, stamp: float) -> Store:
    return Store(base)


def store() -> Store:
    base = _store_dir()
    stamp = os.path.getmtime(f"{base}/manifest.json")
    return _store_cached(base, stamp)


# ----------------------------------------------------------------- the query


class PreparedQuery:
    """The selected query genes and effects, fixed once and reused for every partition.

    The same selection is scored against every reference, so a compound is never compared
    against a differently-chosen version of the disease program.
    """

    def __init__(self, features: list[str], effects: list[float], params: dict[str, Any]):
        pairs: dict[str, float] = {}
        collisions = 0
        for name, value in zip(features, effects):
            if value is None:
                continue
            value = float(value)
            if not np.isfinite(value) or value == 0.0:
                continue
            key = str(name).upper()
            if key in pairs:
                collisions += 1
                if abs(value) <= abs(pairs[key]):
                    continue
            pairs[key] = value
        ordered = sorted(pairs.items(), key=lambda kv: (-abs(kv[1]), kv[0]))
        cap = int(params.get("gene_cap") or GENE_CAP)
        self.duplicate_features = collisions
        self.n_offered = len(features)
        self.genes = [g for g, _v in ordered[:cap]]
        self.effects = np.array([v for _g, v in ordered[:cap]], dtype=np.float64)
        self.n_up = int((self.effects > 0).sum())
        self.n_down = int((self.effects < 0).sum())
        self.params = dict(params)
        payload = json.dumps({"genes": self.genes,
                              "effects": [round(v, 6) for v in self.effects.tolist()],
                              "params": self.params, "method": METHOD_VERSION},
                             sort_keys=True)
        self.hash = hashlib.sha256(payload.encode()).hexdigest()[:16]


def prepare_query(features: list[str], effects: list[float],
                  **params: Any) -> PreparedQuery:
    return PreparedQuery(features, effects, params)


# ----------------------------------------------------------------- the statistic


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Average ranks down axis 0, ties sharing the mean of the positions they occupy.

    scipy.stats.rankdata computes the same thing. It is written out here so the tie rule
    is visible in the engine that depends on it, and so a whole block of signatures is
    ranked in one pass instead of one column at a time. Equality with scipy is asserted on
    random and heavily tied input by checks/test_drug_matching.py.

    A run of equal values occupies positions [first, last]; every member takes the mean of
    those positions. Both ends are found by an accumulate along the sorted axis, so no
    Python loop runs over the signatures.
    """
    values = np.asarray(values, dtype=np.float64)
    single = values.ndim == 1
    if single:
        values = values[:, None]
    n = values.shape[0]
    if n == 0:
        return values[:, 0] if single else values
    order = np.argsort(values, axis=0, kind="stable")
    ordered = np.take_along_axis(values, order, axis=0)
    position = np.arange(n, dtype=np.float64)[:, None]

    is_first = np.empty(ordered.shape, dtype=bool)
    is_first[0, :] = True
    is_first[1:, :] = ordered[1:, :] != ordered[:-1, :]
    first = np.maximum.accumulate(np.where(is_first, position, -1.0), axis=0)

    is_last = np.empty(ordered.shape, dtype=bool)
    is_last[-1, :] = True
    is_last[:-1, :] = ordered[1:, :] != ordered[:-1, :]
    last = np.minimum.accumulate(
        np.where(is_last, position, float(n))[::-1, :], axis=0)[::-1, :]

    ranks = (first + last) / 2.0 + 1.0
    out = np.empty_like(ranks)
    np.put_along_axis(out, order, ranks, axis=0)
    return out[:, 0] if single else out


def _spearman_block(query: np.ndarray, block: np.ndarray) -> np.ndarray:
    """Spearman of the query against each column, on that column's own finite genes.

    Ranks are recomputed on each column's shared set. A column whose shared set differs
    gets its own query ranking, which is what Spearman means; ranking once over a wider
    universe and subsetting is a different statistic.

    Columns are grouped by their pattern of present values so that all columns sharing a
    pattern -- every LINCS signature shares one, since Level 5 reports all its genes -- are
    ranked in one pass. A constant query or a constant reference yields NaN, which the
    caller turns into a typed exclusion rather than a score.
    """
    from scipy import stats

    finite = np.isfinite(block)
    out = np.full(block.shape[1], np.nan)
    shared = finite.sum(axis=0)
    packed = np.packbits(finite, axis=0)
    _, group_of = np.unique(packed, axis=1, return_inverse=True)
    group_of = np.asarray(group_of).ravel()
    order = np.argsort(group_of, kind="stable")
    bounds = np.flatnonzero(np.diff(group_of[order])) + 1
    for take in np.split(order, bounds):
        mask = finite[:, take[0]]
        n = int(mask.sum())
        if n < 2:
            continue
        qr = stats.rankdata(query[mask])
        qc = qr - qr.mean()
        qn = float(np.sqrt((qc ** 2).sum()))
        if qn <= 0:
            continue                      # constant query on this shared set
        for lo in range(0, take.size, 4000):
            part = take[lo:lo + 4000]
            cr = stats.rankdata(block[np.ix_(mask, part)], axis=0)
            cc = cr - cr.mean(axis=0)
            cn = np.sqrt((cc ** 2).sum(axis=0))
            with np.errstate(invalid="ignore", divide="ignore"):
                rho = (qc @ cc) / (qn * cn)
            rho[cn <= 0] = np.nan         # constant reference on this shared set
            out[part] = np.clip(rho, -1.0, 1.0)
    return out, shared


def _direction_counts(query: np.ndarray, block: np.ndarray) -> dict[str, np.ndarray]:
    """How the shared genes actually moved, counted rather than inferred from rho.

    A Spearman of -1 does not mean the genes moved oppositely: for D = [1..20] and
    P = 100 - D every gene is positive in both vectors and rho is still -1. So the four
    concordance cells are counted directly. A drug effect of exactly zero is NEUTRAL and
    enters none of them; the old scorer called it "opposite direction" for every up gene.
    """
    finite = np.isfinite(block)
    up = (query > 0)[:, None]
    down = (query < 0)[:, None]
    d_up = finite & (block > 0)
    d_down = finite & (block < 0)
    d_zero = finite & (block == 0)
    counts = {
        "query_up_drug_down": (up & d_down).sum(axis=0),
        "query_down_drug_up": (down & d_up).sum(axis=0),
        "query_up_drug_up": (up & d_up).sum(axis=0),
        "query_down_drug_down": (down & d_down).sum(axis=0),
        "drug_zero": d_zero.sum(axis=0),
    }
    counts["n_opposite_direction"] = (counts["query_up_drug_down"]
                                      + counts["query_down_drug_up"])
    counts["n_same_direction"] = (counts["query_up_drug_up"]
                                  + counts["query_down_drug_down"])
    counts["n_directionally_evaluable"] = (counts["n_opposite_direction"]
                                           + counts["n_same_direction"])
    return counts


def _midrank_percentiles(scores: np.ndarray) -> np.ndarray:
    """The share of the reference group scoring above this row, ties counted as half.

    pct = 100 * (greater + 0.5 * equal) / (n - 1), excluding the row itself. With every
    row tied this is 50 for all of them; the old `(values > rho).mean()` gave 0, which
    reads as "nothing opposes more strongly" when in truth nothing differs.
    """
    n = scores.size
    if n < 2:
        return np.full(n, np.nan)
    order = np.argsort(scores, kind="stable")
    ranked = scores[order]
    greater = n - np.searchsorted(ranked, scores, side="right")
    equal = (np.searchsorted(ranked, scores, side="right")
             - np.searchsorted(ranked, scores, side="left"))
    return 100.0 * (greater + 0.5 * (equal - 1)) / (n - 1)


# ----------------------------------------------------------------- the screen


TREATMENT_TYPES = ("drug treatment",)


def score_library(prepared: PreparedQuery, *, include_types: tuple[str, ...] = (),
                  gene_space: str = "all", chunk: int = 4000) -> dict[str, Any]:
    """Score every eligible source signature. One row each, nothing combined.

    THE RESULT IS COLUMNAR, NOT A LIST OF DICTS. A screen covers 312,438 signatures, and
    building a dictionary per signature cost 3.6 GiB of resident memory -- about 12 KB a
    row -- which two concurrent requests would not survive. The scores and counts are kept
    as numpy arrays, about 8.7 MB for the whole screen, and a row is built only for the
    page being shown or the line being written to an export.

    `include_types` names the signature types eligible for this screen; the default is
    treatments only, so a CTRP sensitivity association cannot enter a reversal ranking.
    Excluded partitions are reported by name and count rather than quietly skipped.
    """
    kinds = tuple(include_types) or TREATMENT_TYPES
    if gene_space not in GENE_SPACES:
        raise ValueError(f"gene_space must be one of {GENE_SPACES}, not {gene_space!r}")
    st = store()
    excluded: dict[str, int] = {}
    coverage: list[dict[str, Any]] = []
    part_ids: list[str] = []
    row_ids: list[int] = []
    acc = {name: [] for name in
           ("rho", "n_shared", "n_opposite", "n_same", "n_evaluable", "n_zero")}

    for pid, partition in st.partitions.items():
        if partition.signature_type not in kinds:
            excluded[f"{pid}:{partition.signature_type or 'unstated'}"] = (
                partition.n_signatures)
            continue
        gene_rows, keep, dropped = [], [], 0
        for slot, gene in enumerate(prepared.genes):
            i = partition.index.get(gene)
            if i is None:
                continue
            if (gene_space == "measured" and partition.is_landmark is not None
                    and not int(partition.is_landmark[i])):
                dropped += 1
                continue
            gene_rows.append(i)
            keep.append(slot)
        coverage.append({"partition_id": pid, "source": partition.source,
                         "gene_space": gene_space,
                         "n_query_genes_inferred_only": dropped,
                         "release": partition.release,
                         "n_signatures": partition.n_signatures,
                         "n_query_genes_in_axis": len(gene_rows),
                         "gene_axis_id": partition.gene_axis_id,
                         "effect_type": partition.effect_type,
                         "effect_scale": partition.scale,
                         "pvalue_available": partition.pvalue_available})
        if len(gene_rows) < MIN_OVERLAP:
            excluded[f"{pid}:too_few_query_genes_in_axis"] = partition.n_signatures
            continue
        gene_rows = np.asarray(gene_rows)
        query = prepared.effects[np.asarray(keep)]
        for lo in range(0, partition.n_signatures, chunk):
            hi = min(lo + chunk, partition.n_signatures)
            block = partition.rows(gene_rows, lo, hi)
            rho, shared = _spearman_block(query, block)
            counts = _direction_counts(query, block)
            good = np.isfinite(rho) & (shared >= MIN_OVERLAP)
            short = int(((shared < MIN_OVERLAP)).sum())
            if short:
                excluded["below_minimum_overlap"] = (
                    excluded.get("below_minimum_overlap", 0) + short)
            constant = int((~np.isfinite(rho) & (shared >= MIN_OVERLAP)).sum())
            if constant:
                excluded["constant_query_or_reference"] = (
                    excluded.get("constant_query_or_reference", 0) + constant)
            take = np.where(good)[0]
            if not take.size:
                continue
            part_ids.extend([pid] * take.size)
            row_ids.extend((take + lo).tolist())
            acc["rho"].append(rho[take])
            acc["n_shared"].append(shared[take])
            acc["n_opposite"].append(counts["n_opposite_direction"][take])
            acc["n_same"].append(counts["n_same_direction"][take])
            acc["n_evaluable"].append(counts["n_directionally_evaluable"][take])
            acc["n_zero"].append(counts["drug_zero"][take])
            del block, counts

    n = len(row_ids)
    columns = {name: (np.concatenate(parts) if parts else np.zeros(0))
               for name, parts in acc.items()}
    if n:
        percentile = _midrank_percentiles(columns["rho"])
        # FULL PRECISION, THEN THE KEY. The old scorer sorted a four-decimal rounding, so
        # scores differing in the fifth place tied and their order came from list position.
        keys = [f"{part_ids[i]}#{row_ids[i]:09d}" for i in range(n)]
        order = sorted(range(n), key=lambda i: (columns["rho"][i], keys[i]))
        order = np.asarray(order)
    else:
        percentile = np.zeros(0)
        order = np.zeros(0, dtype=int)

    result_id = hashlib.sha256(json.dumps({
        "query": prepared.hash, "kinds": list(kinds), "store": st.version,
        "method": METHOD_VERSION, "min_overlap": MIN_OVERLAP,
        "gene_space": gene_space,
    }, sort_keys=True).encode()).hexdigest()[:16]

    rho = columns["rho"]
    return {
        "result_id": result_id,
        "gene_space": gene_space,
        "gene_space_note": (
            "every gene the signature reports, measured and inferred"
            if gene_space == "all" else
            "only genes the assay measured; for LINCS Level 5 that is the 978 landmarks, "
            "which is the space the iLINCS connectivity screen uses"),
        "method_version": METHOD_VERSION,
        "store_version": st.version,
        "query_hash": prepared.hash,
        "result_unit": "source_signature",
        "local_effect_aggregation": "none",
        "local_score_aggregation": "none",
        "inference_mode": "not_available",
        "calibration_status": "not_run",
        "statistic": ("signed Spearman between the query effects and one source "
                      "signature's effects, on the genes both report"),
        "sign": "negative opposes the query program, positive mimics it",
        "n_query_features_offered": prepared.n_offered,
        "n_query_genes_selected": len(prepared.genes),
        "n_query_up": prepared.n_up,
        "n_query_down": prepared.n_down,
        "duplicate_query_features_resolved": prepared.duplicate_features,
        "eligible_signature_types": list(kinds),
        "n_scored": n,
        "n_opposing": int((rho < 0).sum()),
        "n_mimicking": int((rho > 0).sum()),
        "n_neutral": int((rho == 0).sum()),
        "excluded": excluded,
        "coverage": coverage,
        "_columns": columns,
        "_percentile": percentile,
        "_partition": part_ids,
        "_row": row_ids,
        "_order": order,
    }


def row_at(result: dict[str, Any], i: int) -> dict[str, Any]:
    """Materialise one result row. Called for a page or an export line, not for a screen."""
    st = store()
    pid = result["_partition"][i]
    record = st.record(pid, int(result["_row"][i]))
    columns = result["_columns"]
    evaluable = int(columns["n_evaluable"][i])
    pct = result["_percentile"][i] if result["_percentile"].size else float("nan")
    out = {
        "signature_key": record["signature_key"],
        "source_signature_id": record["source_signature_id"],
        "partition_id": pid, "row_index": int(result["_row"][i]),
        "rho": float(columns["rho"][i]),
        "n_shared": int(columns["n_shared"][i]),
        "n_opposite_direction": int(columns["n_opposite"][i]),
        "n_same_direction": int(columns["n_same"][i]),
        "n_directionally_evaluable": evaluable,
        "drug_zero": int(columns["n_zero"][i]),
        "opposite_direction_fraction":
            (float(columns["n_opposite"][i]) / evaluable if evaluable else None),
        "signature_reversal_percentile": (None if not np.isfinite(pct) else float(pct)),
    }
    for field, key in (("compound", "compound_name"), ("compound_id", "compound_id"),
                       ("pubchem_cid", "pubchem_cid"), ("library", "library"),
                       ("source", "source"), ("source_release", "source_release"),
                       ("signature_type", "signature_type"),
                       ("effect_type", "effect_type"),
                       ("perturbation_data_level", "perturbation_data_level"),
                       ("upstream_replicate_aggregation",
                        "upstream_replicate_aggregation"),
                       ("experimental_species", "experimental_species"),
                       ("identifier_species", "identifier_species"),
                       ("tissue", "tissue"), ("cell_line", "cell_line"),
                       ("dose", "dose_original"), ("time", "time_original"),
                       ("source_url", "source_url")):
        out[field] = record.get(key, "")
    return out


def page(result: dict[str, Any], tail: str = "opposing", offset: int = 0,
         limit: int = 25, compound: str = "") -> tuple[list[dict[str, Any]], int, int]:
    """One page of a screen, how many rows the tail holds, and how many the search matched.

    `compound` is a DISPLAY-ONLY search. It changes which rows are shown and nothing else:
    no score, no percentile and no denominator moves, and every row still carries the rank
    it holds in the whole tail. A named drug has to be findable -- saracatinib has 373
    signatures in this store, and a positive control cannot be checked by paging through
    165,000 rows -- but finding it must not quietly redefine the screen.
    """
    order = result["_order"]
    rho = result["_columns"]["rho"]
    if tail == "opposing":
        chosen = order[rho[order] < 0]
    elif tail == "mimicking":
        chosen = order[rho[order] > 0][::-1]
    else:
        chosen = order
    n_tail = int(chosen.size)
    ranks = np.arange(1, n_tail + 1)
    needle = str(compound or "").strip().lower()
    if needle:
        st = store()
        keep_i, keep_rank = [], []
        for position, i in enumerate(chosen):
            record = st.record(result["_partition"][int(i)],
                               int(result["_row"][int(i)]))
            if (needle in str(record.get("compound_name", "")).lower()
                    or needle in str(record.get("compound_id", "")).lower()
                    or needle in str(record.get("source_signature_id", "")).lower()):
                keep_i.append(int(i))
                keep_rank.append(int(ranks[position]))
        chosen, ranks = np.asarray(keep_i, dtype=int), np.asarray(keep_rank, dtype=int)
    lo = max(0, int(offset))
    span = slice(lo, lo + max(1, min(int(limit), 500)))
    rows = []
    for i, rank in zip(chosen[span], ranks[span]):
        row = row_at(result, int(i))
        row["rank_in_tail"] = int(rank)
        rows.append(row)
    return rows, n_tail, int(chosen.size)


def iter_rows(result: dict[str, Any]):
    """Every scored row in rank order, one at a time, for a full export."""
    for i in result["_order"]:
        yield row_at(result, int(i))


def scored_genes(prepared: PreparedQuery, signature_key: str) -> dict[str, Any]:
    """The exact genes that produced one signature's score, with both sides' statistics.

    Reads ONE signature. The old export read every reference column and then selected one.
    """
    st = store()
    found = st.resolve(signature_key)
    if found is None:
        return {"error": f"{signature_key} is not a signature of this store."}
    pid, column = found
    partition = st.partitions[pid]
    record = partition.record(column)
    gene_rows, keep = [], []
    for slot, gene in enumerate(prepared.genes):
        i = partition.index.get(gene)
        if i is not None:
            gene_rows.append(i)
            keep.append(slot)
    if not gene_rows:
        return {"error": f"no query gene is on the {partition.id} gene axis."}
    gene_rows = np.asarray(gene_rows)
    effects = partition.column(gene_rows, column)
    pvals = partition.pvalues(gene_rows, column)
    query = prepared.effects[np.asarray(keep)]
    finite = np.isfinite(effects)
    from scipy import stats
    q_rank = stats.rankdata(query[finite]) if finite.any() else np.zeros(0)
    d_rank = stats.rankdata(effects[finite]) if finite.any() else np.zeros(0)
    rank_at = {int(np.where(finite)[0][i]): i for i in range(int(finite.sum()))}

    rows = []
    for slot, gene_row in enumerate(gene_rows):
        drug = float(effects[slot])
        used = bool(np.isfinite(drug))
        agreement = ("neutral: the drug effect is zero" if used and drug == 0 else
                     "same direction" if used and (query[slot] > 0) == (drug > 0) else
                     "opposite direction" if used else "not measured")
        p = None
        if pvals is not None and np.isfinite(pvals[slot]):
            p = float(pvals[slot])
        at = rank_at.get(slot)
        rows.append({
            "canonical_gene_id": partition.genes[int(gene_row)],
            "query_effect": float(query[slot]),
            "drug_effect_decoded": drug if used else None,
            "drug_p": p,
            "drug_p_status": ("available" if p is not None else
                              "not published by this source"
                              if not partition.pvalue_available else "missing"),
            "query_rank_used": float(q_rank[at]) if at is not None else None,
            "drug_rank_used": float(d_rank[at]) if at is not None else None,
            "direction_agreement": agreement,
            "used_in_score": used,
            "effect_encoding": str(partition.values.dtype),
            "effect_scale": partition.scale,
            # THE FIELD DESCRIBES THIS SIGNATURE'S VALUE, not the assay in general. A
            # gene the signature does not report is neither measured nor inferred, and
            # calling it "measured" while its effect column is blank contradicts itself.
            "measured_or_inferred": (
                "not reported by this signature" if not used
                else "measured" if partition.is_landmark is None
                else ("measured" if int(partition.is_landmark[int(gene_row)])
                      else "inferred by the source")),
        })
    n_used = sum(1 for r in rows if r["used_in_score"])
    return {
        "signature_key": signature_key,
        "record": record,
        "partition_id": partition.id,
        "effect_type": partition.effect_type,
        "n_query_genes_selected": len(prepared.genes),
        "n_on_axis": len(gene_rows),
        "n_used_in_score": n_used,
        "rows": rows,
    }


# ------------------------------------------------------- one result, shared by all views

_RESULTS: "OrderedDict[str, dict[str, Any]]" = None      # type: ignore[assignment]
_QUERIES: "OrderedDict[str, PreparedQuery]" = None       # type: ignore[assignment]
_CACHE_SIZE = 8


def _caches():
    """The table, the TSV and the gene download must describe ONE screen.

    The old page rebuilt the query from whatever the controls said at the moment a link
    was clicked, so changing a filter after a screen produced a download from a different
    query than the rows on screen. Every view now quotes a result_id and reads it here.
    """
    global _RESULTS, _QUERIES                                        # noqa: PLW0603
    if _RESULTS is None:
        from collections import OrderedDict
        _RESULTS = OrderedDict()
        _QUERIES = OrderedDict()
    return _RESULTS, _QUERIES


# A SCREEN SURVIVES THE PROCESS THAT BUILT IT.
#
# Scoring 323,661 signatures costs 25 to 29 s, measured on 2026-09-20, and the
# answer is then served in 47 ms from `_caches()`. That cache dies with the worker,
# so a restart, a deploy and every additional uvicorn worker each made one reader
# wait half a minute. The finished screen is 12.7 MB on disk and reads back in
# 0.031 s, a factor of 932, with every score identical.
#
# The file name is the result id, which already carries the query hash, the store
# version, the method version and the gene space, so a rebuilt store cannot be
# answered from an old file. Nothing here may break a request: a cache that cannot
# be read or written is logged and the screen is computed as before.
_DISK_CACHE_FILES = 64


def _disk_cache_dir():
    from pathlib import Path

    path = Path(config.state_dir()) / "cache" / "drug_screen"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _disk_load(probe: str) -> dict[str, Any] | None:
    import pickle

    try:
        path = _disk_cache_dir() / f"{probe}.pkl"
        if not path.is_file():
            return None
        with path.open("rb") as fh:
            result = pickle.load(fh)
        os.utime(path, None)                 # newest-used, for the trim below
    except Exception as exc:                 # noqa: BLE001 - never fail a request
        logger.warning("drug screen cache unreadable for %s: %s", probe, exc)
        return None
    if not isinstance(result, dict) or result.get("result_id") != probe:
        logger.warning("drug screen cache holds the wrong result for %s", probe)
        return None
    return result


def _disk_save(probe: str, result: dict[str, Any]) -> None:
    import pickle

    try:
        directory = _disk_cache_dir()
        tmp = directory / f".{probe}.{os.getpid()}.tmp"
        with tmp.open("wb") as fh:
            pickle.dump(result, fh, protocol=5)
        os.replace(tmp, directory / f"{probe}.pkl")
        files = sorted(directory.glob("*.pkl"), key=lambda p: p.stat().st_mtime)
        for old in files[:-_DISK_CACHE_FILES]:
            old.unlink(missing_ok=True)
    except Exception as exc:                 # noqa: BLE001 - never fail a request
        logger.warning("drug screen cache not written for %s: %s", probe, exc)


def screen(contrast: str, cell_state: str, modality: str = "rna",
           min_abs_log2fc: float = 0.0, max_fdr: float = 0.05,
           gene_cap: int = GENE_CAP, gene_space: str = "all",
           include_types: tuple[str, ...] = ()) -> dict[str, Any]:
    """Prepare the query from the atlas and score every eligible source signature."""
    from . import multimodal as mm

    names, folds, _counts = mm.program_vector(
        contrast, cell_state, modality=modality,
        min_abs_log2fc=float(min_abs_log2fc), max_fdr=float(max_fdr))
    prepared = prepare_query(
        names, folds, contrast=contrast, cell_state=cell_state, modality=modality,
        min_abs_log2fc=float(min_abs_log2fc), max_fdr=float(max_fdr),
        gene_cap=int(gene_cap), gene_space=gene_space)
    results, queries = _caches()
    probe = hashlib.sha256(json.dumps(
        {"query": prepared.hash, "kinds": list(include_types or TREATMENT_TYPES),
         "store": store().version, "method": METHOD_VERSION,
         "min_overlap": MIN_OVERLAP, "gene_space": gene_space},
        sort_keys=True).encode()).hexdigest()[:16]
    if probe in results:
        results.move_to_end(probe)
        return results[probe]
    stored = _disk_load(probe)
    if stored is not None:
        results[probe] = stored
        queries[probe] = prepared
        while len(results) > _CACHE_SIZE:
            old, _ = results.popitem(last=False)
            queries.pop(old, None)
        return stored
    result = score_library(prepared, include_types=include_types,
                           gene_space=gene_space)
    result["query"] = {"contrast": contrast, "cell_state": cell_state,
                       "modality": modality, "min_abs_log2fc": float(min_abs_log2fc),
                       "max_fdr": float(max_fdr), "gene_cap": int(gene_cap),
                       "gene_space": gene_space}
    results[result["result_id"]] = result
    queries[result["result_id"]] = prepared
    _disk_save(result["result_id"], result)
    while len(results) > _CACHE_SIZE:
        old, _ = results.popitem(last=False)
        queries.pop(old, None)
    return result


def result_for(result_id: str) -> tuple[dict[str, Any] | None, PreparedQuery | None]:
    results, queries = _caches()
    return results.get(result_id), queries.get(result_id)


def provenance() -> dict[str, Any]:
    """What the deployed store holds, read from the store rather than written here."""
    st = store()
    by_partition = []
    for spec in st.manifest["partitions"]:
        pid = spec["partition_id"]
        partition = st.partitions[pid]
        # A SAMPLE, NOT A SCAN. Counting compounds and cell models over every record would
        # read the whole metadata file on each request; the manifest carries the totals
        # the build computed and this reports those.
        by_partition.append({
            "partition_id": pid,
            "source": spec.get("source", ""),
            "release": spec.get("release", ""),
            "n_signatures": partition.n_signatures,
            "n_genes": spec.get("n_genes"),
            "signature_type": spec.get("signature_type", ""),
            "effect_type": spec.get("effect_type", ""),
            "effect_encoding": spec.get("effect_encoding", ""),
            "effect_scale": spec.get("effect_scale"),
            "pvalue_available": bool(spec.get("pvalue_available")),
            "eligible_for_treatment_scoring":
                spec.get("signature_type") in TREATMENT_TYPES,
            "n_compounds": spec.get("n_compounds"),
            "species": spec.get("species", []),
            "cell_lines": spec.get("n_cell_lines"),
            "tissues": spec.get("tissues", []),
            "source_url": spec.get("source_url", ""),
        })
    return {
        "store_version": st.version,
        "method_version": METHOD_VERSION,
        "result_unit": st.manifest.get("result_unit", "source_signature"),
        "local_effect_aggregation": st.manifest.get("local_effect_aggregation", "none"),
        "inference_mode": "not_available",
        "n_signatures": st.n_signatures,
        "n_signatures_eligible": sum(p["n_signatures"] for p in by_partition
                                     if p["eligible_for_treatment_scoring"]),
        "n_compounds_named": st.manifest.get("n_compounds_named"),
        "partitions": by_partition,
        "note": st.manifest.get("note", ""),
        "not_used": ["DrugBank"],
    }
