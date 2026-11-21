"""scALABLE input preparation; MarkerFinder's statistical/scaling gates stay native.

The inferred-log route was explicitly authorized by Nathan on 2026-11-20.
It reconstructs processed abundances, never raw counts or missing measurements.
"""
from copy import deepcopy
from types import SimpleNamespace

import anndata as ad
import numpy as np
from scipy import sparse

from . import markerFinder as native


class RecoveryUnavailable(native.MarkerFinderInputError):
    """No sufficiently supported source; retain the other completed analyses."""


def _blocks(matrix):
    if sparse.issparse(matrix):
        for start in range(0, matrix.nnz, 4_000_000):
            yield matrix.data[start:start + 4_000_000]
    else:
        step = max(1, 4_000_000 // max(1, matrix.shape[1]))
        for start in range(0, matrix.shape[0], step):
            yield matrix[start:start + step]


def _validate_values(matrix, workspace):
    for values in _blocks(matrix):
        if not np.isfinite(values).all() or (values < 0).any():
            raise RecoveryUnavailable("RNA contains negative or nonfinite values; no log inverse was inferred.")
        _release_block(matrix, values, workspace)


def _release_block(matrix, values, workspace, *, write=False):
    if workspace is None:
        return
    if sparse.issparse(matrix):
        start = (values.ctypes.data - matrix.data.ctypes.data) // matrix.dtype.itemsize
        workspace.release_range(matrix.data, start, start + values.size, write=write)
    else:
        workspace.release_pages()


def _diagnose(matrix, workspace):
    if workspace is not None:
        obj = SimpleNamespace(X=matrix, layers={}, n_obs=matrix.shape[0], n_vars=matrix.shape[1])
        matrix = workspace.scale_probe(obj).X
    return native.detect_input_scaling(matrix)


def _working(adata, matrix):
    result = ad.AnnData(X=matrix, obs=adata.obs.copy(), var=adata.var.copy())
    if "lineage_order" in adata.uns:
        result.uns["lineage_order"] = deepcopy(adata.uns["lineage_order"])
    return result


def prepare_marker_input(adata):
    """Return (hook input, supported hook options, provenance) after X rejection.

    Every working matrix has exactly the supplied RNA cell and gene axes. Raw
    matrices may contain additional genes, but must resolve every supplied gene
    uniquely. No intersection, alias guessing or count rounding is permitted.
    """
    workspace = getattr(adata, "_matrix_workspace", None)
    inspected = []
    candidates = [(f"layers[{key!r}]", matrix, key) for key, matrix in adata.layers.items()]
    candidates.sort(key=lambda item: (item[2] != "counts", item[2]))
    for source, matrix, layer in candidates:
        try:
            report = _diagnose(matrix, workspace)
            _validate_values(matrix, workspace)
        except native.MarkerFinderInputError as error:
            inspected.append({"source": source, "reason": str(error)})
            continue
        inspected.append({"source": source, "diagnostics": report})
        if report["status"] not in {"counts", "counts_like", "linear_normalized"}:
            continue
        return adata, {"layer": layer, "scale_data": True}, {
            "method": "verified_linear_source", "source": source,
            "diagnostics": report, "inspected_sources": inspected,
        }
    if adata.raw is not None:
        raw = adata.raw
        if raw.var_names.is_unique and adata.var_names.is_unique:
            positions = raw.var_names.get_indexer(adata.var_names)
            if (positions >= 0).all():
                matrix = (workspace.select(raw.X, cols=positions) if workspace is not None
                          else raw.X[:, positions])
                try:
                    report = _diagnose(matrix, workspace)
                    _validate_values(matrix, workspace)
                except native.MarkerFinderInputError as error:
                    inspected.append({"source": "raw.X", "reason": str(error)})
                    report = {"status": "unusable"}
                inspected.append({"source": "raw.X", "diagnostics": report})
                if report["status"] in {"counts", "counts_like", "linear_normalized"}:
                    return _working(adata, matrix), {"scale_data": True}, {
                        "method": "verified_linear_source", "source": "raw.X",
                        "diagnostics": report, "inspected_sources": inspected,
                    }
            else:
                inspected.append({"source": "raw.X", "reason": "Incomplete supplied RNA gene coverage",
                                  "missing_genes": adata.var_names[positions < 0].tolist()})
        else:
            inspected.append({"source": "raw.X", "reason": "Nonunique gene identifiers"})

    _validate_values(adata.X, workspace)
    report = _diagnose(adata.X, workspace)
    if report["status"] in {"counts", "counts_like", "linear_normalized"}:
        return adata, {"scale_data": True}, {
            "method": "verified_linear_source", "source": "X", "diagnostics": report,
            "inspected_sources": inspected,
        }
    log_info = adata.uns.get("log1p")
    recorded = isinstance(log_info, dict) and "base" in log_info
    if recorded:
        try:
            base = np.e if log_info["base"] is None else float(log_info["base"])
        except (TypeError, ValueError) as error:
            raise RecoveryUnavailable("Recorded log base is invalid; uploaded RNA was retained.") from error
        if not np.isfinite(base) or base <= 1:
            raise RecoveryUnavailable("Recorded log base is invalid; uploaded RNA was retained.")
        evidence = "recorded log1p metadata"
    else:
        winner = report.get("best_transform")
        other = "natural-log" if winner == "log2" else "log2"
        cv = report.get("cv_cell_total", np.inf)
        supported = (winner in {"log2", "natural-log"} and 0 <= cv <= .10
                     and report.get(f"cv_cell_total[{other}]", 0) >= max(2 * cv, .01)
                     and report.get("cv_cell_total[identity]", 0) >= max(1.5 * cv, .01))
        if not supported:
            raise RecoveryUnavailable("No verified counts or unambiguous log-scale evidence; uploaded RNA and completed alignment were retained.")
        base = 2.0 if winner == "log2" else float(np.e)
        evidence = "inferred from separated inverse-log cell-total diagnostics"

    # A separate working copy; no new counts layer or changes to source X.
    work = (workspace.copy_matrix(adata.X, dtype=np.float64) if workspace is not None
            else adata.X.astype(np.float64, copy=True))
    for values in _blocks(work):
        with np.errstate(over="ignore", invalid="ignore"):
            np.multiply(values, np.log(base), out=values)
            np.expm1(values, out=values)
        if not np.isfinite(values).all():
            raise RecoveryUnavailable("Recorded/inferred log inverse overflowed; source RNA was retained.")
        _release_block(work, values, workspace, write=True)
    normalized, normalization = native.scale_expression_matrix(work, workspace=workspace)
    # The public hook also runs this strict gate. Check here to diagnose a failed
    # recovery before it reaches statistics; do not disable or relax the gate.
    check = _diagnose(normalized, workspace)
    if check["status"] != "ok":
        raise RecoveryUnavailable("Recovered working expression failed native depth-normalization validation.")
    return _working(adata, normalized), {}, {
        "method": "processed_abundance_log_inverse", "source": "X", "log_base": base,
        "scale_evidence": evidence, "diagnostics": report, "normalization": normalization,
        "validation": check, "inspected_sources": inspected,
        "original_expression_retained": True, "raw_counts_reconstructed": False,
        "policy_version": "scalable-marker-input-v1",
    }
