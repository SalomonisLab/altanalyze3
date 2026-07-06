"""
fast_infer_affinity
===================

Pure-NumPy (BLAS-backed) reimplementation of the forward pass of mhcflurry's
class1 *pan-allele* affinity ensemble, producing output identical (to within
floating-point tolerance) to::

    p.predict_affinity(peptides=<9/10-mers>, alleles={'s': [allele]},
                       include_affinity_percentile=False, verbose=0).affinity

but without any per-call TensorFlow / PyTorch overhead.

The public entry point is :func:`predict_affinity_fast`.

How the reference pipeline works (mhcflurry 2.2.1, PyTorch backend)
------------------------------------------------------------------
For a single allele ``predict_affinity`` reduces to
``affinity_predictor.predict(peptides, allele=allele)`` which:

1. Encodes peptides with ``left_pad_centered_right_pad`` / max_length=15 /
   BLOSUM62 -> ``(n, 45, 21)`` float, flattened to ``(n, 945)``.
2. Encodes the allele pseudosequence (37 aa) with BLOSUM62 -> ``(37, 21)`` ->
   flattened ``(777,)``; broadcast to every peptide.
3. The loaded ensemble has been *optimized* (``optimize()``) into a single
   ``MergedClass1NeuralNetwork`` that concatenates the outputs of its N (=10)
   sub-networks. Each sub-network:
       merged = concat(peptide_flat[945], allele_vec[777])  -> 1722
       -> a stack of Dense(+tanh) layers (feedforward or DenseNet
          "with-skip-connections") -> Dense(1) -> sigmoid
   Dropout layers are identity at inference.
4. Each sub-network output ``o_i`` (a regression target in [0,1]) is mapped to
   nM affinity by ``to_ic50(o) = 50000 ** (1 - o)`` (float64).
5. The N per-model nM affinities are combined with the ensemble "mean"
   centrality *in log space*: ``affinity = exp(mean_i(log(ic50_i)))`` (i.e. the
   geometric mean).

This module extracts the Dense/output weights of every sub-network once
(cached on the predictor), reuses mhcflurry's own numpy encoding helpers for
peptides and the allele pseudosequence, and evaluates steps 1-5 with numpy
matmuls.

Notes
-----
* All matmuls / activations run in float32 to mirror the PyTorch network
  (whose parameters and intermediate tensors are float32); the ic50 transform
  and ensemble geometric mean run in float64 exactly as mhcflurry does.
* This module does NOT import torch. Weights are pulled out of the already-built
  torch modules via ``.detach().cpu().numpy()`` at cache-build time only.
"""

import ctypes
import sys
import weakref

import numpy as np

from mhcflurry.encodable_sequences import EncodableSequences
from mhcflurry.allele_encoding import AlleleEncoding
from mhcflurry.regression_target import to_ic50
from mhcflurry.ensemble_centrality import CENTRALITY_MEASURES


# Keras/PyTorch BatchNorm epsilon (only used if a model has BatchNorm layers;
# the standard pan-allele models do not).
_BATCH_NORM_EPSILON = 1e-3

# Chunk size for the forward pass (does not affect results; rows are
# independent). Bounds peak memory for large peptide batches.
_FORWARD_CHUNK = 32768

# Cache of extracted weights, keyed (weakly) on the affinity predictor object.
_CACHE = weakref.WeakKeyDictionary()


# ----------------------------------------------------------------------------
# BLAS backend: single-precision GEMM.
#
# numpy's bundled OpenBLAS on Apple Silicon uses generic armv8 kernels and is
# ~8x slower than the platform's Accelerate (AMX) BLAS -- which is what the
# PyTorch backend uses. To be genuinely faster than mhcflurry we call Apple
# Accelerate's ``cblas_sgemm`` directly via ctypes (still operating on plain
# numpy float32 arrays -- no torch / tensorflow involved). On any platform where
# Accelerate is unavailable we transparently fall back to numpy's ``@``.
# ----------------------------------------------------------------------------
_CBLAS_ROW_MAJOR = 101
_CBLAS_NO_TRANS = 111
_CBLAS_TRANS = 112


def _load_accelerate():
    """Load Apple Accelerate and return (cblas_sgemm, vvtanhf) or (None, None)."""
    if sys.platform != "darwin":
        return None, None
    lib = None
    for path in (
        "/System/Library/Frameworks/Accelerate.framework/Accelerate",
        "Accelerate",
    ):
        try:
            lib = ctypes.CDLL(path)
            break
        except OSError:
            continue
    if lib is None:
        return None, None
    try:
        sgemm = lib.cblas_sgemm
        sgemm.restype = None
        sgemm.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_int,   # order, transA, transB
            ctypes.c_int, ctypes.c_int, ctypes.c_int,   # M, N, K
            ctypes.c_float,                              # alpha
            ctypes.c_void_p, ctypes.c_int,               # A, lda
            ctypes.c_void_p, ctypes.c_int,               # B, ldb
            ctypes.c_float,                              # beta
            ctypes.c_void_p, ctypes.c_int,               # C, ldc
        ]
    except AttributeError:
        sgemm = None
    try:
        vvtanhf = lib.vvtanhf   # vForce vectorized tanh (float32)
        vvtanhf.restype = None
        vvtanhf.argtypes = [
            ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
    except AttributeError:
        vvtanhf = None
    return sgemm, vvtanhf


_ACCEL_SGEMM, _ACCEL_VVTANHF = _load_accelerate()
_INT32_MAX = 2 ** 31 - 1


def _gemm_xwt(x, W):
    """Return float32 ``x @ W.T`` where x is (n, k) and W is (out, k) (torch
    Linear weight layout). Uses Accelerate when available, else numpy."""
    x = np.ascontiguousarray(x, dtype=np.float32)
    n, k = x.shape
    out = W.shape[0]
    if _ACCEL_SGEMM is None or n == 0:
        return x @ W.T
    W = np.ascontiguousarray(W, dtype=np.float32)
    C = np.empty((n, out), dtype=np.float32)
    # C(n,out) = x(n,k) @ W(out,k)^T   -> transB = Trans, ldb = k
    _ACCEL_SGEMM(
        _CBLAS_ROW_MAJOR, _CBLAS_NO_TRANS, _CBLAS_TRANS,
        n, out, k, 1.0,
        x.ctypes.data_as(ctypes.c_void_p), k,
        W.ctypes.data_as(ctypes.c_void_p), k,
        0.0,
        C.ctypes.data_as(ctypes.c_void_p), out,
    )
    return C


# ----------------------------------------------------------------------------
# Activations (numpy, float32) mirroring torch.tanh / torch.sigmoid
# ----------------------------------------------------------------------------
def _tanh(x):
    if _ACCEL_VVTANHF is not None and x.size <= _INT32_MAX:
        x = np.ascontiguousarray(x, dtype=np.float32)
        y = np.empty_like(x)
        n = ctypes.c_int(x.size)
        _ACCEL_VVTANHF(
            y.ctypes.data_as(ctypes.c_void_p),
            x.ctypes.data_as(ctypes.c_void_p),
            ctypes.byref(n))
        return y
    return np.tanh(x)


def _sigmoid(x):
    # Numerically stable logistic sigmoid.
    out = np.empty_like(x)
    pos = x >= 0
    neg = ~pos
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    ex = np.exp(x[neg])
    out[neg] = ex / (1.0 + ex)
    return out


def _get_activation(name):
    if not name or name == "linear":
        return None
    name = name.lower()
    if name == "tanh":
        return _tanh
    if name == "sigmoid":
        return _sigmoid
    if name == "relu":
        return lambda x: np.maximum(x, 0.0)
    raise ValueError("Unknown activation: %s" % name)


def _linear(x, W, b):
    """Dense layer: torch stores weight as (out, in); y = x @ W.T + b."""
    y = _gemm_xwt(x, W)
    y += b
    return y


def _batch_norm(x, bn):
    gamma, beta, mean, var = bn
    return (x - mean) / np.sqrt(var + _BATCH_NORM_EPSILON) * gamma + beta


def _locally_connected(x, lc):
    """
    Faithful numpy port of mhcflurry.pytorch_layers.LocallyConnected1D.forward.

    x: (batch, in_len, in_ch)
    weight: (out_len, out_ch, kernel_size * in_ch)   [kernel-major, then channel]
    bias:   (out_len, out_ch)
    returns (batch, out_len, out_ch)
    """
    W = lc["weight"]
    b = lc["bias"]
    k = lc["kernel_size"]
    in_ch = lc["in_channels"]
    out_len = lc["output_length"]
    batch = x.shape[0]
    # Build patches with kernel positions first then channels, matching the
    # torch code: unfold -> permute(0,1,3,2) -> reshape(batch, out_len, k*in_ch)
    patches = np.stack(
        [x[:, j:j + k, :] for j in range(out_len)], axis=1
    )  # (batch, out_len, k, in_ch)
    patches = patches.reshape(batch, out_len, k * in_ch)
    # out[b,o,f] = sum_i patches[b,o,i] * W[o,f,i] + bias[o,f]
    # -> batched matmul over out_len positions.
    p = np.transpose(patches, (1, 0, 2))          # (out_len, batch, k*in_ch)
    Wt = np.transpose(W, (0, 2, 1))               # (out_len, k*in_ch, out_ch)
    out = np.matmul(p, Wt)                         # (out_len, batch, out_ch)
    out = np.transpose(out, (1, 0, 2)) + b         # (batch, out_len, out_ch)
    act = _get_activation(lc["activation"])
    if act is not None:
        out = act(out)
    return out


# ----------------------------------------------------------------------------
# Weight extraction
# ----------------------------------------------------------------------------
def _np32(t):
    return t.detach().cpu().numpy().astype(np.float32)


def _extract_subnet(sub):
    """Extract all inference-relevant weights of one Class1NeuralNetworkModel."""
    info = {
        "topology": sub.topology,
        "merge_method": sub.peptide_allele_merge_method,
        "merge_activation": sub.peptide_allele_merge_activation,
        "activation": sub.activation_name,
        "output_activation": sub.output_activation_name,
        "has_allele": sub.has_allele,
    }

    info["lc_layers"] = [
        {
            "weight": _np32(lc.weight),
            "bias": _np32(lc.bias),
            "kernel_size": lc.kernel_size,
            "in_channels": lc.in_channels,
            "output_length": lc.output_length,
            "activation": lc.activation_name,
        }
        for lc in sub.lc_layers
    ]

    info["peptide_dense"] = [
        (_np32(l.weight), _np32(l.bias)) for l in sub.peptide_dense_layers
    ]

    if sub.batch_norm_early is not None:
        bn = sub.batch_norm_early
        info["batch_norm_early"] = (
            _np32(bn.weight), _np32(bn.bias),
            bn.running_mean.detach().cpu().numpy().astype(np.float32),
            bn.running_var.detach().cpu().numpy().astype(np.float32),
        )
    else:
        info["batch_norm_early"] = None

    info["allele_dense"] = [
        (_np32(l.weight), _np32(l.bias)) for l in sub.allele_dense_layers
    ]

    dense = []
    for i, l in enumerate(sub.dense_layers):
        bn = sub.batch_norms[i]
        bn_w = None
        if bn is not None:
            bn_w = (
                _np32(bn.weight), _np32(bn.bias),
                bn.running_mean.detach().cpu().numpy().astype(np.float32),
                bn.running_var.detach().cpu().numpy().astype(np.float32),
            )
        dense.append((_np32(l.weight), _np32(l.bias), bn_w))
    info["dense"] = dense

    info["output"] = (_np32(sub.output_layer.weight), _np32(sub.output_layer.bias))
    return info


def _resolve_affinity_predictor(predictor):
    """Accept either a Class1PresentationPredictor or a Class1AffinityPredictor."""
    if hasattr(predictor, "affinity_predictor"):
        return predictor.affinity_predictor
    return predictor


def _build_cache(ap):
    """Extract and cache sub-network weights and metadata for an affinity predictor."""
    if not ap.class1_pan_allele_models:
        raise ValueError("Predictor has no class1_pan_allele_models (pan-allele).")
    if len(ap.class1_pan_allele_models) != 1:
        # Standard optimized predictor has exactly one (merged) entry. If it
        # were unoptimized we would still handle it (list of separate models).
        pass

    model = ap.class1_pan_allele_models[0]
    net = model.network()

    # Collect the list of sub-networks (Class1NeuralNetworkModel) whose outputs
    # are combined by the ensemble geometric mean.
    subnets = []
    if hasattr(net, "networks"):          # MergedClass1NeuralNetwork
        subnets = list(net.networks)
    else:                                  # single, unmerged network
        subnets = [net]
        for extra in ap.class1_pan_allele_models[1:]:
            subnets.append(extra.network())

    subnet_weights = [_extract_subnet(s) for s in subnets]

    cache = {
        "model": model,                    # for peptides_to_network_input / allele input
        "subnets": subnet_weights,
        "allele_vec_cache": {},            # allele name -> (777,) float32 vector
    }
    _CACHE[ap] = cache
    return cache


def _get_cache(ap):
    cache = _CACHE.get(ap)
    if cache is None:
        cache = _build_cache(ap)
    return cache


def _allele_vector(ap, cache, allele):
    """Return the flattened BLOSUM62 allele-pseudosequence vector, exactly as
    used by mhcflurry's predict path (canonicalize -> AlleleEncoding -> compact)."""
    if allele in cache["allele_vec_cache"]:
        return cache["allele_vec_cache"][allele]

    model = cache["model"]
    normalized = ap.canonicalize_allele_name(allele)
    allele_encoding = AlleleEncoding(
        [normalized], borrow_from=ap.master_allele_encoding).compact()
    indices, reps = model.allele_encoding_to_network_input(allele_encoding)
    reshaped = reps.reshape((reps.shape[0], -1)).astype(np.float32)
    vec = reshaped[np.asarray(indices)][0]     # (embed_dim,) float32
    cache["allele_vec_cache"][allele] = vec
    return vec


# ----------------------------------------------------------------------------
# Forward pass
# ----------------------------------------------------------------------------
def _compute_merged(sub, peptide_enc, allele_mat):
    """Peptide branch + allele branch + merge (+ merge activation).

    Returns the merged feature matrix (batch, merged_size) that feeds the main
    dense stack. Mirrors the first half of Class1NeuralNetworkModel.forward.
    """
    act = _get_activation(sub["activation"])

    # ---- peptide branch ----
    x = peptide_enc
    for lc in sub["lc_layers"]:
        x = _locally_connected(x, lc)
    x = x.reshape(x.shape[0], -1)
    for (W, b) in sub["peptide_dense"]:
        x = _linear(x, W, b)
        if act is not None:
            x = act(x)
    if sub["batch_norm_early"] is not None:
        x = _batch_norm(x, sub["batch_norm_early"])

    # ---- allele branch + merge ----
    if sub["has_allele"]:
        a = allele_mat
        for (W, b) in sub["allele_dense"]:
            a = _linear(a, W, b)
            if act is not None:
                a = act(a)
        if sub["merge_method"] == "concatenate":
            x = np.concatenate([x, a], axis=-1)
        elif sub["merge_method"] == "multiply":
            x = x * a
        else:
            raise ValueError("Unknown merge method: %s" % sub["merge_method"])
        merge_act = _get_activation(sub["merge_activation"])
        if merge_act is not None:
            x = merge_act(x)
    return np.ascontiguousarray(x, dtype=np.float32)


def _forward_from_merged(sub, merged):
    """Main dense stack (feedforward / DenseNet skip) + output. Mirrors the
    second half of Class1NeuralNetworkModel.forward. Returns (batch,) float32."""
    act = _get_activation(sub["activation"])
    x = merged
    merged_input = merged
    prev_outputs = []
    topology = sub["topology"]
    for i, (W, b, bn) in enumerate(sub["dense"]):
        if topology == "with-skip-connections" and i > 0:
            if i == 1:
                x = np.concatenate([merged_input, prev_outputs[-1]], axis=-1)
            else:
                x = np.concatenate([prev_outputs[-2], prev_outputs[-1]], axis=-1)
        x = _linear(x, W, b)
        if act is not None:
            x = act(x)
        if bn is not None:
            x = _batch_norm(x, bn)
        # dropout is identity at inference
        prev_outputs.append(x)

    Wo, bo = sub["output"]
    out = _linear(x, Wo, bo)
    out_act = _get_activation(sub["output_activation"])
    if out_act is not None:
        out = out_act(out)
    return out[:, 0]


def _forward_subnet(sub, peptide_enc, allele_mat):
    """Full faithful forward pass for one sub-network."""
    merged = _compute_merged(sub, peptide_enc, allele_mat)
    return _forward_from_merged(sub, merged)


def _is_simple_merge(sub):
    """True if the sub-network's peptide/allele preprocessing reduces to
    concat(flatten(peptide), allele_vector), i.e. no locally-connected /
    peptide-dense / early-batchnorm / allele-dense layers and a plain
    concatenate merge with no merge activation. In that case the merged input
    is identical across all such sub-networks and can be computed once."""
    return (
        not sub["lc_layers"]
        and not sub["peptide_dense"]
        and sub["batch_norm_early"] is None
        and not sub["allele_dense"]
        and sub["has_allele"]
        and sub["merge_method"] == "concatenate"
        and not sub["merge_activation"]
    )


def _predict_ic50_matrix(cache, peptide_enc, allele_vec):
    """Return (n, num_subnets) float64 nM affinities from all sub-networks."""
    n = peptide_enc.shape[0]
    subnets = cache["subnets"]
    num = len(subnets)
    ic50 = np.empty((n, num), dtype=np.float64)

    # Fast path: when every sub-network shares the trivial peptide/allele
    # preprocessing (the case for the standard pan-allele ensemble), the merged
    # input is identical across sub-networks, so compute it once per chunk.
    shared_merge = all(_is_simple_merge(s) for s in subnets)

    for start in range(0, n, _FORWARD_CHUNK):
        end = min(start + _FORWARD_CHUNK, n)
        m = end - start
        pep = np.ascontiguousarray(peptide_enc[start:end], dtype=np.float32)

        if shared_merge:
            peptide_flat = pep.reshape(m, -1)
            allele_block = np.broadcast_to(allele_vec, (m, allele_vec.shape[0]))
            merged = np.concatenate([peptide_flat, allele_block], axis=1)
            merged = np.ascontiguousarray(merged, dtype=np.float32)
            for j, sub in enumerate(subnets):
                out = _forward_from_merged(sub, merged)
                ic50[start:end, j] = to_ic50(out.astype(np.float64))
        else:
            allele_mat = np.ascontiguousarray(
                np.broadcast_to(allele_vec, (m, allele_vec.shape[0])),
                dtype=np.float32)
            for j, sub in enumerate(subnets):
                out = _forward_subnet(sub, pep, allele_mat)
                ic50[start:end, j] = to_ic50(out.astype(np.float64))
    return ic50


# ----------------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------------
def predict_affinity_fast(predictor, peptides, allele):
    """
    Fast pure-numpy replacement for ``predictor.predict_affinity`` restricted to
    a single allele.

    Parameters
    ----------
    predictor : Class1PresentationPredictor or Class1AffinityPredictor
        A loaded mhcflurry predictor. Weights are extracted and cached on first
        use.
    peptides : sequence of str
        Peptide sequences (e.g. 8-15mers; validation done on 9/10mers).
    allele : str
        A single MHC class I allele name, e.g. "HLA-A*02:01".

    Returns
    -------
    numpy.ndarray of float64, shape (len(peptides),)
        Predicted binding affinities in nM, aligned to ``peptides`` and
        identical (to within float tolerance) to
        ``predictor.predict_affinity(peptides, alleles={'s': [allele]},
        include_affinity_percentile=False).affinity``.
    """
    ap = _resolve_affinity_predictor(predictor)
    cache = _get_cache(ap)
    model = cache["model"]

    encodable = EncodableSequences.create(list(peptides))
    if len(encodable.sequences) == 0:
        return np.zeros(0, dtype=np.float64)

    # Reuse mhcflurry's own numpy peptide encoder.
    peptide_enc = np.asarray(
        model.peptides_to_network_input(encodable), dtype=np.float32)

    allele_vec = _allele_vector(ap, cache, allele)

    ic50 = _predict_ic50_matrix(cache, peptide_enc, allele_vec)

    # Ensemble combination: geometric mean in log space, exactly as
    # Class1AffinityPredictor.predict_to_dataframe (centrality "mean").
    logs = np.log(ic50)
    log_centers = CENTRALITY_MEASURES["mean"](logs)
    affinities = np.exp(log_centers)
    return affinities
