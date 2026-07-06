"""
fast_infer_processing
=====================

Pure-NumPy (BLAS-backed) re-implementation of the forward pass of mhcflurry's
antigen-processing neural network ensemble
(``Class1ProcessingPredictor`` -> ``Class1ProcessingModel``).

Motivation
----------
In SNAF, ``Class1PresentationPredictor.predict_processing`` accounts for the
large majority of the MHC-binding runtime.  mhcflurry 2.2.x runs the network
with PyTorch (one small CNN per peptide, ensemble of 8).  The forward pass is
tiny arithmetically but the PyTorch/eager overhead dominates.  Running the whole
batch at once as a handful of NumPy ``matmul`` calls (which dispatch to BLAS)
gives an order-of-magnitude speedup while producing *numerically identical*
output (float32 arithmetic; differences at the level of float32 rounding).

Public API
----------
``predict_processing_fast(predictor, peptides)`` -> numpy.ndarray

Returns the same array as::

    predictor.predict_processing(peptides, verbose=0)

where ``predictor`` is an mhcflurry ``Class1PresentationPredictor`` (it uses the
``processing_predictor_without_flanks`` ensemble, which is what SNAF hits).  A
``Class1ProcessingPredictor`` may also be passed directly.

The network weights are extracted from the passed predictor and cached on first
call, so repeated calls only pay for encoding + the NumPy forward pass.

This module reuses mhcflurry's own peptide encoding (``FlankingEncoding`` /
``network_input``) verbatim -- encoding is not the bottleneck and reusing it
guarantees bit-identical inputs.  Only the neural-network forward pass is
re-implemented here.

Notes on numerical equivalence
------------------------------
* All matmuls are done in float32, matching PyTorch's CPU float32 conv/linear.
* Each model's scalar output is cast to float64 and the 8-model ensemble is
  averaged in float64 -- exactly as ``Class1ProcessingPredictor`` does.
* Conv1d uses cross-correlation with symmetric 'same' padding (all production
  kernel sizes are odd), matching PyTorch's ``padding='same'``.
"""

import numpy as np

try:  # NumPy >= 1.20
    from numpy.lib.stride_tricks import sliding_window_view as _sliding_window_view
except Exception:  # pragma: no cover
    _sliding_window_view = None


# Module-level cache of extracted weights, keyed by id() of the processing
# predictor (the Class1ProcessingPredictor ensemble object).
_WEIGHT_CACHE = {}


# --------------------------------------------------------------------------- #
# Weight extraction
# --------------------------------------------------------------------------- #
def _extract_model_params(net):
    """
    Extract NumPy weights + config from a single Class1ProcessingNeuralNetwork.

    Returns a dict of float32 arrays / scalars describing the forward pass.
    """
    model = net.network()  # the torch nn.Module (Class1ProcessingModel)
    hp = net.hyperparameters

    def npf(t):
        return t.detach().cpu().numpy().astype(np.float32)

    conv_w = npf(model.conv1.weight)            # (out_ch, in_ch, k)
    conv_b = npf(model.conv1.bias)              # (out_ch,)
    k = conv_w.shape[2]
    # PyTorch padding='same', stride 1, dilation 1: total pad = k - 1.
    # All production kernels are odd -> symmetric pad on both sides.
    total_pad = k - 1
    pad_left = total_pad // 2
    pad_right = total_pad - pad_left

    # Flatten conv weight to (in_ch*k, out_ch) with in_ch-major, k-minor order
    # to match the im2col patch layout below.
    out_ch, in_ch, _ = conv_w.shape
    conv_w_flat = conv_w.reshape(out_ch, in_ch * k).T.copy()  # (in_ch*k, out_ch)

    def post_layers(module_list):
        layers = []
        for conv in module_list:
            w = npf(conv.weight)   # (out, in, 1)
            b = npf(conv.bias)     # (out,)
            w2 = w[:, :, 0].T.copy()  # (in, out)  -> x@(in,out)
            layers.append((w2, b))
        return layers

    n_post = post_layers(model.n_flank_post_convs)
    c_post = post_layers(model.c_flank_post_convs)

    out_w = npf(model.output_layer.weight)[0].copy()  # (4,)
    out_b = npf(model.output_layer.bias)[0].copy()     # scalar (as 0-d)

    activation = str(hp["convolutional_activation"])

    return {
        "conv_w_flat": conv_w_flat,     # (in_ch*k, out_ch)
        "conv_b": conv_b,               # (out_ch,)
        "k": int(k),
        "in_ch": int(in_ch),
        "out_ch": int(out_ch),
        "pad_left": int(pad_left),
        "pad_right": int(pad_right),
        "act": activation,
        "n_post": n_post,
        "c_post": c_post,
        "out_w": out_w,                 # (4,)
        "out_b": np.float32(out_b),
        "n_flank_length": int(hp["n_flank_length"]),
        "c_flank_length": int(hp["c_flank_length"]),
        "peptide_max_length": int(hp["peptide_max_length"]),
    }


def _get_processing_ensemble(predictor):
    """
    Resolve the Class1ProcessingPredictor (ensemble) from whatever the caller
    passed.

    Accepts a Class1PresentationPredictor (uses its
    ``processing_predictor_without_flanks``) or a Class1ProcessingPredictor
    directly.
    """
    proc = getattr(predictor, "processing_predictor_without_flanks", None)
    if proc is not None:
        return proc
    if hasattr(predictor, "models"):
        return predictor
    raise ValueError(
        "Passed predictor is neither a Class1PresentationPredictor nor a "
        "Class1ProcessingPredictor.")


def _get_cached_params(proc):
    key = id(proc)
    cached = _WEIGHT_CACHE.get(key)
    if cached is None:
        cached = [_extract_model_params(net) for net in proc.models]
        _WEIGHT_CACHE[key] = cached
    return cached


# --------------------------------------------------------------------------- #
# NumPy forward pass primitives
# --------------------------------------------------------------------------- #
def _activation(x, kind):
    if kind == "relu":
        return np.maximum(x, np.float32(0.0))
    if kind == "tanh":
        return np.tanh(x)
    if kind == "sigmoid":
        return np.float32(1.0) / (np.float32(1.0) + np.exp(-x))
    # mhcflurry default fallback is tanh
    return np.tanh(x)


def _im2col_windows(xp, k, out_L=None):
    """
    xp : (batch, L_pad, in_ch) float32
    out_L : if given, only build the first ``out_L`` output positions.
    returns (batch, out_L, in_ch*k) float32 with in_ch-major, k-minor layout.
    """
    if _sliding_window_view is not None:
        # sliding along axis=1 -> (batch, full_out_L, in_ch, k)
        sw = _sliding_window_view(xp, window_shape=k, axis=1)
        if out_L is not None and out_L < sw.shape[1]:
            sw = sw[:, :out_L, :, :]
        batch, ol, in_ch, kk = sw.shape
        return np.ascontiguousarray(sw).reshape(batch, ol, in_ch * kk)
    # Fallback (older numpy): manual gather.
    batch, L_pad, in_ch = xp.shape
    full_out_L = L_pad - k + 1
    ol = full_out_L if out_L is None else min(out_L, full_out_L)
    cols = np.empty((batch, ol, in_ch, k), dtype=np.float32)
    for j in range(k):
        cols[:, :, :, j] = xp[:, j:j + ol, :]
    return cols.reshape(batch, ol, in_ch * k)


def _conv1d_same(x, m, out_L_needed=None):
    """
    Main convolution with 'same' padding.

    x : (batch, L, in_ch) float32
    out_L_needed : if given, compute only the first ``out_L_needed`` output
        positions.  The forward pass only ever reads positions 0..max(plen)-1
        (everything beyond is masked out / never indexed), so restricting the
        output length to max(plen) is exactly equivalent and saves work.
    returns (batch, out_L, out_ch) float32 (pre-activation).
    """
    batch, L, in_ch = x.shape
    k = m["k"]
    pad_l = m["pad_left"]
    pad_r = m["pad_right"]

    if pad_l == 0 and pad_r == 0:
        xp = x
    else:
        xp = np.zeros((batch, L + pad_l + pad_r, in_ch), dtype=np.float32)
        xp[:, pad_l:pad_l + L, :] = x

    cols = _im2col_windows(xp, k, out_L=out_L_needed)  # (batch, out_L, in_ch*k)
    out_L = cols.shape[1]
    # Collapse batch/position axes into ONE large BLAS gemm
    #   (batch*out_L, in_ch*k) @ (in_ch*k, out_ch)
    # rather than `batch` tiny batched gemms.  The 3D `cols @ weight` form makes
    # NumPy issue one gemm per batch element, which is ~30x slower here.
    cols2d = cols.reshape(batch * out_L, in_ch * k)
    out = cols2d @ m["conv_w_flat"]        # (batch*out_L, out_ch)
    out += m["conv_b"]                     # broadcast over (out_ch,)
    return out.reshape(batch, out_L, m["out_ch"])


def _post_convs(conv_result, layers, act):
    """
    Apply the per-flank 1x1 post-conv stack.

    conv_result : (batch, L, filters) float32
    layers      : list of (W (in,out), b (out,))
    Returns (batch, L) float32  -- the final single-channel per-position output.
    """
    batch, L, C = conv_result.shape
    # Reshape to 2D so each 1x1 conv is a single large gemm
    #   (batch*L, C) @ (C, out)
    # instead of a per-position batched gemm.
    x = conv_result.reshape(batch * L, C)
    n = len(layers)
    for i, (w2, b2) in enumerate(layers):
        x = x @ w2 + b2                    # (batch*L, out)
        if i < n - 1:
            x = _activation(x, act)
        else:
            x = np.tanh(x)                 # final post-conv layer is always tanh
    return x.reshape(batch, L)             # final channel width is 1


def _neg_masked_max(single, starts, ends, L):
    """
    Replicates mhcflurry's masked max-pool:
        masked = (single + 1) * mask
        result = -( max_over_positions(masked) - 1 )

    single : (batch, L) float32
    starts : scalar int (inclusive lower bound)
    ends   : (batch,) int array (exclusive upper bound)
    """
    pos = np.arange(L, dtype=np.int64)[None, :]           # (1, L)
    mask = (pos >= starts) & (pos < ends[:, None])        # (batch, L) bool
    shifted = single + np.float32(1.0)
    masked = shifted * mask.astype(np.float32)
    mx = masked.max(axis=1) - np.float32(1.0)
    return (-mx).astype(np.float32)


def _forward_model(m, seq, plen, chunk_size=20000):
    """
    Full NumPy forward pass for one ensemble member.

    seq  : (N, L, in_ch) float32   (mhcflurry-encoded sequence array)
    plen : (N,) int                (peptide lengths)
    Returns (N,) float64 processing scores.
    """
    N, L, _ = seq.shape
    nfl = m["n_flank_length"]
    cfl = m["c_flank_length"]
    out = np.empty(N, dtype=np.float64)

    for s in range(0, N, chunk_size):
        e = min(s + chunk_size, N)
        x = seq[s:e]                       # (b, L, in_ch)
        pl = plen[s:e].astype(np.int64)    # (b,)
        b = x.shape[0]

        # Only positions 0..max(plen)-1 (= n_flank_length + max(plen) - 1 in the
        # general case) are ever read downstream; positions past that are always
        # masked out or never indexed.  Restricting the conv output length there
        # is exactly equivalent and cuts the dominant matmul proportionally.
        out_L = int(nfl + pl.max())        # positions 0..nfl+max(plen)-1
        out_L = min(out_L, L)

        conv = _conv1d_same(x, m, out_L_needed=out_L)   # (b, out_L, filters)
        conv = _activation(conv, m["act"])

        n_single = _post_convs(conv, m["n_post"], m["act"])  # (b, out_L)
        c_single = _post_convs(conv, m["c_post"], m["act"])  # (b, out_L)

        rows = np.arange(b)

        # --- n_flank features ---
        n_cleaved = n_single[:, nfl]                                       # (b,)
        n_maxpool = _neg_masked_max(n_single, nfl + 1, nfl + pl, out_L)    # (b,)

        # --- c_flank features ---
        c_idx = nfl + pl - 1
        c_cleaved = c_single[rows, c_idx]                                 # (b,)
        c_maxpool = _neg_masked_max(c_single, nfl, nfl + pl - 1, out_L)    # (b,)

        feats = np.stack(
            [n_cleaved, n_maxpool, c_cleaved, c_maxpool], axis=1
        ).astype(np.float32)               # (b, 4)

        z = feats @ m["out_w"] + m["out_b"]                # (b,)
        y = np.float32(1.0) / (np.float32(1.0) + np.exp(-z))
        out[s:e] = y.astype(np.float64)

    return out


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def predict_processing_fast(predictor, peptides, chunk_size=20000):
    """
    Fast NumPy equivalent of ``predictor.predict_processing(peptides, verbose=0)``.

    Parameters
    ----------
    predictor : mhcflurry Class1PresentationPredictor (or Class1ProcessingPredictor)
        The processing weights are extracted and cached from this object on the
        first call.
    peptides : list of str
        Peptide sequences (no flanks; N-/C-flanks are treated as empty, matching
        ``predict_processing`` with ``n_flanks=None``).
    chunk_size : int
        Row-batch size for the NumPy forward pass (bounds peak memory; does not
        affect the numerical result).

    Returns
    -------
    numpy.ndarray
        Processing scores, one per peptide (float64), ensemble-averaged over the
        8 networks -- identical to mhcflurry's output.
    """
    from mhcflurry.flanking_encoding import FlankingEncoding

    peptides = list(peptides)
    if len(peptides) == 0:
        return np.array([], dtype=float)

    proc = _get_processing_ensemble(predictor)
    params = _get_cached_params(proc)

    # Reuse mhcflurry's encoding verbatim (cached inside the FlankingEncoding
    # instance, so it is computed once and shared across all ensemble members).
    empty = [""] * len(peptides)
    sequences = FlankingEncoding(peptides=peptides, n_flanks=empty, c_flanks=empty)

    # Encode ONCE per distinct encoding config and reuse across the whole
    # ensemble.  All production networks share peptide_max_length / n_flank /
    # c_flank, so this normally encodes a single time.  (Encoding is cheap --
    # ~0.4 ms per 100 peptides -- but re-`asarray`-ing it 8x is pure waste.)
    enc_cache = {}

    def _encoded_for(net):
        hp = net.hyperparameters
        key = (
            hp["amino_acid_encoding"],
            hp["peptide_max_length"],
            hp["n_flank_length"],
            hp["c_flank_length"],
        )
        cached = enc_cache.get(key)
        if cached is None:
            x_dict = net.network_input(sequences, throw=True)
            cached = (
                np.asarray(x_dict["sequence"], dtype=np.float32),
                np.asarray(x_dict["peptide_length"]),
            )
            enc_cache[key] = cached
        return cached

    acc = None
    for net, m in zip(proc.models, params):
        seq, plen = _encoded_for(net)
        scores = _forward_model(m, seq, plen, chunk_size=chunk_size)
        acc = scores if acc is None else (acc + scores)

    return acc / len(params)
