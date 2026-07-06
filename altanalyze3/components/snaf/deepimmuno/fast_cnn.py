'''
Pure-NumPy (BLAS-backed) reimplementation of the DeepImmuno ``seperateCNN`` forward pass.

This module reproduces, bit-for-bit within float32 rounding, the immunogenicity scores
produced by ``deepimmuno.file_process`` -- but WITHOUT running TensorFlow at inference
time. TensorFlow / Keras is used exactly once, on the first call, only to read the trained
checkpoint and hand us the raw weight arrays; every subsequent prediction is a handful of
``numpy`` matmuls (Conv2D, fused BatchNorm, ReLU, MaxPool, Dense) that run on BLAS.
Dropout is inactive at inference (identity) and is therefore omitted.

Two exact optimisations make it faster than the Keras predict it replaces:
  * every conv either spans the full input width or has width 1, so each Conv2D is a sum
    of height-shifted BLAS matmuls (no large im2col patch tensor is ever materialised);
  * the two branches are pure per-sample functions with inference BatchNorm (no batch
    coupling), so branch outputs are computed once per UNIQUE peptide / allele and gathered
    -- bit-identical to per-row scoring, and dramatically cheaper when alleles repeat.

The public entry point mirrors ``deepimmuno.file_process``:

    >>> from fast_cnn import predict_immunogenicity_fast
    >>> out = predict_immunogenicity_fast(df)   # df: 2 columns [peptide, HLA]
    >>> out['immunogenicity']                   # sigmoid outputs, one per row

Weights and the PCA / HLA encoding tables are loaded once and cached for the process.
This file does not modify deepimmuno.py or any other existing file.
'''

import os
import numpy as np
import pandas as pd

# Import the vendored deepimmuno module for the (unchanged) encoding functions and the
# Keras model builder. Works both as a package member and when the folder is on sys.path.
try:                          # imported as part of the ``deepimmuno`` package
    from . import deepimmuno as _di
except ImportError:           # imported flat (folder on sys.path)
    import deepimmuno as _di


_BASE = os.path.dirname(os.path.abspath(__file__))
_STATE = None   # cached (_Weights, after_pca, hla_dic, dic_inventory) for this process


# --------------------------------------------------------------------------- #
#  Weight extraction (Keras used ONCE here, then never again)
# --------------------------------------------------------------------------- #
# Each of the four Conv2D kernels has a UNIQUE shape, so we identify convs by shape
# rather than by Keras layer name (names carry a process-global counter and are not
# stable if another seperateCNN was already built in this interpreter).
_CONV_ROLE_BY_KSHAPE = {
    (2, 12, 1, 16): 'pep_conv0',   # peptide branch, first conv
    (2, 1, 16, 32): 'pep_conv1',   # peptide branch, second conv
    (15, 12, 1, 16): 'hla_conv0',  # HLA branch, first conv
    (9, 1, 16, 32): 'hla_conv1',   # HLA branch, second conv
}


def _producer_layer(layer):
    '''Return the layer feeding `layer`'s first input (its graph predecessor).'''
    inbound = layer.inbound_nodes[0].inbound_layers
    if isinstance(inbound, (list, tuple)):
        inbound = inbound[0]
    return inbound


def _extract_weights():
    '''Build the Keras model, load the checkpoint, and pull every layer's arrays as
    float32 -- keyed by role via structure/connectivity, not by (unstable) layer names.

    BatchNormalization is folded into the exact (scale, offset) form that TensorFlow's
    fused inference kernel itself uses:
        scale  = gamma / sqrt(moving_var + eps)
        offset = beta  - moving_mean * scale
        y      = x * scale + offset
    which is algebraically and numerically identical to Keras BN inference. Each BN is
    mapped to its branch by looking up the Conv2D that directly feeds it.'''
    os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
    os.environ.setdefault('TF_USE_LEGACY_KERAS', '1')

    model = _di.seperateCNN()
    model.load_weights(os.path.join(_BASE, 'models', 'cnn_model_331_3_7', ''))

    W = {}
    conv_role_by_name = {}
    bn_layers = []
    for layer in model.layers:
        cls = layer.__class__.__name__
        if cls == 'Conv2D':
            kernel, bias = layer.get_weights()
            role = _CONV_ROLE_BY_KSHAPE[tuple(kernel.shape)]
            conv_role_by_name[layer.name] = role
            W[role] = (kernel.astype(np.float32), bias.astype(np.float32))
        elif cls == 'Dense':
            w, b = layer.get_weights()
            role = 'dense0' if w.shape[1] == 128 else 'dense1'  # (256,128) vs (128,1)
            W[role] = (w.astype(np.float32), b.astype(np.float32))
        elif cls == 'BatchNormalization':
            bn_layers.append(layer)

    for layer in bn_layers:
        conv_role = conv_role_by_name[_producer_layer(layer).name]  # BN's input is a Conv2D
        bn_role = conv_role.replace('conv', 'bn')                   # pep_conv0 -> pep_bn0
        gamma, beta, mean, var = layer.get_weights()
        eps = np.float32(layer.epsilon)
        scale = (gamma / np.sqrt(var + eps)).astype(np.float32)
        offset = (beta - mean * scale).astype(np.float32)
        W[bn_role] = (scale, offset)

    return W


def _load_state():
    global _STATE
    if _STATE is not None:
        return _STATE
    after_pca = np.loadtxt(os.path.join(_BASE, 'data', 'after_pca.txt'))
    hla = pd.read_csv(os.path.join(_BASE, 'data', 'hla2paratopeTable_aligned.txt'), sep='\t')
    hla_dic = _di.hla_df_to_dic(hla)
    dic_inventory = _di.dict_inventory(list(hla_dic.keys()))
    W = _extract_weights()
    _STATE = (W, after_pca, hla_dic, dic_inventory)
    return _STATE


# --------------------------------------------------------------------------- #
#  NumPy forward-pass primitives (channels-last, valid padding, float32)
# --------------------------------------------------------------------------- #
def _conv2d_valid(x, kernel, bias):
    '''2-D cross-correlation with 'valid' padding, matching Keras Conv2D.
    x:      (B, H, W, Cin)
    kernel: (kh, kw, Cin, Cout)   bias: (Cout,)

    Every conv in seperateCNN either spans the full input width (kw == W, so the output
    width Wo == 1) or has kw == 1 -- i.e. the kernel never slides along the width axis.
    In that case the whole (kw, Cin) block is a fixed inner-product dimension and the conv
    reduces to a sum of `kh` height-shifted BLAS matmuls, which avoids materialising the
    large im2col patch tensor (the single biggest cost). A general im2col fallback covers
    the Wo > 1 case for completeness.'''
    B, H, Wd, Cin = x.shape
    kh, kw, _, Cout = kernel.shape
    Ho, Wo = H - kh + 1, Wd - kw + 1

    if Wo == 1:                       # width fully spanned (kw == Wd): no sliding on width
        xf = x[:, :, :kw, :].reshape(B, H, kw * Cin)     # (B, H, kw*Cin), C-order [w, cin]
        Kf = kernel.reshape(kh, kw * Cin, Cout)          # matches xf inner ordering
        out = xf[:, 0:Ho, :] @ Kf[0]
        for di in range(1, kh):
            out += xf[:, di:di + Ho, :] @ Kf[di]
        out += bias
        return out.reshape(B, Ho, 1, Cout)

    # general fallback (not exercised by this network): im2col + one matmul
    patches = np.empty((B, Ho, Wo, kh, kw, Cin), dtype=np.float32)
    for di in range(kh):
        for dj in range(kw):
            patches[:, :, :, di, dj, :] = x[:, di:di + Ho, dj:dj + Wo, :]
    patches = patches.reshape(B * Ho * Wo, kh * kw * Cin)
    out = patches @ kernel.reshape(kh * kw * Cin, Cout) + bias
    return out.reshape(B, Ho, Wo, Cout)


def _batchnorm(x, scale_offset):
    scale, offset = scale_offset
    return x * scale + offset            # broadcast over last (channel) axis


def _relu(x):
    return np.maximum(x, np.float32(0))


def _maxpool_2x1(x):
    '''MaxPool2D(pool=(2,1), strides=(2,1)), valid padding: max over adjacent height
    pairs; odd remainder rows are dropped (Keras 'valid' behaviour).'''
    B, H, Wd, C = x.shape
    Ho = H // 2
    return x[:, :Ho * 2].reshape(B, Ho, 2, Wd, C).max(axis=2)


def _peptide_branch(input1, W):
    '''Peptide branch: (B,10,12,1) -> (B,128) flattened embedding.'''
    x = np.ascontiguousarray(input1, dtype=np.float32)
    x = _relu(_batchnorm(_conv2d_valid(x, *W['pep_conv0']), W['pep_bn0']))
    x = _relu(_batchnorm(_conv2d_valid(x, *W['pep_conv1']), W['pep_bn1']))
    x = _maxpool_2x1(x)
    return x.reshape(x.shape[0], -1)          # Flatten (C-order) -> 128


def _hla_branch(input2, W):
    '''HLA branch: (B,46,12,1) -> (B,128) flattened embedding.'''
    y = np.ascontiguousarray(input2, dtype=np.float32)
    y = _relu(_batchnorm(_conv2d_valid(y, *W['hla_conv0']), W['hla_bn0']))
    y = _maxpool_2x1(y)
    y = _relu(_batchnorm(_conv2d_valid(y, *W['hla_conv1']), W['hla_bn1']))
    y = _maxpool_2x1(y)
    return y.reshape(y.shape[0], -1)          # Flatten (C-order) -> 128


def _head(pep_emb, hla_emb, W):
    '''concatenate([peptide, HLA]) -> Dense(128,relu) -> Dropout(id) -> Dense(1,sigmoid).'''
    z = np.concatenate([pep_emb, hla_emb], axis=1)   # (B,256), peptide first (seperateCNN order)
    w0, b0 = W['dense0']
    z = _relu(z @ w0 + b0)
    w1, b1 = W['dense1']
    z = (z @ w1 + b1).astype(np.float32)
    return (np.float32(1.0) / (np.float32(1.0) + np.exp(-z))).ravel()


def _forward(input1, input2, W):
    '''Full per-row batch forward pass (no de-duplication). input1: (B,10,12,1),
    input2: (B,46,12,1). Returns sigmoid immunogenicity of shape (B,). This is the
    apples-to-apples equivalent of Keras `model.predict([input1, input2])`.'''
    return _head(_peptide_branch(input1, W), _hla_branch(input2, W), W)


# --------------------------------------------------------------------------- #
#  Public API
# --------------------------------------------------------------------------- #
def _encode_peptides(peptides, after_pca):
    out = np.empty((len(peptides), 10, 12, 1), dtype=np.float64)
    for k, p in enumerate(peptides):
        out[k] = _di.peptide_data_aaindex(p, after_pca)
    return out


def _encode_hlas(hlas, after_pca, hla_dic, dic_inventory):
    out = np.empty((len(hlas), 46, 12, 1), dtype=np.float64)
    for k, h in enumerate(hlas):
        out[k] = _di.hla_data_aaindex(hla_dic, h, after_pca, dic_inventory)
    return out


def predict_immunogenicity_fast(df_input, chunk_size=100000):
    '''Drop-in replacement for deepimmuno.file_process.

    df_input : DataFrame whose first two columns are (peptide, HLA). As in file_process,
               the columns are renamed to ['peptide', 'HLA'] and an 'immunogenicity'
               column (the sigmoid output) is added; the same DataFrame object is returned.

    Both CNN branches are pure per-sample functions -- and BatchNormalization at inference
    uses fixed moving statistics (no coupling across the batch) -- so the peptide branch
    output depends only on the peptide string and the HLA branch output only on the allele.
    We therefore run each branch once per UNIQUE peptide / allele and gather the results,
    which is bit-for-bit identical to scoring every row independently. On real workloads
    (a peptide scored against a handful of shared alleles) this is dramatically faster.

    chunk_size : rows scored per head block (bounds peak memory only).
    '''
    W, after_pca, hla_dic, dic_inventory = _load_state()

    ori = df_input
    ori.columns = ['peptide', 'HLA']
    peps = ori['peptide'].to_numpy()
    hlas = ori['HLA'].to_numpy()

    # unique peptides / alleles + inverse index back to every row
    uniq_pep, pep_inv = np.unique(peps, return_inverse=True)
    uniq_hla, hla_inv = np.unique(hlas, return_inverse=True)

    # HLA branch: always a tiny number of distinct alleles
    hla_emb = _hla_branch(_encode_hlas(uniq_hla, after_pca, hla_dic, dic_inventory), W)

    # peptide branch: encode + run in chunks over unique peptides (bounds memory)
    pep_emb = np.empty((len(uniq_pep), 128), dtype=np.float32)
    for start in range(0, len(uniq_pep), chunk_size):
        stop = min(start + chunk_size, len(uniq_pep))
        pep_emb[start:stop] = _peptide_branch(
            _encode_peptides(uniq_pep[start:stop], after_pca), W)

    # gather per-row embeddings + run the (cheap) head, chunked over rows
    n = len(peps)
    scores = np.empty(n, dtype=np.float32)
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        scores[start:stop] = _head(pep_emb[pep_inv[start:stop]],
                                    hla_emb[hla_inv[start:stop]], W)

    ori['immunogenicity'] = scores.astype(float)
    return ori
