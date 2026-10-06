"""Write per-cell predictions incrementally in an ordinary, portable H5AD."""
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
from types import SimpleNamespace

from ..imputed_scale import normalize_scale, prediction_encoding, inverse_predictions, log_base_value


class PredictionWriter:
    def __init__(self, path, obs, var, uns, obsm, *, compression='lzf', expression_scale='linear', log_base=None):
        uns = dict(uns)
        self.expression_scale = normalize_scale(expression_scale)
        if 'expression_scale' in uns and normalize_scale(uns['expression_scale']) != self.expression_scale:
            raise ValueError('Conflicting prediction encodings in writer and metadata')
        uns['expression_scale'] = self.expression_scale
        if log_base is not None:
            if 'log_base' in uns and not np.isclose(log_base_value(uns['log_base'], default=np.e),
                                                   log_base_value(log_base, default=np.e)):
                raise ValueError('Conflicting declared log bases in writer and metadata')
            uns['log_base'] = log_base
        self.encoding = (None if self.expression_scale == 'log2' and 'log_pseudocount' not in uns
                         else prediction_encoding(SimpleNamespace(uns=uns)))
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        shape = (len(obs), len(var))
        skeleton = ad.AnnData(X=sp.csr_matrix(shape, dtype=np.float32),
                             obs=obs.copy(), var=var.copy(), uns=uns, obsm=obsm)
        skeleton.write_h5ad(self.path, compression=compression)
        self.handle = h5py.File(self.path, 'r+')
        del self.handle['X']
        self.x = self.handle.create_dataset('X', shape=shape, dtype=np.float32,
                                            chunks=(min(256, max(1, shape[0])), min(512, max(1, shape[1]))),
                                            compression=compression)
        self.x.attrs.update({'encoding-type': 'array', 'encoding-version': '0.2.0'})
        self.counts = None
        self.log_base = log_base
        self.position = 0

    def _linear(self, values):
        linear = inverse_predictions(values, self.encoding).astype(np.float32)
        if not np.isfinite(linear).all() or (self.encoding[1] is not None and self.encoding[2] == 0 and np.any(linear <= 0)):
            raise ValueError('Abundance export overflow/underflow; do not clip or fill predictions')
        return linear

    def append(self, values):
        values = np.asarray(values, dtype=np.float32)
        if self.expression_scale == 'native_relative_log2':
            if not np.isfinite(values).all():
                raise ValueError('Nonfinite native-log2 predictions; do not silently fill values')
        else:
            values = np.nan_to_num(values, copy=True)
        end = self.position + len(values)
        if values.shape[1] != self.x.shape[1] or end > self.x.shape[0]:
            raise ValueError('Prediction block does not fit the output matrix')
        self.x[self.position:end] = values
        if self.encoding is not None and self.counts is None and self.expression_scale != 'linear':
            self.counts = self.handle.create_dataset('layers/counts', shape=self.x.shape,
                                                     dtype=np.float32, chunks=self.x.chunks,
                                                     compression=self.x.compression)
            self.counts.attrs.update({'encoding-type': 'array', 'encoding-version': '0.2.0'})
            for start in range(0, self.position, 256):
                stop = min(start + 256, self.position)
                self.counts[start:stop] = self._linear(self.x[start:stop])
        if self.counts is not None:
            self.counts[self.position:end] = self._linear(values)
        self.position = end

    def close(self):
        try:
            if self.position != self.x.shape[0]:
                raise ValueError('Prediction output is incomplete')
            if self.counts is None and self.encoding is not None:
                # Linear predictions (including signed scores) and their layer are identical.
                # A hard link saves a second full disk matrix; ordinary H5AD
                # readers still see both X and layers/counts.
                self.handle['layers/counts'] = self.x
        finally:
            self.handle.close()

    def abort(self):
        self.handle.close()
        self.path.unlink(missing_ok=True)
