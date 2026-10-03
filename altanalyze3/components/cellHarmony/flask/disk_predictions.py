"""Write per-cell predictions incrementally in an ordinary, portable H5AD."""
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp


class PredictionWriter:
    def __init__(self, path, obs, var, uns, obsm, *, compression='lzf', expression_scale='linear', log_base=None):
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
        self.expression_scale, self.log_base = expression_scale, log_base
        self.position = 0

    def _linear(self, values):
        if self.expression_scale == 'log2':
            return np.maximum(np.exp2(values.astype(np.float64)) - 1, 0).astype(np.float32)
        if self.expression_scale == 'log1p':
            base = self.log_base
            base = np.e if base is None or str(base).lower() in ('e', 'ln', 'natural') else float(base)
            values64 = values.astype(np.float64)
            converted = np.expm1(values64) if np.isclose(base, np.e) else np.power(base, values64) - 1
            return np.maximum(converted, 0).astype(np.float32)
        return np.maximum(values, 0)

    def append(self, values):
        values = np.nan_to_num(np.asarray(values, dtype=np.float32), copy=True)
        end = self.position + len(values)
        if values.shape[1] != self.x.shape[1] or end > self.x.shape[0]:
            raise ValueError('Prediction block does not fit the output matrix')
        self.x[self.position:end] = values
        if self.counts is None and (self.expression_scale != 'linear' or np.any(values < 0)):
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
            if self.counts is None:
                # Nonnegative linear predictions and their counts are identical.
                # A hard link saves a second full disk matrix; ordinary H5AD
                # readers still see both X and layers/counts.
                self.handle['layers/counts'] = self.x
        finally:
            self.handle.close()

    def abort(self):
        self.handle.close()
        self.path.unlink(missing_ok=True)
