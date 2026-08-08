"""Import shim that lets the pristine upstream sctriangulate 0.13.0 import on this Mac.

Upstream `sctriangulate/preprocessing.py` does a module-level `import umap`.
`umap/__init__.py` imports `umap.parametric_umap`, which imports TensorFlow, which
fails against numpy 1.23.5 with `AttributeError: module 'numpy' has no attribute 'dtypes'`.

This shim installs a lazy `umap` package object that resolves public attributes from
`umap.umap_` only when they are touched. It changes no numeric behaviour; it only
defers the TensorFlow import that upstream never needs.

Import this module BEFORE importing sctriangulate.
"""

import sys
import types
import importlib
import importlib.util


class _LazyUmap(types.ModuleType):
    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        real = importlib.import_module('umap.umap_')
        return getattr(real, name)


def install():
    if 'umap' in sys.modules:
        return
    spec = importlib.util.find_spec('umap')
    if spec is None:
        return
    mod = _LazyUmap('umap')
    mod.__spec__ = spec
    mod.__loader__ = spec.loader
    mod.__path__ = list(spec.submodule_search_locations or [])
    sys.modules['umap'] = mod


install()
