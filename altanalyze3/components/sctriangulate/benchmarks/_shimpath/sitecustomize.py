"""Installed on PYTHONPATH for the upstream benchmark runs only.

Python imports `sitecustomize` automatically at interpreter startup, including
in multiprocessing 'spawn' children. That is the only hook that reaches those
children, and upstream scTriangulate needs it: a spawned worker re-imports
`sctriangulate`, which does `import umap`, which imports TensorFlow, which
fails under numpy 1.23.5. Without this the workers die during bootstrap and the
parent's pool.join() waits forever.

This exists so the upstream-vs-optimized timing compares the parallel DESIGNS
rather than just re-measuring the import failure.
"""
import sys, types, importlib, importlib.util


class _LazyUmap(types.ModuleType):
    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return getattr(importlib.import_module('umap.umap_'), name)


if 'umap' not in sys.modules:
    _spec = importlib.util.find_spec('umap')
    if _spec is not None:
        _m = _LazyUmap('umap')
        _m.__spec__ = _spec
        _m.__loader__ = _spec.loader
        _m.__path__ = list(_spec.submodule_search_locations or [])
        sys.modules['umap'] = _m
