'''
Deferred loading of scTriangulate's heavyweight optional dependencies.

Upstream scTriangulate 0.13.0 imports ``umap``, ``scrublet`` and ``squidpy`` at
module scope. On this Mac that makes the package unimportable and slow:

* ``import umap`` pulls in ``umap.parametric_umap``, which imports TensorFlow.
  TensorFlow built against a newer numpy raises
  ``AttributeError: module 'numpy' has no attribute 'dtypes'`` under numpy
  1.23.5, so ``import sctriangulate`` fails outright before any analysis runs.
* ``scrublet`` is only needed when ``predict_doublet=True``; the default path
  fills a constant 0.5 and never touches it.
* ``squidpy`` is only needed by the spatial helpers.

Each of the three is used in exactly one place. Loading them on first use keeps
every default code path working with none of them installed, and cuts package
import time. Nothing numeric changes.

If a deferred import fails, the error names the package and the feature that
wanted it, instead of surfacing as a TensorFlow traceback at import time.
'''

import importlib


class _Deferred:
    '''Attribute proxy that imports its target module on first attribute access.'''

    def __init__(self, module_name, needed_for):
        self._module_name = module_name
        self._needed_for = needed_for
        self._module = None

    def _load(self):
        if self._module is None:
            try:
                self._module = importlib.import_module(self._module_name)
            except Exception as e:
                raise ImportError(
                    "scTriangulate needs '{}' for {}, but importing it failed: {}: {}".format(
                        self._module_name, self._needed_for, type(e).__name__, e)) from e
        return self._module

    def __getattr__(self, name):
        return getattr(self._load(), name)


umap = _Deferred('umap.umap_', 'umap_dual_view_save / reference-mapping embeddings')
scrublet = _Deferred('scrublet', 'doublet_predict(predict_doublet=True)')
squidpy = _Deferred('squidpy', 'the spatial module')
