"""Numerical runtime settings applied before an isolated worker imports NumPy."""
from __future__ import annotations

import sys
from collections.abc import Mapping


def worker_environment(environ: Mapping[str, str], *, platform: str | None = None) -> dict[str, str]:
    env = dict(environ)
    if (sys.platform if platform is None else platform) == "darwin":
        # NumPy and SciPy can load separate OpenBLAS libraries on macOS. We
        # observed SciPy's LU initialization deadlock in blas_thread_init after
        # the graph stage initialized the other thread pool. Set this before
        # package imports; threadpoolctl after import cannot repair that lock.
        # Only BLAS changes: Numba/OpenMP and analysis settings are retained.
        for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
            env[key] = "1"
    return env
