"""Fast receptor-ligand communication scoring grounded by receiver response."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .api import FastCommParams, FastCommResult, run_fastcomm

__all__ = ["FastCommParams", "FastCommResult", "run_fastcomm"]


def __getattr__(name: str):
    if name in __all__:
        from . import api

        return getattr(api, name)
    raise AttributeError(name)
