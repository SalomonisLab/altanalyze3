"""Unsupervised clustering utilities."""

from .hopach import HopachResult, distancematrix, hopach, labelstomss

__all__ = [
    "HopachResult",
    "ICGS3Config",
    "ICGS3Result",
    "distancematrix",
    "hopach",
    "labelstomss",
    "run_icgs3",
]


def __getattr__(name):
    if name in {"ICGS3Config", "ICGS3Result", "run_icgs3"}:
        from .ICGS import ICGS3Config, ICGS3Result, run_icgs3

        values = {
            "ICGS3Config": ICGS3Config,
            "ICGS3Result": ICGS3Result,
            "run_icgs3": run_icgs3,
        }
        return values[name]
    raise AttributeError(name)
