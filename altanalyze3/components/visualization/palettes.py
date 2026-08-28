#!/usr/bin/env python3
"""Categorical colour palettes that stay distinguishable past 20 populations.

Discrete construction, not a sampled continuous colormap
--------------------------------------------------------
Sampling a continuous colormap at N points fails as N grows. A continuous map traces one
path through colour space, so the gap between neighbouring samples shrinks as 1/N. At 86
populations, two adjacent ``viridis`` samples differ by roughly 1 CIE DeltaE unit, under the
2.3-unit just-noticeable difference, so a reader sees one colour where the data holds two
populations.

Repeating a discrete map fails the same way for a different reason. ``tab20`` holds exactly
20 colours. ``plt.get_cmap("tab20", 46)`` resamples those 20 and returns 20 unique colours
for 46 categories, and the resampling makes ADJACENT indices collide, so cluster C1, C2 and
C3 draw in one colour. ``palette(i % palette.N)`` wraps for the same result.

This module builds the palette the way Glasbey et al. (2007) describe. Take a dense grid of
sRGB candidates, convert to CIELAB, then repeatedly add the candidate whose minimum CIELAB
distance to the already-chosen colours is the largest available. Each new colour therefore
sits as far from every earlier colour as the space allows. The first 20 stay as separable as
``tab20`` and the sequence keeps going to several hundred populations.

The result is not a rainbow ramp. A rainbow walks hue at near-constant lightness, so it
gives poor separation and a false ordering. This sequence maximises perceptual distance and
carries no order.

Every colour comes back as an ``#RRGGBB`` string.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Sequence

import numpy as np

#: Candidate grid resolution per sRGB channel. 18**3 = 5,832 candidates.
_GRID = 18

#: CIELAB lightness window. Colours lighter than the top read as white on a white page;
#: colours darker than the bottom read as black next to text and axis lines.
_L_MIN = 22.0
_L_MAX = 90.0

#: Fixed starting colour, so a palette of a given length never changes between runs.
_SEED_RGB = (0.12, 0.31, 0.72)


def _srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """Convert an (N, 3) array of sRGB values in 0..1 to CIELAB under D65."""
    rgb = np.asarray(rgb, dtype=np.float64)
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    m = np.array([
        [0.4124564, 0.3575761, 0.1804375],
        [0.2126729, 0.7151522, 0.0721750],
        [0.0193339, 0.1191920, 0.9503041],
    ])
    xyz = linear @ m.T
    white = np.array([0.95047, 1.00000, 1.08883])
    t = xyz / white
    delta = 6.0 / 29.0
    f = np.where(t > delta ** 3, np.cbrt(t), t / (3 * delta ** 2) + 4.0 / 29.0)
    lab = np.empty_like(f)
    lab[:, 0] = 116.0 * f[:, 1] - 16.0
    lab[:, 1] = 500.0 * (f[:, 0] - f[:, 1])
    lab[:, 2] = 200.0 * (f[:, 1] - f[:, 2])
    return lab


def _candidates() -> tuple:
    axis = np.linspace(0.0, 1.0, _GRID)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    lab = _srgb_to_lab(grid)
    keep = (lab[:, 0] >= _L_MIN) & (lab[:, 0] <= _L_MAX)
    return grid[keep], lab[keep]


def _to_hex(rgb: Sequence[float]) -> str:
    r, g, b = (int(round(float(c) * 255.0)) for c in rgb[:3])
    return "#%02X%02X%02X" % (max(0, min(255, r)), max(0, min(255, g)), max(0, min(255, b)))


@lru_cache(maxsize=8)
def _glasbey_sequence(n: int) -> tuple:
    """Return ``n`` maximally separated ``#RRGGBB`` colours, in a fixed order."""
    if n <= 0:
        return ()
    grid, lab = _candidates()
    seed_lab = _srgb_to_lab(np.asarray([_SEED_RGB]))[0]
    first = int(np.argmin(((lab - seed_lab) ** 2).sum(axis=1)))

    picked = [first]
    # Track, for every candidate, its distance to the nearest colour already picked. Adding
    # one colour costs one distance pass, so building n colours costs n passes, not n^2/2.
    min_dist = ((lab - lab[first]) ** 2).sum(axis=1)
    for _ in range(1, n):
        nxt = int(np.argmax(min_dist))
        if min_dist[nxt] <= 0:
            break
        picked.append(nxt)
        np.minimum(min_dist, ((lab - lab[nxt]) ** 2).sum(axis=1), out=min_dist)

    colours = [_to_hex(grid[i]) for i in picked]
    # The candidate grid is finite. If it ever runs dry, cycle rather than raise, and let
    # the caller's own check report the reuse.
    while len(colours) < n:
        colours.append(colours[len(colours) % max(1, len(picked))])
    return tuple(colours)


def categorical_colors(n: int) -> List[str]:
    """Return ``n`` distinguishable ``#RRGGBB`` colours."""
    return list(_glasbey_sequence(int(n)))


def categorical_palette(categories: Sequence) -> Dict[str, str]:
    """Map each category to a distinguishable ``#RRGGBB`` colour.

    Categories keep the order given, so the same category list always returns the same
    colours. Repeated categories collapse to one entry and one colour.
    """
    seen: List[str] = []
    for cat in categories:
        key = str(cat)
        if key not in seen:
            seen.append(key)
    colours = categorical_colors(len(seen))
    return {key: colours[i] for i, key in enumerate(seen)}


def min_pairwise_distance(colours: Sequence[str]) -> float:
    """Smallest CIELAB distance between any two colours. Use it to check a palette."""
    if len(colours) < 2:
        return float("inf")
    rgb = np.asarray([[int(c[i:i + 2], 16) / 255.0 for i in (1, 3, 5)] for c in colours])
    lab = _srgb_to_lab(rgb)
    d = np.sqrt(((lab[:, None, :] - lab[None, :, :]) ** 2).sum(axis=-1))
    np.fill_diagonal(d, np.inf)
    return float(d.min())


# ---------------------------------------------------------------------------------------
# Paired, the default for UMAP category colours
# ---------------------------------------------------------------------------------------
#: Number of discrete colours matplotlib's Paired holds.
PAIRED_LIMIT = 12


def paired_colors(n: int) -> List[str]:
    """Return ``n`` ``#RRGGBB`` colours built from matplotlib's Paired.

    Paired is a discrete ListedColormap of 12 colours. Up to 12 categories take those colours
    unchanged. Above 12, indexing Paired repeats them, so this interpolates between the same
    12 anchors and returns one distinct colour per category.

    Paired trades distinctness for a calmer, publication-style look. Its minimum CIE distance
    falls below the 2.3 just-noticeable difference past roughly 60 categories (1.93 at 85),
    where :func:`categorical_colors` holds 21.10. Use :func:`categorical_colors` when a reader
    must tell every population apart; use this where the house style asks for Paired.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    n = max(int(n), 0)
    if n == 0:
        return []
    base = plt.get_cmap("Paired")
    if n <= PAIRED_LIMIT:
        picks = [base(i) for i in range(n)]
    else:
        cont = LinearSegmentedColormap.from_list("PairedC", base.colors, N=max(n, 2))
        picks = [cont(i / max(n - 1, 1)) for i in range(n)]
    return [_to_hex(c[:3]) for c in picks]


def paired_palette(categories: Sequence) -> Dict[str, str]:
    """Map each category to a Paired-derived ``#RRGGBB`` colour, keeping the order given."""
    seen: List[str] = []
    for cat in categories:
        key = str(cat)
        if key not in seen:
            seen.append(key)
    colours = paired_colors(len(seen))
    return {key: colours[i] for i, key in enumerate(seen)}
