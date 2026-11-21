"""Vector fishplot rendering; no inference or VAF conversion.

Adapted from Xuan's 10-17-25_fishplot/fishplot_python.py Matplotlib
fill_between rendering (source SHA256
9fe0533595cc55d3c989a45065a384cdfc7ddf6732f6f96fb9817185c7e5684f).
Its VAF values, assumed hierarchy, CCF conversion, clipping and
renormalization are deliberately not part of this renderer.
"""
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from .integrity import require

PALETTE = ['#e41a1c', '#377eb8', '#08519c', '#984ea3', '#4daf4a',
           '#f781bf', '#ff7f00', '#bcbd22', '#a65628', '#00a6a6',
           '#756bb1', '#636363', '#66c2a5', '#e6ab02', '#8c510a']


def fish_coordinates(F, nodes, parents, points=121):
    """Exact sample knots plus a separately flagged, synthetic ancestral lead-in.

    In the lead-in, each diagnosis clone's mass moves down its *displayed*
    ancestral path. This is illustration, not an estimate of historical burden.
    Clinical intervals use convex interpolation, avoiding spline overshoot.
    """
    F = np.asarray(F, float)
    require(F.shape[0] == len(nodes) and F.shape[1] > 0 and
            np.isfinite(F).all() and (F >= -1e-9).all() and
            np.allclose(F.sum(0), 100), 'Invalid exclusive clone frequencies')
    idx = {n: i for i, n in enumerate(nodes)}
    paths = {}
    for n in nodes:
        if n <= 0:
            continue
        path = [n]
        while path[-1] != 0:
            path.append(int(parents[path[-1]]))
        paths[n] = [a for a in path[::-1] if a in idx]
    steps = max([len(p)-1 for p in paths.values()] + [1])
    # One smooth expansion per lineage, with staggered illustrative onsets.
    # Restarting an easing curve at several synthetic knots creates scallops.
    u = np.linspace(0, 1, max(points*2, 241))
    def activation(level):
        start = (level-1)/(steps+1)
        t = np.clip((u-start)/(1-start), 0, 1)
        return t*t*(3-2*t)
    growth = activation(1)
    synthetic = np.zeros((len(nodes), len(u)))
    synthetic[idx[0]] = 100*(1-growth)+growth*F[idx[0], 0]
    synthetic[idx[-1]] = growth*F[idx[-1], 0]
    for n, path in paths.items():
        for level, target in enumerate(path[1:], 1):
            share = activation(level)
            if level < len(path)-1:
                share = share-activation(level+1)
            synthetic[idx[target]] += share*F[idx[n], 0]
    require(np.allclose(synthetic[:, -1], F[:, 0]), 'Lead-in changed diagnosis')
    xx = list(-1.25+1.25*u)
    values = [synthetic]
    for k in range(F.shape[1]-1):
        u = np.linspace(0, 1, points)
        w = (1-np.cos(np.pi*u))/2
        sl = slice(1, None)
        xx.extend((k+u)[sl])
        values.append((F[:, k, None]*(1-w)+F[:, k+1, None]*w)[:, sl])
    xx, V = np.asarray(xx), np.concatenate(values, axis=1)
    require(np.allclose(V.sum(0), 100) and (V >= -1e-8).all(), 'Fish mass loss')
    for k in range(F.shape[1]):
        require(np.allclose(V[:, np.flatnonzero(np.isclose(xx, k))[0]], F[:, k]),
                'Changed measured frequency')
    return xx, V


def draw_clone_fish(ax, F, nodes, parents, colors, stage_labels):
    xx, V = fish_coordinates(F, nodes, parents)
    exclusive = dict(zip(nodes, V))
    children = {n: [i for i, p in enumerate(parents) if i and p == n]
                for n in range(len(parents))}
    lineage = {}
    def mass(n):
        lineage[n] = exclusive.get(n, np.zeros_like(xx)) + sum(
            (mass(ch) for ch in children[n]), np.zeros_like(xx))
        return lineage[n]
    mass(0)
    U = exclusive[-1]
    require(np.allclose(lineage[0]+U, 100), 'Tree mass loss')
    # Vector polygons keep the reference cream/gold/brown background editable.
    cmap = LinearSegmentedColormap.from_list('reference_WT', ['#fbe8c8', '#f6cf70', '#c88b52'])
    cuts = np.linspace(xx[0], xx[-1], 257)
    quads = [[(a, 0), (b, 0), (b, 100), (a, 100)] for a, b in zip(cuts[:-1], cuts[1:])]
    ax.add_collection(PolyCollection(quads, facecolors=cmap(np.linspace(0, 1, 256)),
                                     edgecolors='none'))
    ax.fill_between(xx, 0, U/2, color=colors[-1], linewidth=0)
    ax.fill_between(xx, 100-U/2, 100, color=colors[-1], linewidth=0)
    def band(n, lo):
        if np.max(lineage[n]) > 1e-9:
            # Zero-mass latent nodes must not draw artificial lines through a clone.
            ax.fill_between(xx, lo, lo+lineage[n], color=colors.get(n, '#dddddd'), linewidth=0)
        childlo = lo+exclusive.get(n, np.zeros_like(xx))/2
        for ch in children[n]:
            band(ch, childlo)
            childlo = childlo+lineage[ch]
        require((childlo <= lo+lineage[n]+1e-7).all(), 'Child outside ancestor')
    lo = U/2+exclusive[0]/2
    for child in children[0]:
        band(child, lo)
        lo = lo+lineage[child]
    for k in range(F.shape[1]):
        ax.axvline(k, color='#aaaaaa', lw=.55, ls='--' if k == 0 else '-', alpha=.8)
    ax.set_xlim(xx[0], xx[-1]); ax.set_ylim(0, 100)
    ax.set_xticks(np.arange(F.shape[1]), stage_labels)
    ax.xaxis.tick_top(); ax.tick_params(axis='x', length=0, labelsize=12, pad=9)
    ax.set_yticks([0, 50, 100]); ax.set_ylabel('% of all profiled cells', fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(-.625, 103, 'Prediagnosis\n(schematic)', ha='center', va='bottom', fontsize=10)
    return xx, V
