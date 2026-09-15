import io, sys
P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

def sub(old, new, count=1):
    global src
    n = src.count(old)
    assert n == count, f"expected {count} match, found {n} for:\n{old[:200]}"
    src = src.replace(old, new)

# --- 1. _group_axis: numpy array is not a truth value -------------------------
sub(
"""    elif column == cluster_key and cache.get("populations"):
        groups = [str(s) for s in cache["populations"]]""",
"""    elif column == cluster_key and len(cache.get("populations", [])):
        # `populations` is a numpy array. `array or []` calls bool() on it and
        # raises "truth value of an array ... is ambiguous", which took down
        # every DotPlot and CombPlot request.
        groups = list(dict.fromkeys(str(s) for s in cache["populations"]))""")

# --- 2. _state_colors: same truthiness bug, and the wrong ordering ------------
sub(
'''def _state_colors(cache: Dict[str, Any], states: List[str]) -> List[str]:
    """The colour each cell state is drawn in, falling back to a neutral grey."""
    adata = cache["adata"]
    cluster_key = cache["cluster_key"]
    stored = adata.uns.get(f"{cluster_key}_colors")
    order = [str(s) for s in cache.get("populations") or []]
    lookup = {}
    if stored is not None and len(stored) == len(order):
        lookup = {s: str(c) for s, c in zip(order, stored)}
    return [lookup.get(s, "#BBBBBB") for s in states]''',
'''def _state_colors(cache: Dict[str, Any], states: List[str]) -> List[str]:
    """The colour each cell state is drawn in, falling back to a neutral grey.

    `adata.uns[f"{cluster_key}_colors"]` holds one colour per cell-state
    CATEGORY. Two defects lived in the previous version. It read
    `cache["populations"]`, which is one label per CELL, so the length test
    never matched and every state came back grey. And it wrote
    `cache.get("populations") or []`, which calls bool() on a numpy array and
    raises, so both gene-set figures returned HTTP 500 on every request.
    """
    adata = cache["adata"]
    cluster_key = cache["cluster_key"]
    stored = adata.uns.get(f"{cluster_key}_colors")
    order = _state_order(cache)
    lookup = {}
    if stored is not None and len(stored) == len(order):
        lookup = {s: str(c) for s, c in zip(order, stored)}
    return [lookup.get(s, "#BBBBBB") for s in states]


def _state_order(cache: Dict[str, Any]) -> List[str]:
    """The cell-state categories, in the dataset's own order.

    The categorical order is the centroid order every other figure reads, so a
    dataset that stores the column as a plain string column keeps first-seen
    order rather than being re-sorted alphabetically here.
    """
    adata = cache["adata"]
    cluster_key = str(cache["cluster_key"])
    series = adata.obs[cluster_key] if cluster_key in adata.obs.columns else None
    if series is not None and str(series.dtype) == "category":
        return [str(c) for c in series.cat.categories]
    return list(dict.fromkeys(str(s) for s in cache.get("populations", [])))''')

open(P, "w").write(src)
print("patch1 applied")
