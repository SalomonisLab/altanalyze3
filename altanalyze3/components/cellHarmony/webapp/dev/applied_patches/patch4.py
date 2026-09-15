P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

old = '''    adata = cache["adata"]
    cluster_key = cache["cluster_key"]
    stored = adata.uns.get(f"{cluster_key}_colors")
    order = _state_order(cache)
    lookup = {}
    if stored is not None and len(stored) == len(order):
        lookup = {s: str(c) for s, c in zip(order, stored)}
    return [lookup.get(s, "#BBBBBB") for s in states]'''

new = '''    adata = cache["adata"]
    cluster_key = cache["cluster_key"]
    stored = adata.uns.get(f"{cluster_key}_colors")
    order = _state_order(cache)
    lookup = {}
    if stored is not None and len(stored) == len(order):
        lookup = {s: str(c) for s, c in zip(order, stored)}
    else:
        # cellHarmony output h5ads carry no `<cluster_key>_colors`, so every bar
        # of the CombPlot came back grey and its state blocks were unreadable.
        # The fallback is the same Paired palette the front end draws cell states
        # with, assigned by position in the dataset's own state order, so a state
        # keeps one colour across every figure.
        lookup = {state: _PAIRED_HEX[index % len(_PAIRED_HEX)]
                  for index, state in enumerate(order)}
    return [lookup.get(s, "#BBBBBB") for s in states]'''
assert src.count(old) == 1
src = src.replace(old, new)

anchor = 'def _state_colors(cache: Dict[str, Any], states: List[str]) -> List[str]:'
palette = '''#: The 12 Paired colours the Explore panel already uses for cell states, as hex.
#: Kept identical to PAIRED_COLOR_STOPS in static/app.js so a state is the same
#: colour whichever side of the app drew it.
_PAIRED_HEX = [
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
]


'''
assert src.count(anchor) == 1
src = src.replace(anchor, palette + anchor)
open(P, "w").write(src)
print("patch4 applied")
