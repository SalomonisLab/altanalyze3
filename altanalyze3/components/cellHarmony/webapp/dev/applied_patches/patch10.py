P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:160]}"
    src = src.replace(old, new)

helpers = '''def _chat_states_in_question(question: str, states: List[str]) -> List[str]:
    """The cell states this sentence actually names, in the order it names them.

    The router matches a state name anywhere in the sentence, so "pre-aceNKP"
    also reports "aceNKP" and a two-state question came back comparing a state
    with itself. A match that sits inside a longer match is dropped here, which
    leaves the states the reader wrote.
    """
    text = str(question or "").lower()
    spans = []
    for state in sorted([str(s) for s in states if s], key=len, reverse=True):
        needle = state.lower()
        start = 0
        while True:
            at = text.find(needle, start)
            if at < 0:
                break
            if not any(begin <= at and at + len(needle) <= end for begin, end, _ in spans):
                spans.append((at, at + len(needle), state))
            start = at + 1
    spans.sort()
    out: List[str] = []
    for _, _, state in spans:
        if state not in out:
            out.append(state)
    return out


#: Words that are also gene symbols in some annotations. The question scan below
#: only runs when the router found no gene, and these would turn an ordinary
#: sentence into a gene lookup.
_CHAT_GENE_STOPWORDS = {
    "and", "are", "can", "cell", "cells", "for", "gene", "genes", "has", "how",
    "impact", "many", "max", "most", "not", "rest", "set", "she", "state",
    "states", "the", "was", "what", "when", "where", "which", "who", "why",
}


def _chat_genes_in_question(question: str, cache: Dict[str, Any]) -> List[str]:
    """Genes this sentence names, matched against the dataset's own gene list.

    The router returns an empty gene list for a plain sentence such as "Where is
    Cdca3 expressed?". Nothing is invented here: a token is only taken when the
    dataset holds a gene of that name.
    """
    index = {}
    for name in cache.get("var_names", []):
        index.setdefault(str(name).lower(), str(name))
    out: List[str] = []
    for token in re.findall(r"[A-Za-z0-9_.\\-]{3,}", str(question or "")):
        key = token.lower()
        if key in _CHAT_GENE_STOPWORDS:
            continue
        name = index.get(key)
        if name and name not in out:
            out.append(name)
    return out


'''
anchor = "def _chat_read_question("
assert src.count(anchor) == 1
src = src.replace(anchor, helpers + anchor)

# --- use both scans in the route ---------------------------------------------
sub('''        states = [s for s, _ in _chat_states_by_size(cache)]
        contrast = _chat_contrast(meta)''',
    '''        states = [s for s, _ in _chat_states_by_size(cache)]
        contrast = _chat_contrast(meta)

        # The router reads the sentence with a keyword matcher that overlaps
        # state names and drops gene names. Both are repaired against the names
        # this job actually holds, and the repaired reading is returned so the
        # answer says which state and which gene it used.
        named = _chat_states_in_question(question, states)
        if named:
            if state not in named:
                state = named[0]
            if state2 and state2 not in named:
                state2 = next((s for s in named if s != state), "")
            if not state2 and len(named) > 1:
                state2 = next((s for s in named if s != state), "")
        if not genes:
            genes = _chat_genes_in_question(question, cache)
        reading = dict(reading)
        reading["cell_state"] = state
        reading["cell_state_2"] = state2
        reading["genes"] = genes
        result["reading"] = reading''')

open(P, "w").write(src)
print("patch10 applied")
