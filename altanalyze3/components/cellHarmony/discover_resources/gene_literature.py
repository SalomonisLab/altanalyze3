"""Native Discover rate-limited PubMed transport; no model-generated sources."""
import json, re, threading, urllib.parse, urllib.request
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TIMEOUT=8.0
_LOCK=threading.Lock()
_LAST=[0.0]

def _fetch(url: str) -> bytes:
    """One E-utilities call, spaced to NCBI's 3 requests a second without an API key.

    The first version sent seven calls back to back and PubMed answered HTTP 429 on the
    third, so the age theme arrived and the other two did not. One retry follows a 429.
    """
    import time
    import urllib.error

    for attempt in (0, 1, 2):
        # One lock across threads: the drug panel looks up several compounds at once,
        # and each thread spacing its own calls would still exceed 3 a second together.
        with _LOCK:
            wait = 0.36 - (time.monotonic() - _LAST[0])
            if wait > 0:
                time.sleep(wait)
            _LAST[0] = time.monotonic()
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as r:
                return r.read()
        except urllib.error.HTTPError as exc:
            # Several site workers share one IP and each has its own lock, so a 429 can
            # still arrive; two retries with a growing wait absorb it.
            if exc.code != 429 or attempt == 2:
                raise
            time.sleep(1.1 * (attempt + 1))
    return b'{}'


def _get(url: str) -> dict:
    return json.loads(_fetch(url).decode())


def pubmed_url(query: str) -> str:
    return "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.parse.quote(query)


def phrase_terms(*phrases: str) -> str:
    """An OR of quoted title/abstract phrases, for a theme built from the page's context."""
    seen, out = set(), []
    for p in phrases:
        p = re.sub(r"[\"()\[\]]", " ", p or "").strip()
        if len(p) >= 3 and p.lower() not in seen:
            seen.add(p.lower())
            out.append(f'"{p}"[tiab]')
    return " OR ".join(out)
