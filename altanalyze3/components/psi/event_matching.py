"""Exact, oriented event-arm joins; annotation only, with no event collapsing.

Sub-exons, introns and novel-coordinate suffixes are literal identities. Sharing
one endpoint or an exon block is insufficient. An event's first arm is its PSI
numerator; a second-arm match reverses its directional interpretation.
"""
import re

_GENE = re.compile(r"ENS[A-Z]*G\d+(?:\.\d+)?:")


def event_junctions(label):
    """Return ordered gene-qualified endpoint pairs, retaining all literal arms.

    Accepts a leading display symbol and repeated or omitted same-gene prefixes.
    Explicit cross-gene junction endpoints remain distinct. No region rounding,
    coordinate stripping, genomic direction reversal or pair deduplication occurs.
    """
    text = str(label).strip()
    first = _GENE.search(text)
    if first is None:
        raise ValueError(f"Event has no explicit Ensembl gene identity: {label}")
    text = text[first.start():]
    default_gene = text.split(":", 1)[0]
    result = []
    for arm in text.split("|"):
        endpoints = arm.split("-")
        if len(endpoints) != 2:
            raise ValueError(f"Expected two junction endpoints: {arm}")
        resolved = []
        gene = default_gene
        for endpoint in endpoints:
            if ":" in endpoint:
                gene, region = endpoint.split(":", 1)
            else:
                region = endpoint
            if not _GENE.fullmatch(gene + ":") or not re.fullmatch(
                r"[EIU]\d+(?:\.\d+)?(?:_\d+)?", region
            ):
                raise ValueError(f"Unresolved junction endpoint: {endpoint}")
            resolved.append(gene + ":" + region)
        result.append(tuple(resolved))
    if len(result) != 2:
        raise ValueError(f"Expected two event arms: {label}")
    return tuple(result)


def junction_matches(left, right):
    """All shared-arm evidence, with a sign multiplier relative to left arm 1.

    The factor is +1 for 1:1 or 2:2, -1 for 1:2 or 2:1. A single matching
    junction suffices. Duplicate/reversed event records remain the caller's rows.
    """
    return [dict(left_arm=i + 1, right_arm=j + 1,
                 junction="-".join(a), direction_multiplier=1 if i == j else -1)
            for i, a in enumerate(left) for j, b in enumerate(right) if a == b]
