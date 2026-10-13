"""Refresh the versioned, offline LungMAP biology digest (explicit maintainer command).

No dataset is read here. Alias joins use the published vocabulary mapping, never
fuzzy matching. CellCards text is retained as brief source excerpts, with its
species qualifiers and citation URLs; it is background, not job evidence.
"""
import csv
import hashlib
import io
import json
import re
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import urlopen

SOURCES = {
    "vocabulary": "https://www.lungmap.net/vocabulary/",
    "mapping": "https://www.lungmap.net/research/vocabulary/downloads/cell_state_mapping",
    "cards": "https://www.lungmap.net/research/cell-cards/cellcards/cell_cards.json",
}

class PlainText(HTMLParser):
    def __init__(self):
        super().__init__(); self.parts = []; self.references = []
    def handle_starttag(self, tag, attrs):
        if tag in {"br", "p", "li"}: self.parts.append(" ")
        if tag == "a":
            url = dict(attrs).get("href", "")
            if url.startswith(("https://doi.org/", "https://pubmed.ncbi.nlm.nih.gov/", "https://www.nature.com/")) and " " not in url:
                self.references.append(url)
    def handle_data(self, data): self.parts.append(data)


def build(raw):
    cards = []
    for card in json.loads(raw["cards"]):
        identity = re.search(r"LMCC\d+", card.get("ID", ""))
        if not identity: continue
        facts, references = {}, set()
        for field in ("Function", "Developmental_origin", "Heterogeneity_and_cellular_states", "Regenerative_potential"):
            parser = PlainText(); parser.feed(" ".join(str(s) for s in card.get(field, []) if s))
            words = " ".join(parser.parts).split()
            if words: facts[field] = {"excerpt": " ".join(words[:65]), "truncated": len(words)>65,
                                      "species_scope": "as stated in excerpt; do not transfer mouse findings to human"}
            references.update(parser.references)
        publications = []
        for entry in card.get("publication_markers", []):
            url = entry.get("dataset_publication")
            if url and url.startswith("https://"):
                item = {key: entry.get(key, "") for key in ("dataset_citation", "dataset_publication", "dataset_assay", "species", "dataset_cellname")}
                if item not in publications: publications.append(item)
        names = [str(v) for key in ("Label", "display_name", "Other names") for v in card.get(key, []) if v]
        cards.append({"id": identity[0], "names": names, "facts": facts,
                      "literature": sorted(references), "marker_publications": publications,
                      "human_present":card.get("human_present"), "mouse_present":card.get("mouse_present"),
                      "url": "https://www.lungmap.net/research/cell-cards/?cell_cards_id="+identity[0]})
    mapping = list(csv.DictReader(io.StringIO(raw["mapping"].decode('utf-8-sig', errors='replace')), delimiter='\t'))
    disease_section = raw['vocabulary'].decode().split('id="diseases"',1)[-1].split('</section>',1)[0]
    diseases=[]
    for row in re.findall(r'<tr[^>]*>(.*?)</tr>',disease_section,re.S):
        fields=re.findall(r'<td[^>]*>(.*?)</td>',row,re.S)
        if len(fields)<3: continue
        text=[]
        for field in fields[:3]:
            parser=PlainText(); parser.feed(field);text.append(' '.join(' '.join(parser.parts).split()))
        # Exclude the fourth column: these are reference atlas counts, not job data.
        diseases.append(dict(zip(('term','name','ontology'),text)))
    return {"schema_version":1, "retrieved_at":datetime.now(timezone.utc).isoformat(),
            "sources": {k:{"url":url, "sha256":hashlib.sha256(raw[k]).hexdigest()} for k,url in SOURCES.items()},
            "scope":"Lung reference background. Brief source excerpts, not a systematic literature review. Species qualifiers must be preserved. Missing aliases are reported; no state is discarded.",
            "cell_states": mapping, "cards":cards, "diseases":diseases}

if __name__ == '__main__':
    raw = {key:urlopen(url, timeout=45).read() for key,url in SOURCES.items()}
    target = Path(__file__).with_name('lungmap.json')
    target.write_text(json.dumps(build(raw), indent=2, ensure_ascii=False)+'\n')
    print(target)
