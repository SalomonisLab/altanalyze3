# Stored cell-state biology

`lungmap.json` is a versioned offline digest of the LungMAP controlled vocabulary,
its published cell-state mapping, and CellCards. Source URLs, retrieval date and
SHA-256 checksums are embedded. It contains brief factual source excerpts with
species qualifiers, linked literature and published marker citations. It is not
a systematic literature review. Vocabulary reference atlas counts are excluded;
only the uploaded job supplies dataset counts.

The summary uses exact published aliases (allowing a terminal cluster suffix and
`Cell(s)` suffix). Ambiguous/unmatched optional background entries remain explicit;
no dataset state, cell, sample or model feature is discarded or reannotated.
Lung-specific background is not applied to a non-lung reference.

Refresh explicitly from the repository root:

```sh
python altanalyze3/components/cellHarmony/resources/biology/build_lungmap.py
```

Review the resulting diff before committing. Job requests never scrape these sites.
