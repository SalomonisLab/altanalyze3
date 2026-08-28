# Validation of the scALABLE-viewer chat

Measured 2026-08-26 against the running viewer on port 8062, dataset
`COPD-metacells` (123,076 cells, 39 cell states, 178 donors).

## Result

| suite | asked | passed | failed | median | slowest |
|---|---|---|---|---|---|
| protocol examples | 55 | **55** | 0 | 0.030 s | 1.250 s |
| paraphrase variants | 68 | **68** | 0 | 0.034 s | 1.329 s |
| contrast labels | 15 | **15** | 0 | — | — |

The variant suite scores every question four ways: as written, reworded,
with the gene or cell state swapped, and with the clauses reordered.

| phrasing | passed |
|---|---|
| canonical | 17 of 17 |
| paraphrase | 17 of 17 |
| gene-swap | 17 of 17 |
| reorder | 17 of 17 |

Reproduce the two HTTP suites:

    cd /Users/saljh8/Dropbox/LungMAP/refactored_website
    .venv/bin/python tests/test_viewer_chat_live.py
    .venv/bin/python tests/test_viewer_chat_variants.py

Both need the viewer on 8062 and the LungMAP site on 8001.

Reproduce the contrast-label suite, which reads the bundle and needs no server.
The site venv lacks pandas, so run it under the viewer's own interpreter:

    /opt/homebrew/opt/python@3.11/bin/python3.11 tests/test_contrast_labels.py

## What a passing answer means

The chat returns one of four outcomes, and never a neighbouring analysis:

| status | count of 68 | meaning |
|---|---|---|
| `answered` | 59 | a table of real rows from this dataset |
| `use_existing_view` | 8 | the analysis exists in another tab, and the answer names it |
| `not_covered` | 1 | the analysis exists but does not cover the cell state asked about |

The suite fails any question that returns a different protocol's answer. One
question carried the original defect: "Show me the transcriptional targets of
RUNX1 in alveolar macrophages" returned the marker genes of alveolar
macrophages.

15 of the 17 protocols now answer inline. `regulatory_driver` and
`communication_rewiring` still route to an existing tab rather than computing
an answer, and both name the tab they send the user to.

## The model never produces a number

The router chooses which of 17 protocols to run and copies slot values from
lists the viewer sends. The viewer then runs that analysis against the bundle.
Every figure and every statistic comes from the data.

In the 68-question variant sweep the router called the model **0 times**: the
keyword table decided all 68. The model stays as a fallback for phrasings the
table misses, and the router caps it at 6 seconds.

## Timing, and what made it fast

Median 0.034 s across 68 questions, slowest 1.329 s, and 6 of 68 exceed one
second. The six are the per-donor correlation protocols, which read the
expression matrix.

Four changes account for the speed. Each was a correctness fix first and a
latency fix second, because a missing cue meant falling through to the model:

| change | before | after |
|---|---|---|
| cues for "differ", "separate them" | 20.7 s | 0.01 s |
| cues for "only some", "all donors" | 15.1 s | 0.01 s |
| skip the model when the sentence names nothing in the dataset | 11.2 s | 0.01 s |
| decide an absent modality before the model runs | 11.5 s | 0.01 s |

A fifth change bounds the worst case: the router caps the model fallback at 6 s,
so a phrasing nobody anticipated degrades to a fast clarification instead of a
twenty-second wait.

## Defects this validation found

Routing tests that called the classifier directly passed 44 of 44 while the
viewer was still wrong, because the viewer runs in a different process and was
wired to a retired set of four intents. Only the running viewer, asked over
HTTP, exposed the six defects below:

1. **The viewer was never connected to the protocol router.** Its chat still
   dispatched on `markers`, `expression`, `differential`, `compare`. A protocol
   name it did not recognise fell past every branch and the nearest one
   answered.
2. **The intent endpoint dropped `covariates` and `modalities`,** so every
   severity, stratification and composition question could not fill its slot.
3. **No cue matched "transcriptional targets".** The sentence scored zero, fell
   through to the model, and the model chose `contrast_specificity`.
4. **The marker table covers 14 of 39 cell states.** DC1, DC2, NK cells and Mast
   cells returned an empty table. Markers are now computed as a one-versus-rest
   contrast on the per-state means where the table is silent, and each row
   carries a `source` field saying which produced it.
5. **The COPD differential covers 18 of 39 cell states.** Asking about Mast
   cells returned an empty table, which reads as "nothing changes here". That is
   a different and wrong claim, so the answer now names the gap and lists the
   states the contrast does cover.
6. **Two protocols were missing.** An earlier pass dropped `cell_identity` and
   `state_comparison` as too simple. Users ask both first of any atlas, so both
   are back.
7. **`bundle_meta._norm_label` confused three GOLD labels.** The function
   deleted every non-alphanumeric character, so `GOLD_I_II`, `GOLD I, II` and
   `GOLD III` all read as `goldiii`. `_contrast_group_field` then matched
   GOLD III as the control group of `GOLD_IV_vs_GOLD_I_II`. The bug reached
   past the chat: `build_differential_block` wrote
   `group2_samples: ["GOLD III"]` into the Differential tab's config for that
   comparison, so the tab named donors the differential never compared.
   `_label_tokens` now keeps separators as word boundaries, which separates
   ('gold','i','ii') from ('gold','iii') and still equates `no_cancer` with
   `no cancer`. The fix changes 2 of the bundle's 8 contrasts and leaves the
   other 6 identical. `tests/test_contrast_labels.py` fails against the old
   key, and restoring the old matcher brings both halves of the bug back.
8. **The frequency figure re-guessed which group was the case.** The renderer
   tested the first group's name against `control|non|healthy|normal|never`.
   The test matched nothing in "GOLD I, II" against "GOLD IV", so the milder
   stage took the red bar. The server now orders the groups reference first,
   from the contrast's own `control_label`, and the renderer trusts that order.

## Scope and limits

- **One dataset.** `COPD-metacells` only. The suites cover no other bundle.
- **RNA only.** This bundle carries no imputed layer, so questions about ADT,
  lipid, metabolite or GRN return `unsupported` naming the layers that exist.
  The suites route those questions but never check their answers, because the
  bundle holds no such layer.
- **No browser.** All three suites drive code, not pixels. The two HTTP suites
  call the endpoint the chat panel calls, and the label suite reads the bundle.
  Nobody clicked the rendered panel, so a figure that draws wrongly from a
  correct payload would still pass.
- **The variant set is 68 questions.** It is not a sample of real user
  language; it is four phrasings each of seventeen questions, written by the
  same person who wrote the cues. A pass here shows the router is not brittle to
  wording. It does not show that the seventeen protocols are the right
  seventeen.
- **The frequency figure draws 14 cell states.** `most_affected_state` ranks up
  to 39 and the figure annotates how many of them it shows.
- **`most_affected_state` answers with two measures.** The table ranks states by
  how many genes clear FDR 0.05, and the bars show the same states' abundance in
  each group. A state can change in expression, in abundance, or in both, and
  the two confound each other: a depleted state looks changed in any pooled
  expression comparison.
- **Composition substitutes a categorical variable for a numeric one.** The
  GOLD IV composition question routes the slot to `gold_ordinal`, which is
  numeric. The executor falls back to `Group`, whose levels name the GOLD
  stages, and the answer says which variable it used.

## Files

- `/Users/saljh8/Dropbox/LungMAP/refactored_website/app/lungmap/assistant/viewer_protocols.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/app/lungmap/assistant/viewer.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/app/lungmap/web/data_routes.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/tests/test_viewer_chat_live.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/tests/test_viewer_chat_variants.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/tests/test_viewer_protocols.py`
- `/Users/saljh8/Dropbox/LungMAP/refactored_website/tests/test_contrast_labels.py`
- `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/scalable_app.py`
- `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/bundle_meta.py`
- `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/static/app.js`
- `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/templates/index.html`
