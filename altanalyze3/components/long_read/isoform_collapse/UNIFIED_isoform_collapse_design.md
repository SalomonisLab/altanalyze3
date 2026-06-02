# Unified Cross-Sample Long-Read Isoform Collapse — Authoritative Design

**Status:** APPROVED SPEC. Supersedes `DESIGN_cross_sample_collapse.md` (Design A) and `canonical_isoform_collapse_design.md` (Design B). Replaces the untested `collapse_transcript_associations.py`.
**Date:** 2026-05-30
**Scope:** Collapse raw long reads into a standardized, cross-dataset isoform catalog and re-key each sample's h5ad, at **≤200 samples × ≤~20M reads/sample**, without loading all samples' reads at once and without dropping rare-but-real novel isoforms. Structures observed exactly once across all samples are dropped.

> This document is the merge of two prior designs after a 5-way adversarial stress-test. Every decision below cites verified `file:line` evidence. Where the two prior docs contradicted each other, the contradiction and its resolution are stated explicitly.

---

## 0. The single source of truth (read this first)

> **`gff-output/candidate_to_final.tsv.gz` is the ONE authoritative collapse record.** It stores every collapse decision with full context. **`gff-output/final_isoforms.tsv.gz`** is the authoritative *winner* table (one row per final isoform). **Everything else is DERIVED** from these two — including the legacy `isoform_links.txt`, which becomes a deterministic, hash-verified VIEW. **No consumer may recompute a collapse decision independently of `candidate_to_final`.**

This rule exists because the stress-test found a real divergence risk: today the collapse decision lives in `gff_process` and the matrix builder `exportConsensusIsoformMatrix` *trusts* `isoform_links.txt` with zero validation, silently skipping any isoform not in the map (`if isoform in isoform_to_meta_isoform`, [isoform_ratios.py:92](../isoform_ratios.py#L92)). Two independent truths → silent count loss. One authoritative artifact + pure projections eliminates this class of bug.

---

## 1. Why the current path doesn't scale (verified)

Isoform consolidation runs through `gff_process.consolidateLongReadGFFs(...)` ([gff_process.py:670](../gff_process.py#L670)), called from `export_isoform_h5ad` ([isoform_automate.py:345](../isoform_automate.py#L345)). Three verified failure modes:

1. **All-reads-in-RAM.** Every transcript from every sample GFF accumulates into `junction_db[(gene, structure)] -> [(file, info), ...]` with one tuple **per read** ([gff_process.py:705-809](../gff_process.py#L705-L809)). Memory grows with total read count.
2. **O(k²) per-gene collapse** via Python substring containment (`if shorter_isoform in isoform`, [gff_process.py:585](../gff_process.py#L585)). This is also **incorrect**: string containment matches `E12|E23` inside `E1|E12|E23|E4` — a token-boundary bug (see §5).
3. **Dense matrix build.** The downstream matrix builder densifies via `cell_data = adata.X.toarray()` ([isoform_ratios.py:86](../isoform_ratios.py#L86)) then loops per isoform ([isoform_ratios.py:89-120](../isoform_ratios.py#L89-L120)) — ~10s of GB at 200 samples.

The **junction/splicing path is untouched** by this design; it runs first (`pre_process_samples` → junctions; `combine_processed_samples` → isoforms, [isoform_automate.py:598-618](../isoform_automate.py#L598-L618)).

**Structural insight:** the number of *distinct exon-structures* is bounded by biology, not read count. Collapse per-sample first; the cross-sample merge then operates on `n_samples × structures_per_gene` rows.

---

## 2. Output contract (non-negotiable)

`exportConsensusIsoformMatrix` reads `isoform_links.txt` and builds `source_isoform -> "gene:ref_isoform"` ([isoform_ratios.py:62-78](../isoform_ratios.py#L62-L78)). The **5-column order is hard-unpacked** at [isoform_ratios.py:66](../isoform_ratios.py#L66) and filtered by `gff_source` against **both** `source_gff` and `ref_gff` ([isoform_ratios.py:67,71](../isoform_ratios.py#L67-L71)):

```
gene <tab> ref_isoform <tab> ref_gff <tab> source_isoform <tab> source_gff
```

- `source_isoform` = the **bare molecule ID** as it appears in a sample h5ad `var_name` *after* stripping the `gene:` prefix ([isoform_structure_extract.py:454](../../bam/isoform_structure_extract.py#L454), `isoform_id = "GENE:molecule_id"`).
- A single-transcript winner writes a row with **empty source fields**: `gene, ref, ref_gff, '', ''` ([gff_process.py:965](../gff_process.py#L965)).
- Multi-member winners fan out one row per `(source_isoform, source_gff)` ([gff_process.py:968-984](../gff_process.py#L968-L984)).

`gff_translate` additionally needs `combined.gff` + `transcript_associations.txt` ([isoform_automate.py:348](../isoform_automate.py#L348)).

**Tier 2 must emit byte-compatible `isoform_links.txt`, `transcript_associations.txt`, `isoform_annotations.txt`, `combined.gff`** as DERIVED views. Then `export_isoform_h5ad` needs only its call site swapped; nothing downstream changes.

Also preserved verbatim from `gff_process`: the AltAnalyze structure string `E1.1|E2.2|E2.4|` (`exonAnnotate`+`exon_str`, [gff_process.py:772-788](../gff_process.py#L772-L788)); terminal-novel-coord stripping ([gff_process.py:420-431](../gff_process.py#L420-L431)); Ensembl/RefSeq priority ([gff_process.py:950-955](../gff_process.py#L950-L955)).

---

## 3. Architecture

```
  BAM ──isoform_structure_extract (UNCHANGED)──▶ sample.gff.gz + sample.h5ad
        (raw reads, 1 transcript/read;            (var = "GENE:molecule_id",
         chromosome-partitioned, position-order)   X = int32 CSR per-cell counts)
                              │
   ── TIER 1 (NEW) collapse_sample.py — parallel across samples ─────────────
      reuse exonAnnotate/exon_str; streaming hash agg keyed (gene,strand,S);
      per-CHROMOSOME flush; contiguous-subsequence containment; singleton flag
      out (beside GFF):  sample.isoform_candidates.tsv.gz
                         sample.molecule_to_candidate.tsv.gz
                              │  ≤~0.5M rows/sample
   ── TIER 2 (NEW) merge_catalog.py — k-way sorted-stream merge, gene-by-gene ─
      reference-first seeding; global total_count>1 filter; contiguous collapse;
      stable novel IDs.  AUTHORITATIVE:
                         gff-output/candidate_to_final.tsv.gz
                         gff-output/final_isoforms.tsv.gz
      DERIVED VIEWS (byte-compatible):
                         gff-output/isoform_links.txt   (+ self-check hash)
                         gff-output/transcript_associations.txt
                         gff-output/isoform_annotations.txt
                         gff-output/combined.gff
                              │
   ── TIER 3 (NEW) re-key per-sample h5ad ───────────────────────────────────
      v1: existing exportConsensusIsoformMatrix consumes isoform_links.txt (A)
      v2: exportConsensusIsoformMatrix_sparse — X_meta = X_raw @ G  (B, flagged)
                         out: <sample>-isoform.h5ad (var = gene:ref_isoform)
```

---

## 4. Tier 1 — per-sample candidate reduction (`collapse_sample.py`)

### 4.1 Reuse validated annotation
Refactor `exonAnnotate`, `exon_str`, `importEnsemblGenes` into importable helpers (**no behavior change**) shared by old and new paths. Per transcript → `(gene, strand, S=trimmed structure, S_raw=structure with terminal coords)`.

### 4.2 Collapse = streaming hash aggregation, **flushed per chromosome**
Collapse key = trimmed structure string `S`. Aggregate:

```
struct_db[(gene, strand, S)] = {count, rep_molecule_id, rep_S_raw, is_ref, ref_id,
                                first_start_min/max, last_end_min/max, cell_count}
molecule_to_candidate  → streamed to gzip TSV (never held in RAM)
```

**Memory correction (resolved blocker).** The raw GFF is **chromosome-partitioned and position-ordered, NOT gene-sorted** — `isoform_structure_extract` fetches per chromosome in coordinate order and `_combine_chunk_gffs` concatenates chunks ([isoform_structure_extract.py:674](../../bam/isoform_structure_extract.py#L674), [:803](../../bam/isoform_structure_extract.py#L803)). A dict with no boundary flush therefore holds **all distinct structures for at least a whole chromosome**, not "one gene's." **Decision: flush `struct_db` at each chromosome boundary** (natural seam, since input is chromosome-clustered). Peak RAM ≈ reference tables + one chromosome's distinct structures (~20–30K rows, ~100–150 MB). Trans-chromosome trans-spliced reads (already a special case in legacy code) go to a dedicated overflow bucket flushed at end, never per-chromosome.

### 4.3 Containment collapse — **CONTIGUOUS subsequence only** (resolved blocker)
> **The two prior designs contradicted each other here.** Design A said "ordered subsequence"; Design B said "contiguous subsequence." **Contiguous is correct.**

A 5′/3′-truncated read shares a **contiguous run** of interior tokens with its parent. An exon-skip / intron-retention / alt-donor / alt-acceptor isoform is a *non-contiguous* subsequence and is a **biologically distinct isoform that must NOT fold**. Worked example:

```
parent:        E1.1|E2.1|E3.1|E4.1|E5.1
E2.1|E3.1|E4.1            → contiguous run  ⇒ FOLD (terminal truncation)
E1.1|E3.1|E4.1|E5.1       → E2.1 skipped, non-contiguous ⇒ DO NOT FOLD (exon skip)
```

Contiguity is **self-sufficient**: it subsumes the "junction-absent" anti-overcollapse guard, because a skip can never be a contiguous run. The legacy string test (`shorter in longer`, [gff_process.py:585](../gff_process.py#L585)) is **retired** — it matches across token boundaries.

Shared utility used identically in Tier 1 and Tier 2:
```python
def tokens(s): return [t for t in s.split('|') if t]
def is_contiguous_subsequence(short, long):
    if len(short) > len(long): return False
    for i in range(len(long) - len(short) + 1):
        if long[i:i+len(short)] == short: return True
    return False
```

**Token-posting-list optimization is a PRE-FILTER, not the test.** Within a gene, index token → structures containing it; intersect the posting lists of `s`'s tokens to get candidate parents (typically 0–3), then run the full `is_contiguous_subsequence` on each survivor. O(k·candidates), candidates ≪ k.

**Parent priority** (encodes "reference first, then longest, then best expression"): sort candidate parents by `(is_ref desc, len desc, count desc)`; a fragment folds into the highest-priority contiguous parent. **Never fold a reference isoform as a fragment.**

### 4.4 Singletons
Structures with `count==1` that did not fold are **kept** (flagged `singleton=True`) — a singleton here may be abundant in another sample. The global "seen-once-across-all" drop happens only in Tier 2.

### 4.5 Outputs (gzipped sorted TSV; Parquet optional via `ISOFORM_USE_PARQUET=1`)
`sample.isoform_candidates.tsv.gz`: `sample, gene, strand, candidate_id, structure_key, raw_structure_key, first_exon_start_min, first_exon_start_max, last_exon_end_min, last_exon_end_max, junction_count, read_count, cell_count, is_ensembl, ensembl_transcript_id, best_exemplar_molecule, singleton`
`sample.molecule_to_candidate.tsv.gz`: `sample, molecule_id, barcode, gene, candidate_id, structure_key, read_count`

Sorted by `(gene, strand, structure_key)` so Tier 2's k-way merge needs no re-sort. `transcript_associations.txt` still emitted as plain TSV for `gff_translate`.

### 4.6 Parallelism / resumability
One sample = one independent job; reuse the SLURM/LSF/PBS-aware worker pool from `parallel_extract_isoform_structures` ([isoform_structure_extract.py:46-54](../../bam/isoform_structure_extract.py#L46-L54)). Reference loaded once per worker (read-only). Skip-if-exists per sample unless `--force`.

---

## 5. Tier 2 — cross-sample canonical resolution (`merge_catalog.py`)

### 5.1 Reference injected first, wins ties
Parse reference GFF(s) through the same annotation helper → reference candidates `is_ref=True, ref_id=ENST*/NM_*`. Seed each gene so equal-or-contained observed structures collapse onto the reference ID ([gff_process.py:950-955](../gff_process.py#L950-L955)).

### 5.2 K-way sorted-stream merge, **one gene at a time** (gzipped TSV default)
Use the existing sort-merge pattern (`_merge_sorted_runs` in `collapse_transcript_associations.py`). **SQLite is NOT the default** — gene-by-gene processing already bounds memory, so the extra dependency buys nothing; offer `--use-sqlite` only as an escape hatch for >10M-candidate cohorts. For each gene:

```
aggregate candidates by (strand, structure):
    total_count = Σ read_count over samples ; sample_count = #samples ; is_ref = any ; ref_id = priority ENST/NM
drop structures with total_count <= 1                      # global singleton rule
winners, fold_map = contiguous_collapse(structs, priority=(is_ref, len, total_count))   # §4.3 test, global counts
assign final_isoform_id: ref_id if present else f"{gene}:NOVEL_{stable_hash(structure_key)}"
emit candidate_to_final rows + final_isoforms rows
```
RAM ceiling = one gene's cross-sample structure set (tiny). Total ≈ few hundred MB regardless of sample count. Gene-block parallelizable.

### 5.3 Stable novel IDs
`gene:NOVEL_<hash(trimmed structure_key)>` — rerun-stable and subset-stable (legacy used whatever molecule happened to be exemplar = unstable). Golden-parity tests compare **structure sets**, not raw novel ID strings.

### 5.4 AUTHORITATIVE outputs
**`candidate_to_final.tsv.gz`** (the one source of truth — schema includes the fields needed to reconstruct the legacy map and to remap h5ads):
```
sample, source_sample, gene, strand, candidate_id, candidate_structure_key,
source_isoform, final_isoform_id, final_structure_key,
collapse_reason, collapse_score, candidate_total_reads, final_total_reads
```
- `source_isoform` = bare molecule ID (joins to h5ad `var_name` after `:`-split).
- `collapse_reason ∈ {EXACT_MATCH, CONTAINMENT, REFERENCE_PRIORITY, EXPRESSION_WINNER, SINGLETON_WINNER, DROPPED_SINGLETON}`. **Single-transcript winners are recorded explicitly as `SINGLETON_WINNER`** (resolved high finding) — not relegated to a lossy empty-source row.

**`final_isoforms.tsv.gz`** (one row per winner): `final_isoform_id, gene, strand, structure_key, source, ensembl_transcript_id, total_reads, sample_count, cell_count, candidate_count, is_novel, is_reference, representative_sample, representative_candidate_id, first_exon_start_min/max, last_exon_end_min/max`.

### 5.5 DERIVED views + self-check
`isoform_links.txt` is generated by a single deterministic transform over `candidate_to_final`: for each row with non-empty `source_isoform`, emit `gene, final_isoform_id, ref_gff, source_isoform, source_sample`; for `SINGLETON_WINNER` rows emit the empty-source form `gene, final_isoform_id, ref_gff, '', ''`. **Self-check:** after writing, re-derive in memory, hash, compare to the file; warn on mismatch. `transcript_associations.txt`, `isoform_annotations.txt`, `combined.gff` (using `rep_S_raw` terminal coords) likewise derived.

---

## 6. Tier 3 — re-key each sample's h5ad

**v1 (Option A, ship first, zero downstream change):** Tier 2 emits `isoform_links.txt`; existing `exportConsensusIsoformMatrix` ([isoform_ratios.py:17](../isoform_ratios.py#L17)) works unchanged.

**v2 (Option B, sparse, behind `use_sparse=False` flag):** new `exportConsensusIsoformMatrix_sparse`, replacing the `toarray()` densification ([isoform_ratios.py:86-120](../isoform_ratios.py#L86-L120)):

```python
# 1. recover bare molecule IDs (BLOCKER: var_names are "GENE:molecule_id")
mol = [v.split(':', 1)[1] if ':' in v else v for v in adata.var_names]
# 2. map molecule -> final column via candidate_to_final (scoped to this sample's source_gff)
final_ids = sorted(set(map_for_sample.values()))                 # full catalog order, identical across samples
col_of   = {f: j for j, f in enumerate(final_ids)}
rows, cols = [], []
for i, m in enumerate(mol):
    f = map_for_sample.get(m)
    if f is None:                                                # singleton-dropped molecule
        continue                                                 # report n_dropped
    rows.append(i); cols.append(col_of[f])
G = coo_matrix((np.ones(len(rows), np.int64), (rows, cols)),
               shape=(adata.n_vars, len(final_ids))).tocsr()
X_meta = (adata.X.astype(np.int64)) @ G                          # int64 upcast — no overflow; no densification
meta = ad.AnnData(X=X_meta, obs=adata.obs.copy(),                # obs/cluster metadata preserved by copy
                  var=pd.DataFrame(index=final_ids))
```
Handles every stress-test blocker: `:`-split alignment; int64 upcast (input is int32, [isoform_structure_extract.py:483](../../bam/isoform_structure_extract.py#L483)); `obs.copy()` preserves cluster assignment from `calculate_barcode_match_percentage`; full-catalog `final_ids` gives identical `var` across samples; unmapped singletons drop to no column with a reported count. Enable for cohorts >50 samples after pseudobulk-ratio parity vs Option A passes.

---

## 7. Integration into `isoform_automate`

- **Tier 1** runs in `combine_processed_samples` *after* `barcode_sample_dict` is built and *after* the junction path completes (junctions must finish first — they're independent and must not change). Add `reduce_sample_exemplars(sample_dict, ensembl_exon_dir)` (parallel, missing-samples-only).
- **Tier 2** replaces the consolidate call: swap `gff_process.consolidateLongReadGFFs(gff_files, ensembl_exon_dir, mode='Ensembl')` ([isoform_automate.py:345](../isoform_automate.py#L345)) → `merge_catalog.build(gff_files, ensembl_exon_dir, reference_gff, force=deleteGFF)`. It emits the 4 contract artifacts so `gff_translate` ([isoform_automate.py:348](../isoform_automate.py#L348)) and the per-sample loop are unchanged.
- **Tier 3** = existing `export_isoform_matrix`/`exportConsensusIsoformMatrix` (Option A) until Option B is validated.
- **Resumability:** Tier-2 skip only if ALL N sample exemplars exist AND the reference GFF is unchanged (mtime/hash) — a stale skip with a changed reference would silently emit an old catalog. `--force` wired to all three tiers.

---

## 8. Sparse / efficiency opportunities (plan item #4)
1. Tier-3 sparse matmul `X_raw @ G` (§6) — the headline win; eliminates `toarray()`.
2. Dictionary-encode `structure`/`candidate_id` columns → integer codes → build `G` directly from codes.
3. Token-posting bitsets in containment (§4.3/§5.2): `s_mask & S_mask == s_mask` pre-filter, then verify contiguity.
4. Reuse the viewer's gene-index cache (`gene_indexes_v2/*.npz`, [isoform_structure_view.py:725-761](../../visualization/isoform_structure_view.py#L725-L761)) to slice sample h5ads by gene.

---

## 9. Resource budget (200 samples × 20M reads)
| stage | per-unit RAM | parallelism | scales with |
|---|---|---|---|
| Tier 1 | ref tables + **one chromosome's** structures (~100–150 MB) | N samples in parallel | reads (CPU only, distributed) |
| Tier 2 | one gene at a time (few hundred MB total) | gene-block optional | N small TSVs (I/O) |
| Tier 3 | one sample sparse h5ad | N samples in parallel | nnz (sparse) |

No stage holds all raw reads in memory. Reductions: ~20M→≤0.5M (Tier 1); 200×0.5M → ≤~1M after `total_count>1` (Tier 2).

---

## 10. Ordered implementation milestones
1. **Shared utilities.** Refactor `exonAnnotate`/`exon_str`/`importEnsemblGenes` into importable helpers (no behavior change); add `tokens` + `is_contiguous_subsequence`. Unit-test the 5 separation cases (exon-skip, alt-donor, alt-acceptor, intron-retention, terminal-truncation). *(Areas 3+4.)*
2. **Tier-1 `collapse_sample.py`** — per-chromosome flush, gzip-TSV candidate + molecule_to_candidate, singleton flag, resumability. *(Area 3.)*
3. **Tier-2 `merge_catalog.py`** — k-way merge, reference-first, `total_count>1` filter, contiguous collapse, stable novel IDs. Authoritative `candidate_to_final` + `final_isoforms`; derive byte-compatible `isoform_links.txt` (+ self-check hash) + the other 3 contract artifacts. *(Areas 1+5.)*
4. **Integration swap** in `combine_processed_samples`: add `reduce_sample_exemplars`, replace consolidate call, wire `--force`. Golden-gene parity vs legacy (reference rows exact, novel by structure set). *(Area 5.)*
5. **Tier-3 sparse** `exportConsensusIsoformMatrix_sparse` (`use_sparse=False` default): `:`-split, COO→CSR `G`, int64 upcast, `obs.copy()`, full-catalog var. Validate pseudobulk-ratio parity; enable for >50-sample cohorts. *(Area 2.)*

---

## 11. Validation plan
1. Exact Ensembl match wins over observed candidate ID.
2. Terminal-truncated read folds into stronger full candidate (contiguous).
3. Alt internal exon / alt donor / alt acceptor / intron-retention isoforms stay **separate** (non-contiguous).
4. Rare novel structure with ≥2 total reads retained; structure seen once across all samples dropped.
5. Isoform abundant in only one sample retained; both samples' molecules map to the same `final_isoform_id`.
6. Molecule→final remap preserves total counts per sample; sparse remap preserves total counts per barcode (Option B == Option A within rounding).
7. Final isoform IDs deterministic across reruns and across sample subsets.
8. Junction/splicing outputs unchanged.
9. Derived `isoform_links.txt` self-check hash matches; legacy `exportConsensusIsoformMatrix` produces identical matrix from the derived view as from a legacy run on the same small input.

---

## 12. Recommended defaults
```
min_total_reads = 2          # drop structures seen once across ALL samples
collapse_terminal_truncations = true
collapse_internal_alternatives = false      # contiguous-only enforces this
prefer_ensembl = true
backend = gzipped_sorted_tsv  # pyarrow + sqlite optional escape hatches
use_sparse_remap = false      # enable >50 samples after parity check
force = false                 # resumable; skip-if-exists per stage
```
