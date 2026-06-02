# Cross-Sample Long-Read Isoform Collapse — Design

**Status:** Proposed (supersedes the untested `collapse_transcript_associations.py`)
**Scope:** Collapse raw long reads into a standardized, cross-dataset isoform catalog and re-key each sample's h5ad, at the scale of **≤200 samples × ≤~20M reads/sample**, without dropping rare-but-real novel isoforms and without loading all samples' reads at once.
**Author:** design draft, 2026-05-30

---

## 1. The problem with the current path, stated precisely

Isoform consolidation today runs through one function:

`gff_process.consolidateLongReadGFFs(gff_files, ensembl_exon_dir, mode='Ensembl')`
([gff_process.py:670](../gff_process.py#L670)), called from
`export_isoform_h5ad` ([isoform_automate.py:345](../isoform_automate.py#L345)).

It does three things that are individually correct but collectively don't scale:

1. **All-reads-in-RAM.** Every transcript record from *every* sample GFF (plus the Ensembl/GENCODE reference GFFs) is parsed and accumulated into in-memory dicts: `gene_db[gene] -> [isoform_string,...]` and `junction_db[(gene, isoform_string)] -> [(file, info),...]` ([gff_process.py:705-809](../gff_process.py#L705-L809)). With 200 samples each contributing millions of *raw, uncollapsed* reads, `junction_db` alone holds one entry per `(gene, distinct_structure)` and one `(file, info)` tuple **per read**. This is the dominant memory cost and it grows with total read count, not with the (much smaller) number of distinct structures.

2. **O(n²) per-gene substring collapse.** `collapseIsoforms` sorts a gene's isoform strings by length and, for each, scans all shorter ones testing Python `in` substring containment ([gff_process.py:577-587](../gff_process.py#L577-L587)). For a gene with *k* distinct structures this is O(k² · L). Across 200 samples *k* per gene can be large, and the work is wasted because most of those *k* structures are identical across samples.

3. **One molecule = one var column, deferred to per-sample matrix build.** Each sample's per-read h5ad (from `isoform_structure_extract.py`) has one `var` column **per molecule ID** (`var_names = "GENE:molecule_id"`, [isoform_structure_extract.py:454](../../bam/isoform_structure_extract.py#L454)). The molecule→meta-isoform reduction only happens later, in `exportConsensusIsoformMatrix` ([isoform_ratios.py:54-97](../isoform_ratios.py#L54-L97)), via the `isoform_links.txt` map. So the expensive structural collapse and the per-molecule expansion are coupled to a single giant in-memory pass.

The splicing/junction path is **not affected** and must not change — it already works on validated junctions and runs *before* isoform collapse (`pre_process_samples` → junctions; `combine_processed_samples` → isoforms, [isoform_automate.py:598-618](../isoform_automate.py#L598-L618)).

### The structural insight that makes this tractable

> The number of *distinct exon-structures* is bounded by biology, not by read count. A gene has tens to low-hundreds of real isoforms. 20M reads collapse to that. **If we collapse per-sample first, the cross-sample merge operates on `n_samples × n_structures_per_gene` rows — megabytes, not terabytes.**

That is the entire strategy: **two-tier collapse**. Tier 1 is embarrassingly parallel per sample and reduces ~20M reads → ≤~0.5M exemplars. Tier 2 merges 200 small exemplar tables → one catalog (≤~1M, then filtered). Tier 3 re-keys each sample's h5ad against the catalog. No step ever holds all raw reads.

---

## 2. Output contract we must preserve (non-negotiable)

The downstream consumer `exportConsensusIsoformMatrix` reads **`isoform_links.txt`** and builds, for one sample's GFF source, a map `source_isoform -> "gene:ref_isoform"` ([isoform_ratios.py:62-78](../isoform_ratios.py#L62-L78)). The file schema is exactly:

```
gene <tab> ref_isoform <tab> ref_gff <tab> source_isoform <tab> source_gff
```
([gff_process.py:931, 965-986](../gff_process.py#L931-L986))

- `source_isoform` is the molecule/transcript ID that appears as a `var_name` (after the `gene:` prefix is stripped) in a sample's per-read h5ad.
- `ref_isoform` is the **representative** (collapsed) ID — Ensembl `ENST*`/RefSeq `NM_/XM_/XR_` preferred, else a novel exemplar ID.
- `source_gff` / `ref_gff` are the basenames (sans `.g*`) used to scope a sample, matched against `gff_source` ([isoform_ratios.py:67,71](../isoform_ratios.py#L67-L71)).

**Whatever the new pipeline does internally, its final emitted artifact must be a byte-compatible `isoform_links.txt`** (plus `transcript_associations.txt`, `isoform_annotations.txt`, `combined.gff` for `gff_translate`, [isoform_automate.py:348](../isoform_automate.py#L348)). Then *zero* downstream code changes are required. This is the seam we build to.

Also preserved verbatim from `gff_process`:
- The **AltAnalyze structure string** `E1.1|E2.2|E2.4|` built by `exonAnnotate` + `exon_str` ([gff_process.py:772-788](../gff_process.py#L772-L788)).
- **Terminal-coordinate stripping** of novel 5′/3′ sites before the structure becomes a collapse key (`exon_str` drops a leading/trailing token containing `_`; `_strip_terminal_coords` does the same for clustering, [gff_process.py:420-431](../gff_process.py#L420-L431)).
- **Ensembl/RefSeq priority** when choosing a group's representative ([gff_process.py:950-955](../gff_process.py#L950-L955), [gff_process.py:631](../gff_process.py#L631)).

---

## 3. Architecture

```
                 (already exists, unchanged)
  BAM  ──isoform_structure_extract──▶  sample.gff.gz   +   sample.h5ad
                                          │  (raw reads,        (var = GENE:molecule_id,
                                          │   1 transcript/read) X = per-cell read counts)
                                          ▼
 ┌─────────────────────────────────────────────────────────────────────┐
 │ TIER 1  (NEW)  per-sample exemplar collapse — parallel across samples │
 │   collapse_sample.py                                                  │
 │   in:  sample.gff.gz, ensembl_exon_dir                                │
 │   out (beside the GFF):                                               │
 │     sample.exemplars.parquet   gene,strand,structure,count,           │
 │                                rep_molecule_id, is_ref, ref_id        │
 │     sample.mol2exemplar.parquet  molecule_id -> structure  (or to     │
 │                                  the dropped-singleton sentinel)      │
 └─────────────────────────────────────────────────────────────────────┘
                                          │  ≤~0.5M rows/sample
                                          ▼
 ┌─────────────────────────────────────────────────────────────────────┐
 │ TIER 2  (NEW)  cross-sample structure merge — one pass, low memory    │
 │   merge_catalog.py                                                    │
 │   in:  all sample.exemplars.parquet  +  Ensembl/GENCODE reference     │
 │   out: gff-output/catalog.parquet      (the meta-isoform table)       │
 │        gff-output/isoform_links.txt    (CONTRACT artifact)            │
 │        gff-output/transcript_associations.txt, isoform_annotations.txt│
 │        gff-output/combined.gff                                        │
 └─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
 ┌─────────────────────────────────────────────────────────────────────┐
 │ TIER 3  (NEW, thin)  re-key per-sample h5ad to catalog IDs            │
 │   relabel_sample.py  (or: just let exportConsensusIsoformMatrix run,  │
 │   since isoform_links.txt now carries the full molecule->ref map)     │
 │   out (beside the GFF):  sample_isoform.h5ad  (var = gene:ref_isoform)│
 └─────────────────────────────────────────────────────────────────────┘
```

The cleanest integration: **Tier 2 emits a real `isoform_links.txt`**, so Tier 3 is *already implemented* by the existing `exportConsensusIsoformMatrix`. We only replace `consolidateLongReadGFFs` *for the multi-sample isoform case* and leave the entire matrix-building tail intact.

---

## 4. Tier 1 — per-sample exemplar collapse (`collapse_sample.py`)

**Goal:** Turn one sample's raw-read GFF into a tiny table of distinct structures with read counts, prioritizing reference models, while respecting expression when deciding which long structure "wins."

### 4.1 Reuse the validated annotation, not the brute-force collapse

We **reuse** the GFF→structure machinery from `gff_process` (so annotation rules stay byte-identical) but bypass the O(n²) `collapseIsoforms`:

- Parse the GFF with the existing reader loop ([gff_process.py:828-873](../gff_process.py#L828-L873)).
- For each transcript call `exonAnnotate(...)` → `exon_str(...)` to get the **trimmed structure string** `S` and the **raw structure** (with terminal coords) `S_raw` ([gff_process.py:786-788](../gff_process.py#L786-L788)).
- These two functions are the only parts of `consolidateLongReadGFFs` we depend on; refactor them into importable helpers (`annotate_transcript(chr, strand, exons, tid) -> (gene, S, S_raw, genes)`) so both the old and new paths share one definition. No behavior change.

### 4.2 Collapse = streaming hash aggregation (O(reads), tiny memory)

The collapse key is the **trimmed structure string** `S` (terminal coords already stripped). Aggregation is a dict keyed by `(gene, strand, S)`:

```
struct_db[(gene, strand, S)] = {
    count           : int                     # reads with exactly this structure
    rep_molecule_id : str                     # a representative molecule (longest S_raw, then ENST/NM)
    rep_S_raw       : str                      # retained terminal coords of the representative
    is_ref          : bool                     # rep_molecule_id is ENST*/NM_/XM_/XR_
}
mol2struct (streamed to disk, NOT held) : molecule_id -> (gene,strand,S)
```

Why this is enough and why it's cheap:

- **Identical structures collapse for free** — they hash to the same key. This already eliminates the vast majority of redundancy (most reads of a gene share a handful of structures). No similarity math, no O(n²).
- **Memory is bounded by distinct structures in *this one sample*** (tens of thousands per genome, not millions of reads). `count` is an int; we never keep the per-read `(file, info)` lists that bloat the legacy `junction_db`.
- `mol2struct` is **append-only to a Parquet/CSV writer** as we stream reads — never materialized in RAM. This is the molecule→interim-exemplar relation the plan calls for, stored beside the input GFF.

### 4.3 Containment collapse (the "shorter folds into longer" step) — bounded and expression-aware

Equal-structure folding above does not yet fold a **fragment** (a short read whose junction chain is a *subsequence* of a longer structure) into its parent. The legacy code does this with full O(k²) substring matching. We replace it with a **prefix/suffix-indexed, expression-ordered, bounded** containment pass *within a gene* — and critically only over the now-small set of distinct structures, not over reads:

1. Within a gene, list distinct structures sorted by **(is_ref desc, length desc, count desc)**. This directly encodes the plan's requirement: *reference models first, then longest, then best expression evidence* — so a long structure only becomes a "winner"/parent if it is reference or sufficiently supported, never just because one stray long read existed.
2. Build a token index: map each **junction token** (`E2.2`, `I3.1`, …) → the set of structures containing it. A short structure `s` can only be contained in a long structure `S` if every token of `s` appears in `S`; intersect the token-posting lists of `s`'s tokens to get *candidate* parents (typically 0–3), then verify `s`'s token list is an ordered subsequence of `S`'s. This turns O(k²) into ≈ O(k · candidates), candidates ≪ k.
3. A fragment folds into the **highest-priority** parent it is contained in (parents already sorted by the rule above). Its `count` is added to the parent's `count`; its `mol2struct` rows are rewritten to point at the parent's `S` in a second cheap pass over the (small) fold map.
4. **Anti-overcollapse guard** (addresses the stated risk that "the longest reads may not be the ones the shorter collapse into" and that subsequence clustering over-merges): a fragment is *not* folded into a longer parent when (a) the fragment itself is a reference isoform (`is_ref`), or (b) the fragment carries a junction token (e.g. an `I*` intron-retention or a novel `E*_coord` interior site) **absent** from the parent — i.e. only true ordered-subsequences fold, never near-misses. This is stricter than `cluster_by_subsequence`'s `threshold=0.4` fuzzy match ([isoform_structure_view.py:3298](../../visualization/isoform_structure_view.py#L3298)) and is intentional: at Tier 1 we only merge structures that are *literally* compatible. Fuzzy biological clustering, if ever wanted, belongs in the viewer, not in quantification.

**Singletons:** structures with `count == 1` that did **not** fold into any parent are flagged. They are *kept* in the per-sample exemplar file (a structure seen once *here* may be abundant in another sample — plan requirement #2) but marked `singleton=True` so Tier 2 can apply the global rule "drop structures seen exactly once *across all samples*."

### 4.4 Tier 1 outputs (beside the input GFF)

| file | columns | ~size |
|---|---|---|
| `sample.exemplars.parquet` | `gene, strand, structure, count, rep_molecule_id, rep_S_raw, is_ref, ref_id, singleton` | ≤~0.5M rows |
| `sample.mol2exemplar.parquet` | `molecule_id, structure` (gene/strand recoverable from structure's gene) | ~n reads kept, columnar+dict-encoded ⇒ small |

Parquet (via pyarrow) is chosen over TSV for the intermediates because: dictionary encoding crushes the highly repetitive `structure`/`gene` columns; column pruning lets Tier 2 read only `gene,strand,structure,count,is_ref` for the merge without touching `mol2exemplar`; and it's typed (no re-parsing ints). `transcript_associations.txt` is **still emitted as TSV** for `gff_translate` compatibility. If pyarrow is undesirable as a dependency, fall back to gzipped TSV sorted by `(gene,strand,structure)` — the merge algorithm in §5 only needs sorted streams.

### 4.5 Parallelism & memory for Tier 1

- One sample = one independent job ⇒ trivially parallel (the existing `parallel_extract_isoform_structures` already establishes the SLURM/LSF/PBS-aware worker-count pattern, [isoform_structure_extract.py:46-54](../../bam/isoform_structure_extract.py#L46-L54)); reuse that helper for the sample-level pool.
- The Ensembl exon reference (`importEnsemblGenes`, [gff_process.py:1039](../gff_process.py#L1039)) is loaded **once per worker process** and is read-only; it's the only large shared structure.
- Peak RAM per worker ≈ reference tables + one gene's worth of distinct structures. Independent of read count.

---

## 5. Tier 2 — cross-sample structure merge (`merge_catalog.py`)

**Input:** N `sample.exemplars.parquet` + the reference GFFs (Ensembl, GENCODE). **Output:** the catalog and the contract artifacts.

### 5.1 Reference structures are injected first and win ties

Parse the reference GFF(s) through the *same* `annotate_transcript` helper to produce reference exemplars with `is_ref=True` and `ref_id = ENST*/NM_*`. Seed the catalog with these so that any sample structure equal to (or a fragment of) a reference collapses onto the reference ID — exactly the existing Ensembl-priority behavior ([gff_process.py:950-955](../gff_process.py#L950-L955)), now applied globally and only once.

### 5.2 Merge by sorted streaming (no global hash table of reads)

Because each sample file is small and we can sort each by `(gene, strand, structure)`, do a **k-way merge** (heap over N+refs sorted iterators — the pattern already in `collapse_transcript_associations._merge_sorted_runs`, [collapse_transcript_associations.py](../isoform-collapse/collapse_transcript_associations.py)). Walk all rows for one `gene` at a time:

```
for gene in merged_stream:                       # all rows for one gene held at once — small
    structs = aggregate rows by (strand, structure):
                 total_count = Σ count over samples
                 sample_count = #samples with this structure
                 is_ref = any(is_ref)
                 ref_id = first ENST/NM by priority, else None
    # global singleton rule: drop structures with total_count <= 1
    structs = [s for s in structs if s.total_count > 1]   # plan: ignore seen-once-across-all
    # containment collapse over the *gene's* distinct structures (same bounded,
    # token-indexed, expression-ordered algorithm as §4.3, now on global counts)
    winners, fold_map = collapse_containment(structs, priority=(is_ref, len, total_count))
    assign representative id per winner:
        ref_id if present else f"{gene}:NOVEL_{stable_hash(structure)}"
    emit catalog rows + isoform_links rows + annotations
```

Holding **one gene at a time** is the memory ceiling, and a single gene's cross-sample structure set is tiny. Total Tier-2 RAM ≈ a few hundred MB regardless of sample count. Wall-clock is dominated by reading/sorting N small Parquet files, parallelizable by chromosome/gene-block if needed.

### 5.3 Representative-ID assignment & stability

- Reference present → `ref_id` (e.g. `ENST00000655252.1`). Matches legacy output.
- Novel winner → `gene:NOVEL_<hash>` where `<hash>` is a stable hash of the **trimmed structure string** (not of read content), so the same novel isoform gets the same ID across reruns and across the inevitable re-processing of subsets of samples. (Legacy novel IDs were whatever molecule happened to be the exemplar — unstable across runs. This is a strict improvement and is backward-tolerable because downstream only treats the string as an opaque key.)
- `rep_S_raw` (terminal coordinates of the most-supported representative read) is retained in the catalog per the plan, so a final GFF/annotation can show real 5′/3′ ends even though they weren't used as collapse keys.

### 5.4 Tier 2 outputs

| file | purpose | schema |
|---|---|---|
| `gff-output/isoform_links.txt` | **CONTRACT** for `exportConsensusIsoformMatrix` | `gene, ref_isoform, ref_gff, source_isoform, source_gff` |
| `gff-output/transcript_associations.txt` | input to `gff_translate` | `gene, strand, structure, transcript_id, source_gff` |
| `gff-output/isoform_annotations.txt` | known/novel + members | `gene, ref_isoform, ref_gff, members(csv), known|novel` |
| `gff-output/combined.gff` | for translation/visualization | standard GFF, one record per winner using `rep_S_raw` |
| `gff-output/catalog.parquet` | machine-readable master table | `gene,strand,structure,ref_isoform,is_ref,total_count,sample_count,n_members` |

The first four are **drop-in replacements** for what `consolidateLongReadGFFs` writes today; building them from the catalog means `export_isoform_h5ad` needs only its call site swapped (`consolidateLongReadGFFs` → `merge_catalog.build`), nothing downstream.

### 5.5 Emitting `isoform_links.txt` correctly

For each winner with representative `ref_isoform` and member structures `{m}`, and for each member structure's **per-sample** appearances, write one row per `(sample source_gff, source_isoform)`. The `source_isoform` values are the molecule/exemplar IDs from `sample.mol2exemplar.parquet` for that structure — i.e. Tier 2 joins each winning/member structure back to the molecule IDs that map to it, scoped by sample (`source_gff`). This is the same fan-out the legacy loop does over `isoform_pair_map` ([gff_process.py:974-984](../gff_process.py#L974-L984)), but driven by the compact mol2exemplar tables instead of an in-RAM `junction_db`.

> Optimization: if `exportConsensusIsoformMatrix` is updated to consume `catalog.parquet` + per-sample `mol2exemplar.parquet` directly (vectorized join of the sample h5ad's `var_names` against a `molecule_id -> gene:ref_isoform` Series), we skip materializing the potentially large `isoform_links.txt` entirely. See §7.

---

## 6. Tier 3 — re-key each sample's h5ad

Two equivalent options; pick by how much downstream we want to touch:

**Option A (zero downstream change):** Tier 2 writes `isoform_links.txt`. The existing `exportConsensusIsoformMatrix` ([isoform_ratios.py:17](../isoform_ratios.py#L17)) already reads it, builds `source_isoform -> gene:ref_isoform`, and sums per-cell read counts from the sample h5ad into the meta-isoform matrix ([isoform_ratios.py:89-97](../isoform_ratios.py#L89-L97)). **Done.** The only inefficiency is its `cell_data = adata.X.toarray()` densification ([isoform_ratios.py:86](../isoform_ratios.py#L86)) — see §7.

**Option B (faster, slightly more code):** a thin `relabel_sample.py` that loads `sample.h5ad`, builds a `var_name -> ref_isoform` vector from `mol2exemplar`+catalog, and does a **sparse column-sum reduction** (group var columns by ref id) entirely in CSR space:
`M_meta = M_raw @ G`, where `G` is a sparse `(n_molecules × n_meta)` 0/1 grouping matrix. This is the sparse-matrix efficiency opportunity the plan asks for (#4): no densification, one sparse matmul per sample, output written beside the GFF as `sample_isoform.h5ad`.

Recommended: ship **A** for immediate correctness/compatibility, then offer **B** behind a flag once validated against A (pseudobulk cluster ratios should match within rounding).

---

## 7. Sparse / vectorization opportunities (plan item #4)

1. **Per-sample re-key as a sparse matmul** (§6 Option B): `X_raw (cells×molecules) · G (molecules×meta)` → `X_meta (cells×meta)`. Replaces the Python per-isoform loop + `toarray()` in `exportConsensusIsoformMatrix` ([isoform_ratios.py:84-120](../isoform_ratios.py#L84-L120)) which currently densifies and iterates columns. Big win at 200 samples.
2. **Tier-1 mol2exemplar as categorical codes.** Store `structure` as a dictionary-encoded (integer-coded) column; the molecule→structure map becomes an int array, and the grouping matrix `G` is built from `scipy.sparse` codes directly.
3. **Token postings as bitsets** in the containment step (§4.3/§5.2): represent each gene's structures as bitmasks over the gene's token vocabulary; "is `s` a subset of `S`?" becomes `s_mask & S_mask == s_mask`, then verify order. Vectorizable across candidate parents.
4. **Reuse the gene-index cache** from the viewer (`gene_indexes_v2/*.gene_index.npz`, [isoform_structure_view.py:725-761](../../visualization/isoform_structure_view.py#L725-L761)) so Tier 3 can slice a sample h5ad by gene without scanning all vars.

---

## 8. What we keep, replace, and add

| Component | Action |
|---|---|
| `isoform_structure_extract.py` (BAM→raw GFF/h5ad) | **Keep** unchanged. |
| Junction/splicing path (`pre_process_samples`) | **Keep** unchanged; runs first. |
| `exonAnnotate` / `exon_str` / `importEnsemblGenes` | **Refactor** into shared importable helpers (no behavior change). |
| `consolidateLongReadGFFs` brute-force collapse | **Replace** for the multi-sample isoform case with Tier 1 + Tier 2. Keep the single-GFF/`gene_id` debug paths. |
| `collapse_transcript_associations.py` | **Replace.** It had the right skeleton (disk sort, per-sample counts, high-vs-low confidence) but: tokenizes on bare `\|` with no exon awareness, uses fuzzy Jaccard/LCS that can chimerize/over-merge, isn't wired into `isoform_automate`, and never injects reference priority. Salvage: its external sort-merge (`_write_sorted_runs`/`_merge_sorted_runs`) is reused for Tier 2's k-way merge. |
| `exportConsensusIsoformMatrix` | **Keep** (Option A) or **add** sparse fast path (Option B). |
| `export_isoform_h5ad` call site ([isoform_automate.py:345](../isoform_automate.py#L345)) | **One-line swap** to `merge_catalog.build(...)`. |

---

## 9. Resource budget (back-of-envelope, 200 samples × 20M reads)

| stage | per-unit RAM | parallelism | dominant cost |
|---|---|---|---|
| Tier 1 | ref tables (~1–3 GB) + 1 gene's structures | N samples in parallel (cap by cores/SLURM) | streaming GFF parse, O(reads) |
| Tier 2 | few hundred MB (one gene at a time) | sort I/O; gene-block parallel optional | reading/sorting N small Parquets |
| Tier 3 | one sample h5ad sparse | N samples in parallel | one sparse matmul/sample |

No stage scales with *total* reads in memory; Tier 1's CPU is the only thing that does, and it's distributed. Target reductions hold: ~20M→≤0.5M (Tier 1), 200×0.5M→≤~1M after the `total_count>1` filter (Tier 2).

---

## 10. Correctness / validation plan

1. **Golden-gene parity:** for a handful of genes (incl. a reference-rich one like a well-annotated TF and a novel-heavy one), run legacy `consolidateLongReadGFFs(mode='Ensembl')` vs new pipeline on the *same* small set of GFFs; assert identical `isoform_links.txt` rows for reference winners (novel IDs will differ by the stable-hash naming — compare structure sets instead).
2. **Pseudobulk invariant:** collapsed per-cluster isoform TPM/ratio (the `*_pseudo_cluster_*` outputs) must match between Option A and Option B within rounding.
3. **No-drop check:** assert every structure with `total_count > 1` in any sample appears in the catalog (rare-but-real isoforms preserved); assert every `total_count == 1` structure is the *only* thing dropped.
4. **Singleton-rescue check:** a structure that is a singleton in sample X but abundant in sample Y survives and both samples' molecules map to the same `ref_isoform`.

---

## 11. Open questions for review

1. **pyarrow dependency** acceptable for the intermediates, or stick to gzipped sorted TSV? (Design works either way; Parquet is faster/smaller.)
2. **Novel-ID scheme** — `gene:NOVEL_<hash-of-structure>` good, or do you want a running integer per gene (`gene:PB.N`) for readability? Hash is rerun-stable; integer is prettier but order-dependent.
3. **Fragment-fold strictness** (§4.3 guard) — confirm we want *literal ordered-subsequence only* (no fuzzy merge) at quantification time. I believe yes, given the over-clustering concern.
4. Ship **Option A first** (zero downstream change), add sparse **Option B** behind a flag after validation — agreed?
