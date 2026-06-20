# iso2function — Isoform-resolved molecular-interaction & functional annotation

**Component:** `altanalyze3/components/iso2function/`
**Status:** scaffolding / v1 in progress (started 2026-06-17)
**Owner:** Salomonis lab (saljh8)

---

## 1. Vision

Long-read modules in AltAnalyze3 (`sclr`, `sclr-isoforms`, `sclr-diff`) produce a catalog of
**known and novel isoforms** per cohort, each represented by a reproducible **Ensembl91 exon/intron
structure string** and a predicted protein (NMD status, intron retention, length). What they do *not*
yet provide is **molecular function**: which proteins an isoform binds (PPI), which DNA elements it
binds (PDI), whether it activates transcription (M1H), and whether it forms condensates / changes
localization.

`iso2function` closes that gap by attaching **definitively measured, isoform-resolved interaction and
function data** to observed isoforms, so that an isoform switch or expression difference between
conditions can be interpreted as a concrete, testable change in a molecular-interaction network — and
visualized (Cytoscape) and tested for pathway enrichment (a new, isoform-focused GSEA).

**Long-term goal:** automate the analysis → visualization → interpretation of *known-function*
isoforms, integrated into the AltAnalyze pipeline and the ISV web interface. Functional *inference*
for unannotated isoforms is explicitly **out of scope for v1** (a later phase).

---

## 2. The identifier model (foundation)

The single most important design decision. Isoform identity is anchored on the **structure string**:
the pipe-delimited Ensembl91 exon/intron token string, e.g.

```
E1.2|E1.3|I1.1|E2.1|E2.2|E4.1|E5.1|E6.1|E6.2|E7.2|E7.3|E8.1|E9.1|E10.1|I10.2|E10.2|I10.1|E11.1
```

(this example = AEBP2 / ENST00000266508). It is the **canonical, cross-project key** because it is
defined against a *fixed* reference annotation (`Hs_Ensembl_exon.txt`, Ensembl91). It is reproducible
across any dataset/cohort/run that uses the same reference — unlike:

- `final_isoform_id` — a *per-cohort* collapse hash (e.g. `84860836.Y2982-E3`); project-specific.
- clone IDs — three different namespaces: TFIso GFF `SYMBOL|i/n|well` (`AEBP2|1/3|11E07`), paper2
  atlas `SYMBOL-N` (`AEBP2-2`) + numeric ORFeome `orf_id`, paper1 `SYMBOL_number` + Entrez + GenBank.
- `ENST` / `ENSP` / UniProt — only exist for *known* (reference) isoforms.

**Rule:** the structure string is the join key for matching isoforms *across* datasets. All other IDs
are crosswalked *to* it, and every crosswalk hop is **verified** (sequence identity), never guessed.

---

## 3. Data sources

### v1 (primary) — paper2: Mol Cell 2025 TFIso atlas — `~/Downloads/TFIso/paper2/`
Isoform-resolved, keyed by `clone_id` (`SYMBOL-N`) + numeric `orf_id`.

| Modality | File (inside zip) | Shape | Key columns |
|---|---|---|---|
| PDI (eY1H) | `mmc3.zip:SuppTable_eY1HResults.txt` | 171 iso × 186 DNA baits, Boolean | `gene_symbol`, `clone_id`, `<bait_id>`... |
| PDI bait seqs | `mmc3.zip:SuppTable_DNABaits.txt` | 186 | `bait_id`, `seq` |
| PDI validation | `mmc3.zip:SuppTable_PDI_validation.txt` | 113 | `gene_symbol`, `clone_id`, `Bait`, `Y1H_result`, reps, `Log2(FC)` |
| PPI (Y2H) | `mmc5.zip:SuppTable_PairwiseY2HResults.txt` | 9562 | `ad_clone_id`, `ad_gene_symbol`, `ad_orf_id`, `db_gene_symbol`, `db_orf_id`, `Y2H_result`, `db_gene_category`, `db_gene_cofactor_type` |
| PPI (N2H) | `mmc5.zip:SuppTable_N2HResults.txt` | 765 | `clone_id`, `gene_symbol_tf`, `gene_symbol_partner`, `test_orf_ida/b`, `source`, `score_pair`, `log2 NLR` |
| Activation (M1H) | `mmc4.zip:Data_S3.tsv` | 622 | `clone_id`, `gene_symbol`, `M1H_rep1..3`, `M1H_mean` |
| Paralog divergence | `mmc7.zip:Data_S6.tsv` | 2139 | `gene_symbol_a/b`, `clone_id_a/b`, `category`, `aa_seq_pct_identity`, `PDI_Jaccard_d`, `PPI_Jaccard_d`, `activation_abs_fold_change_log2` |
| Condensate / localization | `mmc8.zip:Data_S7.tsv` | 189 | `clone_id`, `cloned_reference_or_alternative`, condensate/localization × 2 cell lines × 2 reps |
| **Functional categories (per switch)** | `mmc9.zip:Data_S8.tsv` | 447 | `gene_symbol`, `reference_isoform`, `alternative_isoform`, `DBD_pct_lost_in_alt`, `PDI_category`, `PPI_category`, `M1H_activation_category`, `localization_category`, `alt_iso_classification` — **the authors' own high-confidence switch consequence calls** |
| Isoform GSEA (alt vs ref) | `mmc11.zip:Data_S10.tsv` | 1156 | `alt_iso`, `ref_iso` (clone-keyed), `msigdb_term`, `gsea_nes`, `gsea_qval` (from Joung TF-ORF overexpression) |
| (context) TCGA | `mmc10.zip` | — | tumor isoform usage (later) |
| (niche) PBM | `mmc6.zip` | 32896×2 | CREB1/TBX5 8-mer DNA specificity (not ingested in v1) |

(Headers verified 2026-06-17 by streaming each member; an early inventory had swapped the S3/S6/S7 labels.)

### Interaction-calling thresholds (binary-only, paper-defined — verified verbatim from STAR Methods)

Source: Lambourne et al., *Mol Cell* 2025, `mmc12.pdf`. Policy (user-directed): **use only binary calls; do
not invent cutoffs; if a numeric cutoff is absent from the metadata, do not threshold.**

| Assay | Cutoff in metadata? | Rule applied |
|---|---|---|
| **Y2H** (PPI, primary) | yes | authors' `Y2H_result` Boolean (*"growth score ≥ 2 … higher than auto-activator test"*, p.e7) — used directly as the called PPI edge |
| **eY1H** (PDI, primary) | manual | matrix True/False = *"manually analyzed by three independent researchers"* (p.e6) — used directly as the called PDI edge |
| **M1H** (activation) | yes | `M1H_mean ≥ 1` → activator, `≤ −1` → repressor (p.e12) — applied as a per-clone `m1h_call` |
| **N2H** (PPI, orthogonal) | **NO** | no numeric `log2 NLR` cutoff in text (dashed line only, Fig S2I) → **NOT thresholded; excluded from called edges**, scores kept for reference only |
| **luciferase** (PDI validation) | **NO** | no numeric `Log2(FC)` cutoff in text (Fig S2K) → not thresholded; the table's own `Y1H_result` Boolean is retained as the binary call |

Switch-consequence interpretation uses the authors' own `Data_S8` categories (PDI/PPI/M1H/localization)
in preference to re-deriving gained/lost from raw edges — these are the high-confidence published calls.
| (context) TCGA | `mmc10.zip` | — | tumor isoform usage (later) |

### Bridge data (already produced)
- `~/Claude_Project/TFiso_AltAnalyze3_annotation/tfiso_to_final_assignment.protein.tsv` — 1108 TFIso
  GFF clones (`SYMBOL|i/n|well`) → `final_isoform_id` (MDS-AML observed) + `tfiso_structure` +
  `final_structure` + match_type + protein consequences. **Provides the GFF-clone → structure →
  observed-isoform legs of the crosswalk (already validated, 77% assigned).**
- `~/Dropbox/Revio/ENCODE/TFiso/c_6k_unique_acc_aligns.gff` — TFIso ORFeome clones aligned to genome
  (source of the GFF clone IDs; translatable to protein via `isoform_translation`).

### Reference
- `…/Hs/Ensembl91/Hs_Ensembl_exon.txt` (structure-string vocabulary), `…/Hs_Ensembl-annotations.txt`
  (symbol↔ENSG), GENCODE v45 GFF, `genome.fa`. (cluster `/data/...` ↔ local `/Volumes/...`).

### Deferred (later phases — NOT v1)
- paper1 (NIHMS754821, 2016): 2028 Y2H PPIs + **linear-motif annotations** (mechanism of gain/loss).
- Public gene-level DBs: BioGRID / IntAct / STRING. **Caveat:** gene-level only → cannot be
  definitively assigned to one isoform without sequence/domain inference (deferred). Use for
  *context* edges only, clearly labeled as gene-level.
- Functional **inference** for unannotated isoforms (domain/motif presence → predicted PPI/PDI).

---

## 4. Architecture

Four layers, each a sub-package. Layers below the crosswalk depend only on tidy, validated tables.

```
iso2function/
  config.py            # canonical input paths + /Volumes<->/data resolver
  cli.py               # `altanalyze3 sclr-iso2func` subcommand (ingest|crosswalk|network|enrich|all)
  ingest/              # LAYER 1: parse raw supp tables -> tidy long-form interaction tables
    paper2.py
  crosswalk/           # LAYER 2: reconcile clone_id <-> structure string <-> observed iso / ENST/ENSP
    structure_key.py   #   structure-string normalization + match (exact/contained_in/...)
    clone_map.py       #   atlas clone_id -> structure (gene+isoform-index join: resolve_clone_to_structure)
  associate.py         # PRODUCT: structure-keyed isoform_function.tsv + switch_function.tsv
  network/             # LAYER 3: assemble interaction graph + Cytoscape export + isoform-switch diff
    build.py
    cytoscape.py
  enrichment/          # LAYER 4: isoform-focused gene-set enrichment over interaction partners/targets
    isoform_gsea.py
  resources/           # small, TRACKED bundled lookups (e.g. cofactor categories) — no large data
  data/                # GITIGNORED: extracted raw tables + parsed tidy outputs (regenerable)
  external/            # GITIGNORED: any downloaded reference snapshots
  artifacts/           # GITIGNORED: scratch / figures / per-run outputs
  tests/
```

### Integration points in altanalyze3 (reuse, don't reinvent)
- `components/long_read/isoform_translation.py:extract_cds_and_protein` — DNA→protein for sequence-
  verified crosswalk (clone & observed isoform protein sequences).
- `components/long_read/gff_process.py` — structure-string annotation of a GFF (same path the prior
  TFIso annotation used; honor `ISOFORM_STRUCT_COORDS=0` for GFFs with start>end exon edge cases).
- `components/goelite/` — pathway/GO enrichment engine for Layer 4.
- `utilities/parser.py` — register the `sclr-iso2func` subcommand (mirror existing `sclr-*` blocks).
- ISV web (`isv_web` / cellHarmony webapp) — later: serve the network + functional annotations.

---

## 5. Crosswalk strategy (the make-or-break layer)

Goal: for every observed isoform (and every paper2 clone) establish a **verified** link to its
interaction data. Namespaces are bridged *through the structure string*, with sequence identity as
the arbiter. No ID-string match is trusted without confirmation.

```
paper2 clone_id (AEBP2-2, orf_id) ──hint──┐
                                          ▼
   gene_symbol + protein-sequence identity  ⇄  TFIso GFF clone (AEBP2|2/3|05F03)
                                          │   [translate both via isoform_translation;
                                          │    require high % AA identity, same gene]
                                          ▼
                          structure string  (tfiso_structure)        ← CANONICAL KEY
                                          │   [from tfiso_to_final_assignment; match_type]
                                          ▼
        MDS-AML observed isoform (final_isoform_id) / ENST / ENSP / protein consequences
```

- **Leg A (clone_id ⇄ GFF clone):** the open problem. Both are the *same* TFIso ORFeome, but in
  different ID schemes. Establish via shared `gene_symbol` + **protein-sequence identity** between the
  paper2 clone (needs its AA sequence — derive from ORFeome accession/CDS or, if absent, from the GFF
  clone's own translation) and the GFF clone. ID-number alignment (`-2` ↔ `|2/3|`) is a *hint to
  verify*, never an assumption. Output a confidence-scored map; report unresolved clones explicitly.
- **Leg B (GFF clone ⇄ structure):** already done — `tfiso_to_final_assignment.protein.tsv`.
- **Leg C (structure ⇄ observed/ENST):** already done — `match_type` + `final_*` columns.

Deliverable: `crosswalk.tsv` keyed on structure string, columns = {structure, gene, all clone IDs in
every namespace, ENST/ENSP if known, final_isoform_id, match_type, sequence-identity evidence,
confidence}. This table is the contract every downstream layer consumes.

**Honesty requirements (critical-thinking skill):**
- Preserve the **tested-vs-untested** distinction. A Y2H/Y1H `False` means *tested, not detected* — it
  is evidence of absence only where the pair was assayed. Never render "no edge" as "interaction
  absent" when the pair was simply not in the assay matrix.
- Report crosswalk coverage and unresolved fraction; downstream claims are scoped to resolved isoforms.
- Distinguish **known** (ENST-backed) from **novel** observed isoforms; novel isoforms inherit
  function only via a verified structure match, with match_type carried through.

---

## 6. Network & isoform-switch interpretation (Layer 3)

- **Graph model:** nodes = {isoform (structure-keyed), interactor protein, DNA bait/target gene};
  edges = {PPI (Y2H/N2H), PDI (eY1H), activation (M1H, node attribute), condensate (node attribute)}.
  Edge attributes: assay, direction (tested/detected), score (NLR / log2FC), evidence count.
- **Isoform-switch differential:** given two conditions (or two cell states) with a switch in the
  dominant isoform of a gene (or an expression difference), compute the **differential interactome** —
  partners/targets *gained* and *lost* between the up- and down-isoform — using only measured edges.
  This is descriptive (measured), not inferred.
- **Cytoscape export:** Cytoscape.js JSON + a Cytoscape Desktop `.cyjs`/SIF + style file; node/edge
  attribute tables. Colors follow lab convention (RGB hex only; no named colors; cyan→yellow ramp for
  continuous; single 2-point straight edges in any static vector figure).
- Consumes `crosswalk.tsv` + the cohort's isoform pseudobulk / `sclr-diff` outputs.

---

## 7. Isoform-focused gene-set enrichment (Layer 4)

A new GSEA variant that operates on **isoform-level** evidence rather than gene-level:
- Gene sets = pathways/GO (via `goelite`) + interactor sets (the PPI partners / PDI target genes of an
  isoform from the atlas).
- Test: are the partners/targets *gained or lost* across an isoform switch enriched for a pathway?
  i.e. does the switch concentrate its interaction change on a coherent functional module?
- Background must be the **assayed** space (tested partners/targets), not the whole genome, to respect
  the tested-vs-untested distinction. Report effect size + FDR; non-parametric where appropriate.

---

## 8. Phased task plan

- [x] **P0 — Recon & scaffold.** Map altanalyze3, inventory inputs, fix identifier model, create
      component skeleton, GOALS.md, .gitignore, memory. *(this turn)*
- [ ] **P1 — Ingest (Layer 1).** Parse all paper2 tables from the zips into tidy long-form TSVs in
      `data/`. Validate row/edge counts against the inventory (QC). *Runnable + tested first.*
- [x] **P2 — Crosswalk (Layer 2).** DONE. Legs B+C from `tfiso_to_final_assignment`. Leg A resolved by
      the **gene + isoform-index join** (`resolve_clone_to_structure`): the atlas `SYMBOL-N` and the
      bridge `SYMBOL|i/n|well` are the same paper clone labels, so atlas `N` == bridge `i` maps each
      interaction clone to its standardized structure key (user-confirmed; the `|i/n|well` label is just
      the paper's annotation and is dropped downstream). **724/745 atlas clones resolved to a structure**
      (520 with a known ENST, 643 distinct structures); the 21 unresolved (19 genes absent from the
      bridge, 2 atlas indices beyond the bridge isoform count) are reported, never force-matched.
      GSE253638 / sequence verification deemed unnecessary by the user.
- [x] **P2b — Associate (product layer, `associate.py`).** Structure-keyed deliverables:
      `isoform_function.tsv` (per isoform: gene, ENST/observed, NMD, M1H call, detected PPI partners +
      PDI targets with tested counts) and `switch_function.tsv` (447 ref→alt switches, 434 fully mapped,
      carrying the authors' `Data_S8` PDI/PPI/M1H/localization categories + DBD loss). The structure
      string is the key; the clone label is provenance only.
- [ ] **P3 — Network (Layer 3).** Graph builder + Cytoscape export + isoform-switch differential.
      Demo on a TF with a real isoform switch in MDS-AML-KINNEX-4.
- [x] **P4 — Enrichment (Layer 4, `enrichment/switch_enrichment.py`).** Per-switch differential
      interactions (`switch_differential_interactions.tsv`: PPI/PDI gained/lost in co-tested space) +
      the authors' isoform GSEA re-keyed to structures (`isoform_gsea_by_structure.tsv`, 1083/1156
      mapped). **Validation:** the independent binary gained/lost counts match the authors' `Data_S8`
      PPI categories by direction (loss→lost≫gained, rewire→both, no-change→~0) — confirms the
      crosswalk + calls. `enrich_partner_categories()` runs the hypergeometric layer over interactor
      classes (TF/cofactor) needing no external sets. Remaining: wire GO-Elite pathway gene sets.
- [x] **P5 — CLI + integration.** `sclr-iso2func {ingest|crosswalk|associate|network|enrich|all}`
      registered in `utilities/parser.py` (import is module-load-safe; global CLI verified to still load);
      dispatches end-to-end through the real `altanalyze3` entry point. Tests under `tests/` (13/13).
- [x] **P7 — Paper1 PPIs (Yang 2016, measured, incl. non-TF).** DONE. Ingest (`ingest/paper1.py`):
      `paper1_ppi` (2026; 1043+/763- binary), `paper1_isoforms` (1035), `paper1_orfs` (1423 ORF seqs),
      `paper1_linear_motifs` (151), `paper1_domain_domain` (60). Crosswalk (`crosswalk/paper1_map.py`):
      translate ORFs -> match GENCODE -> **698 ENST-mapped (603 exact + 95 clone-variant)** ->
      **631 -> Ensembl91 structure + genomic splice coords** (references 389/506 = 77%). Product
      (`associate.build_paper1_function`): `paper1_isoform_function.tsv` (source=paper1, evidence=measured).
      **Crosswalk-correctness fixes (2026-06-18):** (1) gap-aware identity (`difflib`) replacing the
      buggy ungapped N-terminal scoring; (2) accumulate ALL ENSTs per protein so the build-gap ENST that
      has an Ensembl91 structure is found (`reference_structures.resolve_structure`) — this raised
      structure coverage 545 -> 631; (3) clone-variant recovery at gap-aware identity >=99%.
      ~725 unmatched are genuinely NOVEL alternative isoforms (no reference protein) — they would need
      splice-aware genome alignment of their ORFs to derive a de-novo Ensembl91 structure (next step for
      full coverage; not force-mapped to a canonical structure, which would be wrong).
- [x] **P8 — UniProt isoform -> ENST(by sequence) -> domain-contact -> PPI.** DONE (`uniprot.py`).
      4187 human UniProt isoforms (1232 canonical / 2955 alt) sequence-matched to GENCODE ENST (77%
      exact) -> Ensembl91 structure + genomic coords (`uniprot_isoform_map.tsv`). Domain-contact
      features (PPI/DNA interface counts+spans, domains, zinc fingers) + gene-level UniProt-curated
      interaction_partners attached; per-isoform PPI impact inferred from PPI-interface RETENTION vs the
      gene's canonical (`uniprot_isoform_function.tsv`). `evidence=domain_inferred` — kept separate from
      measured. Adds non-TF coverage.
- [x] **Shared infra (`crosswalk/sequence_map.py` + `crosswalk/reference_structures.py`).** A
      protein-sequence -> (ENST) matcher (exact + same-gene fuzzy) over the GENCODE pc_translations
      FASTA; ENST -> Ensembl91 structure via `ENST_reference_structures.tsv` (244k); structure ->
      genomic splice-site coordinates via `structure_coords.sqlite` (decoded with
      `bam.coord_store.unpack_exons`). All reused from existing Daedalus/MDS-AML outputs — no download.
- [x] **Portable dual annotation (user requirement).** Every crosswalked isoform now carries BOTH its
      Ensembl91 E/I structure string (canonical key for this build) AND its genomic splice-site
      coordinates (`chrom:strand:donor-acceptor,...`) so a different annotation build can re-derive the
      structure from the raw junctions.
- [ ] **P6 — ISV web integration** *(later)*. Serve network + functional annotations in the viewer.
- [ ] **Deferred:** public gene-level DBs (BioGRID/IntAct/STRING) as labeled context edges.

### Parallelization opportunities
- P1 modality parsers are independent → fan out one agent per supp table.
- P2 Leg-A sequence matching is embarrassingly parallel per gene.
- P3 isoform-switch differentials are independent per gene → fan out.

---

## 9. Risks & open questions

1. **Leg A (clone_id ⇄ GFF clone) feasibility** — top risk, now scoped: a verifiable bridge requires
   the atlas-clone sequences, which are **not in the supplement** (deposited at GEO GSE253638). Options:
   (a) fetch GSE253638 and sequence-verify; (b) user/author-provided clone↔transcript map;
   (c) clone-level-only delivery until then. **No guessed ID alignment will be shipped** (the
   `id_index_hint` is recorded but never treated as verified).
2. **paper2 clone protein sequences** — CONFIRMED absent from the provided supp tables; present only in
   the deposited consensus FASTA / TFIso1.0 collection (GSE253638). paper1 has full ORF sequences but a
   different ID scheme (`SYMBOL_N`) with no crosswalk.
2b. **Thresholds** — RESOLVED: Y2H/eY1H/M1H cutoffs verified verbatim and applied; N2H/luciferase have
   no numeric cutoff in metadata → excluded from calls (binary-only). See the thresholds table above.
3. **Reference-isoform PPI annotation gap** (user's concern) — for the TFIso atlas this is solved
   (isoform-resolved). For *general* DBs it is not; hence those are deferred and labeled gene-level.
4. **Coverage** — only ~171–2139 clones are phenotyped per modality; the cohort has millions of
   isoforms. Functional annotation is sparse by construction; always report coverage, never imply
   genome-wide function.

---

## 10. Constraints (hard rules)

- **No git/gh commands** — the user owns all git operations and pushes. Claude never creates GitHub
  repos or pushes.
- **No mass/dangerous deletions** of folders.
- **Never modify input files** (`~/Downloads/TFIso/*`, Dropbox GFF, Ensembl/MDS references are
  read-only). Extract zips into `data/` (gitignored), never into the source folders.
- **Artifacts gitignored**: `data/`, `external/`, `artifacts/` (and the repo-wide `components/*/artifacts/`).
- **Python over numpy/scipy/sklearn only** — no R, no rpy2, no compiled-language backends.
- **No guessing identifiers/coordinates** — every ID/coordinate looked up or sequence-verified; say
  "not verified" otherwise.
- Lab output conventions: CSV (matrices as TXT), editable PDF (Arial, fonttype 42), RGB-hex colors,
  single 2-point straight lines, comprehensive README per analysis.
