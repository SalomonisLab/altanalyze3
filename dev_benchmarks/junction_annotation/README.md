# Junction annotation speed work

I profiled and then optimized the AltAnalyze3 splice-junction annotation code on 2026-08-08 and
2026-08-09. I proved on 29,643 junctions that `legacy` mode reproduces the prior output exactly.

On 2026-08-09 I validated `corrected` mode on the real Cal27P5 test data and made it the default.
See
`/Users/saljh8/Dropbox/Revio/BoneMarrow/test_env/altanalyze3/tests/annotation_speed_benchmark/README.md`.
Pass `--novel-gene-mode legacy` to reproduce an older result.

All code edits live inside altanalyze3 and I validated every one of them inside altanalyze3.

## Overview

`annotate_junctions` scanned all 68,008 genes for every junction whose two splice sites both miss
the reference. The scan cost 3.40 ms per call on the plus strand and 14.78 ms on the minus strand,
and the code called it twice per junction. I replaced the scan with a precomputed sweep index that
answers the same question by one `bisect`, at 0.24 us per call.

## Method

I measured against the real human reference file. No junction BED file and no aggregated h5ad exist
on this Mac, so I built the junction tables from real reference coordinates. The reference
dictionaries and the gene scans therefore run at real data volumes. The junction mix is synthetic,
and I state each mix below.

- Reference file:
  `/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt`
  (1,329,979 rows, modified 2025-07-17)
- Interpreter: `/opt/homebrew/opt/python@3.11/bin/python3.11` (pandas 1.5.2, anndata 0.10.9)
- Frozen copy of the original function, taken before I edited anything:
  `/Users/saljh8/Documents/GitHub/altanalyze3/dev_benchmarks/junction_annotation/annotate_reference_frozen.py`

## What I changed

| File | Change |
|---|---|
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/aggregate/gene_index.py` | New. Builds the sweep index. Serves both modes. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/aggregate/annotate.py` | Queries the index. Iterates columns, not `iterrows`. Reuses the duplicated call. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/aggregate/main.py` | Passes `novel_gene_mode` and logs it. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/utilities/parser.py` | Adds `--novel-gene-mode {corrected,legacy}`, default `corrected`. |
| `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/long_read/gff_process.py` | Declares `gene_chr` global, clears it in `clearEnsemblCache`, and detects whether line 1 is a header. |

I did not change the return signature of `importEnsemblGenes`. Every existing caller still unpacks
three values. The 3 long-read tests still pass.

## Speed

At 20,000 junctions with 15% of them fully novel:

| | Total | Annotation loop | Reference load |
|---|---|---|---|
| Original | 48.86 s | 46.67 s | 2.19 s |
| New, `legacy` | 2.46 s | 0.27 s | 2.19 s |

The loop runs 173x faster. Projected to 500,000 junctions at the same 15% novel fraction, the loop
drops from 19.4 min to 6.8 s. I assume that 15% fraction. I did not measure it, because no junction
BED file exists on this Mac.

Component measurements:

| Step | Before | After |
|---|---|---|
| Gene lookup, plus strand | 3.40 ms | 0.24 us |
| Gene lookup, minus strand | 14.78 ms | 0.24 us |
| Row iteration, 200,000 rows | 2.85 s | 0.019 s |
| Index build | n/a | 0.09 s |

## Proof that `legacy` does not move the output

`validate_identity.py` compares the frozen original against the new `legacy` mode, element by
element, on 1,140 junctions across 19 branches at 60 junctions per branch.

**1,140 of 1,140 annotations match. 0 mismatches.**

The branches cover both strands and hit: both splice sites known, donor novel, acceptor novel, both
novel, a position inside an intron within 50 nt of each edge, a position mid-intron, a position
upstream of the first exon, a position downstream of the last exon, and scaffold contigs.

`scale_compare.py` repeats the comparison on 20,000 junctions. **20,000 of 20,000 match.**

On the real Cal27P5 test BED files, `legacy` matches the frozen original and the stored 2025-07-26
output on **8,503 of 8,503** junctions, including every count column.

Combined, 29,643 of 29,643 annotations match.

## What `corrected` changes, and why you may want it

The old scan carries two faults, and I kept both under `legacy`.

1. The scan never tests the chromosome. A chr2 position returned ENSG00000081913, which the
   reference file places on chr18.
2. The scan can never return a minus-strand gene. The loader reverses minus-strand coordinates, so
   the range test reads `high <= pos <= low`. That test holds for 0 of 32,995 minus-strand genes.
   Every minus-strand fully novel junction therefore reads `None:None-None`.

A third fault sits in `annotate.py`. The old code called the gene scan twice and passed `start` both
times. Line 64 handles the junction end, so it wants `end`. `corrected` passes `end`. `legacy`
reuses the first result, which keeps the output identical and halves the work.

`corrected` also spans each gene from its lowest to its highest coordinate. That span differs from
`geneData[gene][-1][1]` for 376 of 35,013 plus-strand genes, whose exon blocks overlap.

On the 1,140-junction validation set, `corrected` changes 680 junctions:

| Count | Change |
|---|---|
| 322 | Both modes assign a gene, but a different gene or exon |
| 265 | `legacy` returns None, `corrected` assigns a gene |
| 93 | `legacy` assigns a gene, `corrected` returns None |

Example, where `legacy` names a gene on the wrong contig:

```
legacy    LRG_410:I19.1_107310-I19.1_107310=chrLRG_766:107310-108056
corrected LRG_766:E82.1_107310-E83.1_108056=chrLRG_766:107310-108056
```

Those 680 of 1,140 come from a junction set that deliberately over-samples novel positions. The
fraction in a real dataset will be far lower, because most junctions match the reference at both
splice sites and never reach this code. I have not measured the real fraction.

## How to run

```
# corrected mode, the default since 2026-08-09
python3 -m altanalyze3 aggregate --juncounts <dir> --output <prefix> --ref <ref.bed>

# reproduce a pre-2026-08 result
python3 -m altanalyze3 aggregate --juncounts <dir> --output <prefix> --ref <ref.bed> --novel-gene-mode legacy
```

## Files

- `profile_annotate.py` times `annotate_junctions` per junction class and prints a cProfile table.
- `verify_semantics.py` reports what the old scan can and cannot match.
- `check_claims.py` counts unreachable genes and overlapping exon blocks.
- `scale_test.py` measures row iteration and dense export.
- `prototype_index.py` checks the sweep index against the old scan on 6,200 probe positions.
- `validate_identity.py` proves branch-by-branch identity and reports the `corrected` difference.
- `scale_compare.py` compares old against new at 20,000 junctions.
- `annotate_reference_frozen.py` holds the original function.

All files sit in `/Users/saljh8/Documents/GitHub/altanalyze3/dev_benchmarks/junction_annotation/`.

## What I did not do

- I did not speed up `importEnsemblGenes`. I profiled it on 2026-08-09: the file read takes 0.069 s
  and the Python loop takes 1.798 s. A vectorized parse still builds about 2 million dictionary
  entries, so it saves roughly 1 s, or 2x. The long-read pipeline shares the function.
- I did not touch `findNovelSpliceSite`. The long-read collapse pipeline shares it, and it costs
  roughly 100 us per novel splice site.
- I did not touch `export_dense_matrix`. I measured 0.49 s for 200,000 junctions by 20 samples.
- I did not optimize `collect_all_uids` or `collect_sample_counts`. Both read every BED file, so the
  pair reads each file twice.
- `altanalyze3/components/aggregate/test_aggregate.py` fails to import. It calls
  `annotate_junctions_and_export`, which the frozen original never defined. The failure predates my
  edits.
