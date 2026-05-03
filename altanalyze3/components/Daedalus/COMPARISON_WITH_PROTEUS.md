# Daedalus versus Proteus

This document compares the two local projects as they exist now in the
repository. The comparison is grounded in the current documentation and current
available outputs, not on planned future work.

## Scope

- `Daedalus`
  - benchmark substrate, weak-label corpus, leakage-resistant baselines, and
    reference-conditioned tabular/multitask baseline models
- `Proteus`
  - multimodal reference-conditioned architecture using frozen Orthrus, ESM2,
    and disorder encoders, plus a surfaceome validation path

## Current project state

| Project | Current state | Evidence |
|---|---|---|
| `Daedalus` | implemented through Phase A and Phase B with frozen splits and completed held-out baseline reports | `Daedalus/phase_a/FEASIBILITY.md`, `Daedalus/phase_b/STATUS.md`, `Daedalus/phase_b/checkpoints/reference_delta_*` |
| `Proteus` | architecture, methods, datasets, and training workflow documented; supervised fine-tuning not yet completed | `Proteus/README.md`, `Proteus/PROTEUS_METHODS.md` |

## Pros and Cons

| Project | Pros | Cons |
|---|---|---|
| `Daedalus` | Has real held-out benchmark numbers now; benchmark leakage has been actively checked and reduced; uses a simpler stack that is feasible on one H100; already integrates UniProt, HPA, BioGRID, ClinVar, APPRIS, MANE, and sequence-derived topology features; current docs describe how to run Phase A and B end to end | Current models are still relatively simple compared with the intended end-state; no frozen RNA/protein foundation encoders yet; mouse sequence-derived topology features are not yet fully rebuilt; external TF assay validation is planned but not yet wired in |
| `Proteus` | Stronger architectural ambition; explicitly uses Orthrus for RNA and ESM2 for protein; includes disorder modeling and a surfaceome validation module; better positioned conceptually for a Nature Methods narrative if the supervised results land; includes synthetic pretraining design | No completed supervised benchmark outputs yet; benchmark table is still `TBD`; depends on heavier infrastructure and more moving parts; current claims are partly forward-looking rather than backed by current held-out reports; surfaceome/topology logic still has documented simplifications |

## Strengths and Weaknesses

| Project | Strengths | Weaknesses |
|---|---|---|
| `Daedalus` | Strong data engineering discipline; reproducible Phase A/Phase B workflow; actual measurable feasibility on current hardware; clean benchmark framing for `global`, `membrane`, `surface`, `kinase`, and `transcription_factor` tasks | Method novelty is currently more in benchmark construction and reference-conditioned formulation than in model architecture; without the next model stage it risks looking like a strong baseline project rather than the final methods contribution |
| `Proteus` | Better model-level novelty; clearer multimodal fusion story; more explicit reference-conditioned delta formulation; stronger direct path to transcript-plus-protein reasoning | Scientific risk is higher because the implemented benchmark evidence is still missing; more vulnerable to encoder/dependency friction; current validation language is ahead of the completed results |

## Initial benchmarking and validation state

### Daedalus

Daedalus has completed held-out benchmark reports on gene-disjoint splits.

#### Phase B sklearn logistic baseline

| Task | AP | AUROC |
|---|---:|---:|
| `global` | 0.810107 | 0.845348 |
| `membrane` | 0.809486 | 0.843945 |
| `surface` | 0.802420 | 0.835832 |
| `kinase` | 0.853772 | 0.874458 |
| `transcription_factor` | 0.827434 | 0.842662 |

#### Phase B torch multitask baseline

| Task | AP | AUROC |
|---|---:|---:|
| `global` | 0.840682 | 0.882048 |
| `membrane` | 0.839828 | 0.875073 |
| `surface` | 0.837000 | 0.871491 |
| `kinase` | 0.866310 | 0.903348 |
| `transcription_factor` | 0.858565 | 0.884078 |

Interpretation:

- `Daedalus` already shows that the current weak-label tasks are learnable
- topology/surface performance became credible only after adding topology,
  signal, glyco, cysteine, PPI, and PTM context
- the torch multitask baseline is consistently stronger than the sklearn
  baseline

### Proteus

Proteus does not currently provide completed supervised benchmark outputs in its
documentation. Its methods document is explicit about that:

> Proteus v1.0 has not yet completed supervised fine-tuning. The benchmarks below describe the evaluation protocol and baselines.

The current benchmark section therefore contains:

- evaluation protocol
- proposed baselines
- an expected benchmark table with `TBD` values

What Proteus does report now:

- synthetic pretraining event counts
  - total synthetic events: `1,960,432`
  - `domain_disrupted = 214,122`
  - `nmd_triggered = 1,143,796`
  - `tm_lost = 34,142`
- structural feature coverage checks
  - structural features non-zero: `290,255 / 290,255`
  - disorder features non-zero: `289,697 / 290,255`

Interpretation:

- `Proteus` has a reasonable pretraining/feature-quality story
- `Proteus` does not yet have a completed held-out supervised performance story
- so it is currently stronger as an architectural proposal than as a proven
  benchmarked system

## Practical conclusion

If the goal is to decide which project is currently more credible as an
implemented result:

- `Daedalus` is ahead

If the goal is to decide which project currently has the more ambitious model
story:

- `Proteus` is ahead

The sensible development path is not to choose one and discard the other.
The sensible path is:

1. keep `Daedalus` as the benchmark and reproducibility substrate
2. treat `Proteus` as the richer model layer to be evaluated on that substrate
3. only elevate `Proteus` claims once it has held-out results comparable to
   `Daedalus`

## Recommended positioning

For internal decision-making:

- `Daedalus` = reliable benchmark/data/validation backbone
- `Proteus` = higher-risk, higher-upside model program

For manuscript strategy:

- do not present `Proteus` as the stronger result until supervised fine-tuning
  and held-out evaluation are complete
- use `Daedalus` metrics as the current implemented baseline state
- use `Proteus` as the next-step multimodal model layer on top of that baseline
