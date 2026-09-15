# How to use scALABLE-viewer

scALABLE-viewer explores a single-cell atlas whose analyses already exist. You upload nothing
and launch nothing. The dataset opens by itself, and one viewer can hold several datasets.

scALABLE, the analysis tool, is a different program. It takes raw data, aligns it and runs the
comparisons; the viewer serves the finished result through the same interface. The controls of
`Explore`, `Differential` and `Chat` are the ones described in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/HOW_TO_USE.md`.
This guide covers what the viewer adds or changes. The methods behind every number are in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/README.md`.

## The four tabs

| Tab | Content |
| --- | --- |
| `Study` | title, abstract, assay, organism, technology, cell count, links to the raw data and the paper; lists of tools, downloads and samples |
| `Explore` | two plot panels over the same cells |
| `Differential` | the precomputed group comparisons |
| `Chat` | questions in plain language, answered from the dataset |

Explore opens first. When the viewer serves more than one dataset, a dataset menu switches
between them.

## Study

The Study tab reads the LungMAP record of the study the bundle belongs to. When no record is
configured, the tab says so instead of showing another study's record. A metacell dataset states
its metacell count and how many donors each metacell mixes.

## Explore

The plot types match the analysis tool: `UMAP cell types`, `UMAP broad`, `Cell frequency`,
`UMAP`, `Violin`, `DotPlot`, `CombPlot`, `MarkerHeatmap`, `MarkerNetwork` and
`Cell communication`. Three things differ.

| Difference | Behaviour |
| --- | --- |
| Modality menu | lists the modalities the bundle carries, for example RNA, ADT, lipid, GRN; a modality predicted per cell state, not per cell, does not appear here, because every cell of a state would show one value |
| Violin covariate | the violin can group by any categorical covariate of the dataset, not only by cell state |
| CombPlot bands | dark bands under the bars mark each donor column's disease status, group, sex and smoking status when the dataset records them |

`Select gene` accepts the names a reader knows, for example a lipid's common name, and the
viewer resolves them to the stored feature. `Download PDF` saves an editable vector PDF.

## Differential

A `Precomputed comparison` menu lists the contrasts the bundle carries for the chosen modality.
Choosing one loads its tables; nothing recomputes. The `Modality` menu shows only modalities
that have at least one contrast, and switching modality keeps the same comparison when it exists
on the new modality.

The views are `Summary`, `Heatmap`, `Volcano`, `GO Terms`, `Network` and `Gene Detail`, with
`Filter by gene` as in the analysis tool. GO terms and networks exist for RNA contrasts only.

The GO Terms plot colours each term by tier: `Representative` terms are the ones GO-Elite kept
after removing redundant parents and children, `Significant` terms passed the FDR cut of 0.05 but
folded into a broader term, and other terms failed the cut. A term can be significant without
being representative.

A comparison lists only the cell states it tested. A state needs at least 2 samples on each
side; a state with fewer is absent, not reported at low confidence. The comparison table gives
the sample count per side.

## Chat

Type a question or click an example. The viewer sends the question to an assistant service that
picks one of 17 named analyses and fills in the gene, cell state, comparison or covariate it
names. The viewer then computes the answer from the bundle's own statistics. The assistant never
sees a value and never states a number.

| Ask about | Example | Answer |
| --- | --- | --- |
| a cell state's markers | What are the best marker genes of AT2 cells? | a table and a DotPlot |
| two cell states | What distinguishes AT1 from AT2 cells? | a table and a DotPlot |
| a gene | Where is SFTPC expressed? | the top states by mean and a DotPlot |
| a comparison in a state | Which genes change in COPD versus control in AT2 cells? | the top rows and a volcano |
| a clinical variable | Which genes track FEV1 in AT2 cells? | a per-donor correlation table and a gradient plot |
| co-expression | Which genes co-vary with SFTPC across donors? | the top partners and a CombPlot |
| composition | Which cell states shift in GOLD IV versus GOLD I, II? | per-donor fractions and a frequency plot |
| donor heterogeneity | Do all COPD donors show the AT2 signature? | a per-donor score and a signature plot |
| the most affected state | Which cell type is most affected in COPD versus control? | states ranked by significant genes |
| pathways | Which pathways change in AT2 cells? | GO-Elite terms as bars |
| regulators | Show me the transcriptional targets of a regulator in AT2 cells | a TF-to-target network, or TF activity bars |

Below the answer, `Table` and `Plot` toggle the view, and up to four follow-up questions appear.

| Reply | Meaning |
| --- | --- |
| a list of states or comparisons | the question named none of them; pick one |
| `not covered` | the comparison exists but did not test that state; the reply lists the states it did test |
| `unsupported` | the question asks for a modality the bundle lacks; the reply names the ones it has |
| HTTP 503 | the assistant service is down |

On 2026-08-26, 55 of 55 protocol questions and 68 of 68 paraphrased questions on the
COPD-metacells dataset returned the intended analysis; the record is in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/VALIDATION.md`.
The suites covered one RNA-only dataset and drove the endpoint, not the rendered page.

## What the numbers mean

A cell in this atlas may be a metacell: several real cells summed after matching on donor
metadata, which protects participant privacy and reduces noise. Differential tests run on
pseudobulks, one per sample and cell state, so each point in a comparison is a sample, not a
cell. Imputed modalities are model predictions from the RNA, not measurements; the analysis
tool's README names each model and its held-out accuracy.

## Things that surprise people

| Symptom | Cause |
| --- | --- |
| a cell state is missing from Differential | the comparison had fewer than 2 samples on one side for that state |
| a modality is missing from Explore but present in Differential | its values exist per cell state only |
| MarkerHeatmap takes a moment | it loads a large matrix into an embedded viewer from the Broad Institute's server |
| a gene shows nothing | the dataset did not detect it; check the suggestion list |
| a chat answer says the assistant is unavailable | the service on port 8001 is down |

## Getting a figure out

Use `Download PDF` rather than a screenshot. The PDF holds vector shapes and editable text at
any size. With `Filter by gene` on, the button saves the filtered figure you see.

## Where the data came from

The Study tab names the source study, the raw-data accession and the publication. The
`README.md` beside this file describes how to build, validate and serve a bundle.
