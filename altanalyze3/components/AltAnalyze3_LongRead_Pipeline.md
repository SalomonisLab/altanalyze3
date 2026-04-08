# AltAnalyze3 Long-Read Single-Cell Pipeline
## One Workflow for Isoform and Splicing Biology

### Overview
AltAnalyze3 provides a consolidated, end-to-end long-read single-cell workflow for population-level splicing and isoform analysis. The pipeline starts from scKINNEX BAM files and carries data through quantification, annotation, differential analysis, and visualization, with optional integration of short-read/Ribo-Seq evidence and variant-based clonal analysis. It is implemented in open-source Python and supports Dockerized deployment for reproducibility.

### Major Functional Modules
1. **scKINNEX BAM extraction and quantification**
   - Extracts long-read evidence from BAM files.
   - Produces junction-level and gene-level quantification matrices suitable for downstream single-cell analyses.

2. **Cell-state assignment on gene-level data**
   - Uses **cellHarmony** for reference-guided cell-state alignment.
   - Supports ambient RNA correction before/within assignment workflows.
   - Enables robust transfer of biological labels for downstream isoform/splicing comparisons.

3. **Isoform model consolidation**
   - Performs GFF isoform collapsing to define non-redundant transcript models.
   - Applies exon annotation with published **AltAnalyze/Ensembl** annotation logic.

4. **Pseudobulk splicing quantification**
   - Aggregates cells by biological and/or technical replicate groups.
   - Computes replicate-aware splicing metrics for robust condition-level inference.

5. **Differential isoform and splicing analysis**
   - Runs isoform differential expression across groups.
   - Runs differential splicing analyses at junction/exon/isoform-resolution (context-dependent).
   - Supports replicate-aware testing and effect-size-based ranking.

6. **Isoform functional interpretation**
   - Annotates isoform protein-coding potential.
   - Supports differential isoform interpretation in coding/non-coding context.

7. **Visualization and interpretation**
   - Visualizes isoform structures.
   - Supports isoform clustering-based views.
   - Supports pairwise junction-restriction views for targeted event interrogation.

8. **Cross-modal evidence integration**
   - Integrates splicing evidence from bulk short-read RNA-Seq, bulk long-read RNA-Seq, and Ribo-Seq on shared transcript models.
   - Improves confidence in isoform/event support across assays.

9. **Optional clonal analysis**
   - Extracts genomic variants for clonal structure analysis when requested.

10. **SNAF integration**
   - Outputs can be transferred to SNAF.
   - Supports domain-preservation interpretation.
   - Supports extracellular/transmembrane isoform prediction for **ExNeoEpitope** discovery.

### Inputs and Outputs (High Level)
- **Inputs:** scKINNEX BAMs, references/annotations, sample metadata, optional bulk/Ribo-Seq data, optional variant settings.
- **Core outputs:** gene/junction quantifications, collapsed isoform models, annotated transcript features, pseudobulk matrices, differential isoform/splicing results, visualization-ready files, optional variant/clonal outputs, SNAF-compatible exports.

### Design Goals
- Single integrated workflow rather than disconnected tools.
- Reproducible and portable execution via Python + Docker.
- Multi-resolution biology: cell-state assignment, splicing, isoform function, and translational relevance in one system.
