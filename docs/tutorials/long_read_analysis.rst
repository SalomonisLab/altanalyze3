
Long Read Single-Cell Analysis Tutorial
======================

This tutorial will guide you through the automated combined analysis of multiple long read single-cell datasets. Automation includes defaults for GFF merging and reference transcript selection, protein translation and statistical filtering. These steps can be further customized (alternative options, statistical thresholds) through separate calls to the underlying python functions.

Prerequisites
-------------
**Requirements**:

- Python 3.11 or higher
- AltAnalyze3 installed via pip (`pip install altanalyze3`)
- 10x Genomics long read matrices (.mtx) and associated GFF files

**Sample Metadata**:
Ensure that your samples are properly annotated in a metadata file. A sample metadata file may look like this:

.. code-block:: text

   uid     gff                   matrix              library       reverse    groups
   D001    /Diag1-1/D001.gff     /Diag1-1/sciso      D001-HSC      TRUE       Diagnosis
   D001    /Diag1-2/D001.gff     /Diag1-2/sciso      D001-MPP      TRUE       Diagnosis
   D002    /Diag2/D002.gff       /Diag2/sciso        D002-HSPC     TRUE       Diagnosis
   D003    /Diag3/D003.gff       /Diag3/sciso        D003-HSPC     TRUE       Diagnosis
   D004    /Relapse1/D004.gff    /Relapse1/sciso     D004-HSPC     TRUE       Relapse
   D005    /Relapse2/D005.gff    /Relapse2/sciso     D005-HSPC     TRUE       Relapse
   D006    /Relapse3/D006.gff    /Relapse3/sciso     D006-HSPC     FALSE      Relapse

Note: Multiple sequencing runs or libraries will be combined with the same uid. If the cell barcodes are reverse complemented from the barcode-to-cluster relationships (two column file), enter reverse as True.

**Install Dependencies**:
Use the following command to install dependencies:

.. code-block:: bash

   pip install altanalyze3

   curl -O https://altanalyze.org/isoform/Hs.zip

   unzip Hs.zip


Step-by-Step Preprocessing
--------------------------
**Prepare Metadata and Cluster Files**:
You need metadata and barcode-cluster files for cluster-guided analyses. Extract database files from the Hs.zip file. Example:

.. code-block:: text

   /path/to/metadata.txt
   /path/to/barcode_to_clusters.txt

   /path/to/gencode.annotation.gff3
   /path/to/Hs_Ensembl-annotations.txt
   /path/to/Hs_Ensembl_exon.txt
   /path/to/genome.fa

**Run Preprocessing Script**:
In your Python environment or script, run:
   
.. code-block:: python

   import altanalyze3.components.long_read.isoform_matrix as iso
   import altanalyze3.components.long_read.isoform_automate as isoa

   metadata_file = "/path/to/metadata.txt"
   ensembl_exon_dir = "/path/to/Hs_Ensembl_exon.txt"
   barcode_cluster_dirs = ["/path/to/barcode_to_clusters.txt"]

   sample_dict = isoa.import_metadata(metadata_file)
   isoa.pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir)

Note: ensembl_exon_dir, gene_symbol_file, genome_fasta and gencode_gff must pre-downloaded from the Hs.zip above.

**Combining Processed Samples**:
Once preprocessed, combine them using:

.. code-block:: python

   import altanalyze3.components.long_read.comparisons as comp
   gencode_gff = "/path/to/gencode.annotation.gff3"
   genome_fasta = "/path/to/genome.fa"

   isoa.combine_processed_samples(
      metadata_file,
      barcode_cluster_dirs,
      ensembl_exon_dir,
      gencode_gff,
      genome_fasta
   )

**Compute and Annotate Differential Splicing Events and Isoforms**:
Once preprocessed, combine them using:

.. code-block:: python

   gene_symbol_file = "/path/to/Hs_Ensembl-annotations.txt"

   # Import all cell clusters in order or replace with a list of select cluster(s)
   cluster_order = iso.return_cluster_order(barcode_cluster_dirs)

   # Differential analyses to perform
   analyses = ['junction', 'isoform', 'isoform-ratio']

   condition1 = 'Diagnosis'
   condition2 = 'Relapse'
   conditions = [(condition1, condition2)]

   comp.compute_differentials(
      sample_dict,
      conditions,
      cluster_order,
      gene_symbol_file,
      analyses=analyses
   )

**Expected Outputs**:

- *gff_output* - Directory of isoform exon structure and isoform mappings
- *sample.h5ad* - Anndata for each sample with consensus isoform or junctions IDs
- *protein_sequences.fasta*  - protein sequence for consensus isoforms
- *protein_summary.txt*  - isoform NMD prediction
- *isoform_combined_pseudo_cluster_tpm.txt* - Cluster-level pseudobulks TPMs
- *junction_combined_pseudo_cluster_counts.txt* - Junction, intron & 3' end counts
- *protein_summary.txt*  - isoform NMD prediction
- *psi_combined_pseudo_cluster_counts.txt* - PSI for junctions in >2 cluster pseudobulks
- *junction_combined_pseudo_cluster_counts.txt* - Junction, intron & 3' end counts
- *dPSI-events.txt* - Pairwise group Mann-Whitney U differential PSI events
- *dPSI-cluster/covariate* - Pairwise group Mann-Whitney U differential PSI events
- *diff-cluster/covariate-isoform* - Pairwise group Mann-Whitney U differential isoform log2 TPM
- *diff-cluster/covariate-ratio* - Pairwise group Mann-Whitney U differential isoform/gene ratios

**Verify Output**:
Ensure that the processed outputs include files with differential splicing, isoform, and ratio data in the current working directory.

Next Steps
----------
After preprocessing, you are ready to inspect your results in a spreadsheet editor, **Perform Secondary Analyses** or **Visualize Results**. See the relevant tutorials for these steps.

Support
-------
For issues, please refer to our GitHub repository:  
https://github.com/SalomonisLab/altanalyze3
