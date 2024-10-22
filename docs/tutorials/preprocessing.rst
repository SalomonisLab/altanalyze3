
Preprocessing Tutorial
======================

This tutorial will guide you through the preprocessing steps required for running AltAnalyze3 workflows, specifically for long-read and single-cell RNA-Seq data integration. Follow these steps to prepare your data properly for further analyses.

Prerequisites
-------------
1. **Required Tools**:
   - Python 3.11 or higher
   - AltAnalyze3 installed via pip (`pip install altanalyze3`)
   - Access to 10x Genomics data files or compatible GFF/FASTA files.

2. **Sample Metadata**:
   Ensure that your samples are properly annotated in a metadata file. A sample metadata file may look like this:

   ```
   Sample_ID    Condition    Cell_Type
   AM72         Healthy      HSC-1
   BF71         Diseased     MSC-Fibroblast-1
   ```

3. **Install Dependencies**:
   Use the following command to install dependencies:
   ```
   pip install -r docs/requirements.txt
   ```

Step-by-Step Preprocessing
--------------------------
1. **Prepare Metadata and Cluster Files**:
   You need metadata and barcode-cluster files for proper clustering. Example:
   ```
   /path/to/metadata.txt
   /path/to/cellHarmony-labels-young.txt
   ```

2. **Run Preprocessing Script**:
   In your Python environment, run:
   ```python
   import altanalyze3.components.long_read.isoform_automate as isoa

   metadata_file = "/path/to/metadata.txt"
   barcode_cluster_dirs = ["/path/to/cellHarmony-labels-young.txt"]

   isoa.pre_process_samples(metadata_file, barcode_cluster_dirs, "/path/to/ensembl_exon_dir")
   ```

3. **Combining Processed Samples**:
   Once preprocessed, combine them using:
   ```python
   isoa.combine_processed_samples(
       metadata_file,
       barcode_cluster_dirs,
       "/path/to/ensembl_exon_dir",
       "/path/to/gencode.gff",
       "/path/to/genome.fa"
   )
   ```

4. **Verify Output**:
   Ensure that the processed outputs include files with combined splicing, isoform, and ratio data.

Next Steps
----------
After preprocessing, you are ready to run **differential analyses** or perform **ICGS splicing analysis**. See the relevant tutorials for these steps.

Support
-------
For issues, please refer to our GitHub repository:  
https://github.com/SalomonisLab/altanalyze3
