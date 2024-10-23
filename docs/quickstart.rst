Quickstart Tutorial
===================

This section demonstrates how to preprocess long-read data and perform differential analyses.

**Example Long-Read Workflow:**

.. code-block:: python

   import altanalyze3.components.long_read.isoform_matrix as iso
   import altanalyze3.components.long_read.isoform_automate as isoa

   metadata_file = "/path/to/metadata.txt"
   ensembl_exon_dir = "/path/to/Hs_Ensembl_exon.txt"
   barcode_cluster_dirs = ["/path/to/barcode_to_clusters.txt"]

   sample_dict = isoa.import_metadata(metadata_file)
   isoa.pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir)

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

   gene_symbol_file = "/path/to/Hs_Ensembl-annotations.txt"
   genome_fasta = "/path/to/genome.fa"

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