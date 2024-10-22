Quickstart Tutorial
===================

This section demonstrates how to preprocess long-read data and perform differential analyses.

**Example Workflow:**

.. code-block:: python

    import altanalyze3.components.long_read.isoform_matrix as iso
    import altanalyze3.components.long_read.isoform_automate as isoa
    import altanalyze3.components.long_read.comparisons as comp

    # Define paths
    metadata_file = "/path/to/sample_metadata.txt"
    ensembl_exon_dir = "/path/to/EnsMart91/Hs_Ensembl_exon.txt"

    # Preprocess samples
    isoa.pre_process_samples(metadata_file, ensembl_exon_dir)

    # Perform differential analysis
    condition1 = 'healthy'
    condition2 = 'AML'
    comp.compute_differentials(metadata_file, condition1, condition2)
