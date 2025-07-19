import pandas as pd
import pysam
import os
import argparse
import logging
from pathlib import Path

def filter_and_index_bed(input_file, output_prefix, keep=False):
    df = pd.read_csv(input_file, sep="\t", header=None)

    # Filter to introns
    df = df[df[1].str.startswith('I')]

    # Create correctly ordered BED dataframe
    bed = pd.DataFrame({
        0: df[2],                           # chrom
        1: df[4],                           # start
        2: df[5],                           # end
        3: df[1],                           # exon (e.g., I1.1)
        4: df[0],                           # gene (e.g., ENSG00000223972)
        5: df[3]                            # strand
    })

    # Sort
    bed = bed.sort_values([0, 1, 2])

    # File names
    bed_file = output_prefix + ".bed"
    bed_gz = bed_file + ".gz"

    # Save BED
    bed.to_csv(bed_file, sep="\t", header=False, index=False)

    # Compress and index
    pysam.tabix_compress(bed_file, bed_gz, force=True)
    pysam.tabix_index(bed_gz, preset="bed", force=True)

    if not keep:
        os.remove(bed_file)

    print(f"Created {bed_gz} and {bed_gz}.tbi")

def build_index(args):
    """
    Build an index from a GTF/GFF file.
    This function is called by the parser.
    """
    try:
        # Ensure output directory exists
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # First, run gene_model/main.py to create the TSV file
        from altanalyze3.components.gene_model.main import main as gene_model_main
        
        # Create a namespace object with the required arguments
        class GeneModelArgs:
            def __init__(self, gtf, gene, outdir):
                self.gtf = gtf
                self.gene = gene
                self.outdir = outdir
        
        gene_model_args = GeneModelArgs(args.gtf, 'all', str(output_dir))
        gene_model_main(gene_model_args)
        
        # Now create the index from the TSV file
        tsv_file = output_dir / "gene_model_all.tsv"
        if not tsv_file.exists():
            raise FileNotFoundError(f"Expected TSV file not found: {tsv_file}")
        
        # Create output prefix for the index
        output_prefix = output_dir / "gene_model"
        
        # Process the TSV file to create the index
        filter_and_index_bed(str(tsv_file), str(output_prefix), keep=False)
        
        logging.info(f"Successfully created index in {output_dir}")
        
    except Exception as e:
        logging.error(f"Failed to create index: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter intron rows and create compressed indexed BED file.")
    parser.add_argument("--input", required=True, help="Input tab-delimited file (e.g., Hs_Ensembl_exon.txt)")
    parser.add_argument("--output", default="hs_ref", help="Output prefix (default: hs_ref)")
    parser.add_argument("--keep", action="store_true", help="Keep uncompressed .bed file")
    args = parser.parse_args()

    filter_and_index_bed(args.input, args.output, args.keep)
