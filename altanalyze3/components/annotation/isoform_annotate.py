import os,sys
import pandas as pd
import glob
from tqdm import tqdm

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ..long_read import gff_process as gff_process
from ..long_read import isoform_translation as isot


def main(gff_source, genome_fasta, ensembl_exon_dir=None, ref_gff_file=None, outdir=None):
    """Assign isoform annotations to the statistical output from comparisons.py.

    ``gff_source``       list of GFF/GTF files describing the isoforms to annotate.
    ``genome_fasta``     genome the isoforms are translated against.
    ``ensembl_exon_dir`` exon reference used to collapse long-read GFFs to unique splice-site
                         combinations. Omit it for bulk, where the GFF is already a reference
                         model and no collapse applies.
    ``ref_gff_file``     reference GFF used to recover first exons. Optional.
    ``outdir``           where the sequences and summaries land. Defaults to the directory the
                         collapse step produced, else the working directory.

    Returns ``(query_gff_file, transcript_associations_file)``.

    The three trailing arguments used to be required, so the only caller in the tree raised
    TypeError. They are optional now, and the bulk path runs with a GFF and a genome alone.
    """
    transcript_associations_file = None
    gff_output_dir = None

    #### Collapse the long-read GFFs and find the longest consensus isoforms. Bulk skips this.
    if ensembl_exon_dir:
        output = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)
        if '.txt' in output:
            transcript_associations_file = output
            gff_output_dir = os.path.dirname(transcript_associations_file)
        else:
            gff_output_dir = output
            transcript_associations_file = os.path.join(gff_output_dir, "transcript_associations.txt")

    #### Pick the GFF to translate: the collapsed one when it exists, else the input as given.
    if gff_output_dir and os.path.isfile(os.path.join(gff_output_dir, "combined.gff")):
        query_gff_file = os.path.join(gff_output_dir, "combined.gff")
    else:
        query_gff_file = gff_source[0] if isinstance(gff_source, (list, tuple)) else gff_source

    outdir = outdir or gff_output_dir or os.getcwd()
    os.makedirs(outdir, exist_ok=True)

    print(f"Translating {query_gff_file}")
    print(f"  reference GFF: {ref_gff_file}")
    print(f"  transcript associations: {transcript_associations_file}")
    print(f"  output directory: {outdir}")

    #### Translate each isoform into protein sequence and check for NMD
    cds_records, transcript_records, protein_records = isot.gff_translate(
        query_gff_file, genome_fasta, ref_gff_file, transcript_associations_file, outdir=outdir
    )

    from Bio import SeqIO
    with open(os.path.join(outdir, "protein_sequences.fasta"), "w") as protein_file:
        SeqIO.write(protein_records, protein_file, "fasta")

    with open(os.path.join(outdir, "transcript_sequences.fasta"), "w") as cds_file:
        SeqIO.write(transcript_records, cds_file, "fasta")

    with open(os.path.join(outdir, "cds_sequences.fasta"), "w") as cds_file:
        SeqIO.write(cds_records, cds_file, "fasta")

    print(f"Translated {len(protein_records)} isoforms to protein")
    return query_gff_file, transcript_associations_file
