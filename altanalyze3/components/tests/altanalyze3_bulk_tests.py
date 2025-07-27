import os, sys
from pathlib import Path
from altanalyze3.utilities.api import run

import altanalyze3.components.long_read.comparisons as comp
import altanalyze3.components.long_read.gff_process as gff_process
from altanalyze3.components.aggregate.main import aggregate

# Config paths
gtf_path = Path("reference/gencode.v45.annotation.gff3")    # GTF input
index_output_dir = Path("hg38")                             # Where to save BED index
bam_dir = Path("BAMs")                                      # Folder with .bam files
output_dir = Path("BEDs")                                   # Output base
output_dir.mkdir(exist_ok=True, parents=True)
ref_bed = index_output_dir / "gene_model.bed.gz"
counts = output_dir / "aggregated_counts"
counts_dir = counts.with_name(counts.stem + '_annotated.tsv')
psi_dir = output_dir / "psi_values.txt"

# Build index from GTF
run("index", gtf=gtf_path, output=index_output_dir, threads=4, cpus=1)

# Run juncount and intcount on all BAMs
for bam_file in bam_dir.glob("*.bam"):
   sample = bam_file.stem.replace(" ", "_")
   junction_out = output_dir / f"{sample}_juncounts"
   intron_out = output_dir / f"{sample}_intcounts"
   run("juncount", bam=bam_file, output=junction_out, threads=2, cpus=1)
   run("intcount", bam=bam_file, ref=ref_bed, output=intron_out, threads=2, cpus=1)

# Aggregate all count files ===
junc_beds = list(output_dir.glob("*_juncounts.bed"))
int_beds = list(output_dir.glob("*_intcounts.bed"))

class Args:
    def __init__(self):
        self.juncounts = junc_beds
        self.intcounts = int_beds
        self.ref = ref_bed
        self.output = counts
        self.threads = 4
        self.cpus = 8
        self.chr = None  # or specify if needed
args = Args()
aggregate(args)

import asyncio
import altanalyze3.components.psi.psi_single as psi
asyncio.run(psi.main(junction_path=counts_dir, query_gene=None, outdir=psi_dir))