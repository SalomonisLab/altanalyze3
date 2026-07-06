"""bam2hla - fast class-I HLA genotyping from a genome BAM via pysam pileup."""
from .bam2hla import type_bam, detect_build, format_optitype, load_signatures

__all__ = ['type_bam', 'detect_build', 'format_optitype', 'load_signatures']
