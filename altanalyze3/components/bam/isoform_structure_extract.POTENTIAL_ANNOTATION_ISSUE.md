# Deferred potential exon-annotation inconsistency: LAT

Flag recorded 2026-10-17 beside `isoform_structure_extract.py`, the BAM-to-read-GFF/molecule-h5ad extraction entry point. **The origin of this discrepancy has not been traced to this script or to BAM parsing. No extraction, collapse, quantification or annotation algorithm is changed by this note.**

The user instructed: “Document this as a potential issue along side the bam to h5ad associated script as a flag and move on for now,” followed by “ignore it.” Investigation and the proposed coordinate-based junction-association change are deferred; continue the existing analysis protocol.

## Exact observed example

KINNEX-5 retained transcript `59726342.MDS5801-postHMA_CD34`, LAT / `ENSG00000213658`, chr16, positive strand. The current staged `combined.gff.gz` contains these 13 exon intervals (GFF, 1-based inclusive):

```text
28985049-28985115
28985418-28985517
28985713-28985740
28985854-28985888
28986135-28986216
28986382-28986446
28986540-28986569
28986657-28986714
28986799-28986893
28989527-28989589
28989774-28989839
28989933-28990019
28990189-28990454
```

Its `transcript_associations.txt.gz` exon string is:

```text
E2.1|E2.2|E2.3|E2.4|E2.5|E2.6|E2.7|E2.8|E2.9|E2.10|I2.1|E3.1|I3.1|E4.1|E4.2|E4.3|E4.4|E4.5|ENSG00000261067:E11.1_28985115|ENSG00000261067:E11.1|ENSG00000261067:E12.1|ENSG00000261067:E13.1|ENSG00000261067:E14.1|E5.1|E5.3|ENSG00000261067:E15.1|ENSG00000261067:E16.1|ENSG00000261067:E17.1|ENSG00000261067:E18.1|ENSG00000261067:E19.1|E9.3
```

The GFF splice sites 28985740 and 28985854 correspond to the counted junction `ENSG00000213658:E3.1-E4.2=chr16:28985740-28985854`. The string instead contains `E3.1|I3.1|E4.1|E4.2`. `annotation/junction_isoform.py::load_transcript_associations` constructs junction keys from adjacent string entries, so this transcript is absent from the `E3.1-E4.2` candidate list.

The same parent exon string was already present before tertiary cleanup. Merged transcript `54396929.EP_SFSF2_AZA` had the adjacent `E3.1|E4.2` pair. This establishes a coordinate/string discrepancy that predates tertiary cleanup. It does not establish a BAM parser defect, incorrect read aggregation, or its genome-wide prevalence.

## Interpretation and deferred investigation

The possible direct consequence is a missing or incorrect transcript/protein annotation beside a splicing differential result, in fields such as `Isoform_1|Length`, `Isoform_2|Length`, `Event_Type` and `Event_Direction`. An effect on collapse or molecule assignments remains unproven. A read-only, one-link diagnostic changed none of the selected transcript pairs in 21 LAT event rows from one statistics table per campaign; this limited comparison neither validates the complete alternative method nor proves the issue harmless.

If revisited, trace this precise read/transcript through its original BAM alignment, per-source GFF and exon string, and each collapse phase before selecting a correction. Do not silently replace the established method with coordinate matching.

Exact extraction and diagnostic records are stored in the KINNEX-5 analysis directory `isoform_collapse_validation_20261015/LAT_annotation_visual_example_20261017/`, including `LAT_RETAINED_ISOFORM_EXACT_RECORDS.json` and `LAT_ANNOTATION_EXAMPLE.json`.
