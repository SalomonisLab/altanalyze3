"""Build the isoform lookup files a bulk run needs, from an Ensembl build.

``junction_isoform.annotate`` associates a differential junction with the isoforms that carry it.
It reads two files that the long-read pipeline normally produces from observed GFFs:

``transcript_associations.txt``  gene, strand, exon sequence, isoform, source
``protein_summary.txt``         gene, transcript, protein length, NMD status, intron retention,
                                longest isoform length

A bulk cohort has no long-read GFF, so those files have to come from the reference. Both exist
inside an AltAnalyze Ensembl build already, in a different shape:

``mRNA-ExonIDs.txt``                     gene, transcript, protein, exon sequence
``Hs_ProteinCoordinates_build_<b>.tab``  per-protein coding blocks, whose highest AA position is
                                         the protein length
``Hs_Ensembl_Protein__<b>.tab``          gene, transcript, protein

Building from the same directory that annotated the junctions keeps exon identifiers consistent.
Mixing builds silently lowers the match rate: on one measured chromosome the examined junction
matched an isoform 58.7% of the time under Ensembl 91, 50.1% under Ensembl 100 and 28.0% under
gencode v45.
"""

import logging
import os
from collections import defaultdict


def load_gene_strands(exon_file):
    """Strand per gene, read from the build's exon model."""
    strands = {}
    with open(exon_file) as handler:
        for line in handler:
            fields = line.split("\t", 4)
            if len(fields) < 4 or not fields[0].startswith("ENS"):
                continue
            strands.setdefault(fields[0], fields[3])
    return strands


def load_protein_lengths(protein_coordinates_tab):
    """Protein length per Ensembl protein.

    The columns are named AA_NT_Start and AA_NT_Stop but hold amino-acid positions, which the
    Ensembl-100 record for ENSP00000349748 confirms: its highest value is 707, the length our own
    translation produced for ENST00000357214, and a domain on that protein ends at residue 538.
    """
    lengths = {}
    with open(protein_coordinates_tab) as handler:
        handler.readline()
        for line in handler:
            fields = line.split("\t", 5)
            if len(fields) < 4:
                continue
            try:
                stop = int(fields[3])
            except ValueError:
                continue
            if stop > lengths.get(fields[0], 0):
                lengths[fields[0]] = stop
    return lengths


def load_transcript_to_protein(mrna_exon_file):
    """Gene and protein per transcript, from ``mRNA-ExonIDs.txt``.

    Its third column holds the Ensembl protein for a coding transcript, ``<transcript>-PEP`` for a
    non-coding one and ``None`` where there is no translation, so only the first kind is kept.

    The build's ``Hs_Ensembl_Protein__<b>.tab`` would be the obvious source, but in the Ensembl 91
    build that file is a bare list of protein identifiers with no transcript column, so it cannot
    supply this mapping.
    """
    mapping = {}
    with open(mrna_exon_file) as handler:
        for line in handler:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2].startswith("ENSP"):
                mapping[fields[1]] = (fields[0], fields[2])
    return mapping


def build_transcript_associations(ensembl_dir, output_file, source_label=None):
    """Write ``transcript_associations.txt`` from ``mRNA-ExonIDs.txt``."""
    mrna_file = os.path.join(ensembl_dir, "mRNA-ExonIDs.txt")
    exon_file = os.path.join(ensembl_dir, "Hs_Ensembl_exon.txt")
    for path in (mrna_file, exon_file):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required reference file not found: {path}")
    source_label = source_label or os.path.basename(os.path.normpath(ensembl_dir))

    strands = load_gene_strands(exon_file)
    written = 0
    missing_strand = 0
    with open(mrna_file) as source, open(output_file, "w") as target:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or not fields[3]:
                continue
            gene, transcript, _protein, exons = fields[0], fields[1], fields[2], fields[3]
            strand = strands.get(gene)
            if strand is None:
                missing_strand += 1
                strand = "+"
            target.write(f"{gene}\t{strand}\t{exons}\t{transcript}\t{source_label}\n")
            written += 1
    logging.info(f"Wrote {written:,} transcript associations to {output_file}")
    if missing_strand:
        logging.warning(f"{missing_strand:,} rows had no strand in the exon model; wrote '+'")
    return written


def build_protein_summary(ensembl_dir, output_file, build="91_38"):
    """Write ``protein_summary.txt`` from the build's protein tables.

    NMD status and intron retention come from translating an observed transcript, which a
    reference build does not record, so both columns read 'Not-NMD' and 'False'. A caller that
    needs real values must translate the isoforms.
    """
    coordinates = os.path.join(ensembl_dir, f"Hs_ProteinCoordinates_build_{build}.tab")
    mrna_file = os.path.join(ensembl_dir, "mRNA-ExonIDs.txt")
    for path in (coordinates, mrna_file):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required reference file not found: {path}")

    lengths = load_protein_lengths(coordinates)
    transcript_to_protein = load_transcript_to_protein(mrna_file)

    longest = defaultdict(int)
    for transcript, (gene, protein) in transcript_to_protein.items():
        length = lengths.get(protein)
        if length and length > longest[gene]:
            longest[gene] = length

    written = 0
    no_length = 0
    with open(output_file, "w") as target:
        target.write("Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\t"
                     "Longest Isoform Length\n")
        for transcript, (gene, protein) in transcript_to_protein.items():
            length = lengths.get(protein)
            if not length:
                no_length += 1
                continue
            target.write(f"{gene}\t{transcript}\t{length}\tNot-NMD\tFalse\t{longest[gene]}\n")
            written += 1
    logging.info(f"Wrote {written:,} protein records to {output_file}")
    if no_length:
        logging.warning(f"{no_length:,} transcripts had a protein with no coordinate record")
    return written
