"""Decide which protein domains a splicing event gains or loses.

The reference is Ensembl 100, the same build the gene model comes from, so protein identifiers and
coordinates line up. ``Hs_ProteinFeatures_build_100_38.tab`` gives every domain an InterPro
accession together with its genomic span, and that span is what makes the comparison possible
without any protein sequence lookup.

A domain counts as present in an isoform when the isoform's coding blocks cover it. Coding blocks
are the isoform's exons clipped to its CDS, so a domain that falls in an intron of one isoform is
absent from that isoform even though it sits between the isoform's first and last coding base.

Task #7 turns the gained and lost lists into the ``ProteinPredictions`` string.
"""

import os
from collections import defaultdict, namedtuple

### ``intervals`` holds the domain's CODING intervals: its genomic span intersected with the coding
### blocks of the protein it was annotated on. The raw span covers the domain's own introns, so a
### multi-exon domain can never be fully covered by any isoform if the raw span is used directly.
Domain = namedtuple("Domain", "protein interpro name description aa_start aa_stop start stop intervals")


def load_protein_coding_blocks(protein_coordinates_tab):
    """Map each Ensembl protein to its coding blocks, from ``Hs_ProteinCoordinates_build_<b>.tab``."""
    blocks = defaultdict(list)
    with open(protein_coordinates_tab) as handler:
        handler.readline()                                        # protienID exonID AA_NT_Start AA_NT_Stop Genomic_Start Genomic_Stop
        for line in handler:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                start, stop = int(f[4]), int(f[5])
            except ValueError:
                continue
            if start > stop:                                      # minus-strand rows run high to low
                start, stop = stop, start
            blocks[f[0]].append((start, stop))
    for protein in blocks:
        blocks[protein].sort()
    return blocks


def intersect(intervals, start, stop):
    """Clip ``[start, stop]`` to ``intervals``."""
    out = []
    for lo, hi in intervals:
        a, b = max(lo, start), min(hi, stop)
        if b >= a:
            out.append((a, b))
    return out

### A domain must have this fraction of its genomic span covered by the isoform's coding blocks
### before it counts as present. 1.0 means every coding base of the domain survives.
DEFAULT_MIN_COVERAGE = 1.0


def load_protein_to_gene(ensembl_protein_tab):
    """Map every Ensembl protein to its gene, from ``Hs_Ensembl_Protein__<build>.tab``."""
    protein_to_gene, protein_to_transcript = {}, {}
    with open(ensembl_protein_tab) as handler:
        header = handler.readline()
        for line in handler:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[2]:
                continue
            gene, transcript, protein = fields[0], fields[1], fields[2]
            protein_to_gene[protein] = gene
            protein_to_transcript[protein] = transcript
    return protein_to_gene, protein_to_transcript


def load_domains(protein_features_tab, protein_to_gene, protein_blocks=None):
    """Group domain features by gene. Rows whose protein has no gene are reported, not dropped."""
    by_gene = defaultdict(list)
    protein_blocks = protein_blocks or {}
    unmapped = 0
    no_blocks = 0
    with open(protein_features_tab) as handler:
        handler.readline()                                        # ID AA_Start AA_Stop Start Stop Name Interpro Description
        for line in handler:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            protein = f[0]
            gene = protein_to_gene.get(protein)
            if gene is None:
                unmapped += 1
                continue
            try:
                aa_start, aa_stop = int(f[1]), int(f[2])
                start, stop = int(f[3]), int(f[4])
            except ValueError:
                continue
            if start > stop:                                      # minus-strand rows are stored high to low
                start, stop = stop, start
            ### Restrict the domain to the coding bases of the protein it was annotated on. Without
            ### this the span still carries that protein's introns.
            reference_blocks = protein_blocks.get(protein)
            if reference_blocks:
                intervals = intersect(reference_blocks, start, stop)
            else:
                no_blocks += 1
                intervals = [(start, stop)]
            if not intervals:
                continue
            by_gene[gene].append(
                Domain(protein, f[6], f[5], f[7], aa_start, aa_stop, start, stop, tuple(intervals))
            )
    return by_gene, unmapped, no_blocks


def _exon_bounds(exon):
    """Accept either ``(start, stop)`` or the ``(chr, start, stop, strand)`` tuple parse_gff emits."""
    if len(exon) == 2:
        start, stop = exon
    elif len(exon) >= 4:
        start, stop = exon[1], exon[2]
    else:
        raise ValueError(f"Cannot read exon bounds from {exon!r}")
    return (start, stop) if start <= stop else (stop, start)


def coding_blocks(exons, cds_min, cds_max):
    """Clip an isoform's exons to its coding range.

    ``exons`` is a list of ``(start, stop)`` or ``(chr, start, stop, strand)``, either strand.
    Returns the coding intervals only, so intronic gaps stay gaps.
    """
    if cds_min is None or cds_max is None or cds_max < cds_min:
        return []
    blocks = []
    for exon in exons:
        lo, hi = _exon_bounds(exon)
        lo, hi = max(lo, cds_min), min(hi, cds_max)
        if hi >= lo:
            blocks.append((lo, hi))
    return sorted(blocks)


def covered_fraction(blocks, start, stop):
    """Fraction of ``[start, stop]`` covered by ``blocks``."""
    span = stop - start + 1
    if span <= 0:
        return 0.0
    covered = 0
    for lo, hi in blocks:
        if hi < start or lo > stop:
            continue
        covered += min(hi, stop) - max(lo, start) + 1
    return covered / span


def domain_coverage(blocks, domain):
    """Fraction of the domain's CODING bases that this isoform still translates."""
    total = sum(hi - lo + 1 for lo, hi in domain.intervals)
    if total <= 0:
        return 0.0
    covered = 0
    for lo, hi in domain.intervals:
        covered += covered_fraction(blocks, lo, hi) * (hi - lo + 1)
    return covered / total


def domain_identity(domain):
    """The key two isoforms are compared on.

    One protein often carries the same domain several times, once per source database, and
    sometimes at several positions. Comparing on identity rather than on the raw row stops the
    same InterPro domain being reported repeatedly. Counting rows and counting identities gives
    different totals, so callers that report a number must say which one they mean.
    """
    return (domain.interpro, domain.name, domain.description)


def domains_present(blocks, domains, min_coverage=DEFAULT_MIN_COVERAGE):
    """Return the domain identities this isoform retains."""
    present = {}
    for domain in domains:
        if domain_coverage(blocks, domain) >= min_coverage:
            present[domain_identity(domain)] = domain
    return present


def records_covered(blocks, domains, min_coverage=DEFAULT_MIN_COVERAGE):
    """Count individual feature ROWS this isoform covers, not collapsed identities."""
    return sum(1 for d in domains if domain_coverage(blocks, d) >= min_coverage)


def compare_isoforms(blocks_a, blocks_b, domains, min_coverage=DEFAULT_MIN_COVERAGE):
    """Compare two isoform forms of one gene.

    Returns ``(gained, lost, shared)``. ``gained`` holds domains present in form A and absent from
    form B; ``lost`` holds the reverse. The caller decides which form is the included one, so this
    function never guesses direction.
    """
    present_a = domains_present(blocks_a, domains, min_coverage)
    present_b = domains_present(blocks_b, domains, min_coverage)
    gained = [present_a[k] for k in present_a if k not in present_b]
    lost = [present_b[k] for k in present_b if k not in present_a]
    shared = [present_a[k] for k in present_a if k in present_b]
    return gained, lost, shared


def render_feature(domain):
    """Name a feature the way AltAnalyze2 writes it.

    An InterPro domain reads ``Description-IPR000504``. A UniProt feature carries no InterPro
    accession, so it reads as its description alone, for example ``MOD_RES-Phosphoserine``.
    """
    if domain.interpro and domain.interpro.startswith("IPR"):
        return f"{domain.description}-{domain.interpro}"
    ### Structural UniProt features carry no detail, so the source stores 'HELIX-', 'STRAND-',
    ### 'TURN-' and 'DISULFID-' with a trailing separator. 250,539 rows look like this, and
    ### AltAnalyze2 prints them without it.
    return domain.description.rstrip("-")


def terminus_change(sequence_a, sequence_b, flank=5):
    """Classify how two protein sequences differ, following IdentifyAltIsoforms.

    Returns one of ``alt-N-terminus``, ``alt-C-terminus``, ``alt-coding`` or None. AltAnalyze2
    compares the last 5 residues to call an altered C-terminus, so the same flank applies here.
    """
    if not sequence_a or not sequence_b or sequence_a == sequence_b:
        return None
    if sequence_a[:flank] != sequence_b[:flank]:
        return "alt-N-terminus"
    if sequence_a[-flank:] != sequence_b[-flank:]:
        return "alt-C-terminus"
    return "alt-coding"


def _signed(call, direction):
    """Apply the event direction to a gained or lost call.

    ``call`` is '+' when the attribute belongs to the included form. ``direction`` is '+' when the
    examined junction is the inclusion junction. A '-' direction swaps the sign, which is how
    AltAnalyze.py assigns it for exons on the down list.
    """
    if call == "~":
        return "~"
    if direction == "-":
        return "-" if call == "+" else "+"
    return call


def format_protein_predictions(event, direction="+", sequence_a=None, sequence_b=None):
    """Build the AltAnalyze2 ``ProteinPredictions`` string for one event.

    ``event`` comes from :func:`describe_event`. Form A is the included form and form B the
    excluded one. The output holds a protein section and a feature section joined by ``|``, and
    the separator only appears when both sections carry content.
    """
    protein_parts, feature_parts = [], []

    change = terminus_change(sequence_a, sequence_b)
    if change:
        protein_parts.append(f"({_signed('+', direction)}){change}")

    length_a, length_b = event["protein_length_a"], event["protein_length_b"]
    if length_a is not None and length_b is not None:
        ### AltAnalyze2 marks an unchanged length with '~' rather than a direction.
        sign = "~" if length_a == length_b else _signed("+", direction)
        protein_id_a = event.get("protein_id_a") or event["transcript_a"] + "-PEP"
        protein_id_b = event.get("protein_id_b") or event["transcript_b"] + "-PEP"
        protein_parts.append(f"({sign})AA:{length_b}({protein_id_b})->{length_a}({protein_id_a})")

    seen = set()
    for call, features in (("+", event["gained"]), ("-", event["lost"])):
        sign = _signed(call, direction)
        for domain in features:
            label = f"({sign}){render_feature(domain)}"
            if label in seen:
                continue
            seen.add(label)
            feature_parts.append(label)

    protein_section = ", ".join(protein_parts)
    feature_section = ", ".join(sorted(feature_parts))
    if protein_section and feature_section:
        return f"{protein_section}|{feature_section}"
    return protein_section or feature_section


IsoformForm = namedtuple("IsoformForm", "transcript gene exons cds_min cds_max protein_length nmd")


def load_isoform_forms(coding_regions_file, transcripts):
    """Join ``coding_regions.txt`` to the parsed GFF, giving one record per translated isoform.

    ``transcripts`` is what ``isoform_translation.parse_gff`` returns.
    """
    import csv
    forms = {}
    with open(coding_regions_file) as handler:
        for row in csv.DictReader(handler, delimiter="\t"):
            transcript = row["Transcript ID"]
            record = transcripts.get(transcript)
            if record is None:
                continue
            try:
                cds_min, cds_max = int(row["CDS Genomic Min"]), int(row["CDS Genomic Max"])
                protein_length = int(row["Protein Length"])
            except (ValueError, KeyError):
                continue
            forms[transcript] = IsoformForm(
                transcript, row.get("Gene ID", ""), record["exons"],
                cds_min, cds_max, protein_length, row.get("NMD Status", ""),
            )
    return forms


def describe_event(form_a, form_b, reference, min_coverage=DEFAULT_MIN_COVERAGE):
    """Compare two isoform forms of one gene and return everything an event report needs.

    ``form_a`` is the form the event includes, ``form_b`` the form it excludes. Direction is the
    caller's to decide, so nothing here infers it.
    """
    gene = (form_a.gene or form_b.gene).split(".")[0]
    domains = reference.for_gene(gene)
    ### Report the Ensembl protein when the transcript has one. AltAnalyze2 falls back to
    ### '<transcript>-PEP' otherwise, and the BEAT-AML file carries both forms.
    protein_id_a = reference.transcript_to_protein.get(form_a.transcript.split(".")[0])
    protein_id_b = reference.transcript_to_protein.get(form_b.transcript.split(".")[0])
    blocks_a = coding_blocks(form_a.exons, form_a.cds_min, form_a.cds_max)
    blocks_b = coding_blocks(form_b.exons, form_b.cds_min, form_b.cds_max)
    gained, lost, shared = compare_isoforms(blocks_a, blocks_b, domains, min_coverage)
    return {
        "gene": gene,
        "transcript_a": form_a.transcript,
        "transcript_b": form_b.transcript,
        "protein_id_a": protein_id_a,
        "protein_id_b": protein_id_b,
        "protein_length_a": form_a.protein_length,
        "protein_length_b": form_b.protein_length,
        "nmd_a": form_a.nmd,
        "nmd_b": form_b.nmd,
        "gained": gained,
        "lost": lost,
        "shared": shared,
        "domain_records": len(domains),
    }


class DomainReference:
    """Ensembl domain features, ready to query per gene."""

    def __init__(self, ensembl_dir, build="100_38", uniprot_dir=None):
        self.protein_features_tab = os.path.join(ensembl_dir, f"Hs_ProteinFeatures_build_{build}.tab")
        self.ensembl_protein_tab = os.path.join(ensembl_dir, f"Hs_Ensembl_Protein__{build}.tab")
        self.protein_coordinates_tab = os.path.join(ensembl_dir, f"Hs_ProteinCoordinates_build_{build}.tab")
        for path in (self.protein_features_tab, self.ensembl_protein_tab, self.protein_coordinates_tab):
            if not os.path.exists(path):
                raise FileNotFoundError(f"Domain reference file not found: {path}")
        self.protein_to_gene, self.protein_to_transcript = load_protein_to_gene(self.ensembl_protein_tab)
        self.transcript_to_protein = {t: p for p, t in self.protein_to_transcript.items()}
        self.protein_blocks = load_protein_coding_blocks(self.protein_coordinates_tab)
        self.by_gene, self.unmapped_rows, self.rows_without_blocks = load_domains(
            self.protein_features_tab, self.protein_to_gene, self.protein_blocks
        )
        ### AltAnalyze2 reports UniProt features beside InterPro domains, so the BEAT-AML file
        ### carries entries such as MOD_RES-Phosphoserine and TOPO_DOM-Extracellular. The file uses
        ### the same column layout, so the same loader reads it.
        ###
        ### Pass uniprot_dir only when you want those features. The EnsMart100 build of
        ### Hs_FeatureCoordinate.txt drops the last character of every description, so it stores
        ### 'MOD_RES-Phosphoserin' and 'REGION-Pyridoxal phosphate bindin'. The bytes on disk end
        ### that way, and the InterPro descriptions in Hs_ProteinFeatures are unaffected. Enabling
        ### UniProt features therefore makes those names one character short of the BEAT-AML file.
        self.uniprot_rows = 0
        if uniprot_dir:
            feature_file = os.path.join(uniprot_dir, "Hs_FeatureCoordinate.txt")
            if not os.path.exists(feature_file):
                raise FileNotFoundError(f"UniProt feature file not found: {feature_file}")
            extra, extra_unmapped, extra_no_blocks = load_domains(
                feature_file, self.protein_to_gene, self.protein_blocks
            )
            for gene, features in extra.items():
                self.by_gene[gene].extend(features)
                self.uniprot_rows += len(features)
            self.unmapped_rows += extra_unmapped
            self.rows_without_blocks += extra_no_blocks

    def for_gene(self, gene):
        return self.by_gene.get(gene.split(".")[0], [])

    def summary(self):
        total = sum(len(v) for v in self.by_gene.values())
        return (f"{total:,} domain records over {len(self.by_gene):,} genes; "
                f"{self.unmapped_rows:,} rows had no gene mapping; "
                f"{self.rows_without_blocks:,} rows had no coding blocks and kept their raw span")
