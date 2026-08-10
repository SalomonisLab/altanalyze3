"""Classify a splicing event and write the EventAnnotation table.

The target format is the one AltAnalyze2 produces, for example
``/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/BEAT-AML/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt``:
eleven annotation columns followed by one column per sample.

An event is a pair of junctions, the examined one and the background one, each written as
``GENE:exon-exon``. An exon identifier carries a block and a region, so ``E3.2`` is block 3
region 2, and a novel splice site appends its coordinate, as in ``E3.2_169888664``.
"""

import re
from collections import namedtuple

Exon = namedtuple("Exon", "gene kind block region novel raw")

EXON_PATTERN = re.compile(r"^(?P<kind>[EI])(?P<block>\d+)\.(?P<region>\d+)(?:_(?P<novel>\d+))?$")


def parse_exon(token, default_gene=None):
    """Read ``ENSG...:E3.2_169888664`` or ``E3.2`` into its parts. Returns None when unreadable."""
    gene = default_gene
    if ":" in token:
        gene, token = token.rsplit(":", 1)
    match = EXON_PATTERN.match(token)
    if match is None:
        return None
    return Exon(gene, match.group("kind"), int(match.group("block")),
                int(match.group("region")), match.group("novel"), token)


def parse_junction(junction):
    """Split ``GENE:E1.1-E2.1`` or ``GENE:E1.1-GENE2:E3.1`` into its two exons."""
    if ":" not in junction:
        return None, None
    gene, remainder = junction.split(":", 1)
    ### A trans-splicing junction names its second gene, so split on the hyphen that precedes it.
    if ":" in remainder:
        left, right = remainder.split("-", 1)
        return parse_exon(left, gene), parse_exon(right, gene)
    if "-" not in remainder:
        return None, None
    left, right = remainder.split("-", 1)
    return parse_exon(left, gene), parse_exon(right, gene)


### Precision of each call, measured on all 155,863 rows of the BEAT-AML reference file by
### re-deriving the event from that file's own Examined-Junction and Background-Major-Junction
### columns and comparing to its EventAnnotation column:
###
###   alt-3'            95.9%      cassette-exon     77.5%
###   trans-splicing    88.4%      intron-retention  58.8%
###                                alt-5'            53.2%
###                                altPromoter       22.3%
###                                alt5'-alt3'       13.3%
###
### Only the first two are close enough to AltAnalyze2 to report as its EventAnnotation. The rest
### come from a full event model in AltAnalyze2's ExpressionBuilder, not from junction geometry,
### and geometry alone reproduces the whole column on 63.8% of rows. ``confident_only`` therefore
### emits the two reliable calls and leaves the rest blank rather than guessing.
CONFIDENT_EVENTS = frozenset({"alt-3'", "trans-splicing"})


def classify_event(examined, background, last_block=None, confident_only=True):
    """Name the splicing event that separates two junctions.

    ``last_block`` is the gene's highest exon block. It was meant to separate a cassette exon from
    an alternative last exon, but measurement showed it makes the whole column worse, 62.5% against
    63.8%, so it stays off unless a caller asks for it.

    ``confident_only`` keeps only the calls that match AltAnalyze2 closely. Set it False to see
    every call this geometry supports, and read the precision table above before trusting them.

    Returns an event name or an empty string. AltAnalyze2 leaves the column blank often too, on
    29,243 of 155,863 rows in the reference file.
    """
    call = _classify(examined, background, last_block)
    if confident_only and call not in CONFIDENT_EVENTS:
        return ""
    return call


def _classify(examined, background, last_block=None):
    a5, a3 = parse_junction(examined)
    b5, b3 = parse_junction(background)
    if a5 is None or a3 is None or b5 is None or b3 is None:
        return ""

    ### A junction spanning two genes is trans-splicing whatever else it looks like.
    if a5.gene != a3.gene or b5.gene != b3.gene:
        return "trans-splicing"

    ### Any intron block on the examined junction makes this a retention event.
    if a5.kind == "I" or a3.kind == "I":
        return "intron-retention"

    same5 = a5.raw == b5.raw
    same3 = a3.raw == b3.raw

    if same5 and not same3:
        if a3.block == b3.block:
            return "alt-3'"
        if last_block is not None and max(a3.block, b3.block) >= last_block:
            return "alt-C-term"
        return "cassette-exon"

    if same3 and not same5:
        if a5.block == b5.block:
            return "alt-5'"
        if min(a5.block, b5.block) <= 1:
            return "altPromoter"
        return "cassette-exon"

    if not same5 and not same3:
        if a5.block == b5.block and a3.block == b3.block:
            return "alt5'-alt3'"
        return ""

    return ""


### The eleven annotation columns of the AltAnalyze2 EventAnnotation file, in order, followed by
### one column per sample.
ANNOTATION_COLUMNS = [
    "Symbol", "Description", "Examined-Junction", "Background-Major-Junction", "AltExons",
    "ProteinPredictions", "dPSI", "ClusterID", "UID", "Coordinates", "EventAnnotation",
]


### The deployed classifier in annotation/splice_event.py names events differently from the
### AltAnalyze2 EventAnnotation column. Measured against 20,503 classifiable BEAT-AML rows, it
### agrees with that column on 85.2%, against 63.8% for the geometry in this module, and it wins on
### every class that matters: cassette-exon 89.0% (against 77.5%), intron-retention 84.0% (58.8%),
### alt-5' 91.9% (53.2%), alt-C-term 84.0% (about 11%), trans-splicing 90.1% (88.4%). Only alt-3'
### is marginally better here, 95.9% against 92.1%.
###
### Two calls are withheld. altPromoter reaches only 42.2%, and alternative_polyA was wrong on all
### 7 rows it produced.
SPLICE_EVENT_VOCABULARY = {
    "Cassette Exon": "cassette-exon",
    "Mutually-Exclusive Exon": "cassette-exon",
    "Intron Retention": "intron-retention",
    "Alt3' SS": "alt-3'",
    "Alt5' SS": "alt-5'",
    "Trans-Splicing": "trans-splicing",
    "Alt C-Terminal Exon": "alt-C-term",
}
LOW_CONFIDENCE_EVENTS = {"AltPromoter", "AltPolyA"}


def load_event_types(annotated_stats_file):
    """Read Event_Type per event from a stats file annotated by ``junction_isoform.annotate``.

    Keys are the examined junction name, so the EventAnnotation writer can look them up whether or
    not its own identifiers carry coordinates.
    """
    import csv
    events = {}
    with open(annotated_stats_file) as handler:
        for row in csv.DictReader(handler, delimiter="\t"):
            feature = row.get("Feature", "")
            event = SPLICE_EVENT_VOCABULARY.get(row.get("Event_Type", "").strip())
            if not feature or not event:
                continue
            examined = feature.split("|", 1)[0].split("=", 1)[0]
            events[examined] = event
    return events


def load_gene_annotations(annotation_file):
    """Read ``Hs_Ensembl-annotations.txt`` into gene to (symbol, description)."""
    annotations = {}
    with open(annotation_file) as handler:
        for line in handler:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0].startswith("ENS"):
                continue
            annotations[fields[0]] = (fields[1], fields[2] if len(fields) > 2 else "")
    return annotations


def load_junction_coordinates(annotated_file):
    """Map each junction identifier to its coordinates, from the annotated aggregate table.

    Rows there read ``GENE:E1.1-E2.1=chr1:17606-17368``, and the PSI file names the part before
    the equals sign.
    """
    coordinates = {}
    with open(annotated_file) as handler:
        handler.readline()
        for line in handler:
            uid = line.split("\t", 1)[0]
            if "=" not in uid:
                continue
            junction, location = uid.split("=", 1)
            coordinates[junction] = location
    return coordinates


def derive_alt_exons(examined, background):
    """Exons the examined junction uses that the background junction does not.

    ``E6.1-E7.1`` against ``E6.1-E8.2`` gives ``E7.1``, because the two share ``E6.1``.
    """
    gene = examined.split(":", 1)[0]
    a5, a3 = parse_junction(examined)
    b5, b3 = parse_junction(background)
    if a5 is None or a3 is None:
        return ""
    background_raw = {exon.raw for exon in (b5, b3) if exon is not None}
    alt = [exon for exon in (a5, a3) if exon.raw not in background_raw]
    return "|".join(f"{exon.gene or gene}:{exon.raw}" for exon in alt)


def write_event_annotation(psi_file, annotated_file, annotation_file, output_file,
                           protein_predictions=None, confident_events_only=True,
                           event_types=None):
    """Write the EventAnnotation table from a PSI file.

    ``protein_predictions`` optionally maps an examined junction to its ProteinPredictions string,
    as built by ``annotation/domain_impact.format_protein_predictions``. Without it that column
    stays empty, and so do dPSI and ClusterID, which come from a group comparison rather than from
    PSI alone.

    Returns a dictionary of counts so the caller can report how much of each column was filled.
    """
    annotations = load_gene_annotations(annotation_file)
    coordinates = load_junction_coordinates(annotated_file)
    protein_predictions = protein_predictions or {}

    filled = {"rows": 0, "symbol": 0, "coordinates": 0, "event": 0, "protein": 0, "alt_exons": 0}
    with open(psi_file) as source, open(output_file, "w") as target:
        samples = source.readline().rstrip("\n").split("\t")
        ### The PSI writer emits no label over its event column, so every field is a sample.
        target.write("\t".join(ANNOTATION_COLUMNS + samples) + "\n")
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if not fields or "|" not in fields[0]:
                continue
            examined_full, background_full = fields[0].split("|", 1)
            ### PSI identifies a junction by name AND coordinates, because the name alone repeats.
            ### The two go to different columns here, matching the AltAnalyze2 layout.
            examined, _, examined_location = examined_full.partition("=")
            background, _, background_location = background_full.partition("=")
            values = fields[1:]
            gene = examined.split(":", 1)[0]
            symbol, description = annotations.get(gene, ("", ""))
            alt_exons = derive_alt_exons(examined, background)
            ### Prefer the deployed classifier, which agrees with AltAnalyze2 far more often than
            ### the geometry in this module. Fall back to the geometry only where it is absent.
            event = (event_types or {}).get(examined) or classify_event(
                examined, background, confident_only=confident_events_only
            )
            ### Prefer the coordinates the PSI identifier carries; fall back to the annotated table
            ### for an older PSI file that stored names only.
            location = "|".join(x for x in (
                examined_location or coordinates.get(examined, ""),
                background_location or coordinates.get(background, ""),
            ) if x)
            prediction = protein_predictions.get(examined, "")

            filled["rows"] += 1
            filled["symbol"] += bool(symbol)
            filled["coordinates"] += bool(location)
            filled["event"] += bool(event)
            filled["protein"] += bool(prediction)
            filled["alt_exons"] += bool(alt_exons)

            row = [
                symbol, description, examined, background, alt_exons, prediction,
                "", "",                                    # dPSI and ClusterID need a comparison
                f"{symbol}:{examined}|{background}" if symbol else f"{examined}|{background}",
                location, event,
            ] + values
            target.write("\t".join(row) + "\n")
    return filled


def load_last_blocks(exon_file):
    """Highest exon block per gene, so an alternative last exon is not called a cassette exon."""
    last = {}
    with open(exon_file) as handler:
        for line in handler:
            fields = line.split("\t", 3)
            if len(fields) < 3:
                continue
            gene, exon = fields[0], fields[1]
            match = EXON_PATTERN.match(exon)
            if match is None:
                continue
            block = int(match.group("block"))
            if block > last.get(gene, 0):
                last[gene] = block
    return last
