"""Predict the full-length isoform that carries a neojunction, without long-read data.

SNAF-B recovers a full-length isoform by searching a long-read catalog. Where no long-read data
exists, this module builds the isoform from reference annotation alone: it edits every reference
isoform of the gene so that it carries the neojunction, translates each edited isoform, and
nominates one.

The rule, in order:

  1. edit EVERY reference isoform of the gene. The exon count of a reference transcript carries
     no information about which isoform is right, so no isoform is skipped.
  2. enumerate every ATG-to-stop open reading frame of every edited isoform.
  3. keep only the reading frames that carry coding sequence on BOTH sides of the junction. A
     reading frame that starts after the junction, or stops before it, cannot encode a
     junction-derived peptide, so it is a failure case and never a candidate.
  4. drop reading frames whose stop codon triggers nonsense-mediated decay, by the 50 nt rule.
     A degraded transcript presents no epitope.
  5. nominate the longest protein that survives. A transcript that nonsense-mediated decay
     destroys is never assigned: when every spanning reading frame is NMD, the junction gets no
     prediction rather than a doomed one.

An edited isoform is named after the isoform it was edited from, with the junctions that edited
it appended, so its parentage stays visible downstream. The output carries the full mRNA, not
only the coding sequence, because a downstream caller may need the untranslated regions.

A long-read isoform often carries more than one novel junction, and one edit cannot reach it.
`co_junctions` lets the caller pass the other neojunctions the same gene carries. Each is applied
alongside the target junction, and the target junction still has to be spanned.

Benchmark on one internal cohort, held against long-read isoforms: requiring the reading frame to
span the junction raised complete-neoepitope recovery from 52.5% to 60.1%, and co-applying a
second junction of the same gene raised it to 67.8%. The residual failures are isoforms that need
three or more novel junctions.

The caller supplies sequence through a `fetch` callable, so this module needs no genome handle,
no database and no network.
"""
import bisect
import re

from .orf_finder import orf2pep

__all__ = ['apply_junction', 'apply_junctions', 'describe_edit', 'enumerate_orfs',
           'enumerate_orfs_with_offsets', 'offset_of',
           'nmd_by_50nt_rule', 'predict_isoform', 'predict_isoforms', 'PredictedIsoform',
           'write_predictions_tsv', 'PREDICTION_COLUMNS']

_START = re.compile(r'ATG')
_STOP = re.compile(r'(TAA|TGA|TAG)')
_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
NMD_WINDOW_NT = 50
MAX_TRANSCRIPT_NT = 200000
# SNAF's own ORF caller refuses anything shorter: orf_finder.prioritize_orf(min_len=30*3).
# Enumerating every reading frame here would otherwise admit a two-codon "protein", and 7.6% of
# one cohort's predictions were shorter than this before the rule was applied.
MIN_ORF_NT = 30 * 3


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def enumerate_orfs_with_offsets(cdna):
    """Every ORF as (offset, sequence).

    enumerate_orfs returns the sequences alone, which forced the caller to search the mRNA for
    each one to recover its position. That search was 59% of prediction time, and it grows with
    transcript length, so a retained intron of a hundred kilobases paid the most. The offsets are
    known here already.
    """
    starts = [m.start() for m in _START.finditer(cdna)]
    if not starts:
        return []
    stops = [m.start() for m in _STOP.finditer(cdna)]
    out = []
    if not stops:
        for s in starts:
            orf = cdna[s:]
            valid = len(orf) % 3 == 0
            for mo in _STOP.finditer(orf):
                if mo.start() % 3 == 0:
                    valid = not valid
                    break
            if valid:
                out.append((s, orf))
        return out
    by_frame = ([], [], [])
    for p in stops:
        by_frame[p % 3].append(p)
    for s in starts:
        frame = by_frame[s % 3]
        i = bisect.bisect_right(frame, s)
        if i >= len(frame):
            continue
        orf = cdna[s:frame[i]]
        valid = True
        for mo in _STOP.finditer(orf):
            if mo.start() % 3 == 0:
                valid = not valid
                break
        if valid:
            out.append((s, orf))
    return out


def enumerate_orfs(cdna):
    """Every ATG-to-first-in-frame-stop reading frame in a cDNA.

    Same output as orf_finder.transcript2orf, which pairs every start with every stop and so
    costs the product of the two counts. A retained intron produces transcripts where that cost
    becomes prohibitive, so the stop codons are scanned once here and each start jumps to the
    first in-frame stop. The two agreed on every one of 401 real transcripts tested, ORF list and
    translated protein alike.
    """
    return [orf for _off, orf in enumerate_orfs_with_offsets(cdna)]


def apply_junction(chain, lo, hi):
    """Edit one exon chain so that it carries the junction lo->hi.

    Exons left of lo keep their shape, the exon holding lo ends at lo, exons between lo and hi
    disappear, the exon holding hi starts at hi, and exons right of hi keep their shape. When lo
    or hi falls inside an intron the neighbouring exon extends to meet it, which is how a novel
    exon, a retained intron and a novel splice site enter the model.

    chain is a list of (start, end) in ascending genomic order. Returns None when the transcript
    cannot reach both sides of the junction.
    """
    left = []
    for s, e in chain:
        if e <= lo:
            left.append((s, e))
        elif s <= lo < e:
            left.append((s, lo))
    if not left:
        return None
    if left[-1][1] < lo:
        left[-1] = (left[-1][0], lo)
    right = []
    for s, e in chain:
        if s >= hi:
            right.append((s, e))
        elif s < hi <= e:
            right.append((hi, e))
    if not right:
        return None
    if right[0][0] > hi:
        right[0] = (hi, right[0][1])
    edited = left + right
    for i in range(1, len(edited)):
        if edited[i][0] <= edited[i - 1][1]:
            return None
    return edited


def apply_junctions(chain, junctions):
    """Apply several junctions to one chain, in genomic order. None when any one fails."""
    edited = chain
    for lo, hi in sorted(junctions):
        edited = apply_junction(edited, lo, hi)
        if edited is None:
            return None
    return edited


def describe_edit(chain, edited, lo, hi):
    """Name what the junction did to this reference transcript, for the report."""
    if edited == chain:
        return 'unchanged'
    if hi - lo == 1:
        return 'intron retention'
    parts = []
    dropped = [e for e in chain if e[0] > lo and e[1] < hi]
    if dropped:
        parts.append('skip(%d exon%s)' % (len(dropped), '' if len(dropped) == 1 else 's'))
    donor = [e for e in chain if e[0] <= lo <= e[1]]
    if not donor:
        parts.append('donor into intron')
    elif donor[0][1] != lo:
        parts.append('novel donor in exon')
    acceptor = [e for e in chain if e[0] <= hi <= e[1]]
    if not acceptor:
        parts.append('acceptor into intron')
    elif acceptor[0][0] != hi:
        parts.append('novel acceptor in exon')
    return ' + '.join(parts) if parts else 'boundary junction'


def transcript_offsets(chain, minus):
    """Genomic coordinate to transcript offset, 5' to 3'."""
    idx = {}
    n = 0
    for s, e in (reversed(chain) if minus else chain):
        for g in (range(e, s - 1, -1) if minus else range(s, e + 1)):
            idx[g] = n
            n += 1
    return idx


def offset_of(chain, minus, coord):
    """Transcript offset of one genomic coordinate, or None when it is not in the transcript.

    Walks the exons rather than building a per-base map, so a retained intron of a hundred
    kilobases costs the same as a short exon.
    """
    n = 0
    for s, e in (reversed(chain) if minus else chain):
        if s <= coord <= e:
            return n + ((e - coord) if minus else (coord - s))
        n += e - s + 1
    return None


def nmd_by_50nt_rule(orf_start, orf_len, chain, minus, window=NMD_WINDOW_NT):
    """True when the stop codon sits more than `window` nt upstream of the last junction.

    The rule the nonsense-mediated decay literature uses. A single-exon transcript has no
    exon-exon junction downstream of any stop, so it never qualifies.
    """
    exons = list(reversed(chain)) if minus else list(chain)
    if len(exons) < 2:
        return False
    last_junction = sum(e - s + 1 for s, e in exons[:-1])
    return (last_junction - (orf_start + orf_len)) > window


class PredictedIsoform(object):
    """One nominated isoform. `backbone` names the reference transcript it was edited from."""

    __slots__ = ('backbone', 'protein', 'orf', 'chain', 'transcript_len', 'edit', 'nmd',
                 'co_applied', 'strand', 'mrna', 'junction_labels', 'cds_start',
                 'first_exon_boundary_source')

    def __init__(self, backbone, protein, orf, chain, transcript_len, edit, nmd, co_applied,
                 strand, mrna='', junction_labels=(), cds_start=-1):
        self.backbone = backbone
        self.protein = protein
        self.orf = orf                      # coding sequence, nucleotides
        self.chain = chain
        self.transcript_len = transcript_len
        self.edit = edit
        self.nmd = nmd
        self.co_applied = co_applied
        self.strand = strand
        self.mrna = mrna                    # the WHOLE spliced transcript, not only the CDS
        self.junction_labels = tuple(junction_labels)
        self.cds_start = cds_start          # 0-based offset of the CDS within the mRNA
        self.first_exon_boundary_source = 'reference'

    @property
    def isoform_id(self):
        """Parent isoform, then the junctions that edited it.

        An unedited parent keeps its own identifier. An edited one becomes, for example,
        PARENT|J1 or PARENT|J1+J2, so the isoform it descends from stays readable.
        """
        if not self.junction_labels:
            return self.backbone
        return '%s|%s' % (self.backbone, '+'.join(self.junction_labels))

    @property
    def cds_end(self):
        return self.cds_start + len(self.orf) if self.cds_start >= 0 else -1

    def __repr__(self):
        return ('PredictedIsoform(backbone=%r, protein_len=%d, exons=%d, edit=%r, nmd=%r, '
                'co_applied=%d)' % (self.backbone, len(self.protein), len(self.chain),
                                    self.edit, self.nmd, self.co_applied))


PREDICTION_COLUMNS = [
    'predicted_isoform_id', 'parent_isoform_id', 'junctions_applied', 'modification_strategy',
    'strand', 'n_exons', 'exon_chain', 'mrna_length', 'cds_start_in_mrna', 'cds_end_in_mrna',
    'cds_length', 'protein_length', 'nmd', 'junction_encodes_residues', 'mrna_sequence',
    'cds_sequence', 'protein_sequence', 'first_exon_boundary_source',
]


def prediction_row(pred, junction_id=''):
    """One TSV row for a nominated isoform. The mRNA is the whole transcript."""
    return {
        'predicted_isoform_id': pred.isoform_id,
        'parent_isoform_id': pred.backbone,
        'junctions_applied': '+'.join(pred.junction_labels) or junction_id,
        'modification_strategy': pred.edit,
        'strand': pred.strand,
        'n_exons': len(pred.chain),
        'exon_chain': ','.join('%d-%d' % x for x in pred.chain),
        'mrna_length': len(pred.mrna),
        'cds_start_in_mrna': pred.cds_start,
        'cds_end_in_mrna': pred.cds_end,
        'cds_length': len(pred.orf),
        'protein_length': len(pred.protein),
        'nmd': 'yes' if pred.nmd else 'no',
        'junction_encodes_residues': 'yes',
        'mrna_sequence': pred.mrna,
        'cds_sequence': pred.orf,
        'protein_sequence': pred.protein,
        'first_exon_boundary_source': getattr(pred, 'first_exon_boundary_source', 'reference'),
    }


def write_predictions_tsv(path, predictions, extra_columns=None):
    """Write nominated isoforms to a TSV.

    :param predictions: iterable of (junction_id, PredictedIsoform). A None prediction is
        written with empty sequence fields, so no junction disappears from the report.
    :param extra_columns: optional {junction_id: {column: value}} merged into each row.
    """
    extra_columns = extra_columns or {}
    extra_names = []
    for values in extra_columns.values():
        for k in values:
            if k not in extra_names:
                extra_names.append(k)
    cols = ['junction_id'] + PREDICTION_COLUMNS + extra_names
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for junction_id, pred in predictions:
            if pred is None:
                row = {c: '' for c in PREDICTION_COLUMNS}
                row['junction_encodes_residues'] = 'no prediction'
            else:
                row = prediction_row(pred, junction_id)
            row['junction_id'] = junction_id
            row.update(extra_columns.get(junction_id, {}))
            fh.write('\t'.join(str(row.get(c, '')) for c in cols) + '\n')
    return path


def _mrna_of(chain, minus, fetch):
    """The spliced transcript. Rebuilt only for a candidate that survived every filter."""
    seq = ''.join(fetch(s, e) for s, e in chain).upper()
    return reverse_complement(seq) if minus else seq


def predict_isoform(gene_models, lo, hi, fetch, co_junctions=(), max_transcript_nt=None,
                    return_all=False, junction_label='', co_junction_labels=None, cache=None,
                    min_orf_nt=MIN_ORF_NT, evidence=None, ranking_model=None,
                    exon_annotation=None, max_edits=4, max_extension_nt=500,
                    first_exon=False):
    """Nominate the isoform that carries the junction lo->hi.

    :param gene_models: {transcript_id: (strand, [(exon_start, exon_end), ...])} for ONE gene.
        strand is '+'/'-' or '1'/'-1'. Exons in ascending genomic order.
    :param lo: lower genomic coordinate of the junction. For a retained intron, pass the two
        coordinates that flank the exon-intron boundary, one base apart.
    :param hi: upper genomic coordinate of the junction.
    :param fetch: callable(start, end) -> forward-strand sequence, 1-based inclusive. Wrap a
        pysam FastaFile, a local dict, or any other source.
    :param co_junctions: other neojunctions of the same gene, as (lo, hi) pairs. Each is applied
        alongside the target junction, which still has to be spanned. A long-read isoform often
        carries several novel junctions, and one edit cannot reach it.
    :param min_orf_nt: shortest reading frame accepted, defaulting to the 30 codons
        orf_finder.prioritize_orf already requires.
    :param cache: optional dict, reused across the junctions of ONE gene. Co-applying a second
        junction often leaves a backbone unchanged, so the same edited exon chain is built many
        times over; the open reading frames of a chain do not depend on which junction is being
        scored, so they are computed once and keyed by (transcript, chain). Pass a fresh dict per
        gene and the results are identical, only faster.
    :param return_all: return every surviving candidate instead of the nomination.
    :param evidence: optional evidence_isoform.Evidence with per-sample junction counts.
        Enables local construction and jointly supported multi-junction hypotheses. When
        supplied, this replaces legacy co_junctions enumeration and longest-ORF ranking.
        NMD is retained as a risk annotation in this mode, rather than a certainty of decay.
    :param ranking_model: optional fitted model for evidence mode. Without one, use the
        explicit sample-support ranking rule. Scores are not calibrated probabilities.
    :param exon_annotation: optional {'exons': set, 'introns': set} from AltAnalyze exon
        blocks, including historical mRNA annotations omitted from Ensembl transcripts.
    :param max_edits: maximum simultaneous junction edits in evidence mode (default 4).
    :param max_extension_nt: maximum exon extension without a second measured splice
        boundary in evidence mode (default 500; None disables this prior).
    :param first_exon: select the target as the first-exon splice junction. Preserve an
        annotated first-exon 5' boundary, otherwise assume a 250-nt first exon.
    :param junction_label: identifier of the target junction, used to name the edited isoform.
    :param co_junction_labels: {(lo, hi): label} for the co-applied junctions.
    :return: a PredictedIsoform, a list of them when return_all, or None when nothing survives.
        Legacy mode excludes predicted NMD; evidence mode retains that annotation.
    """
    if evidence is not None or first_exon:
        from .evidence_isoform import Evidence, generate_hypotheses, rank_candidates
        if evidence is None:
            evidence = Evidence({})
        candidates = generate_hypotheses(
            gene_models, (lo, hi), fetch, evidence, max_edits=max_edits,
            junction_label=junction_label, labels=co_junction_labels,
            exon_annotation=exon_annotation, min_orf_nt=min_orf_nt,
            max_extension_nt=max_extension_nt, first_exon=first_exon)
        if max_transcript_nt is not None:
            candidates = [c for c in candidates if c.prediction.transcript_len <= max_transcript_nt]
        ranked = rank_candidates(candidates, ranking_model)
        if return_all:
            return [c.prediction for _, c in ranked]
        return ranked[0][1].prediction if ranked else None
    labels = dict(co_junction_labels or {})
    limit = max_transcript_nt or MAX_TRANSCRIPT_NT
    target_label = junction_label or '%d-%d' % (lo, hi)
    combos = [[(lo, hi)]] + [[(lo, hi), tuple(j)] for j in co_junctions
                             if tuple(j) != (lo, hi)]
    survivors = []
    for tx, (strand, chain) in gene_models.items():
        minus = str(strand).strip() in ('-', '-1')
        for combo in combos:
            edited = apply_junctions(chain, combo)
            if not edited:
                continue
            key = (tx, tuple(edited)) if cache is not None else None
            if key is not None and key in cache:
                orfs = cache[key]
                if orfs is None:
                    continue
            else:
                mrna = ''.join(fetch(s, e) for s, e in edited).upper()
                if not mrna or len(mrna) > limit:
                    if key is not None:
                        cache[key] = None
                    continue
                if minus:
                    mrna = reverse_complement(mrna)
                orfs = enumerate_orfs_with_offsets(mrna)
                if key is not None:
                    cache[key] = orfs
            off_lo = offset_of(edited, minus, lo)
            off_hi = offset_of(edited, minus, hi)
            if off_lo is None or off_hi is None:
                continue
            left_of, right_of = sorted((off_lo, off_hi))
            for start, orf in orfs:
                if len(orf) < min_orf_nt:
                    continue
                # coding sequence must exist on both sides of the junction
                if not (start <= left_of and start + len(orf) > right_of):
                    continue
                protein = orf2pep(orf)
                if not protein:
                    continue
                combo_labels = [target_label] + [labels.get(tuple(j), '%d-%d' % tuple(j))
                                                 for j in combo[1:]]
                span = sum(e - s + 1 for s, e in edited)
                survivors.append(PredictedIsoform(
                    backbone=tx, protein=protein, orf=orf, chain=edited,
                    transcript_len=span, edit=describe_edit(chain, edited, lo, hi),
                    nmd=nmd_by_50nt_rule(start, len(orf), edited, minus),
                    co_applied=len(combo) - 1, strand='-' if minus else '+',
                    mrna=_mrna_of(edited, minus, fetch), junction_labels=combo_labels,
                    cds_start=start))
    if return_all:
        return survivors
    # An NMD transcript is destroyed before it can present anything, so it is never assigned.
    clean = [c for c in survivors if not c.nmd]
    if not clean:
        return None
    return max(clean, key=lambda c: len(c.protein))


def predict_isoforms(junctions, models, fetch, evidence_by_gene=None,
                     exon_annotation_by_gene=None, **kwargs):
    """Nominate one isoform per junction.

    :param junctions: iterable of (junction_id, gene_id, lo, hi).
    :param models: {gene_id: {transcript_id: (strand, [(start, end), ...])}}.
    :param fetch: callable(start, end) -> forward-strand sequence.
    :return: {junction_id: PredictedIsoform or None}.

    Junctions of the same gene are offered to each other as co_junctions, so an isoform that
    needs two novel junctions can still be built. Supply {gene: Evidence} through
    evidence_by_gene to use sample-supported inference instead. Missing genes receive
    empty evidence and only reference/local single-junction hypotheses.
    """
    by_gene = {}
    for jid, gene, lo, hi in junctions:
        by_gene.setdefault(gene, []).append((jid, lo, hi))
    out = {}
    for gene, items in by_gene.items():
        gene_models = models.get(gene, {})
        for jid, lo, hi in items:
            others = [(l, h) for j2, l, h in items if j2 != jid]
            labels = {(l, h): j2 for j2, l, h in items if j2 != jid}
            options = dict(kwargs)
            if evidence_by_gene is not None:
                from .evidence_isoform import Evidence
                options['evidence'] = evidence_by_gene.get(gene, Evidence({}))
            if exon_annotation_by_gene is not None:
                options['exon_annotation'] = exon_annotation_by_gene.get(gene)
            out[jid] = predict_isoform(gene_models, lo, hi, fetch, co_junctions=others,
                                       junction_label=jid, co_junction_labels=labels,
                                       **options)
    return out
