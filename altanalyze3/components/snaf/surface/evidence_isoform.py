"""Sample-supported isoform hypotheses, with no long-read input at inference.

Coordinates are 1-based inclusive exon boundaries, in ascending genomic order on
both strands. Co-detection supports a hypothesis; it does not establish molecular
phasing. Annotation supplies transcript ends; selected first exons without an
annotated 5' boundary receive an explicitly assumed 250-nt length.

Unlike the legacy single/pair editor, this generator maintains a beam of compatible
junction combinations with a nonempty *joint* sample intersection. It offers local
exon extensions and paired novel exons explicitly and rejects extensions crossing
annotated exons absent from the backbone. NMD is an annotation available to a ranker.
"""
from dataclasses import dataclass
import math

import numpy as np

from . import predict_isoform as P


def load_exon_annotation(path, genes=None):
    """Read AltAnalyze exon blocks, including historical mRNA exon annotations.

    Intronic block coordinates are converted to flanking exon boundary coordinates.
    These local host introns differ from a transcript's exon-skipping splice gaps.
    """
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 6 or (genes is not None and f[0] not in genes):
                continue
            if not f[1].startswith(('E', 'I')):
                continue
            starts, ends = f[4].split('|'), f[5].split('|')
            if len(starts) != len(ends):
                raise ValueError('Ambiguous exon annotation: ' + f[0] + ':' + f[1])
            d = out.setdefault(f[0], {'exons': set(), 'introns': set()})
            for a, b in zip(starts, ends):
                lo, hi = sorted((int(a), int(b)))
                if f[1].startswith('E'):
                    d['exons'].add((lo, hi))
                else:
                    d['introns'].add((lo - 1, hi + 1))
    return out


def introns(chain):
    return tuple((a[1], b[0]) for a, b in zip(chain, chain[1:]) if b[0] > a[1] + 1)


def compatible(junctions):
    js = sorted(set(junctions))
    return all(a < b for a, b in js) and all(a[1] < b[0] for a, b in zip(js, js[1:]))


def contains(chain, junction):
    a, b = junction
    if b == a + 1:  # exon/intron coverage marker, not a splice gap
        return any(s <= a < b <= e for s, e in chain)
    return junction in introns(chain)


def merge_adjacent(chain):
    out = []
    for s, e in chain:
        if out and out[-1][1] + 1 == s:
            out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return out


def first_exon_seeds(chain, strand, target, gene_models):
    """Select a first exon at the target donor, preserving known 5' boundaries.

    The caller must identify the target as a first-exon junction. Junction counts
    alone do not establish that an exon is terminal. Coordinates are inclusive;
    the 250-nt assumption concerns the entire exon, not the UTR or CDS length.
    """
    lo, hi = target
    if hi <= lo + 1:
        raise ValueError('A first-exon junction must be a splice gap')
    minus = str(strand).strip() in ('-', '-1')
    donor = hi if minus else lo
    boundaries = set()
    for st, exons in gene_models.values():
        if not exons or (str(st).strip() in ('-', '-1')) != minus:
            continue
        s, e = exons[-1] if minus else exons[0]
        if s <= donor <= e:
            boundaries.add(e if minus else s)
    source = 'annotated' if boundaries else 'assumed_250nt'
    if not boundaries:
        boundaries.add(donor + 249 if minus else donor - 249)
    for boundary in sorted(boundaries):
        if boundary < 1:
            continue  # Cannot realize 250 nt at this chromosome edge.
        if minus:
            downstream = [(s, min(e, lo)) for s, e in chain if s <= lo]
            seed = downstream + [(donor, boundary)]
        else:
            downstream = [(max(s, hi), e) for s, e in chain if e >= hi]
            seed = [(boundary, donor)] + downstream
        if downstream:
            yield seed, source


def local_edit(chain, junctions, annotated_exons, max_extension_nt=None, host_introns=()):
    """Reject exon extension across an intervening annotated exon.

    Skipping exons *between the splice sites* is legitimate. Swallowing an exon
    outside that gap while extending the backbone is not supported by that junction.
    Checking only the selected backbone would miss the MPZL1/SLC24A4 failure.
    """
    if not compatible(junctions):
        return None
    edited = P.apply_junctions(chain, junctions)
    if not edited:
        return None
    edited = merge_adjacent(edited)
    if not all(contains(edited, j) for j in junctions):
        return None
    # A newly INSERTED exon can fuse reference exons just as an extension can.
    # Forbid every newly covered complete host intron unless the input explicitly
    # includes its exon/intron retention boundary. Checking extension length alone
    # misses a block inserted between an I3 site and an I6 site (MPZL1).
    for a, b in host_introns:
        if (any(s <= a < b <= e for s, e in edited) and
                not any(s <= a < b <= e for s, e in annotated_exons) and
                (a, a + 1) not in junctions and (b - 1, b) not in junctions):
            return None
    for s, e in edited:
        overlaps = [(a, b) for a, b in chain if a <= e and b >= s]
        if not overlaps:
            continue  # a new exon bounded by two measured splice sites
        left, right = min(a for a, _ in overlaps), max(b for _, b in overlaps)
        added = [(s, left - 1), (right + 1, e)]
        for lo, hi in added:
            if max_extension_nt is not None and hi - lo + 1 > max_extension_nt:
                return None
            if lo <= hi and any(lo <= a <= b <= hi for a, b in annotated_exons):
                return None
    return edited


@dataclass
class Evidence:
    """Counts of junctions of ONE gene in a consistent sample order.

    Missing rows remain distinguishable from measured zero rows. `sample_index`
    restricts construction to a particular sample; None builds cohort hypotheses.
    """
    counts: dict
    min_reads: float = 3
    sample_index: int = None

    def __post_init__(self):
        self.counts = {tuple(k): np.asarray(v, dtype=float) for k, v in self.counts.items()}
        sizes = {v.shape for v in self.counts.values()}
        if len(sizes) > 1 or any(len(s) != 1 for s in sizes):
            raise ValueError('Junction counts must have the same one-dimensional sample axis')
        if any(not np.isfinite(v).all() or (v < 0).any() for v in self.counts.values()):
            raise ValueError('Counts must be finite and nonnegative')
        if self.min_reads <= 0:
            raise ValueError('min_reads must be positive')
        self.n_samples = next(iter(sizes))[0] if sizes else 0
        if self.sample_index is not None and not 0 <= self.sample_index < self.n_samples:
            raise ValueError('sample_index outside count matrix')

    def mask(self, js):
        mask = np.ones(self.n_samples, dtype=bool)
        for j in js:
            if j not in self.counts:
                return np.zeros(self.n_samples, dtype=bool)
            mask &= self.counts[j] >= self.min_reads
        return mask

    def supported(self, js):
        mask = self.mask(js)
        return bool(mask[self.sample_index]) if self.sample_index is not None else bool(mask.any())

    def recurrence(self, js):
        return int(self.mask(js).sum())


def paired_exon(a, b, reference_introns, evidence, max_exon_nt=500, recurrent_samples=3):
    """Two complementary, nonoverlapping arms within the same reference intron.

    Longer exons need recurrent co-detection. The recurrence counts matrix columns,
    not independent patients; callers should collapse technical replicates first.
    """
    a, b = sorted((a, b))
    if (a[1] <= a[0] + 1 or b[1] <= b[0] + 1 or
            not compatible((a, b)) or not evidence.supported((a, b))):
        return False
    length = b[0] - a[1] + 1
    host = any(a[0] <= lo < a[1] <= b[0] < hi <= b[1]
               for lo, hi in reference_introns)
    return host and (length <= max_exon_nt or
                     evidence.recurrence((a, b)) >= recurrent_samples)


def combinations(target, evidence, reference_introns, max_edits=4, beam_width=24,
                 max_neighbors=16, max_exon_nt=500, host_introns=None):
    """Bounded search; every multi-event state is observed jointly in >=1 sample."""
    if not evidence.supported((target,)):
        return [(target,)] if evidence.sample_index is None else []
    others = [j for j in evidence.counts if j != target and j[1] > j[0] + 1
              and compatible((target, j)) and evidence.supported((target, j))]
    hosts = reference_introns if host_introns is None else host_introns
    def utility(j):
        pair = paired_exon(target, j, hosts, evidence, max_exon_nt)
        shared = evidence.recurrence((target, j))
        conditional = shared / max(1, evidence.recurrence((target,)))
        # Prefer a complementary novel-exon arm over distant background junctions.
        return (pair, j not in reference_introns, conditional, shared,
                -min(abs(j[0] - target[1]), abs(target[0] - j[1])))
    others = sorted(others, key=utility, reverse=True)[:max_neighbors]
    states = [(target,)]
    all_states = set(states)
    for _ in range(1, max_edits):
        proposals = set()
        for state in states:
            for j in others:
                new = tuple(sorted(set(state) | {j}))
                if len(new) != len(state) + 1 or not compatible(new):
                    continue
                if evidence.supported(new):
                    proposals.add(new)
        def score(js):
            pairs = sum(paired_exon(target, j, hosts, evidence, max_exon_nt)
                        for j in js if j != target)
            return pairs, evidence.recurrence(js), sum(j not in reference_introns for j in js)
        states = sorted(proposals, key=lambda js: (score(js), js), reverse=True)[:beam_width]
        all_states.update(states)
    return sorted(all_states, key=lambda js: (len(js), js))


FEATURE_NAMES = (
    'log_protein_length', 'relative_protein_length', 'nmd', 'cds_start_fraction',
    'log_exon_count', 'log_max_exon_length', 'log_transcript_length',
    'reference_junction_fraction', 'novel_junction_count', 'observed_junction_fraction',
    'mean_conditional_support', 'min_conditional_support', 'joint_novel_support',
    'log_joint_samples', 'log_mean_reads', 'swallowed_supported_introns',
    'edit_count', 'paired_exon', 'length_delta_fraction', 'target_missing',
)


@dataclass
class Candidate:
    prediction: P.PredictedIsoform
    features: tuple
    joint_samples: tuple
    construction: str


def candidate_features(pred, original, target, evidence, ref_introns, combo, hosts=None):
    js = introns(pred.chain)
    observed = [j for j in js if j in evidence.counts]
    novel = [j for j in js if j not in ref_introns]
    target_mask = evidence.mask((target,))
    n_target = max(1, int(target_mask.sum()))
    support = [evidence.recurrence((target, j)) / n_target for j in js]
    required = tuple(sorted(set([target] + novel + list(combo))))
    joint = evidence.mask(required)
    reads = [float(np.mean(evidence.counts[j][target_mask])) for j in observed
             if target_mask.any()]
    swallowed = sum(evidence.recurrence((target, j)) / n_target
                    for j in ref_introns if j in evidence.counts and
                    any(s <= j[0] < j[1] <= e for s, e in pred.chain))
    paired = any(paired_exon(target, j, ref_introns if hosts is None else hosts, evidence)
                 for j in combo if j != target)
    orig_len = sum(e - s + 1 for s, e in original)
    features = (math.log1p(len(pred.protein)), 0.0, float(pred.nmd),
                pred.cds_start / max(1, pred.transcript_len), math.log1p(len(pred.chain)),
                math.log1p(max(e - s + 1 for s, e in pred.chain)),
                math.log1p(pred.transcript_len),
                sum(j in ref_introns for j in js) / max(1, len(js)), len(novel),
                len(observed) / max(1, len(js)), float(np.mean(support)) if support else 0,
                min(support, default=0), int(joint.sum()) / n_target,
                math.log1p(int(joint.sum())), math.log1p(np.mean(reads)) if reads else 0,
                swallowed, len(combo), float(paired),
                (pred.transcript_len - orig_len) / max(1, orig_len),
                float(target not in evidence.counts))
    return features, tuple(np.flatnonzero(joint).tolist()), paired


def generate_candidates(gene_models, target, fetch, evidence, *, max_edits=4,
                        beam_width=24, max_neighbors=16, orfs_per_chain=3,
                        max_chains=160, junction_label='', labels=None,
                        exon_annotation=None, min_orf_nt=P.MIN_ORF_NT,
                        max_extension_nt=500, first_exon=False):
    """Generate structurally constrained hypotheses, then translate spanning ORFs.

    This API accepts reference models and short-read evidence only. Long-read truth
    belongs in a separate training/evaluation process. Search limits are explicit;
    failure to reconstruct is not evidence that the junction is biologically invalid.
    An unbounded intronic arm may extend a local exon by at most max_extension_nt;
    longer extensions need another measured splice boundary. This is a configurable
    construction prior, not proof that longer biological extensions cannot occur.
    first_exon=True selects the target donor as a first exon. Use annotated first
    exon starts when available, otherwise assume a 250-nt exon in strand direction.
    """
    target = tuple(target)
    if max_edits < 1 or beam_width < 1 or max_chains < 1 or orfs_per_chain < 1:
        raise ValueError('Search limits must be positive')
    labels = labels or {}
    ref_introns = {j for _, ch in gene_models.values() for j in introns(ch)}
    annotated = {ex for _, ch in gene_models.values() for ex in ch}
    annotated.update((exon_annotation or {}).get('exons', ()))
    hosts = (exon_annotation or {}).get('introns', ref_introns)
    combos = combinations(target, evidence, ref_introns, max_edits, beam_width, max_neighbors,
                          host_introns=hosts)
    chains = {}
    for tx, (strand, ch) in sorted(gene_models.items()):
        seeds = list(first_exon_seeds(ch, strand, target, gene_models)) if first_exon else [(ch, 'reference')]
        for seed, boundary_source in seeds:
            for combo in combos:
                ed = local_edit(seed, combo, annotated, max_extension_nt=max_extension_nt,
                                host_introns=hosts)
                if not ed:
                    continue
                if first_exon:
                    first_gap = introns(ed)[-1 if str(strand).strip() in ('-', '-1') else 0]
                    if first_gap != target:
                        continue
                valid_exons = True
                for i, (s, e) in enumerate(ed[1:-1], 1):
                    if not any(a <= e and b >= s for a, b in annotated):
                        left, right = (ed[i - 1][1], s), (e, ed[i + 1][0])
                        if not paired_exon(left, right, hosts, evidence):
                            valid_exons = False
                            break
                if not valid_exons:
                    continue
                # Every additional *novel* gap needs joint evidence, including a novel
                # gap inherited from an earlier edit. Reference gaps remain prior knowledge.
                novel = [j for j in introns(ed) if j not in ref_introns and j != target]
                if novel and not evidence.supported([target] + novel):
                    continue
                key = (strand, tuple(ed))
                chains.setdefault(key, (tx, ch, combo, boundary_source))
    def chain_priority(item):
        (strand, ed), (tx, original, combo, boundary_source) = item
        js = introns(ed)
        tm = max(1, evidence.recurrence((target,)))
        support = sum(evidence.recurrence((target, j)) / tm for j in js) / max(1, len(js))
        paired = any(paired_exon(target, j, hosts, evidence) for j in combo if j != target)
        return (paired, support, -max(e - s + 1 for s, e in ed), tx)
    selected = sorted(chains.items(), key=chain_priority, reverse=True)[:max_chains]
    out = []
    for (strand, ed), (tx, original, combo, boundary_source) in selected:
        if sum(e - s + 1 for s, e in ed) > P.MAX_TRANSCRIPT_NT:
            continue
        fragments = {(s, e): fetch(s, e) for s, e in ed}
        if any(len(fragments[s, e]) != e - s + 1 for s, e in ed):
            continue  # A partial reference fetch cannot define a complete exon.
        if boundary_source == 'assumed_250nt':
            first = ed[-1] if str(strand).strip() in ('-', '-1') else ed[0]
            if len(fragments[first]) != 250:
                continue  # Do not silently shorten an assumption at a contig end.
        # `predict_isoform` RE-EDITS the chain it is given, so it asks for spans that the
        # pre-edit `ed` never contained: on the pediatric AML cohort the first such target,
        # LRFN4 ENSG00000173621 at chr11:66857117-66857118, raised
        # KeyError: (66856647, 66857117) and killed the whole run before any target was
        # written. `fragments` is a cache, so a miss must fall through to the real fetcher
        # rather than decide that the span does not exist. The completeness check above still
        # rejects a partial reference fetch.
        def _frag(s, e, _c=fragments, _f=fetch):
            hit = _c.get((s, e))
            if hit is None:
                hit = _c[(s, e)] = _f(s, e)
            return hit

        preds = P.predict_isoform({tx: (strand, list(ed))}, *target,
                                  _frag,
                                  return_all=True, junction_label=junction_label,
                                  min_orf_nt=min_orf_nt)
        # Several starts ending at the same stop are retained; the learned ranker
        # can choose among them. NMD and non-NMD alternatives both remain available.
        preds.sort(key=lambda p: (-len(p.protein), p.cds_start))
        unique = {}
        for pred in preds:
            unique.setdefault(pred.protein, pred)
        preds = list(unique.values())
        retained = preds[:orfs_per_chain]
        clean = next((p for p in preds if not p.nmd), None)
        if clean is not None and clean not in retained:
            retained.append(clean)
        for pred in retained:
            pred.first_exon_boundary_source = boundary_source
            pred.co_applied = len(combo) - 1
            pred.junction_labels = [junction_label or '%d-%d' % target] + [
                labels.get(j, '%d-%d' % j) for j in combo if j != target]
            pred.edit = P.describe_edit(original, pred.chain, *target)
            features, samples, paired = candidate_features(
                pred, original, target, evidence, ref_introns, combo, hosts)
            out.append(Candidate(pred, features, samples,
                                 'first_exon' if first_exon else
                                 'paired_novel_exon' if paired else 'local_edit'))
    longest = max((len(c.prediction.protein) for c in out), default=1)
    for c in out:
        f = list(c.features)
        f[1] = len(c.prediction.protein) / longest
        c.features = tuple(f)
    return out


def generate_hypotheses(gene_models, target, fetch, evidence, *, max_edits=4, **kwargs):
    """Preserve hypotheses at every search depth when the larger beam is capped."""
    unique = {}
    for depth in range(1, max_edits + 1):
        for c in generate_candidates(gene_models, target, fetch, evidence,
                                     max_edits=depth, **kwargs):
            p = c.prediction
            key = (p.strand, tuple(p.chain), p.cds_start, p.protein)
            unique.setdefault(key, c)
    candidates = list(unique.values())
    longest = max((len(c.prediction.protein) for c in candidates), default=1)
    for c in candidates:
        f = list(c.features)
        f[1] = len(c.prediction.protein) / longest
        c.features = tuple(f)
    return candidates


def rank_candidates(candidates, model=None):
    """Rank using an optional fitted sklearn-compatible model, or an explicit rule.

    Scores are ranking values, not calibrated probabilities. Preserve alternatives
    when short reads cannot resolve full-length phasing or transcript ends.
    """
    if not candidates:
        return []
    if model is None:
        scores = [c.features[1] + c.features[10] + c.features[12]
                  + 0.5 * c.features[17] - c.features[2] - c.features[15]
                  for c in candidates]
    else:
        scores = model.predict(np.asarray([c.features for c in candidates]))
    return sorted(zip(map(float, scores), candidates), key=lambda x: x[0], reverse=True)


class LinearQueryRanker:
    """Ridge ranking trained on within-junction feature/quality differences.

    The evaluator centers features and labels within each candidate set during fit.
    At inference, subtracting the query mean would add the same constant to all
    scores, so no labels or reference quality are needed for prediction.
    """
    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, x, y, sample_weight=None):
        from sklearn.linear_model import Ridge
        x = np.asarray(x)
        self.scale_ = np.maximum(np.std(x, axis=0), 1e-6)
        self.model_ = Ridge(alpha=self.alpha, fit_intercept=False)
        self.model_.fit(x / self.scale_, y, sample_weight=sample_weight)
        return self

    def predict(self, x):
        return self.model_.predict(np.asarray(x) / self.scale_)


def load_ranker(path=None):
    """Load the bundled portable linear model, a JSON model, or a fitted pickle bundle."""
    import json
    import pickle
    from pathlib import Path
    path = Path(path) if path else Path(__file__).with_name('data') / 'isoform_ranker.json'
    if path.suffix == '.json':
        bundle = json.loads(path.read_text())
        if tuple(bundle['feature_names']) != FEATURE_NAMES:
            raise ValueError('Ranker feature schema does not match this predictor')
        return PortableLinearRanker(bundle)
    with path.open('rb') as fh:
        bundle = pickle.load(fh)
    if tuple(bundle['feature_names']) != FEATURE_NAMES:
        raise ValueError('Ranker feature schema does not match this predictor')
    return bundle['model']


class PortableLinearRanker:
    """Frozen coefficients with no scikit-learn pickle/version dependency at inference."""
    def __init__(self, bundle):
        self.scale = np.asarray(bundle['scale'], dtype=float)
        self.coefficients = np.asarray(bundle['coefficients'], dtype=float)
        self.intercept = float(bundle['intercept'])
        expected = (len(FEATURE_NAMES),)
        if (self.scale.shape != expected or self.coefficients.shape != expected or
                not np.isfinite(self.scale).all() or (self.scale <= 0).any() or
                not np.isfinite(self.coefficients).all() or not math.isfinite(self.intercept)):
            raise ValueError('Invalid linear ranker parameters')

    def predict(self, x):
        return (np.asarray(x, dtype=float) / self.scale) @ self.coefficients + self.intercept
