"""Hard-coded rules for cross-chromosome joints (EXPERIMENTAL extra-chr-trans-splicing).

Every rule uses only this study's reads, the reference genome and gene annotation already used by the
pipeline (Ensembl 91 exon model, GENCODE v45), and, for fusion-database membership only, curated
cancer fusion collections. No rule uses external expression or other omics data.

Rule groups
  blacklist      membership that removes a molecule before the GFF/h5ad (alignment or compartment
                 artifacts): mitochondrial genome, unplaced contigs, rRNA and mitochondrial-RNA genes,
                 immunoglobulin/T-cell-receptor orphon genes, duplicated sequence near the partner.
  template switch (TS)
                 homology at the joint (genome and mRNA context) above the 99th percentile of the
                 matching null built from the reference annotation; splice-site precision (registered
                 donor and acceptor with canonical dinucleotides); breakpoint position specificity
                 (Poisson test against uniform breakpoints along both mRNAs).
  fusion vs trans-splicing
                 fusion-database membership (curated cancer collections), donor restriction,
                 reciprocal joints and intron-to-intron (genomic breakpoint) molecules in the same donor,
                 and partner-exon imbalance and flanking-gene expression in joint-positive cells.
"""
import math
import re

MAIN_CHROMS = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
CANONICAL_PAIRS = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
CURATED_FUSION_TAGS = {'Mitelman', 'ChimerKB', 'ChimerPub', 'Cosmic', 'chimerdb_omim', 'chimerdb_pubmed', 'HaasMedCancer'}
RNASEQ_DERIVED_TAGS = {'TCGA_StarF2019', 'CCLE_StarF2019', 'ChimerSeq', 'DEEPEST2019', 'TumorFusionsNAR2018',
                       'Klijn_CellLines', 'GUO2018CR_TCGA', 'YOSHIHARA_TCGA', 'Larsson_TCGA',
                       'GTEx_recurrent_StarF2019', 'Babiceanu_Normal', 'Greger_Normal', 'BodyMap'}
PARALOG_TAGS = {'DGD_PARALOGS', 'HGNC_GENEFAM'}
BLACKLIST_GENE_TYPES = {'rRNA', 'rRNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA'}
ORPHON_RE = re.compile(r'^(IGH|IGK|IGL|TRA|TRB|TRG|TRD)V.*OR\d')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- blacklist membership ----
def blacklist_gene_reason(gene_type, gene_name, chrom):
    """Reason a GENCODE gene is a blacklist member, or None."""
    if chrom == 'chrM':
        return 'mitochondrial_genome'
    if gene_type in BLACKLIST_GENE_TYPES:
        return 'rRNA_or_mitochondrial_RNA_gene'
    if gene_name and ORPHON_RE.match(gene_name):
        return 'IG_TR_orphon'
    return None


class IntervalIndex:
    """Sorted intervals per chromosome; overlap query by linear scan of a bisected window."""
    def __init__(self, rows):
        import bisect
        self._bisect = bisect
        self.by_chrom = {}
        for chrom, start, end, payload in rows:
            self.by_chrom.setdefault(chrom, []).append((start, end, payload))
        self.starts = {}
        self.maxlen = {}
        for chrom, items in self.by_chrom.items():
            items.sort()
            self.starts[chrom] = [s for s, _, _ in items]
            self.maxlen[chrom] = max(e - s for s, e, _ in items) + 1

    def overlaps(self, chrom, start, end):
        items = self.by_chrom.get(chrom)
        if not items:
            return []
        lo = self._bisect.bisect_left(self.starts[chrom], start - self.maxlen[chrom])
        hi = self._bisect.bisect_right(self.starts[chrom], end)
        return [p for s, e, p in items[lo:hi] if s <= end and e >= start]


def segment_blacklist(segments, gene_index):
    """First blacklist reason for a molecule (list of segments with chrom and exons), or None."""
    for s in segments:
        if s['chrom'] == 'chrM':
            return 'blacklist_mitochondrial', 'chrM'
    for s in segments:
        if s['chrom'] not in MAIN_CHROMS:
            return 'blacklist_unplaced_contig', s['chrom']
    for s in segments:
        for a, b in s['exons']:
            hits = gene_index.overlaps(s['chrom'], a, b)
            for reason, gene in hits:
                return 'blacklist_' + reason, gene
    return None


# ---------------------------------------------------------------- junction homology ------
def common_suffix(a, b):
    n = 0
    while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def common_prefix(a, b):
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    return n


def genomic_homology(a_exon_end, a_intron_start, b_intron_end, b_exon_start):
    """Breakpoint slide over genome sequence shared by both partners (all transcript-oriented):
    suffix(A exon end, B intron end) + prefix(A intron start, B exon start)."""
    return common_suffix(a_exon_end, b_intron_end) + common_prefix(a_intron_start, b_exon_start)


def mrna_homology(a_mrna, a_index, b_mrna, b_index, span=20):
    """Breakpoint slide over mRNA sequence shared by both partners, for template switching between
    mature mRNAs. a_index = position of the last A base (donor) in A's mRNA; b_index = position of
    the first B base (acceptor) in B's mRNA. Returns (homology, complete) where complete is False
    when A has no downstream mRNA or B no upstream mRNA (that side then contributes 0)."""
    left = common_suffix(a_mrna[max(0, a_index - span + 1):a_index + 1], b_mrna[max(0, b_index - span):b_index])
    right = common_prefix(a_mrna[a_index + 1:a_index + 1 + span], b_mrna[b_index:b_index + span])
    complete = a_index + 1 < len(a_mrna) and b_index > 0
    return left + right, complete


def percentile_from_counts(counts, q):
    """counts: {value: n}. Smallest value v with cumulative fraction >= q."""
    total = sum(counts.values())
    run = 0
    for v in sorted(counts):
        run += counts[v]
        if run / total >= q:
            return v
    return max(counts)


# ---------------------------------------------------------------- position specificity ---
def poisson_sf(k, mu):
    """P(X >= k) for X ~ Poisson(mu); exact sum for small k."""
    if k <= 0:
        return 1.0
    if mu <= 0:
        return 0.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, k):
        term *= mu / i
        cdf += term
    return max(0.0, 1.0 - cdf) if cdf < 1 - 1e-12 else _poisson_upper_tail(k, mu)


def _poisson_upper_tail(k, mu, terms=200):
    term = math.exp(-mu + k * math.log(mu) - math.lgamma(k + 1))
    total = 0.0
    for i in range(k, k + terms):
        total += term
        term *= mu / (i + 1)
    return total


def position_specificity(observed, pair_molecules, a_positions, b_positions, slack=1.0):
    """Expected molecules at one exact donor/acceptor pair if the pair's molecules broke at uniform
    positions along both mRNAs: pair_molecules / (a_positions * b_positions) * slack. Returns
    (expected, Poisson upper-tail p)."""
    expected = pair_molecules / (max(1, a_positions) * max(1, b_positions)) * slack
    return expected, poisson_sf(observed, expected)


def fisher_less(a, b, c, d):
    """One-sided Fisher exact p that a/(a+b) is lower than c/(c+d) (hypergeometric lower tail)."""
    n1, n2, k = a + b, c + d, a + c
    lo = max(0, k - n2)
    lgc = lambda n, r: math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1)
    denom = lgc(n1 + n2, k)
    return min(1.0, sum(math.exp(lgc(n1, x) + lgc(n2, k - x) - denom) for x in range(lo, a + 1)))


def bh_qvalues(pvalues):
    """Benjamini-Hochberg q-values in input order; None entries stay None and are not counted."""
    idx = [i for i, v in enumerate(pvalues) if v is not None]
    m = len(idx)
    q = [None] * len(pvalues)
    running = 1.0
    for rank, i in zip(range(m, 0, -1), sorted(idx, key=lambda k: pvalues[k], reverse=True)):
        running = min(running, pvalues[i] * m / rank)
        q[i] = running
    return q


def stratified_permutation_lower(num, den, positive, strata, n_perm=10000, seed=0):
    """One-sided cell-level test that positive cells have a lower statistic than negative cells.
    Statistic: sum(num)/sum(den) over positive cells minus the same over negative cells (den None:
    mean of num). Labels are permuted within strata, keeping the positive count of each stratum.
    Returns (observed statistic, p) or (None, None) when a group has no denominator."""
    import numpy as np
    num = np.asarray(num, dtype=float)
    den = np.ones_like(num) if den is None else np.asarray(den, dtype=float)
    pos = np.asarray(positive, dtype=bool)
    strata = np.asarray(strata)

    def stat(lbl):
        a, b = den[lbl].sum(), den[~lbl].sum()
        if a <= 0 or b <= 0:
            return None
        return num[lbl].sum() / a - num[~lbl].sum() / b
    observed = stat(pos)
    if observed is None:
        return None, None
    rng = np.random.default_rng(seed)
    groups = [(g, int(pos[g].sum())) for g in (np.where(strata == s)[0] for s in np.unique(strata))]
    groups = [(g, k) for g, k in groups if k]
    lower, valid = 0, 0
    lbl = np.zeros(len(pos), dtype=bool)
    for _ in range(n_perm):
        lbl[:] = False
        for g, k in groups:
            lbl[rng.choice(g, k, replace=False)] = True
        value = stat(lbl)
        if value is None:
            continue
        valid += 1
        lower += value <= observed
    return observed, (1 + lower) / (1 + valid)


# ---------------------------------------------------------------- splice signals ---------
def canonical_splice(donor_dinucleotide, acceptor_dinucleotide):
    return (donor_dinucleotide, acceptor_dinucleotide) in CANONICAL_PAIRS


def concatemer_signature(context):
    """polyA/polyT run of >= 8 in the inserted bases or the 20 nt each side of the joint."""
    seq = context.replace('|', '').upper()
    return 'AAAAAAAA' in seq or 'TTTTTTTT' in seq


# ---------------------------------------------------------------- fusion database --------
def parse_ctat_line(line):
    """CTAT_HumanFusionLib line -> (geneA, geneB, set(tags)) for fusion pairs, else None."""
    import json
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 2 or '--' not in parts[0]:
        return None
    a, b = parts[0].lstrip('+').split('--', 1)
    tags = set(json.loads(parts[1])) if parts[1].startswith('[') else set()
    return a, b, tags


# ---------------------------------------------------------------- decision table ---------
def classify(j):
    """Final class for one joint (dict of evidence). Order matters; the first match wins.
    1 blacklist membership of every molecule; 2 partner sequence similarity (paralog);
    3 template-switch signature (homology above the null 99th percentile, or a polyA/T run at the
    joint in most molecules); 4 not splice-precise (an end is not a registered exon boundary, or the
    splice dinucleotides are not canonical); 5 gene fusion: curated cancer fusion database, or
    donor-restricted (>= 90% of molecules and >= 2 cells in one donor) with a splice-precise reciprocal
    B->A molecule, a recurrent intron-to-intron breakpoint (>= 2 molecules at one position), partner-exon
    imbalance or a flanking-gene expression decrease in joint-positive cells; 6 trans-splicing: template
    switching ruled out (precise, low homology, position-specific recurrence) and >= 2 donors;
    7 unresolved."""
    if j['blacklisted_all']:
        return 'artifact_blacklist'
    if j['paralog_partner']:
        return 'artifact_paralog_sequence'
    if j['ts_status'] == 'template_switch_likely':
        return 'artifact_template_switch'
    if j['ts_status'] == 'not_splice_precise':
        return 'artifact_not_splice_precise'
    if j['fusion_db_curated'] or (j['donor_restricted'] and (j['reciprocal_precise_same_donor'] or j['intronic_breakpoint_recurrent_same_donor']
                                                             or j['partner_exon_imbalance'] or j['flanking_expression_decrease'])):
        return 'gene_fusion_candidate'
    if j['ts_status'] == 'template_switch_ruled_out' and j['donors'] >= 2:
        return 'trans_splicing_candidate'
    return 'unresolved'


def template_switch_status(precise, homology_flag, concatemer_flag, position_specific, adjacent):
    """template_switch_likely: homology above threshold or concatemer signature.
    not_splice_precise: an end is not a registered boundary or dinucleotides are not canonical.
    template_switch_ruled_out: precise, no homology flag, no concatemer, position-specific recurrence and
    at least one molecule whose alignments meet on the read. Otherwise template_switch_not_excluded."""
    if homology_flag or concatemer_flag:
        return 'template_switch_likely'
    if not precise:
        return 'not_splice_precise'
    if position_specific and adjacent:
        return 'template_switch_ruled_out'
    return 'template_switch_not_excluded'
