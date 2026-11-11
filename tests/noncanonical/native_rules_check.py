"""Unit tests for xchr_rules.py (constructed inputs with known answers). Run: python3 test_xchr_rules.py"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from altanalyze3.components.bam import cross_chromosome_rules as R

results = []


def check(name, cond):
    results.append((name, bool(cond)))


# blacklist membership
check('chrM gene is blacklisted', R.blacklist_gene_reason('protein_coding', 'MT-CO1', 'chrM') == 'mitochondrial_genome')
check('rRNA gene is blacklisted', R.blacklist_gene_reason('rRNA', 'RNA5-8SN4', 'chr21') == 'rRNA_or_mitochondrial_RNA_gene')
check('rRNA pseudogene is blacklisted', R.blacklist_gene_reason('rRNA_pseudogene', 'RNA5SP1', 'chr1') == 'rRNA_or_mitochondrial_RNA_gene')
check('IGHV3OR16-9 orphon is blacklisted', R.blacklist_gene_reason('IG_V_pseudogene', 'IGHV3OR16-9', 'chr16') == 'IG_TR_orphon')
check('TRBV20OR9-2 orphon is blacklisted', R.blacklist_gene_reason('TR_V_gene', 'TRBV20OR9-2', 'chr9') == 'IG_TR_orphon')
check('IGHV3-30 (not orphon) is not blacklisted', R.blacklist_gene_reason('IG_V_gene', 'IGHV3-30', 'chr14') is None)
check('ordinary gene is not blacklisted', R.blacklist_gene_reason('protein_coding', 'BCR', 'chr22') is None)
idx = R.IntervalIndex([('chr21', 8_433_000, 8_433_200, ('rRNA_or_mitochondrial_RNA_gene', 'RNA5-8SN4')),
                       ('chr16', 33_000_000, 33_001_000, ('IG_TR_orphon', 'IGHV3OR16-9'))])
seg = lambda c, ex: dict(chrom=c, exons=ex)
check('segment overlapping rRNA gene', R.segment_blacklist([seg('chr1', [(10, 20)]), seg('chr21', [(8_433_100, 8_433_150)])], idx)
      == ('blacklist_rRNA_or_mitochondrial_RNA_gene', 'RNA5-8SN4'))
check('segment on chrM', R.segment_blacklist([seg('chrM', [(100, 200)]), seg('chr1', [(10, 20)])], idx)[0] == 'blacklist_mitochondrial')
check('segment on unplaced contig', R.segment_blacklist([seg('chrUn_GL000220v1', [(1, 50)]), seg('chr1', [(10, 20)])], idx)[0]
      == 'blacklist_unplaced_contig')
check('segment outside all members', R.segment_blacklist([seg('chr1', [(10, 20)]), seg('chr16', [(32_000_000, 32_000_100)])], idx) is None)
check('interval edge overlap counts', idx.overlaps('chr16', 33_001_000, 33_001_005) == [('IG_TR_orphon', 'IGHV3OR16-9')])

# homology
check('genomic homology 0 when unrelated', R.genomic_homology('AAAAACCCCC', 'GTAAGT', 'TTTTTTTTAG', 'CCCCC') == 0)
check('genomic homology left 3', R.genomic_homology('TTTTTCAG', 'GTAAGT', 'CCCCCCAG', 'TACGT') == 3)
check('genomic homology left+right', R.genomic_homology('ACGTCAG', 'GTAAG', 'TTTCAG', 'GTAAC') == 4 + 4)
a_mrna = 'AAAAAAAAAACCCGGGTTTT'      # donor = index 12 (last C of CCC)
b_mrna = 'GGGCCCGGGTTTTACGTACG'      # acceptor = index 6 (first G after GGGCCC)
h, complete = R.mrna_homology(a_mrna, 12, b_mrna, 6)
check('mRNA homology left CCC + right GGGTTTT', h == 3 + 7 and complete)
h, complete = R.mrna_homology('ACGT', 3, 'TTTT', 0)
check('mRNA homology incomplete at transcript ends', complete is False)
check('percentile from counts', R.percentile_from_counts({0: 50, 1: 30, 2: 15, 5: 5}, 0.99) == 5
      and R.percentile_from_counts({0: 50, 1: 30, 2: 15, 5: 5}, 0.8) == 1)

# Poisson tail and position specificity
check('poisson_sf(0) = 1', R.poisson_sf(0, 2.0) == 1.0)
check('poisson_sf(1, mu) = 1 - exp(-mu)', abs(R.poisson_sf(1, 0.5) - (1 - 2.718281828459045 ** -0.5)) < 1e-12)
check('poisson_sf tiny mu stays positive', 0 < R.poisson_sf(3, 1e-6) < 1e-17)
try:
    from scipy.stats import poisson
    ok = all(abs(R.poisson_sf(k, mu) - poisson.sf(k - 1, mu)) <= 1e-9 * max(1e-300, poisson.sf(k - 1, mu)) + 1e-15
             for k in (1, 2, 3, 5, 10, 30) for mu in (1e-8, 1e-3, 0.5, 3.0, 12.0))
    check('poisson_sf equals scipy.stats.poisson.sf', ok)
except ImportError:
    pass
e, p = R.position_specificity(observed=5, pair_molecules=40, a_positions=1500, b_positions=2000)
check('position specificity: 5 molecules at one site out of 40 is specific', e < 1e-4 and p < 1e-15)
e, p = R.position_specificity(observed=1, pair_molecules=40, a_positions=10, b_positions=10, slack=4)
check('position specificity: 1 of 40 over 100 sites is not specific', p > 0.5)

# one-sided Fisher
check('fisher lower tail equal proportions ~ 0.5-1', R.fisher_less(5, 5, 50, 50) > 0.4)
check('fisher lower tail strong depletion', R.fisher_less(0, 20, 50, 50) < 1e-5)
try:
    from scipy.stats import fisher_exact
    ok = all(abs(R.fisher_less(*t) - fisher_exact([[t[0], t[1]], [t[2], t[3]]], alternative='less')[1]) < 1e-9
             for t in [(0, 20, 50, 50), (3, 7, 30, 10), (5, 5, 50, 50), (0, 0, 4, 4), (12, 3, 2, 40)])
    check('fisher_less equals scipy fisher_exact less', ok)
except ImportError:
    pass

# Benjamini-Hochberg and stratified permutation
q = R.bh_qvalues([0.01, 0.04, 0.03, None, 0.5])
check('BH q-values', [round(v, 4) if v is not None else None for v in q] == [0.04, 0.0533, 0.0533, None, 0.5])
try:
    from statsmodels.stats.multitest import multipletests
    import random as _r
    _rng = _r.Random(1); ps = [_rng.random() ** 3 for _ in range(200)]
    check('BH equals statsmodels fdr_bh', max(abs(a - b) for a, b in zip(R.bh_qvalues(ps), multipletests(ps, method='fdr_bh')[1])) < 1e-12)
except ImportError:
    pass
import random as _random
_g = _random.Random(3)
num = [_g.randint(0, 5) for _ in range(300)]; den = [5] * 300
pos = [i < 6 for i in range(300)]; strata = ['s1' if i % 2 else 's2' for i in range(300)]
stat0, p0 = R.stratified_permutation_lower(num, den, pos, strata, n_perm=2000)
num_dep = [0 if i < 6 else n for i, n in enumerate(num)]
stat1, p1 = R.stratified_permutation_lower(num_dep, den, pos, strata, n_perm=2000)
check('permutation: depleted positives give small p', p1 < 0.01 and stat1 < 0)
check('permutation: p in (0, 1]', 0 < p0 <= 1)
null_ps = []
for rep in range(200):
    rg = _random.Random(100 + rep)
    nm = [rg.randint(0, 5) for _ in range(200)]
    null_ps.append(R.stratified_permutation_lower(nm, [5] * 200, [i < 5 for i in range(200)], ['a'] * 200, n_perm=400, seed=rep)[1])
check('permutation: null calibration, P(p<0.05) between 0.01 and 0.10', 0.01 <= sum(x < 0.05 for x in null_ps) / 200 <= 0.10)
check('permutation: no denominator returns None', R.stratified_permutation_lower([0, 0], [0, 0], [True, False], ['a', 'a'])[1] is None)

# splice signals and concatemers
check('GT-AG canonical', R.canonical_splice('GT', 'AG'))
check('GC-AG canonical', R.canonical_splice('GC', 'AG'))
check('CT-AC not canonical', not R.canonical_splice('CT', 'AC'))
check('polyA concatemer', R.concatemer_signature('ACGTAAAAAAAAAC|ttg|GGCA'))
check('no concatemer', not R.concatemer_signature('ACGTACGTAC||GGCATTGACC'))

# CTAT parsing
check('CTAT fusion line', R.parse_ctat_line('+NUP214--ABL1\t["Mitelman"]\t{}') == ('NUP214', 'ABL1', {'Mitelman'}))
check('CTAT gene line ignored', R.parse_ctat_line('1-8D\t.\t{"HGNC":"x"}') is None)

# decision table
base = dict(blacklisted_all=False, paralog_partner=False, ts_status='template_switch_ruled_out', fusion_db_curated=False,
            donor_restricted=False, reciprocal_precise_same_donor=False, intronic_breakpoint_recurrent_same_donor=False,
            partner_exon_imbalance=False, flanking_expression_decrease=False, donors=3)
check('class trans-splicing', R.classify(base) == 'trans_splicing_candidate')
check('class blacklist first', R.classify(dict(base, blacklisted_all=True, fusion_db_curated=True)) == 'artifact_blacklist')
check('class paralog', R.classify(dict(base, paralog_partner=True)) == 'artifact_paralog_sequence')
check('class template switch', R.classify(dict(base, ts_status='template_switch_likely')) == 'artifact_template_switch')
check('TS likely beats fusion database', R.classify(dict(base, ts_status='template_switch_likely', fusion_db_curated=True)) == 'artifact_template_switch')
check('not precise beats fusion evidence', R.classify(dict(base, ts_status='not_splice_precise', donor_restricted=True,
                                                          reciprocal_precise_same_donor=True)) == 'artifact_not_splice_precise')
check('class fusion by database', R.classify(dict(base, fusion_db_curated=True)) == 'gene_fusion_candidate')
check('class fusion by database, single molecule not excluded', R.classify(dict(base, fusion_db_curated=True, donors=1,
                                                                             ts_status='template_switch_not_excluded')) == 'gene_fusion_candidate')
check('fusion by donor-restricted reciprocal', R.classify(dict(base, donor_restricted=True, reciprocal_precise_same_donor=True)) == 'gene_fusion_candidate')
check('fusion by recurrent intronic breakpoint', R.classify(dict(base, donor_restricted=True, intronic_breakpoint_recurrent_same_donor=True)) == 'gene_fusion_candidate')
check('fusion by flanking decrease', R.classify(dict(base, donor_restricted=True, flanking_expression_decrease=True)) == 'gene_fusion_candidate')
check('fusion by exon imbalance', R.classify(dict(base, donor_restricted=True, partner_exon_imbalance=True)) == 'gene_fusion_candidate')
check('reciprocal without donor restriction is not fusion', R.classify(dict(base, reciprocal_precise_same_donor=True)) == 'trans_splicing_candidate')
check('donor-restricted alone is not fusion', R.classify(dict(base, donor_restricted=True, donors=1)) == 'unresolved')
check('one donor without fusion evidence is unresolved', R.classify(dict(base, donors=1)) == 'unresolved')
check('TS not excluded is unresolved', R.classify(dict(base, ts_status='template_switch_not_excluded')) == 'unresolved')
check('not precise class', R.classify(dict(base, ts_status='not_splice_precise')) == 'artifact_not_splice_precise')
check('TS status homology', R.template_switch_status(True, True, False, True, True) == 'template_switch_likely')
check('TS status concatemer', R.template_switch_status(True, False, True, True, True) == 'template_switch_likely')
check('TS status not precise', R.template_switch_status(False, False, False, True, True) == 'not_splice_precise')
check('TS status ruled out', R.template_switch_status(True, False, False, True, True) == 'template_switch_ruled_out')
check('TS status single molecule not excluded', R.template_switch_status(True, False, False, False, True) == 'template_switch_not_excluded')
check('TS status no adjacent molecule not excluded', R.template_switch_status(True, False, False, True, False) == 'template_switch_not_excluded')

failed = [n for n, ok in results if not ok]
for n, ok in results:
    print(('PASS ' if ok else 'FAIL ') + n)
print(f'{len(results) - len(failed)} of {len(results)} unit tests passed')
sys.exit(1 if failed else 0)
