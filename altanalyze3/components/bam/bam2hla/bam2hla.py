#!/usr/bin/env python3
"""
bam2hla - fast class-I HLA genotyping (HLA-A/-B/-C, 2-field) directly from a
genome-aligned BAM via pysam pileup, without read re-alignment or OptiType.

Approach
--------
A precomputed genomic signature database (build_reference.py) gives, for every
IMGT/HLA 2-field allele, its expected forward-strand base at each polymorphic
exon-2/3/4 position on chr6.  At runtime we:
  1. auto-detect the genome build (hg19 vs hg38) from the BAM,
  2. pileup the BAM at the signature positions (spliced-read aware),
  3. solve, per gene, for the diploid allele pair that best explains the
     observed bases (OptiType-style "explain both haplotypes" objective),
  4. resolve homozygous vs heterozygous from minor-allele read support.

Usage
-----
  python bam2hla.py --bam sample.bam [--build auto|hg19|hg38] \
                    [--out sample_hla.tsv] [--genes A,B,C]

Design mirrors the pysam pileup style of components/bam/variant_extraction.py.
"""

import os, sys, json, gzip, argparse, csv, re
from collections import defaultdict

import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))

# GRCh38 / GRCh37 chr6 lengths - definitive build discriminator
CHR6_LEN = {'hg38': 170805979, 'hg19': 171115067}
_BASES = ('A', 'C', 'G', 'T')
_BIT = {'A': 1, 'C': 2, 'G': 4, 'T': 8}          # base -> bitmask
_FREQ_CSV = next((p for p in (
    os.path.join(HERE, 'data', 'HLA_Allele_frequency_21_populations.csv'),
    os.path.join(HERE, '..', '..', 'snaf', 'HLA_Allele_frequency_21_populations.csv'),
) if os.path.exists(p)), os.path.join(HERE, 'data', 'HLA_Allele_frequency_21_populations.csv'))

_SIG_CACHE = {}
_FREQ_CACHE = None


# ----------------------------------------------------------------------------
# signature loading (+ per-gene numpy mask matrix, cached)
# ----------------------------------------------------------------------------
def load_signatures(build):
    if build in _SIG_CACHE:
        return _SIG_CACHE[build]
    path = os.path.join(HERE, 'build', f'signatures_{build}.json.gz')
    if not os.path.exists(path):
        raise FileNotFoundError(f'no signature DB for build {build}: {path}')
    with gzip.open(path, 'rt') as fh:
        sig = json.load(fh)
    freq = load_frequencies()
    for g, gs in sig['genes'].items():
        alleles = list(gs['allele_bases'].keys())
        n = len(gs['positions'])
        M = np.zeros((len(alleles), n), dtype=np.uint8)   # allele bit per position
        for i, tf in enumerate(alleles):
            s = gs['allele_bases'][tf]
            M[i] = [_BIT.get(ch, 0) for ch in s]
        gs['_alleles'] = alleles
        gs['_mask'] = M
        gs['_freq'] = np.array([freq.get(tf, 0.0) for tf in alleles], dtype=np.float64)
        gs['_pos_index'] = {p: k for k, p in enumerate(gs['positions'])}
    _SIG_CACHE[build] = sig
    return sig


def load_frequencies():
    """Mean allele frequency across the 21-population panel, keyed by 2-field
    allele.  Used only as a tie-break prior among equally-fitting solutions."""
    global _FREQ_CACHE
    if _FREQ_CACHE is not None:
        return _FREQ_CACHE
    freq = {}
    if not os.path.exists(_FREQ_CSV):
        _FREQ_CACHE = freq
        return freq
    with open(_FREQ_CSV, newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.reader(fh))
    for row in rows:
        if len(row) < 3 or not row[0] or row[0] == 'Locus':
            continue
        allele = row[1].strip()
        m = re.match(r'^([A-Z0-9]+)\*(\d+):(\d+)', allele)
        if not m:
            continue
        tf = f'{m.group(1)}*{m.group(2)}:{m.group(3)}'
        vals = []
        for x in row[2:]:
            try:
                vals.append(float(x))
            except ValueError:
                pass
        f = sum(vals) / len(vals) if vals else 0.0
        freq[tf] = max(freq.get(tf, 0.0), f)
    _FREQ_CACHE = freq
    return freq


# ----------------------------------------------------------------------------
# build detection
# ----------------------------------------------------------------------------
def _chrom_for(bam, want_len=None, names=('chr6', '6')):
    """Resolve the chr6 contig name actually used in this BAM."""
    refs = dict(zip(bam.references, bam.lengths))
    for n in names:
        if n in refs:
            return n, refs[n]
    # fall back: match by length
    if want_len:
        for n, L in refs.items():
            if L == want_len:
                return n, L
    return None, None


def detect_build(bam_path):
    """Detect hg19/hg38 from the chr6 contig length, then CONFIRM by checking
    that HLA-A/-B/-C exons carry reads at that build's coordinates (5+ highly
    expressed, constitutively spliced chr6 loci) - the two builds place these
    genes ~300 kb apart so only the correct build shows coverage."""
    bam = pysam.AlignmentFile(bam_path, 'rb')
    refs = dict(zip(bam.references, bam.lengths))
    # primary: chr6 length
    call = None
    for n in ('chr6', '6'):
        if n in refs:
            for b, L in CHR6_LEN.items():
                if refs[n] == L:
                    call = b
            break
    # confirmation via coverage at HLA loci for each build
    evidence = {}
    for b in ('hg38', 'hg19'):
        try:
            sig = load_signatures(b)
        except FileNotFoundError:
            continue
        cn, _ = _chrom_for(bam, CHR6_LEN[b])
        if cn is None:
            evidence[b] = 0; continue
        cov = 0
        for g, gs in sig['genes'].items():
            ps = gs['positions']
            mid = ps[len(ps) // 2]
            try:
                cov += bam.count(cn, mid, mid + 1)
            except ValueError:
                pass
        evidence[b] = cov
    bam.close()
    if call is None and evidence:
        call = max(evidence, key=evidence.get)
    if call and evidence.get(call, 0) == 0 and evidence:
        alt = max(evidence, key=evidence.get)
        if evidence[alt] > 0:
            call = alt
    return call, evidence


# ----------------------------------------------------------------------------
# pileup
# ----------------------------------------------------------------------------
def pileup_positions(bam, chrom, positions, min_baseq=13, min_mapq=0, max_depth=4000):
    """Return {genomic_pos_0based: {A,C,G,T counts}} over the given positions.
    Spliced (N) reads are handled natively by pysam; deletions/introns yield no
    base at that column and are ignored."""
    pset = set(positions)
    lo, hi = min(positions), max(positions) + 1
    counts = {p: {b: 0 for b in _BASES} for p in positions}
    for col in bam.pileup(chrom, lo, hi, truncate=True, stepper='nofilter',
                          min_base_quality=min_baseq, ignore_overlaps=False,
                          max_depth=max_depth):
        p = col.reference_pos
        if p not in pset:
            continue
        cp = counts[p]
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if min_mapq and aln.mapping_quality < min_mapq:
                continue
            b = aln.query_sequence[pr.query_position]
            if b in cp:
                cp[b] += 1
    return counts


# ----------------------------------------------------------------------------
# diploid solver
# ----------------------------------------------------------------------------
def type_gene(gene_sig, counts, min_depth=8, top_k=400, balance_thr=0.15,
              min_minor_count=3, min_positions=10):
    """2-field diploid call for one gene via per-position base-set matching.

    (1) Denoise: at each covered position call the *real* base-set = major base
        plus any minor base with frequency >= balance_thr (and >= min_minor_count
        reads).  This removes sequencing-error reads so het/hom is decided by
        real signal, not noise.
    (2) Fit: choose the allele pair whose predicted base-set per position
        ({e_X, e_Y}) matches the observed real base-sets at the most positions.
        Vectorised over all alleles with numpy bitmasks.
    (3) Tie-break equally-fitting solutions by the population-frequency prior,
        then by allele completeness, then prefer homozygous.
    """
    positions = gene_sig['positions']
    M = gene_sig['_mask']                            # (n_alleles, n_pos) uint8 bits
    alleles = gene_sig['_alleles']
    freq = gene_sig['_freq']

    # ---- (1) covered positions + real base-set bitmask --------------------
    idx = []; real = []; depth = []; major_frac = []
    for k, p in enumerate(positions):
        c = counts.get(p)
        if not c:
            continue
        d = c['A'] + c['C'] + c['G'] + c['T']
        if d < min_depth:
            continue
        srt = sorted(((c[b], b) for b in _BASES), reverse=True)
        mask = _BIT[srt[0][1]]
        if srt[1][0] / d >= balance_thr and srt[1][0] >= min_minor_count:
            mask |= _BIT[srt[1][1]]
        idx.append(k); real.append(mask); depth.append(d)
        major_frac.append(srt[0][0] / d)
    ncov = len(idx)
    if ncov < min_positions:
        return {'gene': gene_sig['gene'], 'call': None,
                'reason': f'insufficient coverage ({ncov} pos)', 'n_positions': ncov}

    cols = np.array(idx)
    realv = np.array(real, dtype=np.uint8)           # (ncov,)
    A = M[:, cols]                                    # (n_alleles, ncov) allele bit
    nc = (A == 0).sum(axis=1)                         # no-call count per allele

    # ---- (2a) candidate pre-filter by per-allele consistency --------------
    consistent = ((A & realv) != 0).sum(axis=1)      # positions where allele base in real
    order = np.argsort(-consistent)
    cmax = consistent[order[0]]
    cand = order[(consistent[order] >= cmax - 1)][:top_k]
    if len(cand) < 2:
        cand = order[:max(2, top_k // 4)]
    Ac = A[cand]                                      # (K, ncov)
    K = len(cand)

    # ---- (2b) pair set-mismatch cost, vectorised --------------------------
    C = np.empty((K, K), dtype=np.int32)             # pair set-mismatch cost
    for a in range(K):
        C[a] = ((Ac[a] | Ac) != realv).sum(axis=1)
    fc = freq[cand]; ncc = nc[cand]
    Fsum = fc[:, None] + fc[None, :]                 # freq prior per pair
    NCsum = (ncc[:, None] + ncc[None, :]).astype(np.int32)
    hom = np.ones((K, K), dtype=np.int32); np.fill_diagonal(hom, 0)
    ai, bi = np.triu_indices(K)                      # a <= b
    # lexicographic argmin: cost, then -freq, then nc_sum, then prefer hom
    keys = (hom[ai, bi], NCsum[ai, bi], -Fsum[ai, bi], C[ai, bi])
    winner = np.lexsort(keys)[0]
    ia, ib = int(ai[winner]), int(bi[winner])
    bcost = int(C[ia, ib])
    i_al, j_al = int(cand[ia]), int(cand[ib])
    X, Y = alleles[i_al], alleles[j_al]

    # order het alleles by single-allele consistency (major first)
    if X != Y and consistent[j_al] > consistent[i_al]:
        X, Y = Y, X; i_al, j_al = j_al, i_al

    # ---- minor-allele support (QC) + explained fraction -------------------
    predX = A[i_al]; predY = A[j_al]
    predpair = predX | predY
    minor = None; ndisc = 0
    if X != Y:
        disc = np.where(predX != predY)[0]
        ndisc = int(len(disc))
        if ndisc:
            # support of the weaker allele's base over discriminating positions
            supX = []; supY = []
            for k in disc:
                p = positions[idx[k]]; c = counts[p]; d = depth[k]
                bx = _rev_bit(int(predX[k])); by = _rev_bit(int(predY[k]))
                if bx: supX.append(c[bx] / d)
                if by: supY.append(c[by] / d)
            minor = round(min(np.mean(supX) if supX else 1.0,
                              np.mean(supY) if supY else 1.0), 3)
    zyg = 'hom' if X == Y else 'het'
    return {
        'gene': gene_sig['gene'], 'call': (X, Y), 'zygosity': zyg,
        'n_positions': ncov, 'mean_depth': round(sum(depth) / ncov, 1),
        'setmismatch': int(bcost), 'n_discriminating': ndisc,
        'minor_support': minor,
        'explained_frac': round(float(np.mean(major_frac)), 4),
    }


def _rev_bit(bit):
    for b, v in _BIT.items():
        if v == bit:
            return b
    return None


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def type_bam(bam_path, build='auto', genes=('HLA-A', 'HLA-B', 'HLA-C'),
             min_depth=8, verbose=True, **kw):
    if build == 'auto':
        build, evidence = detect_build(bam_path)
        if verbose:
            print(f'[build] detected {build}  (HLA coverage evidence: {evidence})')
    if build not in CHR6_LEN:
        raise RuntimeError(f'could not determine genome build for {bam_path}')
    sig = load_signatures(build)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chrom, _ = _chrom_for(bam, CHR6_LEN[build])
    if chrom is None:
        raise RuntimeError('no chr6 contig found in BAM')
    out = {'bam': bam_path, 'build': build, 'genes': {}}
    for g in genes:
        gs = sig['genes'][g]
        counts = pileup_positions(bam, chrom, gs['positions'])
        res = type_gene(gs, counts, min_depth=min_depth, **kw)
        out['genes'][g] = res
        if verbose:
            _print_gene(res)
    bam.close()
    return out


def _print_gene(res):
    g = res['gene'].replace('HLA-', '')
    if res.get('call'):
        a, b = res['call']
        print(f'  HLA-{g}: {a} / {b}   [{res.get("zygosity")}] '
              f'pos={res["n_positions"]} depth={res.get("mean_depth")} '
              f'expl={res.get("explained_frac")} minor={res.get("minor_support")}')
    else:
        print(f'  HLA-{g}: NO CALL ({res.get("reason")})')


def format_optitype(out):
    """Comma-joined 'HLA-A*nn:nn,...' string like the OptiType aggregate."""
    order = ['HLA-A', 'HLA-B', 'HLA-C']
    parts = []
    for g in order:
        r = out['genes'].get(g)
        if r and r.get('call'):
            parts += [f'HLA-{c}' for c in r['call']]
    return ','.join(parts)


def main():
    ap = argparse.ArgumentParser(description='Fast class-I HLA typing from a genome BAM (pysam pileup)')
    ap.add_argument('--bam', required=True)
    ap.add_argument('--build', default='auto', choices=['auto', 'hg19', 'hg38'])
    ap.add_argument('--genes', default='A,B,C')
    ap.add_argument('--min-depth', type=int, default=8)
    ap.add_argument('--out', default=None, help='TSV output (default: stdout summary only)')
    args = ap.parse_args()
    genes = [f'HLA-{x.strip().upper().replace("HLA-","")}' for x in args.genes.split(',')]
    out = type_bam(args.bam, build=args.build, genes=genes, min_depth=args.min_depth)
    print(f'[hla] {format_optitype(out)}')
    if args.out:
        with open(args.out, 'w') as fh:
            fh.write('gene\tallele_1\tallele_2\tzygosity\tn_positions\tmean_depth\texplained_frac\tminor_support\n')
            for g, r in out['genes'].items():
                if r.get('call'):
                    a, b = r['call']
                else:
                    a = b = 'NA'
                fh.write(f'{g}\t{a}\t{b}\t{r.get("zygosity","NA")}\t{r.get("n_positions")}\t'
                         f'{r.get("mean_depth","NA")}\t{r.get("explained_frac","NA")}\t'
                         f'{r.get("minor_support","NA")}\n')
        print(f'[out] {args.out}')


if __name__ == '__main__':
    main()
