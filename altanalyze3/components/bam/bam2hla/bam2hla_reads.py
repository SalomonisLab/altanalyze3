#!/usr/bin/env python3
"""
bam2hla_reads - read-level class-I HLA typing (bam2hla v2), HLA-A/-B/-C at 2 fields.

Why v2
------
bam2hla v1 (bam2hla.py) piles up bases at genome positions. RNA reads from one HLA gene
routinely align to another gene's locus: at HLA-A exon 3 of one benchmark sample, 20,918 of
30,281 reads came from HLA-B/-C. Those reads add a second base at many positions, and the v1
set-mismatch solver then picks a rare allele that carries it. v2 never trusts the locus a
read was aligned to; the BAM only supplies the reads of the MHC window.

Method
------
1. Reference (one-time, cached): every expressed IMGT/HLA allele of HLA-A/-B/-C (null 'N'
   alleles excluded) in its gene's alignment layout (A/B/C_nuc.txt); a panel of the 2-field
   alleles with a population frequency for gene assignment; and the 2-field sequences of the
   other class-I genes and pseudogenes in hla_nuc.fasta as competitors.
2. Reads: primary alignments in the MHC window, collapsed to unique sequences.
3. Placement: for each typed gene, the read's 25-mers vote for its start in the gene layout
   (both read orientations).
4. Gene: the read's minimum mismatch count against each gene's panel (and against every
   competitor sequence) decides; the read is used only for a typed gene that beats every
   other gene by at least `margin` mismatches.
5. Fit: each 2-field allele is the set of bases its expressed members carry at each column, so
   synonymous (3rd-field) variants never count against it. Only exons 2-4 (TYPING_EXONS)
   separate alleles; outside them every allele of the gene accepts any base the gene carries. Cost of an allele pair = sum over
   reads of the read's mismatches to the closer allele, over the gene's polymorphic columns. Homozygous unless the best pair lowers the best
   homozygous cost by `het_min` mismatches per read. Equal costs go to the higher
   population frequency, then the allele name, so the call is the same on every host.
"""

import os, sys, gzip, pickle, hashlib, argparse, time
from collections import Counter, defaultdict

import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
try:                                     # package import (altanalyze3.components.bam.bam2hla)
    from . import bam2hla as V1          # build detection, contig lookup, frequency prior
    from . import imgt_parser as IP
except ImportError:                      # run as a script from this directory
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import bam2hla as V1
    import imgt_parser as IP

K = 25
TYPED = ('A', 'B', 'C')
COMPETITORS = ('E', 'F', 'G', 'H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y')
MHC = {'hg38': (29600000, 33400000), 'hg19': (29570000, 33370000)}
IMGT_DIR = os.path.join(HERE, 'data', 'imgt')
NUC_FASTA = os.path.join(IMGT_DIR, 'hla_nuc.fasta')
ALN = {g: os.path.join(IMGT_DIR, 'alignments', f'{g}_nuc.txt') for g in TYPED}
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b'ACGT'):
    CODE[_b] = _i
UNK, DEL = 6, 5                          # allele-matrix codes besides 0-3
REF_VERSION = 4
TYPING_EXONS = (2, 3, 4)          # exons whose bases decide the allele (2-3 = peptide groove)
STRIDE = 5                        # read 25-mers used at every 5th offset


def is_null(name):
    return name.endswith('N')


# ----------------------------------------------------------------------------
# k-mers
# ----------------------------------------------------------------------------
def _kmer_codes_2d(arr):
    """(n, L) uint8 codes -> (n, L-K+1) int64 25-mer codes; -1 where a window holds a non-ACGT."""
    n, L = arr.shape
    m = L - K + 1
    if m <= 0:
        return np.full((n, 0), -1, dtype=np.int64)
    a = arr.astype(np.int64)
    cb = np.zeros((n, L + 1), dtype=np.int32)
    cb[:, 1:] = np.cumsum(a > 3, axis=1)
    code = np.zeros((n, m), dtype=np.int64)
    for j in range(K):
        code = (code << 2) | (a[:, j:j + m] & 3)
    code[(cb[:, K:K + m] - cb[:, :m]) > 0] = -1
    return code


def _kmer_codes(arr):
    return _kmer_codes_2d(arr[None, :])[0]


def _kmer_offsets(L):
    m = L - K + 1
    if m <= 0:
        return np.zeros(0, dtype=np.int64)
    return np.unique(np.r_[np.arange(0, m, STRIDE), m - 1]).astype(np.int64)


def _kmer_codes_at(arr, offsets):
    """25-mer codes of each read at the given offsets only (n, len(offsets)); -1 if invalid."""
    n, L = arr.shape
    a = arr.astype(np.int64)
    cb = np.zeros((n, L + 1), dtype=np.int32)
    cb[:, 1:] = np.cumsum(a > 3, axis=1)
    code = np.zeros((n, len(offsets)), dtype=np.int64)
    for j in range(K):
        code = (code << 2) | (a[:, offsets + j] & 3)
    code[(cb[:, offsets + K] - cb[:, offsets]) > 0] = -1
    return code


def _read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(seq)
                name, seq = line.split()[1], []
            else:
                seq.append(line.strip().upper())
    if name:
        yield name, ''.join(seq)


# ----------------------------------------------------------------------------
# reference
# ----------------------------------------------------------------------------
def _ref_key():
    h = hashlib.sha1()
    for p in [NUC_FASTA] + [ALN[g] for g in TYPED] + [V1._FREQ_CSV]:
        st = os.stat(p)
        h.update(f'{p}:{st.st_size}:{int(st.st_mtime)}'.encode())
    h.update(f'K={K};v={REF_VERSION}'.encode())
    return h.hexdigest()[:12]


def _group_rows(member_of, n_groups):
    order = np.argsort(member_of, kind='stable')
    bounds = np.searchsorted(member_of[order], np.arange(n_groups + 1))
    for gi in range(n_groups):
        yield gi, order[bounds[gi]:bounds[gi + 1]]


def typing_masks(gref, exons=TYPING_EXONS):
    """Group masks for allele choice: outside `exons` every column accepts any base the gene carries."""
    key = ('_Gtyp', tuple(exons))
    if key not in gref:
        G = gref['G'].copy()
        out = ~np.isin(gref['exon_of_layout'], list(exons))
        G[:, out] = np.bitwise_or.reduce(G[:, out], axis=0)
        gref[key] = G
    return gref[key]


def _mask_onehot(G):
    """(n, L) base-set bitmasks -> (n, 4L) float32 indicator of each allowed base."""
    n, L = G.shape
    O = np.zeros((n, L, 4), dtype=np.float32)
    for c in range(4):
        O[:, :, c] = (G >> c) & 1
    return O.reshape(n, 4 * L)


def build_reference(verbose=True):
    t0 = time.time()
    freq2 = V1.load_frequencies()
    genes = {}
    for g in TYPED:
        aln = IP.parse_nuc_alignment(ALN[g])
        names = list(aln['alleles'])
        ncol = len(aln['alleles'][aln['reference']])
        A = np.frombuffer(''.join(aln['alleles'][n] for n in names).encode(),
                          dtype=np.uint8).reshape(len(names), ncol)
        ref_row = A[names.index(aln['reference'])]
        layout_cols = np.where(ref_row != ord('.'))[0]
        col2lay = np.full(ncol, -1, dtype=np.int64)
        col2lay[layout_cols] = np.arange(len(layout_cols))
        sub = A[:, layout_cols]
        Mall = np.full(sub.shape, UNK, dtype=np.uint8)
        for code, ch in enumerate(b'ACGT'):
            Mall[sub == ch] = code
        Mall[sub == ord('.')] = DEL
        expressed = np.array([not is_null(n) for n in names])
        M = Mall[expressed]
        full_names = [n for n, e in zip(names, expressed) if e]
        tf = [IP.two_field(n) for n in full_names]
        freq = np.array([freq2.get(t, 0.0) for t in tf])
        known = (M < 4).sum(axis=1)
        poly = []
        for k in range(M.shape[1]):
            col = M[:, k]
            kb = np.unique(col[col < 4])
            if len(kb) >= 2 or (len(kb) >= 1 and (col == DEL).any()):
                poly.append(k)
        # gene-assignment panel: the most complete member of each 2-field group with freq > 0
        best = {}
        for i, t in enumerate(tf):
            if freq[i] > 0 and (t not in best or known[i] > known[best[t]]):
                best[t] = i
        panel = np.array(sorted(best.values()), dtype=np.int64)
        # placement 25-mers from every allele (null ones too: they are real sequence)
        pc, pv = [], []
        for i in range(len(names)):
            row = A[i]
            cols = np.where(row != ord('.'))[0]
            codes = _kmer_codes(CODE[row[cols]])
            if len(codes) == 0:
                continue
            lay = col2lay[cols]
            a, b = lay[:len(codes)], lay[K - 1:K - 1 + len(codes)]
            okp = (codes >= 0) & (a >= 0) & (b - a == K - 1)
            pc.append(codes[okp]); pv.append(a[okp])
        pc = np.concatenate(pc); pv = np.concatenate(pv)
        o = np.lexsort((pv, pc)); pc, pv = pc[o], pv[o]
        uc, fi, cnt = np.unique(pc, return_index=True, return_counts=True)
        lo, hi = pv[fi], pv[fi + cnt - 1]
        keep = lo == hi                                   # a 25-mer at two starts is not used
        # 2-field group profile: at each column, the set of bases any expressed member carries
        groups = sorted(set(tf))
        gidx = {t: i for i, t in enumerate(groups)}
        G = np.zeros((len(groups), M.shape[1]), dtype=np.uint8)
        member_of = np.array([gidx[t] for t in tf])
        for c in range(4):
            hit = (M == c)
            for gi_, rows_ in _group_rows(member_of, len(groups)):
                G[gi_] |= (hit[rows_].any(axis=0).astype(np.uint8) << c)
        gfreq = np.array([freq2.get(t, 0.0) for t in groups])
        gpanel = np.where(gfreq > 0)[0]
        exon_of_layout = np.zeros(len(layout_cols), dtype=np.int64)
        for e, (c0, c1) in enumerate(aln['exons'], 1):
            sel_ = col2lay[c0:c1]
            exon_of_layout[sel_[sel_ >= 0]] = e
        genes[g] = {'M': M, 'names': full_names, 'tf': tf, 'freq': freq, 'poly': np.array(poly),
                    'groups': groups, 'G': G, 'gfreq': gfreq, 'gpanel': gpanel, 'exon_of_layout': exon_of_layout,
                    'n_layout': len(layout_cols), 'panel': panel,
                    'kmer_codes': uc[keep], 'kmer_start': lo[keep], 'reference': aln['reference']}
        if verbose:
            print(f'[ref] HLA-{g}: {len(names)} alleles ({len(full_names)} expressed, '
                  f'{len(set(tf))} 2-field), panel {len(panel)}, layout {len(layout_cols)}, '
                  f'{len(poly)} polymorphic columns, {int(keep.sum())} placement 25-mers '
                  f'({time.time() - t0:.0f}s)', flush=True)
    # competitors: one sequence per 2-field group (the longest expressed member)
    rep = {}
    for name, seq in _read_fasta(NUC_FASTA):
        g = name.split('*')[0]
        if g not in COMPETITORS or is_null(name):
            continue
        t = IP.two_field(name)
        if t not in rep or len(seq) > len(rep[t][1]):
            rep[t] = (g, seq)
    comp_names = sorted(rep)
    parts, sid = [], []
    for i, t in enumerate(comp_names):
        s = CODE[np.frombuffer(rep[t][1].encode(), dtype=np.uint8)]
        parts += [s, np.array([4], dtype=np.uint8)]
        sid += [i] * (len(s) + 1)
    concat = np.concatenate(parts); sid = np.array(sid, dtype=np.int64)
    cc = _kmer_codes(concat)
    okc = cc >= 0
    cpos = np.where(okc)[0]; cc = cc[okc]
    o = np.argsort(cc, kind='stable')
    comp = {'names': comp_names, 'gene': [rep[t][0] for t in comp_names], 'concat': concat,
            'sid': sid, 'kmer_codes': cc[o], 'kmer_pos': cpos[o]}
    if verbose:
        print(f'[ref] competitors: {len(comp_names)} 2-field sequences over '
              f'{len(set(comp["gene"]))} genes ({time.time() - t0:.0f}s)', flush=True)
    return {'key': _ref_key(), 'genes': genes, 'comp': comp}


_REF = None


def load_reference(cache_dir=None, verbose=True):
    global _REF
    if _REF is not None:
        return _REF
    cache_dir = cache_dir or os.environ.get('BAM2HLA_CACHE') or os.path.join(HERE, 'build')
    path = os.path.join(cache_dir, f'bam2hla_reads_ref_{_ref_key()}.pkl.gz')
    if os.path.exists(path):
        with gzip.open(path, 'rb') as fh:
            _REF = pickle.load(fh)
        return _REF
    _REF = build_reference(verbose)
    os.makedirs(cache_dir, exist_ok=True)
    tmp = path + f'.part{os.getpid()}'
    with gzip.open(tmp, 'wb', compresslevel=3) as fh:
        pickle.dump(_REF, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)
    return _REF


# ----------------------------------------------------------------------------
# reads
# ----------------------------------------------------------------------------
def collect_reads(bam_path, build):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chrom, _ = V1._chrom_for(bam, V1.CHR6_LEN[build])
    if chrom is None:
        raise RuntimeError('no chr6 contig found in BAM')
    lo, hi = MHC[build]
    cnt = Counter()
    n = 0
    for a in bam.fetch(chrom, lo, hi):
        if a.is_secondary or a.is_supplementary or a.is_qcfail or a.query_sequence is None:
            continue
        n += 1
        cnt[a.query_sequence] += 1
    bam.close()
    return cnt, n


def _encode(seqs):
    """Unique read sequences -> (n, Lmax) uint8 base codes, 4 = N or padding."""
    L = max(len(x) for x in seqs)
    arr = np.full((len(seqs), L), 4, dtype=np.uint8)
    by_len = defaultdict(list)
    for i, x in enumerate(seqs):
        by_len[len(x)].append(i)
    for ln, idx in by_len.items():
        blob = np.frombuffer(''.join(seqs[i] for i in idx).encode(), dtype=np.uint8).reshape(len(idx), ln)
        arr[np.array(idx), :ln] = CODE[blob]
    return arr


def _rc(arr):
    return np.where(arr < 4, 3 - arr, 4).astype(np.uint8)[:, ::-1]


def _place(codes, gref, offsets):
    """Per read: layout start voted by its placement 25-mers (int64) and the vote count."""
    kc = gref['kmer_codes']
    flat = codes.ravel()
    ix = np.searchsorted(kc, flat)
    ix[ix >= len(kc)] = 0
    hit = ((flat >= 0) & (kc[ix] == flat)).reshape(codes.shape)
    st = gref['kmer_start'][ix].reshape(codes.shape) - offsets[None, :]
    big = 1 << 40
    smin = np.where(hit, st, big).min(axis=1)
    smax = np.where(hit, st, -big).max(axis=1)
    nv = hit.sum(axis=1)
    start = np.where(smin == smax, smin, 0)
    votes = np.where(smin == smax, nv, 0)
    for r in np.where((nv >= 2) & (smin != smax))[0]:              # indel or repeat: mode
        v, c = np.unique(st[r][hit[r]], return_counts=True)
        if (c == c.max()).sum() == 1:
            start[r] = v[c.argmax()]; votes[r] = c.max()
    return start, votes


def _onehot_rows(M):
    n, L = M.shape
    O = np.zeros((n, L * 4), dtype=np.float32)
    for c in range(4):
        rr, cc = np.nonzero(M == c)
        O[rr, cc * 4 + c] = 1.0
    return O


def _read_onehot(arr, start, Lg):
    """Dense one-hot (n, 4*Lg) of read bases on the layout, and covered-column count."""
    n, L = arr.shape
    P = start[:, None] + np.arange(L)[None, :]
    ok = (P >= 0) & (P < Lg) & (arr < 4)
    X = np.zeros((n, 4 * Lg), dtype=np.float32)
    rr, cc = np.nonzero(ok)
    X[rr, P[rr, cc] * 4 + arr[rr, cc]] = 1.0
    return X, ok.sum(axis=1)


def _read_sparse(arr, start, Lg):
    """Sparse one-hot (n, 4*Lg) of read bases on the layout, and covered-column count."""
    from scipy import sparse
    n, L = arr.shape
    P = start[:, None] + np.arange(L)[None, :]
    ok = (P >= 0) & (P < Lg) & (arr < 4)
    rr, cc = np.nonzero(ok)
    X = sparse.csr_matrix((np.ones(len(rr), dtype=np.float32), (rr, P[rr, cc] * 4 + arr[rr, cc])),
                          shape=(n, 4 * Lg))
    return X, ok.sum(axis=1)


def _panel_mismatch(arr, start, gref):
    """Fewest mismatching covered columns of each placed read against the gene's panel groups
    (a column no member knows counts as a mismatch)."""
    Op = _mask_onehot(gref['G'][gref['gpanel']])
    X, cov = _read_sparse(arr, start, gref['n_layout'])
    match = np.asarray((X @ Op.T).max(axis=1)).ravel() if X.shape[0] else np.zeros(0)
    return (cov - np.rint(match)).astype(np.int64), cov


def _competitor_mismatch(arr_by_ori, codes_by_ori, offsets, comp, seeds=(0, 15, 30, 45, 60, 75)):
    """Fewest unmatched read bases of each read (either orientation) against any competitor sequence; 10**6 if none."""
    kc, kp, concat, sid = comp['kmer_codes'], comp['kmer_pos'], comp['concat'], comp['sid']
    n, L = arr_by_ori[0].shape
    best = np.full(n, 10 ** 6, dtype=np.int64)
    col_of = {int(o): i for i, o in enumerate(offsets)}
    for arr, codes in zip(arr_by_ori, codes_by_ori):
        rows, starts, refs = [], [], []
        for s in seeds:
            if s not in col_of:
                continue
            c = codes[:, col_of[s]]
            lo = np.searchsorted(kc, c, 'left'); hi = np.searchsorted(kc, c, 'right')
            cntk = np.where(c >= 0, hi - lo, 0)
            if cntk.sum() == 0:
                continue
            r_rep = np.repeat(np.arange(n), cntk)
            offs = np.arange(cntk.sum()) - np.repeat(np.cumsum(cntk) - cntk, cntk)
            gp = kp[np.repeat(lo, cntk) + offs]
            rows.append(r_rep); starts.append(gp - s); refs.append(sid[gp])
        if not rows:
            continue
        rows = np.concatenate(rows); starts = np.concatenate(starts); refs = np.concatenate(refs)
        shift = 1 << 20
        _, first = np.unique(rows * (1 << 32) + (starts + shift), return_index=True)
        rows, starts, refs = rows[first], starts[first], refs[first]
        for c0 in range(0, len(rows), 100000):
            r = rows[c0:c0 + 100000]; st = starts[c0:c0 + 100000]; rs = refs[c0:c0 + 100000]
            pos = st[:, None] + np.arange(L)[None, :]
            posc = np.clip(pos, 0, len(concat) - 1)
            same = (pos >= 0) & (pos < len(concat)) & (sid[posc] == rs[:, None])
            rb = arr[r]; ref_b = concat[posc]
            nb = (rb < 4).sum(axis=1)
            match = (same & (rb < 4) & (rb == ref_b)).sum(axis=1)
            np.minimum.at(best, r, nb - match)               # read bases not matched
    return best


def classify_reads(arr, ref, margin=2, max_rate=0.06):
    """Gene (index into TYPED, -1 unused), layout start and oriented read codes for each read."""
    n, L = arr.shape
    rc = _rc(arr)
    offsets = _kmer_offsets(L)
    codes = (_kmer_codes_at(arr, offsets), _kmer_codes_at(rc, offsets))
    big = 10 ** 6
    gene_mm = np.full((n, len(TYPED)), big, dtype=np.int64)
    gene_cov = np.zeros((n, len(TYPED)), dtype=np.int64)
    gene_start = np.zeros((n, len(TYPED)), dtype=np.int64)
    gene_rc = np.zeros((n, len(TYPED)), dtype=bool)
    for gi, g in enumerate(TYPED):
        gref = ref['genes'][g]
        s0, v0 = _place(codes[0], gref, offsets)
        s1, v1 = _place(codes[1], gref, offsets)
        use_rc = v1 > v0
        st = np.where(use_rc, s1, s0); votes = np.where(use_rc, v1, v0)
        idx = np.where(votes >= 2)[0]
        if len(idx):
            oriented = np.where(use_rc[idx][:, None], rc[idx], arr[idx])
            mm, cov = _panel_mismatch(oriented, st[idx], gref)
            gene_mm[idx, gi] = mm; gene_cov[idx, gi] = cov
        gene_start[:, gi] = st; gene_rc[:, gi] = use_rc
    placed_any = np.where((gene_mm < big).any(axis=1))[0]
    comp_mm = np.full(n, big, dtype=np.int64)
    if len(placed_any):
        comp_mm[placed_any] = _competitor_mismatch((arr[placed_any], rc[placed_any]),
                                                   (codes[0][placed_any], codes[1][placed_any]),
                                                   offsets, ref['comp'])
    rows = np.arange(n)
    best_gi = gene_mm.argmin(axis=1)
    best_mm = gene_mm[rows, best_gi]
    other = np.minimum(np.sort(gene_mm, axis=1)[:, 1], comp_mm)
    cov = gene_cov[rows, best_gi]
    use = (best_mm < big) & (cov >= K) & (best_mm + margin <= other) & (best_mm <= np.maximum(3, max_rate * cov))
    gene = np.where(use, best_gi, -1)
    start = gene_start[rows, best_gi]
    oriented = np.where(gene_rc[rows, best_gi][:, None], rc, arr)
    return gene, start, oriented, {'competitor_best': comp_mm, 'gene_mm': gene_mm, 'gene_start': gene_start,
                                   'gene_rc': gene_rc, 'placed': gene_mm < big, 'rc': rc}


# ----------------------------------------------------------------------------
# fit
# ----------------------------------------------------------------------------
def type_gene_reads(gref, arr, start, weight, het_min=0.10, n_cand=200, tie_tol=0.0, typing_exons=TYPING_EXONS):
    """Diploid 2-field call for one gene from its reads (module docstring, step 5)."""
    from scipy import sparse
    Lg = gref['n_layout']
    if len(arr) == 0:
        return {'call': None, 'reason': 'no reads assigned', 'n_reads': 0}
    is_poly = np.zeros(Lg, dtype=bool); is_poly[gref['poly']] = True
    P = start[:, None] + np.arange(arr.shape[1])[None, :]
    sel = (P >= 0) & (P < Lg) & (arr < 4)
    sel &= is_poly[np.clip(P, 0, Lg - 1)]
    keep = sel.any(axis=1)
    if not keep.any():
        return {'call': None, 'reason': 'no informative reads', 'n_reads': 0}
    arr, P, sel, w = arr[keep], P[keep], sel[keep], weight[keep].astype(np.float64)
    ri, ci = np.nonzero(sel)
    feat = P[ri, ci] * 4 + arr[ri, ci].astype(np.int64)
    X = sparse.csr_matrix((np.ones(len(ri), dtype=np.float32), (ri, feat)), shape=(len(w), 4 * Lg))
    X.sum_duplicates()
    keyb = [X.indices[X.indptr[r]:X.indptr[r + 1]].tobytes() for r in range(X.shape[0])]
    _, first, inv = np.unique(np.array(keyb, dtype=object), return_index=True, return_inverse=True)
    w = np.bincount(inv.ravel(), weights=w)
    X = X[first]
    m = np.asarray(X.sum(axis=1)).ravel()
    W = w.sum()
    Mo = _mask_onehot(typing_masks(gref, typing_exons))
    pile = np.asarray(X.T @ w).ravel()
    hom_cost = (w * m).sum() - Mo @ pile                    # every 2-field group, homozygous
    tf = np.array(gref['groups']); freq = gref['gfreq']
    first_rank = np.lexsort((tf, -freq, hom_cost))
    # residual candidates: the groups that best fit the reads the best homozygous group misses
    x0 = first_rank[0]
    miss = (m - np.asarray(X @ Mo[x0]).ravel()) > 0
    if miss.any():
        pile_r = np.asarray(X[miss].T @ w[miss]).ravel()
        res_rank = np.lexsort((tf, -freq, (w[miss] * m[miss]).sum() - Mo @ pile_r))
    else:
        res_rank = first_rank
    cand = list(dict.fromkeys(list(first_rank[:n_cand // 2]) + list(res_rank[:n_cand // 2])))
    for t in first_rank[n_cand // 2:]:
        if len(cand) >= n_cand:
            break
        if t not in cand:
            cand.append(t)
    cand = np.array(cand)
    mm = (m[:, None] - np.asarray(X @ Mo[cand].T)).astype(np.int32)
    prof, pinv = np.unique(mm, axis=0, return_inverse=True)
    pw = np.bincount(pinv.ravel(), weights=w)
    k = len(cand)
    cost = np.empty((k, k), dtype=np.float64)
    for i in range(k):
        cost[i] = pw @ np.minimum(prof[:, i:i + 1], prof)
    iu = np.triu_indices(k, 1)
    het_c = cost[iu]; hom_c = np.diag(cost).copy()
    best_het = het_c.min() if len(het_c) else np.inf
    best_hom = hom_c.min()
    fr = freq[cand]; nm = tf[cand]

    def pick(pc, ia, ib, best):
        within = np.where(pc <= best + tie_tol * W)[0]
        o = sorted(within, key=lambda t: (pc[t], -(fr[ia[t]] + fr[ib[t]]), nm[ia[t]], nm[ib[t]]))
        return ia[o[0]], ib[o[0]]
    gain = (best_hom - best_het) / W
    if gain >= het_min:
        a, b = pick(het_c, iu[0], iu[1], best_het)
        X1, Y1 = nm[a], nm[b]
        if hom_c[b] < hom_c[a]:
            X1, Y1 = Y1, X1
    else:
        d = np.arange(k)
        a, _ = pick(hom_c, d, d, best_hom)
        X1 = Y1 = nm[a]
    return {'call': (str(X1), str(Y1)), 'zygosity': 'hom' if X1 == Y1 else 'het', 'n_reads': float(W),
            'candidates': [str(x) for x in nm], 'n_profiles': int(len(prof)), 'het_gain_per_read': round(float(gain), 4),
            'mismatch_per_read': round(float(min(best_het, best_hom) / W), 4)}


def _unmatched(arr, start, Lg, Mo_c):
    """Read bases not matched by each candidate group (bases outside the layout count as unmatched)."""
    X, _ = _read_sparse(arr, start, Lg)
    nb = (arr < 4).sum(axis=1)
    return (nb[:, None] - np.rint(np.asarray(X @ Mo_c.T))).astype(np.int32)


def joint_refine(ref, arr, w, info, calls, cands, rounds=3, het_min=0.10, n_cand=60, typing_exons=TYPING_EXONS):
    """Re-fit each gene against every placed read, the other genes' current alleles competing.

    calls: {gene: (X, Y)} from the gene-specific fit. A read's cost under a genotype is its
    fewest unmatched bases over the six alleles and the competitor sequences."""
    big = 10 ** 6
    rc = info['rc']
    per = {}
    for gi, g in enumerate(TYPED):
        gref = ref['genes'][g]
        idx = np.where(info['placed'][:, gi])[0]
        ori = np.where(info['gene_rc'][idx, gi][:, None], rc[idx], arr[idx])
        st = info['gene_start'][idx, gi]
        groups = np.array(gref['groups'])
        Mo = _mask_onehot(typing_masks(gref, typing_exons))
        # candidates: the gene-specific fit's best n_cand, the best on all placed reads, the call
        X, _ = _read_sparse(ori, st, gref['n_layout'])
        pile = np.asarray(X.T @ w[idx].astype(np.float64)).ravel()
        hom = -(Mo @ pile)
        gidx = {t: i for i, t in enumerate(groups)}
        cand = [gidx[t] for t in cands.get(g, [])[:n_cand]]
        for j in np.lexsort((groups, -gref['gfreq'], hom))[:n_cand // 3]:
            if j not in cand:
                cand.append(int(j))
        for a in calls.get(g, ()) or ():
            if gidx[a] not in cand:
                cand.append(gidx[a])
        cand = np.array(cand)
        U = np.full((len(arr), len(cand)), big, dtype=np.int32)
        U[idx] = _unmatched(ori, st, gref['n_layout'], Mo[cand])
        per[g] = {'cand': cand, 'names': groups[cand], 'freq': gref['gfreq'][cand], 'U': U}
    comp = np.minimum(info['competitor_best'], big).astype(np.int32)
    cur = {}
    for g in TYPED:
        nm = list(per[g]['names'])
        c = calls.get(g)
        cur[g] = (nm.index(c[0]), nm.index(c[1])) if c else (0, 0)
    W = float(w.sum())
    stats = {}
    for _ in range(rounds):
        changed = False
        for g in TYPED:
            cap = comp.copy()
            for h in TYPED:
                if h != g:
                    Uh = per[h]['U']; a, b = cur[h]
                    cap = np.minimum(cap, np.minimum(Uh[:, a], Uh[:, b]))
            U = per[g]['U']
            rel = (U.min(axis=1) < cap)                      # reads this gene could explain better
            if not rel.any():
                continue
            prof, inv = np.unique(np.concatenate([np.minimum(U[rel], cap[rel, None]), cap[rel, None]], axis=1),
                                  axis=0, return_inverse=True)
            pw = np.bincount(inv.ravel(), weights=w[rel])
            P = prof[:, :-1]
            k = P.shape[1]
            cost = np.empty((k, k))
            for i in range(k):
                cost[i] = pw @ np.minimum(P[:, i:i + 1], P)
            iu = np.triu_indices(k, 1)
            het_c, hom_c = cost[iu], np.diag(cost).copy()
            fr, nm = per[g]['freq'], per[g]['names']
            Wg = float(pw.sum())
            gain = (hom_c.min() - het_c.min()) / Wg
            if gain >= het_min:
                within = np.where(het_c <= het_c.min())[0]
                t = sorted(within, key=lambda t: (-(fr[iu[0][t]] + fr[iu[1][t]]), nm[iu[0][t]], nm[iu[1][t]]))[0]
                new = (int(iu[0][t]), int(iu[1][t]))
            else:
                within = np.where(hom_c <= hom_c.min())[0]
                t = sorted(within, key=lambda t: (-fr[t], nm[t]))[0]
                new = (int(t), int(t))
            if set(new) != set(cur[g]):
                changed = True
            cur[g] = new
            stats[g] = {'joint_reads': Wg, 'joint_het_gain': round(float(gain), 4)}
        if not changed:
            break
    out = {}
    for g in TYPED:
        nm = per[g]['names']; a, b = cur[g]
        X1, Y1 = str(nm[a]), str(nm[b])
        out[g] = ((X1, Y1), stats.get(g, {}))
    return out


def type_from_counts(cnt, ref, margin=2, joint=True, **kw):
    """Type HLA-A/-B/-C from a Counter of read sequences (the BAM-free core)."""
    seqs = list(cnt)
    out = {'n_unique': len(seqs), 'genes': {}}
    if not seqs:
        for g in TYPED:
            out['genes'][f'HLA-{g}'] = {'call': None, 'reason': 'no reads', 'gene': f'HLA-{g}', 'n_assigned_reads': 0}
        return out
    arr = _encode(seqs)
    w = np.array([cnt[x] for x in seqs], dtype=np.int64)
    gene, start, oriented, info = classify_reads(arr, ref, margin=margin)
    calls, cands = {}, {}
    for gi, g in enumerate(TYPED):
        m = gene == gi
        res = type_gene_reads(ref['genes'][g], oriented[m], start[m], w[m], **kw)
        res['gene'] = f'HLA-{g}'
        res['n_assigned_reads'] = int(w[m].sum())
        out['genes'][f'HLA-{g}'] = res
        if res.get('call'):
            calls[g] = res['call']; cands[g] = res.get('candidates', [])
    if joint and calls:
        ref_calls = joint_refine(ref, arr, w, info, calls, cands, het_min=kw.get('het_min', 0.10),
                                 typing_exons=kw.get('typing_exons', TYPING_EXONS))
        for g, (call, st) in ref_calls.items():
            r = out['genes'][f'HLA-{g}']
            r['gene_specific_call'] = r.get('call')
            r['call'] = call
            r['zygosity'] = 'hom' if call[0] == call[1] else 'het'
            r.update(st)
    out['n_unused_reads'] = int(w[gene < 0].sum())
    return out


def type_bam_reads(bam_path, build='auto', verbose=True, **kw):
    t0 = time.time()
    if build == 'auto':
        build, _ = V1.detect_build(bam_path)
    if build not in MHC:
        raise RuntimeError(f'could not determine genome build for {bam_path}')
    ref = load_reference(verbose=verbose)
    cnt, n_primary = collect_reads(bam_path, build)
    out = type_from_counts(cnt, ref, **kw)
    out.update({'bam': bam_path, 'build': build, 'n_primary_reads': n_primary})
    out['seconds'] = round(time.time() - t0, 1)
    if verbose:
        for g, r in out['genes'].items():
            c = '/'.join(r['call']) if r.get('call') else 'NO CALL ' + r.get('reason', '')
            print(f'  {g}: {c}  reads={r["n_assigned_reads"]} profiles={r.get("n_profiles")} '
                  f'mm/read={r.get("mismatch_per_read")} het_gain={r.get("het_gain_per_read")}')
        print(f'  [{out["seconds"]}s, {out["n_primary_reads"]} primary reads, {out["n_unique"]} unique]')
    return out


def format_optitype(out):
    parts = []
    for g in ('HLA-A', 'HLA-B', 'HLA-C'):
        r = out['genes'].get(g)
        if r and r.get('call'):
            parts += [f'HLA-{c}' for c in r['call']]
    return ','.join(parts)


def main():
    ap = argparse.ArgumentParser(description='Read-level class-I HLA typing from a genome BAM (bam2hla v2)')
    ap.add_argument('--bam', required=True)
    ap.add_argument('--build', default='auto', choices=['auto', 'hg19', 'hg38'])
    ap.add_argument('--het-min', type=float, default=0.10)
    ap.add_argument('--margin', type=int, default=2)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = type_bam_reads(a.bam, build=a.build, het_min=a.het_min, margin=a.margin)
    print(f'[hla] {format_optitype(out)}')
    if a.out:
        with open(a.out, 'w') as fh:
            fh.write('gene\tallele_1\tallele_2\tzygosity\tn_assigned_reads\tn_profiles\tmismatch_per_read\thet_gain_per_read\n')
            for g, r in out['genes'].items():
                x, y = r['call'] if r.get('call') else ('NA', 'NA')
                fh.write(f'{g}\t{x}\t{y}\t{r.get("zygosity","NA")}\t{r.get("n_assigned_reads")}\t'
                         f'{r.get("n_profiles")}\t{r.get("mismatch_per_read")}\t{r.get("het_gain_per_read")}\n')


if __name__ == '__main__':
    main()
