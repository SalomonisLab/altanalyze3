"""
Parser for IPD-IMGT/HLA nucleotide alignment files (alignments/<GENE>_nuc.txt).

The alignment format (one block shown):

     cDNA              1
     AA codon          -24
                       |
     A*01:01:01:01     ATG GCC GTC ATG GCG CCC CGA ACC ... G.|GC
     A*01:01:01:02N    --- --- --- --- --- --- --- --- ... -.|--
     ...

Conventions:
  - The FIRST allele listed is the alignment reference; its row carries the
    literal bases.  In every other row:
        '-'  == identical to the reference base at that column
        ACGT == a substitution
        '.'  == a gap (indel) in this allele at that alignment column
        '*'  == sequence unknown / not determined for this allele
  - Spaces separate codons (purely cosmetic, consistent across rows).
  - '|' marks an exon boundary (same column in every row).
  - cDNA position 1 == the 'A' of the initiating ATG (IMGT cDNA numbering).

This module reconstructs, for every allele, its full gap-aware aligned
sequence (equal length across alleles) plus the exon-boundary column indices,
so downstream code can slice out individual exons (exon 2 / exon 3 = the class-I
antigen-recognition domain that defines 2-field type).
"""

import re
from collections import OrderedDict

_ALLELE_RE = re.compile(r'^[A-Z0-9]+\*\d')


def parse_nuc_alignment(path):
    """Parse an IMGT <GENE>_nuc.txt alignment.

    Returns dict with:
      reference        : name of the reference allele (first listed)
      alleles          : OrderedDict allele_name -> reconstructed aligned seq
                         (chars in {A,C,G,T,.,*}); '-' already resolved to the
                         reference base.  Equal length for every allele.
      exon_boundaries  : sorted list of column indices at which a new exon
                         starts (0-based into the aligned sequence).
      exons            : list of (start,end) column slices, one per exon.
    """
    # ---- pass 1: accumulate per-allele fragments, block by block ----------
    frag = OrderedDict()          # allele -> list of fragment strings
    ref_name = None
    boundary_cols = set()         # column indices (in final coord) of exon starts
    ref_running_len = 0           # length of reference reconstructed so far

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.lstrip().startswith('#'):
                continue
            stripped = line.strip()
            if stripped.startswith('cDNA') or stripped.startswith('AA codon') \
               or stripped.startswith('Please') or stripped.startswith('AA '):
                continue
            # marker line (only '|' and spaces) -> skip; boundaries taken from ref row
            if set(stripped) <= {'|'}:
                continue
            # allele row?
            tok = stripped.split(None, 1)
            name = tok[0]
            if not _ALLELE_RE.match(name):
                continue
            seq_part = tok[1] if len(tok) > 1 else ''
            if ref_name is None:
                ref_name = name
            # for the reference row, record exon-boundary positions from '|'
            if name == ref_name:
                # walk seq_part, counting real columns, noting '|' positions
                col = ref_running_len
                for ch in seq_part:
                    if ch == ' ':
                        continue
                    if ch == '|':
                        boundary_cols.add(col)
                        continue
                    col += 1
                ref_running_len = col
            # store fragment with spaces and '|' removed (keep alignment chars)
            clean = seq_part.replace(' ', '').replace('|', '')
            frag.setdefault(name, []).append(clean)

    if ref_name is None:
        raise ValueError(f'no allele rows found in {path}')

    # ---- pass 2: concatenate fragments; resolve '-' against reference ------
    alleles = OrderedDict()
    for name, parts in frag.items():
        alleles[name] = ''.join(parts)
    ref_seq = alleles[ref_name]
    L = len(ref_seq)

    resolved = OrderedDict()
    for name, seq in alleles.items():
        if len(seq) != L:
            # tolerate truncated trailing rows by padding with '*'
            if len(seq) < L:
                seq = seq + '*' * (L - len(seq))
            else:
                seq = seq[:L]
        if name == ref_name:
            resolved[name] = seq
            continue
        out = []
        rs = ref_seq
        for i, ch in enumerate(seq):
            out.append(rs[i] if ch == '-' else ch)
        resolved[name] = ''.join(out)

    boundaries = sorted(b for b in boundary_cols if 0 < b < L)
    # build exon slices from boundaries (exon starts at 0 and at each boundary)
    starts = [0] + boundaries
    ends = boundaries + [L]
    exons = list(zip(starts, ends))

    return {
        'reference': ref_name,
        'alleles': resolved,
        'exon_boundaries': boundaries,
        'exons': exons,
        'length': L,
    }


def two_field(name):
    """Collapse an allele name to 2-field resolution, e.g. A*02:01:01:02 -> A*02:01.
    Preserves an expression suffix letter (N,L,S,Q) only if it is a full null
    at 2 fields (rare); otherwise dropped for typing purposes."""
    m = re.match(r'^([A-Z0-9]+)\*(\d+):(\d+)', name)
    if not m:
        return name
    return f'{m.group(1)}*{m.group(2)}:{m.group(3)}'


if __name__ == '__main__':
    import sys, os
    base = os.path.join(os.path.dirname(__file__), 'data', 'imgt', 'alignments')
    for gene in ('A', 'B', 'C'):
        p = os.path.join(base, f'{gene}_nuc.txt')
        r = parse_nuc_alignment(p)
        ex = r['exons']
        exlens = [e - s for s, e in ex]
        n2 = len(set(two_field(a) for a in r['alleles']))
        print(f'{gene}_nuc: ref={r["reference"]}  n_alleles={len(r["alleles"])}  '
              f'2field={n2}  aln_len={r["length"]}  n_exons={len(ex)}')
        print(f'   exon lengths (cols): {exlens}')
        # exon 2 and exon 3 (1-based exon 2 = index 1)
        if len(ex) >= 3:
            print(f'   exon2 cols={ex[1]} len={exlens[1]}   exon3 cols={ex[2]} len={exlens[2]}')
