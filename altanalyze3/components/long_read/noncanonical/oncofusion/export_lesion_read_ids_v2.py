"""Export the read IDs behind the 3 documented lesions reported on 2026-11-10 (CBFB::MYH11 in AML14, KMT2A-PTD in the
U2AF1-RAS patient, FLT3-ITD in pAML-rel-IDH2). Ad hoc; validate within altanalyze3.

Source: the per-BAM scan outputs samples/NN_<uid>__<library>/ (NN = line NN+1 of mds_metadata_bam_fixed.txt), so each
read carries the BAM it came from (AML14 and MDS-post-U2AF1-RAS have 2 BAMs with the same uid and library).
The selection rules are the ones that produced the reported counts:
  CBFB::MYH11  gene_pairs holds the CBFB/MYH11 pair (combine.py); exon at the CBFB segment 3' end (max coordinate) and
               the MYH11 segment 5' end (max coordinate), MANE exons, +-15 nt; type A = CBFB e5 > MYH11 e33.
  KMT2A-PTD    combine.py region rule (KMT2A noncolinear, >= 2 KMT2A segments, all inside MANE exons 2-11), then the
               first adjacent KMT2A segment join with donor exon >= 8 and acceptor exon 2 or 3 (MANE, +-15 nt), or
               acceptor alt-e2: the GENCODE v45 exon chr11:118,439,024-118,439,122 in intron 1 (exon 2 of KMT2A-206,
               -209, -220, -234; scripts/GENCODE_v45_KMT2A_all_exons.gff3), +-15 nt.
v2 (2026-11-10) adds the alt-e2 acceptor. v1 used MANE exons only and missed 107 patient reads joined e9/e10 > alt-e2
(delete/export_lesion_read_ids_v1_superseded/).
  FLT3-ITD     combine.py insertion grouping (+-10 nt position, +-6 nt length) for pAML-rel-IDH2; the 21-nt and 102-nt
               groups named in the table.
Reads of the same lesion in other samples are written with row_role other_sample (the table cites them as background).
Usage: export_lesion_read_ids_v2.py <scan_dir> <metadata> <MANE exon gff3 subset> <genes.tsv.gz> <out.tsv>
"""
import collections, csv, glob, gzip, os, re, sys

O, META, EXONS, GENES, OUT = sys.argv[1:6]
assert not os.path.exists(OUT), f'{OUT} exists; refusing to overwrite'
csv.field_size_limit(sys.maxsize)

meta = list(csv.DictReader(open(META), delimiter='\t'))
GSTRAND = {}
with gzip.open(GENES, 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        GSTRAND.setdefault(r['gene_name'], r['strand'])
EX = collections.defaultdict(dict)
for line in open(EXONS):
    m = re.search(r'gene_name=(KMT2A|CBFB|MYH11|FLT3);', line)
    if m and '\texon\t' in line and 'MANE_Select' in line:
        c = line.split('\t'); EX[m.group(1)][int(re.search(r'exon_number=(\d+)', c[8]).group(1))] = (int(c[3]), int(c[4]))
assert all(EX[g] for g in ('KMT2A', 'CBFB', 'MYH11', 'FLT3')), {g: len(EX[g]) for g in EX}
print('MANE exons:', {g: len(EX[g]) for g in sorted(EX)})


def exon_of(g, pos, tol=15):
    for n, (a, b) in EX[g].items():
        if a - tol <= pos <= b + tol:
            return n
    return None


ALT_E2 = 118439024   # start of GENCODE v45 exon chr11:118439024-118439122 (exon 2 of ENST00000527869, ENST00000531904)
SEG = re.compile(r'(\S+?):(\d+)-(\d+)\(([+-])\) q(\d+)-(\d+) mq(\d+) \[([^\]]*)\]')


def parse_segments(s):
    return [dict(chrom=m[0], start=int(m[1]), end=int(m[2]), strand=m[3], q=int(m[4]), genes=m[7].split(',')) for m in SEG.findall(s)]


def transcript_order(segs, a, b):   # copied from combine.py
    ia = next((i for i, s in enumerate(segs) if a in s['genes']), None)
    ib = next((i for i, s in enumerate(segs) if b in s['genes'] and i != ia), None)
    if ia is None or ib is None or a not in GSTRAND or b not in GSTRAND:
        return 'same_alignment_or_unresolved'
    sa, sb = segs[ia]['strand'] == GSTRAND[a], segs[ib]['strand'] == GSTRAND[b]
    if sa != sb:
        return 'mixed_strand'
    first_a = (ia < ib) == sa
    return f'{a}>{b}' if first_a else f'{b}>{a}'


TABLE = {'CBFB::MYH11': {'AML14__AML-14_CITE_GEX'},
         'KMT2A-PTD': {'MDS-pre-U2AF1-RAS__Y2982', 'MDS-post-U2AF1-RAS__Y3554', 'MDS-pre-U2AF1-RAS-MISTRG__Y2982-E3'},
         'FLT3-ITD': {'pAML-rel-IDH2__5591-3-LAM'}}
KX = EX['KMT2A']
K_LO, K_HI = min(KX[n][0] for n in range(2, 12)), max(KX[n][1] for n in range(2, 12))

out, bam_of_uid = [], collections.defaultdict(list)
cbfb_check, kmt2a_region_check = collections.Counter(), collections.Counter()
for folder in sorted(glob.glob(os.path.join(O, 'samples', '[0-9][0-9]_*'))):
    nn = int(os.path.basename(folder).split('_', 1)[0])
    m = meta[nn - 1]
    uid = f"{m['uid']}__{m['library']}"
    assert os.path.basename(folder) == f'{nn:02d}_{uid}', (folder, uid)
    bam_of_uid[uid].append(nn)
    base = dict(sample=m['uid'], library=m['library'], group=m['groups'], bam_index=nn, bam=m['bam'])
    for r in csv.DictReader(open(os.path.join(folder, f'{uid}.fusion_reads.tsv')), delimiter='\t'):
        assert r['uid'] == uid
        segs = parse_segments(r['segments'])
        # CBFB::MYH11
        hit = None
        for x in filter(None, r['gene_pairs'].split(';')):
            ab, kind = x.rsplit(':', 1); a, b = ab.split('>', 1)
            if tuple(sorted((a, b))) == ('CBFB', 'MYH11'):
                hit = (a, b, kind); break
        if hit:
            cbfb_check[(uid, r['read'])] += 1
            sm = SEG.findall(r['segments'])
            c = [x for x in sm if 'CBFB' in x[7].split(',')]; my = [x for x in sm if 'MYH11' in x[7].split(',')]
            ec = exon_of('CBFB', max(int(x[2]) for x in c)) if c and my else None
            em = exon_of('MYH11', max(int(x[2]) for x in my)) if c and my else None
            out.append(dict(lesion='CBFB::MYH11', row_role='lesion_sample' if uid in TABLE['CBFB::MYH11'] else 'other_sample', **base,
                            read_id=r['read'], event=f'CBFB e{ec} > MYH11 e{em}' if c and my else 'both genes in one segment',
                            type_A=(ec, em) == (5, 33), transcript_order=transcript_order(segs, hit[0], hit[1]), read_kind=hit[2],
                            insertion_length='', insertion_pos='', in_frame='', tandem_8mer_fraction='', insertion_seq='',
                            segments=r['segments']))
        # KMT2A-PTD (combine.py region rule, then the PTD-type join)
        if 'KMT2A' in r['same_gene_noncolinear'].split(';'):
            sg = [s for s in segs if 'KMT2A' in s['genes']]
            if len(sg) >= 2 and all(s['chrom'] == 'chr11' and s['start'] <= K_HI and K_LO <= s['end'] for s in sg):
                kmt2a_region_check[(uid, r['read'])] += 1
                ks = [x for x in SEG.findall(r['segments']) if 'KMT2A' in x[7].split(',')]
                joins, first = [], None
                for a, b in zip(ks, ks[1:]):
                    dpos, apos = (int(a[2]), int(b[1])) if a[3] == '+' else (int(b[2]), int(a[1]))
                    j = (exon_of('KMT2A', dpos), exon_of('KMT2A', apos))
                    alt = j[1] is None and abs(apos - ALT_E2) <= 15
                    acc = 'alt-e2' if alt else f'e{j[1]}'
                    joins.append(f'e{j[0]}>{acc}')
                    if first is None and j[0] and j[0] >= 8 and (j[1] in (2, 3) or alt):
                        first = (j[0], acc)
                if first:
                    out.append(dict(lesion='KMT2A-PTD', row_role='lesion_sample' if uid in TABLE['KMT2A-PTD'] else 'other_sample', **base,
                                    read_id=r['read'], event=f'KMT2A e{first[0]} > {first[1]}', type_A='',
                                    transcript_order='all KMT2A joins: ' + ', '.join(joins), read_kind='',
                                    insertion_length='', insertion_pos='', in_frame='', tandem_8mer_fraction='', insertion_seq='',
                                    segments=r['segments']))

# validation against the combined tables written by combine.py
doc = collections.Counter((x['uid'], x['read']) for x in csv.DictReader(open(os.path.join(O, 'AML_MDS_FUSION_DOCUMENTED_READS.tsv')), delimiter='\t')
                          if x['pair'] == 'CBFB::MYH11')
assert doc == cbfb_check, ('CBFB::MYH11 read set differs from AML_MDS_FUSION_DOCUMENTED_READS.tsv', sum(doc.values()), sum(cbfb_check.values()))
reg = collections.Counter((x['uid'], x['read']) for x in csv.DictReader(open(os.path.join(O, 'KMT2A_PTD_and_FLT3_ITD_region_noncolinear_reads.tsv')), delimiter='\t')
                          if x['gene'] == 'KMT2A')
assert reg == kmt2a_region_check, ('KMT2A region read set differs', sum(reg.values()), sum(kmt2a_region_check.values()))
print(f'validation: CBFB::MYH11 {sum(cbfb_check.values())} of {sum(doc.values())} reads match AML_MDS_FUSION_DOCUMENTED_READS.tsv; '
      f'KMT2A region {sum(kmt2a_region_check.values())} of {sum(reg.values())} match KMT2A_PTD_and_FLT3_ITD_region_noncolinear_reads.tsv')

# FLT3-ITD: combine.py grouping on the events table (tandem fraction already computed there)
U = 'pAML-rel-IDH2__5591-3-LAM'
assert len(bam_of_uid[U]) == 1, bam_of_uid[U]
nn = bam_of_uid[U][0]; m = meta[nn - 1]
ev = [e for e in csv.DictReader(open(os.path.join(O, 'FLT3_exon14_15_events_all_reads.tsv')), delimiter='\t') if e['uid'] == U]
own = list(csv.DictReader(open(glob.glob(os.path.join(O, 'samples', f'{nn:02d}_{U}', f'{U}.flt3_events.tsv'))[0]), delimiter='\t'))
assert [(e['read'], e['kind'], e['length'], e['ref_pos']) for e in ev] == [(e['read'], e['kind'], e['length'], e['ref_pos']) for e in own]
for e in ev:
    e['length'] = int(e['length']); e['ref_pos'] = int(e['ref_pos'])
ins = sorted((e for e in ev if e['kind'] == 'insertion'), key=lambda e: (e['ref_pos'], e['length']))
used = [False] * len(ins); groups = []
for i, e in enumerate(ins):
    if used[i]:
        continue
    grp = [e]; used[i] = True
    for j in range(i + 1, len(ins)):
        if not used[j] and abs(ins[j]['ref_pos'] - e['ref_pos']) <= 10 and abs(ins[j]['length'] - e['length']) <= 6:
            grp.append(ins[j]); used[j] = True
    groups.append((collections.Counter(x['length'] for x in grp).most_common(1)[0][0], grp))
print('FLT3 insertion groups in', U, [(L, len({x['read'] for x in g})) for L, g in groups])
for L, grp in groups:
    if L not in (21, 102):
        continue
    for x in grp:
        out.append(dict(lesion='FLT3-ITD', row_role='lesion_sample', sample=m['uid'], library=m['library'], group=m['groups'], bam_index=nn,
                        bam=m['bam'], read_id=x['read'], event=f'{L}-nt insertion group', type_A='', transcript_order='', read_kind='',
                        insertion_length=x['length'], insertion_pos=f"chr13:{x['ref_pos']}", in_frame=x['in_frame'],
                        tandem_8mer_fraction=x['tandem_8mer_fraction'], insertion_seq=x['seq'], segments=''))

with open(OUT, 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)

# counts for the reply and README
cnt = collections.Counter((x['lesion'], x['row_role'], x['sample'], x['library'], x['bam_index']) for x in out)
print('\nrows per lesion, role, BAM:')
for k, n in sorted(cnt.items()):
    print('  ', '\t'.join(map(str, k)), n)
for les in ('CBFB::MYH11', 'KMT2A-PTD', 'FLT3-ITD'):
    rows = [x for x in out if x['lesion'] == les]
    print(f'\n{les}: {len(rows)} rows, lesion_sample {sum(x["row_role"] == "lesion_sample" for x in rows)}, '
          f'unique (bam, read) {len({(x["bam_index"], x["read_id"]) for x in rows})}')
    lr = [x for x in rows if x['row_role'] == 'lesion_sample']
    print('  events in lesion samples:', collections.Counter((x['sample'], x['event']) for x in lr).most_common(8))
    if les == 'CBFB::MYH11':
        print('  type A in AML14:', sum(x['type_A'] is True for x in lr), '| CBFB>MYH11 order all samples:',
              sum(x['transcript_order'] == 'CBFB>MYH11' for x in rows), 'of', len(rows))
        by = collections.defaultdict(set)
        for x in lr:
            by[x['bam_index']].add(x['segments'])
        ids = collections.defaultdict(set)
        for x in lr:
            ids[x['bam_index']].add(x['read_id'])
        bs = sorted(by)
        if len(bs) == 2:
            print(f'  AML14 BAMs {bs}: identical segment strings shared {len(by[bs[0]] & by[bs[1]])}; identical read IDs shared {len(ids[bs[0]] & ids[bs[1]])}')
print('\nwrote', OUT, len(out), 'rows')
