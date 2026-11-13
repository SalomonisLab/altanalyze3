"""Combine the per-BAM AML/MDS fusion and FLT3-ITD scans (2026-11-10). Ad hoc; validate within altanalyze3.

Fusion pairs are annotated with (1) the CTAT human fusion library v0.3.0 table of the workflow (reference/fusion_db.tsv:
tags; curated_cancer = Mitelman, ChimerKB, ChimerPub, COSMIC, chimerdb_omim, chimerdb_pubmed, HaasMedCancer), either gene
order, and (2) the WHO 2022 / ICC 2022 AML and MDS/MPN-defining fusions and recurrent pediatric AML fusions (list below).
No splice-site or read-count requirement. FLT3: insertions >= 3 nt in exons 14-15 are grouped per sample by position
(+-10 nt) and length (+-6 nt); an insertion is tandem when >= 80% of its 8-mers occur in the reference within 300 nt.
"""
import collections, csv, glob, gzip, json, os, re, sys

O = sys.argv[1]
REF = sys.argv[2]           # extra-chr-trans-splicing/reference
GENOME = sys.argv[3]        # genome.fa (Ensembl contig names)
csv.field_size_limit(sys.maxsize)
WHO = {('PML', 'RARA'), ('RUNX1', 'RUNX1T1'), ('CBFB', 'MYH11'), ('BCR', 'ABL1'), ('DEK', 'NUP214'), ('RBM15', 'MRTFA'),
       ('KAT6A', 'CREBBP'), ('KAT6A', 'EP300'), ('KAT6A', 'NCOA2'), ('KAT6B', 'CREBBP'), ('FUS', 'ERG'), ('CBFA2T3', 'GLIS2'),
       ('NPM1', 'MLF1'), ('PICALM', 'MLLT10'), ('MNX1', 'ETV6'), ('NUP98', 'NSD1'), ('NUP98', 'KDM5A'), ('NUP98', 'DDX10'),
       ('NUP98', 'TOP1'), ('NUP98', 'HOXA9'), ('NUP98', 'NSD3'), ('RUNX1', 'MECOM'), ('ETV6', 'MECOM'), ('SET', 'NUP214'),
       ('ZBTB16', 'RARA'), ('NPM1', 'RARA'), ('NUMA1', 'RARA'), ('STAT5B', 'RARA'), ('PRKAR1A', 'RARA'), ('FIP1L1', 'RARA'),
       ('BCOR', 'RARA'), ('PCM1', 'JAK2'), ('ETV6', 'JAK2'), ('BCR', 'JAK2'), ('FIP1L1', 'PDGFRA'), ('ETV6', 'PDGFRB'),
       ('ZMYM2', 'FGFR1'), ('BCR', 'FGFR1'), ('MYB', 'GATA1'), ('ETV6', 'ABL1'), ('CBFA2T3', 'RUNX1')}
# the first gene of each WHO pair is the 5' partner of the fusion transcript
KMT2A_PARTNERS_ANY = True   # every KMT2A partner is a KMT2A rearrangement (WHO/ICC)

GROUP = {}
for r in csv.DictReader(open(sys.argv[4]), delimiter='\t'):   # mds_metadata_bam_fixed.txt
    GROUP[f"{r['uid']}__{r['library']}" if r['library'] else r['uid']] = r['groups']
HEALTHY = {'young', 'aged'}
GSTRAND = {}
with gzip.open(os.path.join(REF, 'genes.tsv.gz'), 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        GSTRAND.setdefault(r['gene_name'], r['strand'])
SEG = re.compile(r'(\S+?):(\d+)-(\d+)\(([+-])\) q(\d+)-(\d+) mq(\d+) \[([^\]]*)\]')


def parse_segments(s):
    return [dict(chrom=m[0], start=int(m[1]), end=int(m[2]), strand=m[3], q=int(m[4]), genes=m[7].split(',')) for m in SEG.findall(s)]


def transcript_order(segs, a, b):
    """5'->3' gene order of the read for genes a and b, from each segment's strand relative to its gene's strand."""
    ia = next((i for i, s in enumerate(segs) if a in s['genes']), None)
    ib = next((i for i, s in enumerate(segs) if b in s['genes'] and i != ia), None)
    if ia is None or ib is None or a not in GSTRAND or b not in GSTRAND:
        return 'same_alignment_or_unresolved'
    sa, sb = segs[ia]['strand'] == GSTRAND[a], segs[ib]['strand'] == GSTRAND[b]
    if sa != sb:
        return 'mixed_strand'
    first_a = (ia < ib) == sa
    return f'{a}>{b}' if first_a else f'{b}>{a}'


db = {}
with open(os.path.join(REF, 'fusion_db.tsv')) as h:
    for r in csv.DictReader(h, delimiter='\t'):
        db[(r['gene_a'], r['gene_b'])] = r


def annotate(a, b):
    tags, cur = set(), False
    for k in ((a, b), (b, a)):
        if k in db:
            tags |= set(db[k]['tags'].split(',')); cur |= db[k]['curated_cancer'] == 'True'
    who = (a, b) in WHO or (b, a) in WHO
    kmt2a = 'KMT2A' in (a, b)
    return ','.join(sorted(t for t in tags if t)), cur, who, kmt2a


samples = sorted(glob.glob(os.path.join(O, 'samples', '*', '*.summary.json')))
summ = [dict(json.load(open(p)), bam_folder=os.path.basename(os.path.dirname(p))) for p in samples]
print('BAMs with a summary', len(summ), '| samples (uid__library)', len({s['uid'] for s in summ}))
reads = []
for p in sorted(glob.glob(os.path.join(O, 'samples', '*', '*.fusion_reads.tsv'))):
    reads += list(csv.DictReader(open(p), delimiter='\t'))
pair_reads = collections.defaultdict(list)
for r in reads:
    seen = set()
    for x in filter(None, r['gene_pairs'].split(';')):
        ab, kind = x.rsplit(':', 1)
        a, b = ab.split('>', 1)
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        pair_reads[key].append((r['uid'], r['read'], transcript_order(parse_segments(r['segments']), a, b), kind))
rows = []
for (a, b), lst in pair_reads.items():
    tags, cur, who, kmt2a = annotate(a, b)
    uids = collections.Counter(u for u, _, _, _ in lst)
    order = collections.Counter(o for _, _, o, _ in lst)
    kinds = collections.Counter(k for _, _, _, k in lst)
    canon = next((f'{x}>{y}' for x, y in WHO if {x, y} == {a, b}), '')
    rows.append(dict(gene_1=a, gene_2=b, reads=len(lst), samples=len(uids),
                     reads_healthy=sum(n for u, n in uids.items() if GROUP.get(u) in HEALTHY),
                     reads_disease=sum(n for u, n in uids.items() if GROUP.get(u) not in HEALTHY),
                     canonical_5p_3p=canon, reads_canonical_order=order.get(canon, 0) if canon else '',
                     per_sample='; '.join(f'{u}:{n}' for u, n in uids.most_common()),
                     transcript_order='; '.join(f'{o}:{n}' for o, n in order.most_common()),
                     read_kind='; '.join(f'{k}:{n}' for k, n in kinds.most_common()),
                     WHO_ICC_AML_MDS_fusion=who, KMT2A_rearrangement=kmt2a, CTAT_curated_cancer=cur, CTAT_tags=tags))
rows.sort(key=lambda r: (not (r['WHO_ICC_AML_MDS_fusion'] or r['KMT2A_rearrangement'] or r['CTAT_curated_cancer']), -r['reads']))
with open(os.path.join(O, 'AML_MDS_FUSION_GENE_PAIRS_all_reads.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['gene_1', 'gene_2', 'reads', 'samples', 'reads_healthy', 'reads_disease', 'canonical_5p_3p', 'reads_canonical_order', 'per_sample', 'transcript_order', 'read_kind', 'WHO_ICC_AML_MDS_fusion', 'KMT2A_rearrangement', 'CTAT_curated_cancer', 'CTAT_tags'], delimiter='\t'); w.writeheader(); w.writerows(rows)
doc = [r for r in rows if r['WHO_ICC_AML_MDS_fusion'] or r['KMT2A_rearrangement'] or r['CTAT_curated_cancer']]
with open(os.path.join(O, 'AML_MDS_FUSION_DOCUMENTED_PAIRS.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['gene_1', 'gene_2', 'reads', 'samples', 'reads_healthy', 'reads_disease', 'canonical_5p_3p', 'reads_canonical_order', 'per_sample', 'transcript_order', 'read_kind', 'WHO_ICC_AML_MDS_fusion', 'KMT2A_rearrangement', 'CTAT_curated_cancer', 'CTAT_tags'], delimiter='\t'); w.writeheader(); w.writerows(doc)
# reads of documented pairs, for inspection
doc_keys = {(r['gene_1'], r['gene_2']) for r in doc}
with open(os.path.join(O, 'AML_MDS_FUSION_DOCUMENTED_READS.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['pair','uid','read','n_segments','segments','gene_pairs','same_gene_noncolinear'], delimiter='\t'); w.writeheader()
    for r in reads:
        ks = {tuple(sorted(x.rsplit(':', 1)[0].split('>', 1))) for x in filter(None, r['gene_pairs'].split(';'))} & doc_keys
        for k in sorted(ks):
            w.writerow(dict(pair='::'.join(k), **r))
# same-gene non-colinear (KMT2A partial tandem duplication, FLT3 large ITD, others)
def mane_exons(gene):
    ex = {}
    with open(sys.argv[5]) as h:   # GENCODE v45 GFF3
        for line in h:
            if f'gene_name={gene};' in line and '\texon\t' in line and 'MANE_Select' in line:
                c = line.split('\t'); n = int(re.search(r'exon_number=(\d+)', c[8]).group(1)); ex[n] = (int(c[3]), int(c[4]))
    return ex
KX, FX = mane_exons('KMT2A'), mane_exons('FLT3')
assert KX and FX
REGION = {'KMT2A': ('chr11', min(KX[n][0] for n in range(2, 12)), max(KX[n][1] for n in range(2, 12)), 'KMT2A exons 2-11 (PTD region)'),
          'FLT3': ('chr13', min(FX[n][0] for n in range(13, 17)), max(FX[n][1] for n in range(13, 17)), 'FLT3 exons 13-16 (ITD region)')}
ncl = collections.defaultdict(collections.Counter)
ncl_region = collections.defaultdict(collections.Counter)
region_reads = []
for r in reads:
    for g in filter(None, r['same_gene_noncolinear'].split(';')):
        ncl[g][r['uid']] += 1
        if g in REGION:
            ch, lo, hi, lab = REGION[g]
            sg = [s for s in parse_segments(r['segments']) if g in s['genes']]
            if len(sg) >= 2 and all(s['chrom'] == ch and s['start'] <= hi and lo <= s['end'] for s in sg):
                ncl_region[g][r['uid']] += 1
                region_reads.append(dict(gene=g, region=lab, uid=r['uid'], group=GROUP.get(r['uid'], ''), read=r['read'], segments=r['segments']))
with open(os.path.join(O, 'KMT2A_PTD_and_FLT3_ITD_region_noncolinear_reads.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['gene', 'region', 'uid', 'group', 'read', 'segments'], delimiter='\t'); w.writeheader(); w.writerows(region_reads)
with open(os.path.join(O, 'SAME_GENE_NONCOLINEAR_reads.tsv'), 'w', newline='') as h:
    w = csv.writer(h, delimiter='\t'); w.writerow(['gene', 'reads_anywhere_in_gene', 'samples', 'reads_in_PTD_or_ITD_region', 'per_sample_in_region', 'per_sample_anywhere'])
    for g, c in sorted(ncl.items(), key=lambda x: -sum(x[1].values())):
        cr = ncl_region.get(g, collections.Counter())
        w.writerow([g, sum(c.values()), len(c), sum(cr.values()) if g in REGION else '', '; '.join(f'{u}:{n}' for u, n in cr.most_common()),
                    '; '.join(f'{u}:{n}' for u, n in c.most_common())])

# ---------------------------------------------------------------- FLT3
fai = {}
for line in open(GENOME + '.fai'):
    p = line.split('\t'); fai[p[0]] = tuple(int(x) for x in p[1:5])
G = open(GENOME, 'rb')


def fetch(ch, s, e):
    ln, off, lb, lby = fai[ch]
    a, b = s - 1, e - 1
    G.seek(off + (a // lb) * lby + a % lb)
    raw = G.read((off + (b // lb) * lby + b % lb) - (off + (a // lb) * lby + a % lb) + 1)
    return raw.replace(b'\n', b'').decode().upper()


ref = fetch('13', 28033887 - 400, 28034214 + 400)
ref0 = 28033887 - 400


def tandem_fraction(seq, pos):
    if len(seq) < 8:
        return 1.0 if seq and seq in ref[max(0, pos - ref0 - 300):pos - ref0 + 300] else 0.0
    win = ref[max(0, pos - ref0 - 300):pos - ref0 + 300]
    ks = [seq[i:i + 8] for i in range(len(seq) - 7)]
    return sum(k in win for k in ks) / len(ks)


ev = []
for p in sorted(glob.glob(os.path.join(O, 'samples', '*', '*.flt3_events.tsv'))):
    ev += list(csv.DictReader(open(p), delimiter='\t'))
for e in ev:
    e['length'] = int(e['length']); e['ref_pos'] = int(e['ref_pos'])
    e['tandem_8mer_fraction'] = round(tandem_fraction(e['seq'], e['ref_pos']), 3) if e['kind'] == 'insertion' else ''
cov = collections.Counter()
for s in summ:
    cov[s['uid']] += s['flt3_reads_covering_exons14_15']  # a sample with 2 BAMs (e.g. 251H first_run + size_selected) sums both
clusters = []
for uid in sorted({e['uid'] for e in ev}):
    ins = sorted((e for e in ev if e['uid'] == uid and e['kind'] == 'insertion'), key=lambda e: (e['ref_pos'], e['length']))
    used = [False] * len(ins)
    for i, e in enumerate(ins):
        if used[i]:
            continue
        grp = [e]; used[i] = True
        for j in range(i + 1, len(ins)):
            if not used[j] and abs(ins[j]['ref_pos'] - e['ref_pos']) <= 10 and abs(ins[j]['length'] - e['length']) <= 6:
                grp.append(ins[j]); used[j] = True
        L = collections.Counter(x['length'] for x in grp).most_common(1)[0][0]
        clusters.append(dict(uid=uid, ref_pos=min(x['ref_pos'] for x in grp), length_mode=L, reads=len({x['read'] for x in grp}),
                             reads_covering=cov.get(uid, ''), fraction=round(len({x['read'] for x in grp}) / cov[uid], 4) if cov.get(uid) else '',
                             in_frame=(L % 3 == 0), tandem_8mer_fraction_median=sorted(x['tandem_8mer_fraction'] for x in grp)[len(grp) // 2],
                             example_seq=max(grp, key=lambda x: x['length'])['seq'][:120]))
clusters.sort(key=lambda c: -c['reads'])
with open(os.path.join(O, 'FLT3_exon14_15_insertions_by_sample.tsv'), 'w', newline='') as h:
    if clusters:
        w = csv.DictWriter(h, fieldnames=list(clusters[0]), delimiter='\t'); w.writeheader(); w.writerows(clusters)
with open(os.path.join(O, 'FLT3_exon14_15_events_all_reads.tsv'), 'w', newline='') as h:
    if ev:
        w = csv.DictWriter(h, fieldnames=list(ev[0]), delimiter='\t'); w.writeheader(); w.writerows(ev)
per = []
for u in sorted({s['uid'] for s in summ}):
    s = dict(flt3_reads_covering_exons14_15=cov[u], bams=sum(1 for x in summ if x['uid'] == u))
    cs = [c for c in clusters if c['uid'] == u]
    itd = [c for c in cs if c['in_frame'] and c['tandem_8mer_fraction_median'] >= 0.8 and c['length_mode'] >= 3]
    sc = [e for e in ev if e['uid'] == u and e['kind'].startswith('softclip')]
    per.append(dict(uid=u, bams=s['bams'], reads_covering_exons14_15=s['flt3_reads_covering_exons14_15'],
                    reads_with_insertion_ge3=len({e['read'] for e in ev if e['uid'] == u and e['kind'] == 'insertion'}),
                    in_frame_tandem_insertion_groups=len(itd),
                    itd_reads=sum(c['reads'] for c in itd),
                    itd_lengths='; '.join(f"{c['length_mode']} nt x{c['reads']}" for c in itd),
                    softclip_ge20_reads=len({e['read'] for e in sc}),
                    flt3_noncolinear_reads_exons13_16=ncl_region.get('FLT3', {}).get(u, 0), group=GROUP.get(u, '')))
with open(os.path.join(O, 'FLT3_ITD_per_sample.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=list(per[0]), delimiter='\t'); w.writeheader(); w.writerows(per)
S = dict(bams=len(summ), samples=len({s['uid'] for s in summ}), reads_with_aml_gene_pair_or_noncolinear=len(reads), gene_pairs=len(rows), documented_pairs=len(doc),
         documented_reads=sum(r['reads'] for r in doc), who_icc_pairs=[f"{r['gene_1']}::{r['gene_2']} {r['reads']} reads" for r in doc if r['WHO_ICC_AML_MDS_fusion']],
         kmt2a_partners=[f"{r['gene_1']}::{r['gene_2']} {r['reads']} reads" for r in doc if r['KMT2A_rearrangement']],
         kmt2a_ptd_region_noncolinear_reads={u: n for u, n in ncl_region.get('KMT2A', {}).items()},
         flt3_itd_region_noncolinear_reads={u: n for u, n in ncl_region.get('FLT3', {}).items()},
         flt3_samples_with_in_frame_tandem_insertion=sum(1 for p in per if p['in_frame_tandem_insertion_groups']))
json.dump(S, open(os.path.join(O, 'SUMMARY.json'), 'w'), indent=1)
print(json.dumps(S, indent=1))
