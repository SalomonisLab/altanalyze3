"""KMT2A exon-typed PTD/circle evidence from the validated v2/v3 lesion rules.

The archived scripts also contain cohort-specific sample labels and fixed FLT3
allele exports. This callable adapter applies their KMT2A geometry to every source
without importing those cohort restrictions. It never establishes RNA topology.
"""
import collections
from pathlib import Path
import re

SEG=re.compile(r'(\S+?):(\d+)-(\d+)\(([+-])\) q(\d+)-(\d+) mq(\d+) \[([^\]]*)\]')


def mane_exons(path=None):
    path=path or Path(__file__).parent/'oncofusion/MANE_exons_CBFB_MYH11_KMT2A_FLT3_gencode_v45.gff3'
    ex={}
    for line in Path(path).open():
        if 'gene_name=KMT2A;' in line and '\texon\t' in line and 'MANE_Select' in line:
            c=line.split('\t');ex[int(re.search(r'exon_number=(\d+)',c[8]).group(1))]=(int(c[3]),int(c[4]))
    if sorted(ex)!=list(range(1,37)):raise ValueError('Incomplete validated KMT2A MANE exon profile')
    return ex


def classify(row,exons=None):
    if 'KMT2A' not in row['same_gene_noncolinear'].split(';'):return None
    ex=exons or mane_exons()
    segments=[m for m in SEG.findall(row['segments']) if 'KMT2A' in m[7].split(',')]
    lo=min(ex[n][0] for n in range(2,12));hi=max(ex[n][1] for n in range(2,12))
    if len(segments)<2 or not all(s[0]=='chr11' and int(s[1])<=hi and lo<=int(s[2]) for s in segments):return None
    def exon_of(pos):
        return next((n for n,(a,b) in ex.items() if a-15<=pos<=b+15),None)
    typed=[];first=None;back=[]
    for a,b in zip(segments,segments[1:]):
        dpos,apos=(int(a[2]),int(b[1])) if a[3]=='+' else (int(b[2]),int(a[1]))
        donor,acceptor=exon_of(dpos),exon_of(apos)
        alt=acceptor is None and abs(apos-118439024)<=15
        accepted=bool(donor and donor>=8 and (acceptor in (2,3) or alt))
        typed.append(dict(chrom=a[0],strand=a[3],donor=dpos,acceptor=apos,donor_exon=donor,acceptor_exon='alt-e2' if alt else acceptor,ptd_type=accepted))
        if first is None and accepted:first=typed[-1]
        if a[0]==b[0] and a[3]==b[3]:
            if a[3]=='+' and int(b[1])<int(a[2]):back.append((int(a[2]),int(b[1])))
            elif a[3]=='-' and int(b[2])>int(a[1]):back.append((int(b[2]),int(a[1])))
    if first is None:return None
    if not back:raise ValueError('Native PTD type without a same-strand backward join')
    blo,bhi=min(a for d,a in back),max(d for d,a in back)
    upstream=max(0,blo-min(int(s[1]) for s in segments));downstream=max(0,max(int(s[2]) for s in segments)-bhi)
    label=('linear: sequence outside the joined block' if upstream>15 or downstream>15 else
           'circle-compatible: 2 or more backward joins' if len(back)>=2 else
           'uninformative: 1 join, all sequence inside the block')
    return dict(first_ptd_type_join=first,all_typed_joins=typed,backward_joins=back,
                native_read_class=label,circle_spliced_length_nt=sum(max(0,min(b,bhi)-max(a,blo)+1) for a,b in ex.values()),
                aligned_kmt2a_query_nt=max(int(s[5]) for s in segments)-min(int(s[4]) for s in segments),
                upstream_outside_nt=upstream,downstream_outside_nt=downstream,
                topology='unresolved; multiple different backward joins are not proof of repeated traversal of one circle')
