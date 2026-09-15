"""Create a tiny indexed BAM/GTF fixture; no human data or reference downloads."""
import argparse
from pathlib import Path
import pysam


def create(outdir):
    p = Path(outdir).resolve(); p.mkdir(parents=True, exist_ok=True)
    gtf = p / 'reference.gtf'
    gtf.write_text('chr1\ttest\tgene\t101\t400\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "TEST";\n'
        'chr1\ttest\ttranscript\t101\t400\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";\n'
        'chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_id "ENSE1"; exon_number "1";\n'
        'chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_id "ENSE2"; exon_number "2";\n')
    for sample, n in [('S1', 20), ('S2', 19)]:
        bam = p / (sample + '.bam')
        with pysam.AlignmentFile(str(bam), 'wb', header={'HD': {'VN': '1.6','SO':'coordinate'}, 'SQ':[{'SN':'chr1','LN':1000}]}) as f:
            for i in range(n):
                read = pysam.AlignedSegment(); read.query_name=f'{sample}_{i}'
                read.query_sequence='A'*100; read.flag=0; read.reference_id=0
                read.reference_start=150; read.mapping_quality=60; read.cigar=[(0,50),(3,100),(0,50)]
                read.query_qualities=pysam.qualitystring_to_array('I'*100)
                read.set_tag('XS','+'); f.write(read)
            if sample == 'S1':
                for start in (185, 285):  # evidence at both intron boundaries
                    for i in range(3):
                        read = pysam.AlignedSegment(); read.query_name=f'{sample}_retained_{start}_{i}'
                        read.query_sequence='A'*40; read.flag=0; read.reference_id=0
                        read.reference_start=start; read.mapping_quality=60; read.cigar=[(0,40)]
                        read.query_qualities=pysam.qualitystring_to_array('I'*40)
                        read.set_tag('XS','+'); f.write(read)
        pysam.index(str(bam))
    (p / 'samples.csv').write_text('id,bam,bai,psm\n'+'\n'.join(f'{s},{p}/{s}.bam,{p}/{s}.bam.bai,' for s in ['S1','S2'])+'\n')
    (p / 'hla.tsv').write_text('sample\thla\nS1\tHLA-A*02:01,HLA-B*07:02,HLA-C*07:02\nS2\tHLA-A*02:01\n')
    (p / 'snaf_db').mkdir(exist_ok=True)
    return p

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--outdir', required=True)
    create(parser.parse_args().outdir)
