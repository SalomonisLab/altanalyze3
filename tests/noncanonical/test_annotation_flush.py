"""Compressed raw annotations must retain every annotated source transcript."""
import gzip
import pytest
from altanalyze3.components.long_read import gff_process,io_utils

@pytest.mark.parametrize('n_transcripts',[2,5000])
def test_raw_and_native_annotation_rosters_survive_compression(tmp_path,monkeypatch,n_transcripts):
    monkeypatch.setenv('ISOFORM_STRUCT_COORDS','0')
    compress=io_utils.compress
    monkeypatch.setattr(io_utils,'compress',lambda path,**kw:compress(path,min_bytes=0,**kw))
    exon=tmp_path/'exons.tsv'
    exon.write_text('gene\texon\tchr\tstrand\tstart\tstop\nG\tE1.1\tchr1\t+\t100\t149\nG\tE2.1\tchr1\t+\t300\t349\n')
    gff=tmp_path/'source.gff'
    with gff.open('w') as out:
        for i in range(n_transcripts):
            attr=f'gene_id "G";transcript_id "read{i}";'
            out.write(f'chr1\tbam\ttranscript\t100\t349\t.\t+\t.\t{attr}\nchr1\tbam\texon\t100\t149\t.\t+\t.\t{attr}\nchr1\tbam\texon\t300\t349\t.\t+\t.\t{attr}\n')
    gff_process.consolidateLongReadGFFs(str(gff),str(exon),mode='collapse')
    rosters=[]
    for name in ('transcript_associations.txt.gz','transcript_associations_raw.txt.gz'):
        with gzip.open(tmp_path/'gff-output'/name,'rt') as handle:
            rows=[line.rstrip('\n').split('\t') for line in handle]
        assert len(rows)==n_transcripts
        rosters.append([(r[0],r[1],r[3],r[4]) for r in rows])
    assert rosters[0]==rosters[1]
    assert [r[2] for r in rosters[0]]==[f'read{i}' for i in range(n_transcripts)]
