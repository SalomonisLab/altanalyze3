from altanalyze3.components.long_read import gff_process


def test_final_cross_chromosome_transcript_matches_nonfinal(tmp_path,monkeypatch):
    exon=tmp_path/'exons.tsv'
    exon.write_text('gene\texon\tchr\tstrand\tstart\tstop\nG\tE1.1\tchr1\t+\t100\t149\nG\tE2.1\tchr1\t+\t300\t349\nH\tE1.1\tchr2\t+\t200\t249\nH\tE2.1\tchr2\t+\t400\t449\n')
    gff=tmp_path/'cross.gff'
    def block(name):
        attr=f'gene_id "G";transcript_id "{name}";'
        return (f'chr1\tbam\ttranscript\t100\t149\t.\t+\t.\t{attr}\n'
                f'chr1\tbam\texon\t100\t149\t.\t+\t.\t{attr}segment "1";\n'
                f'chr2\tbam\texon\t200\t249\t.\t+\t.\t{attr}segment "2";\n')
    gff.write_text(block('first')+block('last'))
    # The annotation output is the contract; no collapse or renderer fixture involved.
    gff_process.consolidateLongReadGFFs(str(gff),str(exon))
    rows=(tmp_path/'gff-output/transcript_associations.txt').read_text().splitlines()
    assert len(rows)==2
    first,last=[r.split('\t') for r in rows]
    assert first[:3]==last[:3]
    assert 'H:' in first[2]
    assert first[3]=='first' and last[3]=='last'
