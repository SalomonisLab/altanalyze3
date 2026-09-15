"""Multi-catalog validation-GTF support for SNAF-B.

The contract these tests pin down:
  1. one path behaves exactly as before, and keeps its plain abspath cache key
  2. splitting one catalog into two files and passing both reproduces the single-file parse
  3. a comma-separated string, a list, and repeated CLI values all mean the same thing
  4. the merged transcript lists stay sorted by first-exon start, which the support query needs
  5. an exon chain present in two catalogs is stored once
"""
import os

import pytest

from altanalyze3.components.snaf.surface import main as M


GTF_A = """chr1\tsrc\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1";transcript_id "tA";
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "G1";transcript_id "tA";
chr1\tsrc\texon\t800\t900\t.\t+\t.\tgene_id "G1";transcript_id "tA";
chr1\tsrc\ttranscript\t50\t400\t.\t+\t.\tgene_id "G1";transcript_id "tB";
chr1\tsrc\texon\t50\t120\t.\t+\t.\tgene_id "G1";transcript_id "tB";
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "G1";transcript_id "tB";
chr2\tsrc\ttranscript\t10\t99\t.\t-\t.\tgene_id "G2";transcript_id "tC";
chr2\tsrc\texon\t10\t40\t.\t-\t.\tgene_id "G2";transcript_id "tC";
chr2\tsrc\texon\t70\t99\t.\t-\t.\tgene_id "G2";transcript_id "tC";
"""

GTF_B = """chr1\tsrc\ttranscript\t500\t1500\t.\t+\t.\tgene_id "G1";transcript_id "tD";
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id "G1";transcript_id "tD";
chr1\tsrc\texon\t1400\t1500\t.\t+\t.\tgene_id "G1";transcript_id "tD";
chr3\tsrc\ttranscript\t1\t80\t.\t+\t.\tgene_id "G3";transcript_id "tE";
chr3\tsrc\texon\t1\t30\t.\t+\t.\tgene_id "G3";transcript_id "tE";
chr3\tsrc\texon\t60\t80\t.\t+\t.\tgene_id "G3";transcript_id "tE";
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def _chains(d):
    """{(chrom,strand): [exon-chain,...]} so comparisons ignore the attribute strings."""
    out = {}
    for chrom, sd in d.items():
        for strand in ('+', '-'):
            out[(chrom, strand)] = [tuple(tx[1:]) for tx in sd.get(strand, [])]
    return out


def _clear_caches():
    M._GTF_PARSE_CACHE.clear()
    M._GTF_STARTS_CACHE.clear()


def test_path_list_forms_agree(tmp_path):
    a = _write(tmp_path, 'a.gtf', GTF_A)
    b = _write(tmp_path, 'b.gtf', GTF_B)
    assert M._gtf_path_list(None) == []
    assert M._gtf_path_list(a) == [a]
    assert M._gtf_path_list([a, b]) == [a, b]
    assert M._gtf_path_list('{},{}'.format(a, b)) == [a, b]
    assert M._gtf_path_list([a, a]) == [a], 'a repeated catalog must collapse to one'


def test_single_path_keeps_its_cache_key(tmp_path):
    a = _write(tmp_path, 'a.gtf', GTF_A)
    assert M._gtf_cache_key(a) == os.path.abspath(a)
    assert M._gtf_cache_key([a]) == os.path.abspath(a)


def test_two_catalogs_merge_and_stay_sorted(tmp_path):
    _clear_caches()
    a = _write(tmp_path, 'a.gtf', GTF_A)
    b = _write(tmp_path, 'b.gtf', GTF_B)
    da = M.process_est_or_long_read_with_id(a)
    db = M.process_est_or_long_read_with_id(b)
    merged = M.process_est_or_long_read_with_id([a, b])

    n_a = sum(len(v['+']) + len(v['-']) for v in da.values())
    n_b = sum(len(v['+']) + len(v['-']) for v in db.values())
    n_m = sum(len(v['+']) + len(v['-']) for v in merged.values())
    assert n_m == n_a + n_b, 'transcript loss: {} + {} != {}'.format(n_a, n_b, n_m)
    assert set(merged) == set(da) | set(db)

    for chrom, sd in merged.items():
        for strand in ('+', '-'):
            starts = [int(tx[1][0]) for tx in sd.get(strand, [])]
            assert starts == sorted(starts), 'unsorted {} {}'.format(chrom, strand)

    key = M._gtf_cache_key([a, b])
    starts = M._GTF_STARTS_CACHE.get(key)
    assert starts is not None
    for chrom, sd in merged.items():
        for strand in ('+', '-'):
            assert starts[chrom][strand] == [int(tx[1][0]) for tx in sd.get(strand, [])]


def test_comma_string_equals_list(tmp_path):
    _clear_caches()
    a = _write(tmp_path, 'a.gtf', GTF_A)
    b = _write(tmp_path, 'b.gtf', GTF_B)
    assert _chains(M.process_est_or_long_read_with_id([a, b])) == \
        _chains(M.process_est_or_long_read_with_id('{},{}'.format(a, b)))


def test_split_then_merge_reproduces_one_file(tmp_path):
    """The invariant that matters: one catalog cut in two, then merged, is the original."""
    _clear_caches()
    whole = _write(tmp_path, 'whole.gtf', GTF_A + GTF_B)
    part1 = _write(tmp_path, 'part1.gtf', GTF_A)
    part2 = _write(tmp_path, 'part2.gtf', GTF_B)
    one = _chains(M.process_est_or_long_read_with_id(whole))
    two = _chains(M.process_est_or_long_read_with_id([part1, part2]))
    assert set(one) == set(two)
    for key in one:
        assert sorted(one[key]) == sorted(two[key]), 'chain mismatch at {}'.format(key)


def test_duplicate_chain_across_catalogs_stored_once(tmp_path):
    _clear_caches()
    a = _write(tmp_path, 'a.gtf', GTF_A)
    a_copy = _write(tmp_path, 'a_copy.gtf', GTF_A)
    single = M.process_est_or_long_read_with_id(a)
    both = M.process_est_or_long_read_with_id([a, a_copy])
    n_single = sum(len(v['+']) + len(v['-']) for v in single.values())
    n_both = sum(len(v['+']) + len(v['-']) for v in both.values())
    assert n_both == n_single, 'duplicate exon chains were not collapsed: {} vs {}'.format(
        n_both, n_single)


def test_last_transcript_of_a_classic_gtf_is_kept(tmp_path):
    """Regression: the streaming parser wrote a transcript only when it met the NEXT
    `transcript` line, so the last transcript of every file was dropped. Splitting a catalog
    then multiplied the loss by the number of files."""
    _clear_caches()
    whole = _write(tmp_path, 'whole.gtf', GTF_A + GTF_B)
    assert M._gtf_layout(whole) == 'transcript_first'
    d = M.process_est_or_long_read_with_id(whole)
    ids = sorted(tx[0].split('transcript_id "')[1].split('"')[0]
                 for sd in d.values() for s in ('+', '-') for tx in sd.get(s, []))
    assert ids == ['tA', 'tB', 'tC', 'tD', 'tE'], ids


def test_missing_catalog_still_raises(tmp_path):
    a = _write(tmp_path, 'a.gtf', GTF_A)
    with pytest.raises(Exception):
        M.process_est_or_long_read_with_id([a, str(tmp_path / 'absent.gtf')])


def test_cli_normalizes_every_argument_form(tmp_path):
    """What `--validation_gtf` accepts on the command line, and what the pipeline receives."""
    from altanalyze3.components.snaf.cli import _validation_gtf_paths as norm
    a = _write(tmp_path, 'a.gtf', GTF_A)
    b = _write(tmp_path, 'b.gtf', GTF_B)
    assert norm(None) is None
    assert norm([]) is None
    assert norm([a]) == a, 'one catalog must stay a plain string, to keep its cache key'
    assert norm(a) == a
    assert norm([a, b]) == [a, b]
    assert norm(['{},{}'.format(a, b)]) == [a, b]
    assert norm([a, a]) == a, 'the same catalog twice is one catalog'
    with pytest.raises(FileNotFoundError):
        norm([a, str(tmp_path / 'absent.gtf')])


def test_cli_parser_accepts_several_catalogs(tmp_path):
    """The argparse layer itself: `--validation_gtf x.gtf y.gtf` must not be a usage error."""
    from altanalyze3.utilities.parser import ArgsParser
    a = _write(tmp_path, 'a.gtf', GTF_A)
    b = _write(tmp_path, 'b.gtf', GTF_B)
    args = ArgsParser(['snaf-b', '--juncounts', a, '--db_dir', str(tmp_path),
                       '--output', str(tmp_path), '--validation_gtf', a, b])
    assert list(args.validation_gtf) == [a, b]
