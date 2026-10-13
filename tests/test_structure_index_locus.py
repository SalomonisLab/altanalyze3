"""StructureIndex.exons: the gene-agnostic fallback must not return another locus's molecule.

Bare molecule ids repeat across libraries. On KINNEX-5, FLI1 (ENSG00000151702, chr11) final isoform
40148595.WF80-N-0805 is not registered under FLI1, and the old fallback returned ENSG00000106052 molecule
40148595 (chr7:27,769,688-27,828,861), which the ISV drew on the chr11 axis."""
import sqlite3

from altanalyze3.components.visualization.isv_web.precompute_gff_structures import StructureIndex

ROWS = [
    # gene, tid, tid_base, chrom, strand, source, exons
    ("ENSG_A", "100", "100", "chr11", "+", "bam", "1000-1100,2000-2100"),          # registered under the gene
    ("ENSG_OTHER", "200", "200", "chr7", "+", "bam", "27769688-27769834,27828628-27828861"),  # off-locus collision
    ("ENSG_READTHRU", "300", "300", "chr11", "+", "bam", "1500-1600,2500-2600"),   # other gene id, same locus
    ("ENSG_ANTI", "400", "400", "chr11", "-", "bam", "1500-1600"),                 # same span, opposite strand
    ("ENSG_X", "500", "500", "chr11", "+", "bam", "1200-1300"),                    # two distinct rows at the locus
    ("ENSG_Y", "500", "500", "chr11", "+", "bam", "1400-1450"),
    ("ENSG_B", "ENST00000999999", "ENST00000999999", "chr11", "+", "HAVANA", "1700-1800"),
]
LOCUS = ("chr11", "+", 900, 3000)


def _index(tmp_path):
    p = tmp_path / "final_structures.db"
    con = sqlite3.connect(p)
    con.execute("CREATE TABLE structures (gene TEXT, tid TEXT, tid_base TEXT, chrom TEXT, strand TEXT, "
                "source TEXT, exons TEXT)")
    con.executemany("INSERT INTO structures VALUES (?,?,?,?,?,?,?)", ROWS)
    con.commit(); con.close()
    return StructureIndex(str(p))


def test_gene_match_unchanged(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "100.LIB1", locus=LOCUS) == [(1000, 1100), (2000, 2100)]
    assert s.exons("ENSG_A", "100.LIB1") == [(1000, 1100), (2000, 2100)]
    assert s.stats["gene_match"] == 2


def test_off_locus_collision_rejected(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "200.WF80-N-0805", locus=LOCUS) is None
    assert s.stats["fallback_rejected"] == 1
    assert "chr7" in s.rejected[0][1]


def test_bare_id_without_locus_rejected(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "200.WF80-N-0805") is None
    assert s.exons("ENSG_A", "300.LIB1") is None


def test_same_locus_other_gene_id_accepted(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "300.LIB1", locus=LOCUS) == [(1500, 1600), (2500, 2600)]
    assert s.stats["fallback_at_locus"] == 1


def test_opposite_strand_rejected(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "400.LIB1", locus=LOCUS) is None


def test_two_distinct_structures_at_locus_ambiguous(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "500.LIB1", locus=LOCUS) is None
    assert s.stats["fallback_ambiguous"] == 1


def test_ensembl_id_fallback_kept(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "ENST00000999999.3") == [(1700, 1800)]
    assert s.exons("ENSG_A", "ENST00000999999", locus=LOCUS) == [(1700, 1800)]


def test_unregistered(tmp_path):
    s = _index(tmp_path)
    assert s.exons("ENSG_A", "999.LIB1", locus=LOCUS) is None
    assert s.stats["not_registered"] == 1


def test_source_scoped_ids_do_not_share_numeric_cache_or_legacy_rows(tmp_path):
    s = _index(tmp_path)
    c = sqlite3.connect(s.path)
    c.executemany('INSERT INTO structures VALUES (?,?,?,?,?,?,?)', [
        ('ENSG_A', '100.LIB_A__source001', '100.LIB_A__source001', 'chr11', '+', 'bam', '3000-3100'),
        ('ENSG_A', '100.LIB_B__source002', '100.LIB_B__source002', 'chr11', '+', 'bam', '4000-4100'),
    ])
    c.commit(); c.close()
    assert s.exons('ENSG_A', '100.LIB_A__source001') == [(3000, 3100)]
    assert s.exons('ENSG_A', '100.LIB_B__source002') == [(4000, 4100)]
    assert s.exons('ENSG_A', '100.LIB_A__source001') == [(3000, 3100)]
    assert s.exons('ENSG_A', '100.LEGACY') == [(1000, 1100), (2000, 2100)]


def test_exact_scoped_id_does_not_fall_back_to_another_genes_bare_id(tmp_path):
    s = _index(tmp_path)
    c = sqlite3.connect(s.path)
    c.execute('INSERT INTO structures VALUES (?,?,?,?,?,?,?)',
              ('ENSG_OTHER', '100.LIB_B', '100.LIB_B', 'chr7', '+', 'bam', '7000-7100'))
    c.commit(); c.close()
    assert s.exons('ENSG_A', '100.LIB_B') is None
    assert s.exons('ENSG_OTHER', '100.LIB_B') == [(7000, 7100)]
