"""Read rows come from the stored genomic blocks, not from parsing the exon string.

gff_process writes the token string with `exon_str`, which drops the first and last token's
coordinate, and `uniqueList`, which drops a repeated token. Two different block layouts can
therefore share one string, so the string cannot always be decoded. gff_process also persists every
structure's exon blocks in <sample>/gff-output/structure_coords.sqlite keyed by (gene, structure).

FLI1 read 17261388 is the case that proves the point: 'E3.2_128694073|E3.2|E3.3|E3.4|
I3.1_128714735|I3.1_128715401' has to become 128,694,073-128,694,276 and 128,714,735-128,715,401,
and I3.1 abuts E3.4, so no rule over the string alone can know the read splices there.
"""
import os
import sqlite3
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.bam import coord_store  # noqa: E402
from altanalyze3.components.visualization import isoform_structure_view as isv  # noqa: E402
from altanalyze3.components.visualization.isv_web import data_api as da  # noqa: E402

G = "ENSG00000151702"
EXONS = {"E3.2": (128693916, 128694088), "E3.3": (128694088, 128694092),
         "E3.4": (128694092, 128694276), "I3.1": (128694277, 128748236),
         "E15.1": (128809157, 128809204), "I15.1": (128809205, 128810458),
         "E16.1": (128810459, 128812416)}
LOOKUP = {(G, k): {"chrom": "chr11", "strand": "+", "start": v[0], "end": v[1],
                   "type": k[0]} for k, v in EXONS.items()}

RAW = "E3.2_128694073|E3.2|E3.3|E3.4|I3.1_128714735|I3.1_128715401".split("|")
GFF = [(128694073, 128694276), (128714735, 128715401)]
BARE_END = "E15.1_128809163|E15.1|E16.1".split("|")
BARE_END_GFF = [(128809163, 128809204), (128810459, 128812416)]


class _Ctx:
    def __init__(self, db):
        self.exon_lookup = LOOKUP
        self.viewer_sample_dict = {"LIB": [{"transcript_associations": os.path.join(db, "ta.txt")}]}


def _store(tmp_path, key, blocks):
    path = tmp_path / "structure_coords.sqlite"
    con = sqlite3.connect(str(path))
    con.execute("CREATE TABLE structures (gene TEXT, structure TEXT, chrom TEXT, strand TEXT, "
                "source TEXT, blob BLOB, PRIMARY KEY (gene, structure))")
    con.execute("INSERT INTO structures VALUES (?,?,?,?,?,?)",
                (G, key, "chr11", "+", "t", coord_store.pack_exons(blocks)))
    con.commit()
    con.close()
    return path


def _blocks(segs):
    out = []
    for s in sorted(segs, key=lambda s: s["start"]):
        if out and s["start"] <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], s["end"]))
        else:
            out.append((s["start"], s["end"]))
    return out


def test_the_store_key_drops_only_the_two_terminal_coordinates():
    assert isv.exon_string_key(RAW) == "E3.2|E3.3|E3.4|I3.1_128714735"


def test_store_blocks_beat_the_string(tmp_path):
    _store(tmp_path, "E3.2|E3.3|E3.4|I3.1_128714735", GFF)
    da._struct_conns.clear()
    segs = da.read_segments(_Ctx(str(tmp_path)), G, RAW, "LIB")
    assert _blocks(segs) == GFF
    assert _blocks(isv.build_isoform_segments(RAW, LOOKUP, G)) != GFF   # the string alone cannot


def test_a_bare_terminal_token_ends_at_its_feature_edge(tmp_path):
    """The stored exemplar stops at 128,811,874; this read's last token is bare, so it runs to
    E16.1's own end."""
    _store(tmp_path, "E15.1|E16.1", [(128809157, 128809204), (128810459, 128811874)])
    da._struct_conns.clear()
    segs = da.read_segments(_Ctx(str(tmp_path)), G, BARE_END, "LIB")
    assert _blocks(segs) == BARE_END_GFF


def test_a_missing_store_falls_back_to_the_parser(tmp_path):
    da._struct_conns.clear()
    segs = da.read_segments(_Ctx(str(tmp_path)), G, RAW, "LIB")
    assert _blocks(segs) == _blocks(isv.build_isoform_segments(RAW, LOOKUP, G))


# ---------------------------------------------------------------------------------------------
# Where the coordinate store sits, in the four layouts this run actually uses. viewer_sample_dict
# points `transcript_associations` at gff-output/gene_indexes_v2/transcripts.db, one directory
# BELOW the store, and a run's own isv/reads/<lib>/gff-output holds no store at all, so the
# library's BAM directory is the last candidate.
class _PathCtx:
    def __init__(self, ta, matrix, metadata_file):
        self.viewer_sample_dict = {"LIB": [{"transcript_associations": ta, "matrix": matrix}]}
        self.metadata_file = metadata_file


def _metadata(tmp_path, bam):
    p = tmp_path / "meta.txt"
    p.write_text("uid\tbam\tlibrary\treverse\tgroups\nU\t%s\tLIB\tTRUE\tg\n" % bam)
    return str(p)


def test_store_found_beside_transcript_associations_parent(tmp_path):
    gff = tmp_path / "2-SLoutput" / "gff-output"
    (gff / "gene_indexes_v2").mkdir(parents=True)
    (gff / "structure_coords.sqlite").write_bytes(b"")
    ctx = _PathCtx(str(gff / "gene_indexes_v2" / "transcripts.db"),
                   str(tmp_path / "2-SLoutput" / "x.h5ad"),
                   _metadata(tmp_path, str(tmp_path / "2-SLoutput" / "x.bam")))
    assert da._struct_db_path(ctx, "LIB") == str(gff / "structure_coords.sqlite")


def test_store_found_beside_the_bam_when_the_run_staged_only_associations(tmp_path):
    staged = tmp_path / "isv" / "reads" / "LIB" / "gff-output"
    (staged / "gene_indexes_v2").mkdir(parents=True)
    src = tmp_path / "outputs" / "gff-output"
    src.mkdir(parents=True)
    (src / "structure_coords.sqlite").write_bytes(b"")
    ctx = _PathCtx(str(staged / "gene_indexes_v2" / "transcripts.db"),
                   str(tmp_path / "isv" / "reads" / "LIB" / "LIB.h5ad"),
                   _metadata(tmp_path, str(tmp_path / "outputs" / "scisoseq.mapped.bam")))
    assert da._struct_db_path(ctx, "LIB") == str(src / "structure_coords.sqlite")


def test_no_store_anywhere_returns_none(tmp_path):
    ctx = _PathCtx(str(tmp_path / "a" / "transcripts.db"), str(tmp_path / "b.h5ad"),
                   _metadata(tmp_path, str(tmp_path / "c.bam")))
    assert da._struct_db_path(ctx, "LIB") is None


def test_store_found_through_a_staged_symlink(tmp_path):
    """A run stages a sample by symlinking its files into isv/reads/<lib>/ while creating a real
    gene_indexes_v2 there, so the store only appears once the links are resolved. All 50 libraries
    of isoform_final_20261013 are staged this way."""
    src = tmp_path / "original" / "gff-output"
    (src).mkdir(parents=True)
    (src / "structure_coords.sqlite").write_bytes(b"")
    (tmp_path / "original" / "lib.gff.gz").write_bytes(b"")
    staged = tmp_path / "isv" / "reads" / "LIB"
    (staged / "gff-output" / "gene_indexes_v2").mkdir(parents=True)
    (staged / "lib.gff.gz").symlink_to(tmp_path / "original" / "lib.gff.gz")
    meta = tmp_path / "meta.txt"
    meta.write_text("uid\tgff\tmatrix\tlibrary\treverse\tgroups\n"
                    "U\t%s\t%s\tLIB\tTRUE\tg\n"
                    % (staged / "lib.gff.gz", staged / "LIB.h5ad"))
    ctx = _PathCtx(str(staged / "gff-output" / "gene_indexes_v2" / "transcripts.db"),
                   str(staged / "LIB.h5ad"), str(meta))
    assert da._struct_db_path(ctx, "LIB") == str(src / "structure_coords.sqlite")


# ---------------------------------------------------------------------------------------------
# A terminal coordinate is applied only while the block stays a real interval.
# RUNX1 (minus strand) read 34716389 ends '...|E116.1_34792542|E116.1_34791917'. Its lowest GFF
# block is 34,791,786-34,791,917, so 34,791,917 is that block's END, not the read's outer edge.
# Writing it as the block's START collapsed it to 1 bp. 2 of 1,529 RUNX1 reads of AML-13 took that
# shape.
RUNX1 = "ENSG00000159216"
RUNX1_LOOKUP = {
    (RUNX1, "E92.1"): {"chrom": "chr21", "strand": "-", "start": 34880557, "end": 34880659, "type": "E"},
    (RUNX1, "E116.1"): {"chrom": "chr21", "strand": "-", "start": 34791786, "end": 34792610, "type": "E"},
}
RUNX1_TOKENS = "E92.1_34880659|E92.1|E116.1|E116.1_34792542|E116.1_34791917".split("|")
RUNX1_STORE = [(34791786, 34791917), (34792542, 34792610), (34880557, 34880659)]


def test_a_terminal_coordinate_never_collapses_its_block():
    got = isv.apply_terminal_bounds(RUNX1_STORE, RUNX1_TOKENS, RUNX1, RUNX1_LOOKUP, "-")
    assert got == RUNX1_STORE
    assert (34791917, 34791917) not in got


def test_a_terminal_coordinate_still_trims_and_extends():
    """LRFN4 35292108 extends its first block to 66,856,538 and trims its last to 66,860,473."""
    lookup = {(G, "U1.1"): {"chrom": "chr11", "strand": "+", "start": 66856400,
                            "end": 66856646, "type": "U"},
              (G, "E4.2"): {"chrom": "chr11", "strand": "+", "start": 66860247,
                            "end": 66860475, "type": "E"}}
    store = [(66856647, 66856657), (66859746, 66860475)]
    toks = "U1.1_66856538|E4.2_66860473".split("|")
    assert isv.apply_terminal_bounds(store, toks, G, lookup, "+") == \
        [(66856538, 66856657), (66859746, 66860473)]
