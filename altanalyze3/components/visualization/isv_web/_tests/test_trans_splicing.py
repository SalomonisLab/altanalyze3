"""Regression checks for multi-locus display; no biological inference is fitted."""
from copy import deepcopy
from types import SimpleNamespace
import unittest

from altanalyze3.components.visualization import isoform_structure_view as isv
from altanalyze3.components.visualization.isv_web import data_api as da
from altanalyze3.components.visualization.isv_web.locus_view import classify_loci, is_trans_splicing_view, coding_segments_for_loci


def context(partner_chrom="chr2", partner_start=100):
    segments = {}
    for gene, chrom, strand, offset in [("A", "chr1", "+", 100), ("B", partner_chrom, "-", partner_start)]:
        segments[gene] = [dict(gene=gene, chrom=chrom, strand=strand, exon_id="E1.1", type="E", start=offset, end=offset + 40),
                          dict(gene=gene, chrom=chrom, strand=strand, exon_id="E2.1", type="E", start=offset + 100, end=offset + 140)]
    return SimpleNamespace(gene_segments=segments,
                           exon_lookup={(g, s["exon_id"]): s for g, parts in segments.items() for s in parts},
                           gene_symbol_dict={"A": "GENE_A", "B": "GENE_B"}, resolve_gene=lambda g: g)


class TransSplicingDisplayTests(unittest.TestCase):
    def test_reference_rows_extend_through_same_chromosome_windows_only(self):
        from altanalyze3.components.visualization.isv_web.isv_web_remote import prepare_locus_coordinates
        for chrom, expected_windows in [("chr1", 2), ("chr2", 1)]:
            ctx = context(partner_chrom=chrom, partner_start=1000)
            segments = [da._seg_json(s) for s in isv.build_isoform_segments(
                ["E1.1", "B:E1.1"], ctx.exon_lookup, "A")]
            row = dict(isoform_id="read", count=7, exon_segments=segments,
                       registered_exon_segments=deepcopy(segments))
            response = da._with_loci(ctx, dict(gene="A", panels=[dict(molecules=[row])]))
            original = deepcopy(response)
            projected, frames = prepare_locus_coordinates(response, .1, True, True)
            self.assertEqual(response, original)
            self.assertEqual(projected["panels"][0]["molecules"][0]["count"], 7)
            self.assertTrue(all(len(f["reference_windows"]) == expected_windows for f in frames))
            if chrom == "chr1":
                self.assertEqual(frames[1]["reference_windows"][0]["model"]["blocks"], [])

    def test_saved_cds_is_split_per_gene_before_axis_projection(self):
        segments = [dict(gene="A", start=100, end=200), dict(gene="B", start=1000, end=1100)]
        models = {"A": dict(chrom="chr1"), "B": dict(chrom="1")}
        self.assertEqual(coding_segments_for_loci(segments, models, 150, 1050, "1"),
                         [dict(gene="A", start=150, end=200), dict(gene="B", start=1000, end=1050)])
        models["B"]["chrom"] = "chr2"
        self.assertIsNone(coding_segments_for_loci(segments, models, 150, 1050, "1"))
        self.assertIsNone(coding_segments_for_loci(segments, models, 150, 1050, None))

    def test_extra_reference_rows_only_for_displayed_trans_splicing(self):
        ordinary = dict(isoform_id="ordinary", locus_status="overlapping_gene_annotations",
                        trans_splicing_candidate=False)
        candidate = dict(isoform_id="candidate", locus_status="distinct_loci_candidate",
                         trans_splicing_candidate=True)
        self.assertFalse(is_trans_splicing_view(dict(has_multiple_loci=True, isoforms=[ordinary])))
        self.assertTrue(is_trans_splicing_view(dict(isoforms=[ordinary, candidate])))
        self.assertFalse(is_trans_splicing_view(dict(isoforms=[], panels=[dict(molecules=[ordinary])])))
        self.assertTrue(is_trans_splicing_view(dict(panels=[dict(molecules=[candidate])])))
        composite = dict(ordinary, locus_joins=[dict(locus_status="distinct_loci_candidate")])
        self.assertTrue(is_trans_splicing_view(dict(isoforms=[composite])))
        self.assertFalse(composite["trans_splicing_candidate"])
        composite["locus_joins"][0]["locus_status"] = "overlapping_gene_annotations"
        self.assertFalse(is_trans_splicing_view(dict(isoforms=[composite])))

    def test_independent_gene_blocks_with_shared_numbers_and_opposite_strands(self):
        ctx = context()
        segments = isv.build_isoform_segments(["E1.1", "E2.1", "B:E2.1", "B:E1.1"], ctx.exon_lookup, "A")
        self.assertEqual([(s["gene"], s["start"], s["end"]) for s in segments],
                         [("A", 100, 140), ("A", 200, 240), ("B", 100, 140), ("B", 200, 240)])

    def test_gene_boundary_cannot_close_primary_block_at_partner_coordinate(self):
        ctx = context(partner_start=10)
        segs = isv.build_isoform_segments(["E1.1", "B:E1.1_20", "B:E1.1_35"], ctx.exon_lookup, "A")
        self.assertEqual((segs[0]["gene"], segs[0]["start"], segs[0]["end"]), ("A", 100, 140))
        self.assertTrue(all(s["start"] < 100 for s in segs if s["gene"] == "B"))

    def test_single_gene_geometry_and_terminal_bounds(self):
        ctx = context()
        a = isv.build_isoform_segments(["E1.1_110", "E1.1_135", "E2.1"], ctx.exon_lookup, "A")
        self.assertEqual([(s["start"], s["end"]) for s in a], [(110, 135), (200, 240)])
        b = isv.build_isoform_segments(["B:E2.1", "B:E1.1"], ctx.exon_lookup, "A")
        self.assertEqual([(s["start"], s["end"]) for s in b], [(100, 140), (200, 240)])

    def test_display_metadata_preserves_every_record_and_count(self):
        ctx = context()
        segs = isv.build_isoform_segments(["E1.1", "B:E2.1", "B:E1.1"], ctx.exon_lookup, "A")
        row = dict(isoform_id="read1", count=7, exon_segments=[da._seg_json(s) for s in segs], expression={"HSC": 7})
        response = dict(gene="A", panels=[dict(condition="HSC", molecules=[row], structure_reads={"read1": 7})])
        original = deepcopy(response)
        result = da._with_loci(ctx, response)
        self.assertTrue(result["has_multiple_loci"])
        self.assertEqual(len(result["panels"][0]["molecules"]), 1)
        self.assertEqual(row["count"], original["panels"][0]["molecules"][0]["count"])
        self.assertEqual(row["expression"], {"HSC": 7})
        self.assertEqual(result["panels"][0]["structure_reads"], {"read1": 7})
        self.assertTrue(row["trans_splicing_candidate"])
        self.assertEqual(row["gene_path"], ["A", "B"])
        self.assertEqual(row["locus_joins"][0]["to_coord"], 240)
        self.assertEqual({s["chrom"] for s in row["exon_segments"]}, {"chr1", "chr2"})

    def test_overlapping_annotations_are_not_trans_splicing_candidates(self):
        ctx = context(partner_chrom="chr1")
        models = [da.gene_model_track(ctx, g) for g in ["A", "B"]]
        self.assertEqual(classify_loci(models), "overlapping_gene_annotations")
        models[1]["extent"] = [240, 300]
        self.assertEqual(classify_loci(models), "overlapping_gene_annotations")
        models[1]["extent"] = [241, 300]
        self.assertEqual(classify_loci(models), "distinct_loci_candidate")
        self.assertEqual(classify_loci([dict(chrom="chrM", extent=[1, 20]),
                                       dict(chrom="MT", extent=[10, 30])]), "overlapping_gene_annotations")

    def test_unresolved_reference_is_not_called_trans_splicing(self):
        self.assertEqual(classify_loci([dict(chrom="chr1", extent=[1, 2]), dict(chrom=None, extent=None)]), "unresolved_loci")

    def test_multi_locus_reads_do_not_use_untagged_coordinate_store(self):
        ctx = context()
        segs = da.read_segments(ctx, "A", ["E1.1", "B:E1.1"], "nonexistent_library")
        self.assertEqual({s["gene"] for s in segs}, {"A", "B"})

    def test_unplaced_partner_is_retained_and_not_called_a_candidate(self):
        ctx = context()
        row = dict(isoform_id="read1", count=3, annotation_tokens=["E1.1", "B:E999.1"],
                   exon_segments=[da._seg_json(s) for s in isv.build_isoform_segments(["E1.1"], ctx.exon_lookup, "A")])
        result = da._with_loci(ctx, dict(gene="A", isoforms=[row]))
        self.assertTrue(result["has_multiple_loci"])
        self.assertEqual(row["count"], 3)
        self.assertEqual(row["unplaced_genes"], ["B"])
        self.assertFalse(row["trans_splicing_candidate"])
        self.assertEqual(row["locus_joins"], [])

    def test_return_to_first_gene_uses_the_correct_run_endpoints(self):
        ctx = context()
        segs = isv.build_isoform_segments(["E1.1", "B:E1.1", "E2.1"], ctx.exon_lookup, "A")
        row = dict(isoform_id="read1", count=1, exon_segments=[da._seg_json(s) for s in segs])
        da._with_loci(ctx, dict(gene="A", panels=[dict(molecules=[row])], isoforms=[row]))
        self.assertEqual(row["gene_path"], ["A", "B", "A"])
        self.assertEqual(row["locus_joins"][0]["from_coord"], 140)
        self.assertEqual(row["locus_joins"][1]["to_coord"], 200)

    def test_native_pileup_projects_shared_registered_models_once(self):
        from altanalyze3.components.visualization.isv_web.isv_web_remote import prepare_locus_coordinates
        ctx = context()
        segs = [da._seg_json(s) for s in isv.build_isoform_segments(["E1.1", "B:E1.1"], ctx.exon_lookup, "A")]
        registered = deepcopy(segs)
        rows = [dict(isoform_id="read" + str(i), final_isoform_id="final", count=1,
                     exon_segments=deepcopy(segs), registered_exon_segments=registered,
                     cds_min=100, cds_max=240) for i in range(2)]
        res = da._with_loci(ctx, dict(gene="A", panels=[dict(n_molecules=2, molecules=rows)]))
        original = deepcopy(res)
        projected, frames = prepare_locus_coordinates(res, .1, True, True)
        self.assertEqual(res, original)
        self.assertEqual(len(frames), 2)
        for row in projected["panels"][0]["molecules"]:
            self.assertEqual(row["count"], 1)
            self.assertEqual(len(row["_registered_locus_merged"]["exons"]), 2)
            self.assertTrue(all(s["b"] > s["a"] for s in row["_registered_locus_merged"]["exons"]))
            self.assertIsNone(row["cds_min"])
        self.assertLess(projected["panels"][0]["molecules"][0]["registered_exons"][0][1],
                        projected["panels"][0]["molecules"][0]["registered_exons"][1][0])

    def test_native_structure_track_requires_the_registered_model(self):
        from altanalyze3.components.visualization.isv_web.isv_web_remote import prepare_locus_coordinates
        ctx = context()
        row = dict(isoform_id="read1", final_isoform_id="final1", count=1,
                   exon_segments=[da._seg_json(s) for s in isv.build_isoform_segments(["E1.1", "B:E1.1"], ctx.exon_lookup, "A")])
        res = da._with_loci(ctx, dict(gene="A", panels=[dict(n_molecules=1, molecules=[row])]))
        with self.assertRaisesRegex(ValueError, "registered structure required"):
            prepare_locus_coordinates(res, .1, True, True)


if __name__ == "__main__":
    unittest.main()
