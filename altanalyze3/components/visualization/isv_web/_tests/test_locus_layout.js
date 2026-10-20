"use strict";
const assert = require("node:assert/strict");
const { prepare, isTransSplicingView, referenceWindows } = require("../static/locus_layout.js");
const row = { trans_splicing_candidate: true, isoform_id: "read1", count: 7, expression: { HSC: 7 }, gene_path: ["A", "B"],
  locus_genes: ["A", "B"], cds_min: 110, cds_max: 240,
  exon_segments: [{ gene: "A", chrom: "chr1", start: 100, end: 140, type: "E" },
    { gene: "B", chrom: "chr2", start: 100, end: 140, type: "E" }],
  locus_joins: [{ from_gene: "A", to_gene: "B", from_coord: 140, to_coord: 100 }] };
const source = { gene: "A", has_multiple_loci: true, gene_model: { gene: "A" },
  gene_models: [{ gene: "A", extent: [100, 240], chrom: "chr1" }, { gene: "B", extent: [100, 240], chrom: "chr2" }],
  isoforms: [row], panels: [{ molecules: [row], structure_reads: { read1: 7 } }] };
const before = JSON.stringify(source);
const result = prepare(source);
assert.equal(JSON.stringify(source), before, "API counts and genomic coordinates must not be mutated");
assert.equal(result.isoforms[0].count, 7);
assert.deepEqual(result.isoforms[0].expression, row.expression);
assert.deepEqual(result.panels[0].structure_reads, { read1: 7 });
assert.notEqual(result.isoforms[0].exon_segments[0].start, result.isoforms[0].exon_segments[1].start,
  "identical chromosome-coordinate numbers must remain in separate locus panels");
assert.equal(result.isoforms[0].exon_segments[1].genomic_start, 100);
assert.equal(result.isoforms[0].exon_segments[1].chrom, "chr2");
assert.equal(result.isoforms[0].cds_min, null, "a single CDS span must not be projected across different genes");
assert.equal(result.gene_model.stack_reference_loci, true);
const ordinarySource = { ...source, isoforms: [{ ...row, trans_splicing_candidate: false, locus_status: "overlapping_gene_annotations" }], panels: [] };
assert.equal(prepare(ordinarySource).gene_model.stack_reference_loci, false, "ordinary multi-gene annotations do not add reference rows");
const compositeSource = { ...ordinarySource, isoforms: [{ ...ordinarySource.isoforms[0],
  locus_joins: [{ locus_status: "distinct_loci_candidate" }] }] };
assert.equal(prepare(compositeSource).gene_model.stack_reference_loci, true,
  "a displayed disjoint join retains separate reference rows when a third annotation overlaps");
assert.equal(compositeSource.isoforms[0].trans_splicing_candidate, false,
  "reference row layout does not reclassify the whole composite model");
assert.equal(isTransSplicingView({ panels: [{ molecules: [{ trans_splicing_candidate: true }] }] }), true);
assert.equal(prepare(result), result, "resize/redraw cannot project coordinates twice");
const normal = { has_multiple_loci: false, isoforms: [row] };
assert.equal(prepare(normal), normal, "single-locus response uses its existing rendering path");
const codingSource = JSON.parse(JSON.stringify(source));
codingSource.gene_models[1].chrom = "1";
Object.assign(codingSource.isoforms[0], { cds_min: 110, cds_max: 130, cds_chrom: "1" });
const coding = prepare(codingSource).isoforms[0];
assert.equal(coding._locus_cds_intervals.length, 2);
assert.equal(coding._locus_cds_intervals[0].a, coding._locus_cds_intervals[1].a,
  "overlapping annotations of the same chromosome coordinate must share an axis");
assert.equal(codingSource.isoforms[0].cds_min, 110, "saved CDS coordinates remain untouched");
codingSource.gene_models[1].chrom = "chr2";
assert.equal(prepare(codingSource).isoforms[0]._locus_cds_intervals, undefined,
  "a single-chromosome CDS cannot imply UTR assignments on another chromosome");
const windows = [{ gene: "A", chrom: "chr6" }, { gene: "B", chrom: "6" }, { gene: "C", chrom: "chr2" }];
assert.deepEqual(referenceWindows(windows, windows[1], true), windows.slice(0, 2));
assert.deepEqual(referenceWindows(windows, windows[1], false), [windows[1]],
  "ordinary reference tracks keep their original extent");
console.log("PASS multi-locus axes, counts, real coordinates, CDS guard, redraw, and single-locus parity");
