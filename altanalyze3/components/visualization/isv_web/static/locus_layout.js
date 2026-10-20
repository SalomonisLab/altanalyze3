/* Independent genomic axes for multi-locus transcripts; input/API data remain untouched. */
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.ISVLocusLayout = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";
  function referenceWindows(loci, locus, stacked) {
    if (!stacked || !locus.chrom) return [locus];
    const key = (s) => String(s || "").toUpperCase().replace(/^CHR/, "").replace(/^MT$/, "M");
    const seen = new Set();
    return loci.filter((window) => {
      if (!window.chrom || key(window.chrom) !== key(locus.chrom)) return false;
      if (!window.display_extent) return true;
      const axis = window.display_extent.join(":");
      if (seen.has(axis)) return false;
      seen.add(axis);
      return true;
    });
  }
  function isTransSplicingView(res) {
    const rows = (res.isoforms || []).concat(res.reference_isoforms || [],
      (res.panels || []).flatMap((p) => p.molecules || []));
    return rows.some((r) => r.trans_splicing_candidate === true || r.locus_status === "distinct_loci_candidate"
      || (r.locus_joins || []).some((j) => j.locus_status === "distinct_loci_candidate"));
  }
  function prepare(res) {
    if (!res.has_multiple_loci || res._locusProjected) return res;
    const rows = (res.isoforms || []).concat(res.reference_isoforms || [], (res.panels || []).flatMap((p) => p.molecules || []));
    const chromosomes = new Map();
    (res.gene_models || []).forEach((model) => {
      const chrom = String(model.chrom || "").toUpperCase().replace(/^CHR/, "").replace(/^MT$/, "M");
      if (!chrom || !model.extent) throw new Error("Verified reference chromosome and extent required");
      if (!chromosomes.has(chrom)) chromosomes.set(chrom, []);
      chromosomes.get(chrom).push(model);
    });
    const groups = [];
    chromosomes.forEach((models) => {
      let current = [], high = -Infinity;
      models.slice().sort((a, b) => Math.min(...a.extent) - Math.min(...b.extent)).forEach((model) => {
        if (current.length && Math.min(...model.extent) > high) { groups.push(current); current = []; high = -Infinity; }
        current.push(model); high = Math.max(high, ...model.extent);
      });
      if (current.length) groups.push(current);
    });
    const axisByGene = new Map();
    groups.forEach((models, i) => {
      let lo = Infinity, hi = -Infinity;
      const include = (c) => { if (Number.isFinite(c)) { lo = Math.min(lo, c); hi = Math.max(hi, c); } };
      const genes = new Set(models.map((m) => m.gene));
      models.forEach((model) => (model.extent || []).forEach(include));
      rows.forEach((r) => (r.exon_segments || []).forEach((s) => {
        if (genes.has(s.gene)) { include(s.start); include(s.end); }
      }));
      const coordinates_available = Number.isFinite(lo);
      if (!coordinates_available) { lo = 0; hi = 1; }
      models.forEach((model) => axisByGene.set(model.gene, { ...model, genomic_extent: [lo, Math.max(lo + 1, hi)],
        display_extent: [i * 1120, i * 1120 + 1000], coordinates_available }));
    });
    const loci = (res.gene_models || []).map((m) => axisByGene.get(m.gene));
    const byGene = new Map(loci.map((l) => [l.gene, l]));
    function project(gene, coord) {
      const l = byGene.get(gene);
      return coord == null || !l ? null : l.display_extent[0]
        + (coord - l.genomic_extent[0]) / (l.genomic_extent[1] - l.genomic_extent[0]) * 1000;
    }
    function row(r) {
      const multi = (r.locus_genes || []).length > 1;
      const gene = (r.locus_genes || [res.gene])[0];
      let coding;
      const chromKey = (s) => String(s || "").toUpperCase().replace(/^CHR/, "").replace(/^MT$/, "M");
      const segments = r.exon_segments || [];
      if (multi && r.cds_min != null && r.cds_max != null && r.cds_chrom
          && segments.every((s) => byGene.has(s.gene) && byGene.get(s.gene).chrom
            && chromKey(byGene.get(s.gene).chrom) === chromKey(r.cds_chrom))) {
        coding = segments.flatMap((s) => {
          const a = Math.max(Math.min(s.start, s.end), r.cds_min);
          const b = Math.min(Math.max(s.start, s.end), r.cds_max);
          return b > a ? [{ gene: s.gene, a: project(s.gene, a), b: project(s.gene, b) }] : [];
        });
      }
      return { ...r,
        exon_segments: (r.exon_segments || []).map((s) => ({ ...s,
          genomic_start: s.start, genomic_end: s.end,
          start: project(s.gene || res.gene, s.start), end: project(s.gene || res.gene, s.end) })),
        // One genomic CDS interval cannot describe coding sequence across loci.
        cds_min: multi ? null : project(gene, r.cds_min),
        cds_max: multi ? null : project(gene, r.cds_max),
        _locus_cds_intervals: coding,
        locus_joins: (r.locus_joins || []).map((j) => ({ ...j,
          display_from: project(j.from_gene, j.from_coord), display_to: project(j.to_gene, j.to_coord) })),
      };
    }
    return { ...res, _locusProjected: true,
      isoforms: res.isoforms ? res.isoforms.map(row) : res.isoforms,
      reference_isoforms: res.reference_isoforms ? res.reference_isoforms.map(row) : res.reference_isoforms,
      panels: res.panels ? res.panels.map((p) => ({ ...p, molecules: p.molecules.map(row) })) : res.panels,
      gene_model: { ...res.gene_model, chrom: "", strand: null, blocks: [],
        extent: [0, groups.length * 1120 - 120], loci, stack_reference_loci: isTransSplicingView(res) },
    };
  }
  return { prepare, isTransSplicingView, referenceWindows };
});
