"""Display metadata for transcripts containing more than one annotated gene.

This does not change counts, read assignments, clustering or transcript models.
Non-overlapping reference loci identify candidates, not a molecular mechanism.
"""
from itertools import combinations


def locus_axis_groups(models):
    """Share an axis for overlapping annotations; order disjoint axes genomically.

    Gene-reference rows remain separate. A genomic base must not be moved to
    another horizontal position just because its gene label changes.
    """
    chromosomes = {}
    for model in models:
        chrom = str(model.get("chrom", "")).upper().removeprefix("CHR")
        chrom = "M" if chrom == "MT" else chrom
        if not chrom or not model.get("extent"):
            raise ValueError("Verified reference chromosome and extent required")
        chromosomes.setdefault(chrom, []).append(model)
    groups = []
    for members in chromosomes.values():
        current, high = [], None
        for model in sorted(members, key=lambda m: min(m["extent"])):
            lo, hi = min(model["extent"]), max(model["extent"])
            if current and lo > high:
                groups.append(current)
                current = []
            current.append(model["gene"])
            high = max(high, hi) if current[:-1] else hi
        if current:
            groups.append(current)
    return groups


def coding_segments_for_loci(segments, models, cds_min, cds_max, cds_chrom):
    """Intersect a saved CDS with gene-qualified exons on its verified chromosome.

    A scalar genomic interval is insufficient for transcripts on different
    chromosomes. Return None for unavailable coverage, rather than calling
    unsupported portions UTR. Empty coverage and missing coverage are distinct.
    """
    if cds_min is None or cds_max is None or not cds_chrom:
        return None
    def chrom_key(value):
        key = str(value).upper().removeprefix("CHR")
        return "M" if key == "MT" else key
    if any(not models.get(s["gene"], {}).get("chrom") or
           chrom_key(models[s["gene"]]["chrom"]) != chrom_key(cds_chrom)
           for s in segments):
        return None
    out = []
    for s in segments:
        a, b = max(min(s["start"], s["end"]), cds_min), min(max(s["start"], s["end"]), cds_max)
        if b > a:
            out.append(dict(gene=s["gene"], start=a, end=b))
    return out


def is_trans_splicing_view(response):
    """Extra reference rows require a displayed candidate, not just multiple gene labels."""
    rows = list(response.get("isoforms", [])) + list(response.get("reference_isoforms", []))
    for panel in response.get("panels", []):
        rows.extend(panel.get("molecules", []))
    return any(row.get("trans_splicing_candidate") is True
               or row.get("locus_status") == "distinct_loci_candidate"
               or any(join.get("locus_status") == "distinct_loci_candidate"
                      for join in row.get("locus_joins", [])) for row in rows)


def annotate_response(ctx, response, model_track):
    records = list(response.get("isoforms", []))
    for panel in response.get("panels", []):
        records.extend(panel.get("molecules", []))
    # A composed response may reference one Python row object from several
    # panels. Decorate it once without removing any response row or count.
    records = list({id(row): row for row in records}.values())
    primary = response["gene"]
    genes = dict.fromkeys([primary])
    for row in records:
        path = []
        runs = []
        for segment in row.get("exon_segments", []):
            gene = segment.get("gene") or primary
            segment["gene"] = gene
            genes[gene] = None
            if not path or path[-1] != gene:
                path.append(gene)
                runs.append([])
            runs[-1].append(segment)
        token_path = []
        for token in row.get("annotation_tokens", []):
            gene = token.split(":", 1)[0] if ":" in token else primary
            genes[gene] = None
            if not token_path or token_path[-1] != gene:
                token_path.append(gene)
        row["gene_path"] = token_path or path
        row["unplaced_genes"] = [g for g in dict.fromkeys(row["gene_path"])
                                 if not any(s.get("gene") == g and s.get("start") is not None
                                            for s in row.get("exon_segments", []))]
        row["_locus_runs"] = runs

    models = {gene: {**model_track(ctx, gene),
                     "symbol": ctx.gene_symbol_dict.get(gene, gene)} for gene in genes}
    for row in records:
        for segment in row.get("exon_segments", []):
            model = models[segment["gene"]]
            segment["chrom"] = segment.get("chrom") or model.get("chrom")
            segment["strand"] = segment.get("strand") or model.get("strand")
        participating = list(dict.fromkeys(row["gene_path"]))
        row["locus_genes"] = participating
        row["locus_status"] = ("unresolved_loci" if row["unplaced_genes"] else
                               classify_loci([models[g] for g in participating]))
        row["trans_splicing_candidate"] = row["locus_status"] == "distinct_loci_candidate"
        row["locus_joins"] = []
        runs = row.pop("_locus_runs")
        if row["unplaced_genes"]:
            runs = []  # Do not draw an apparent direct join across an unplaced gene.
        for left_run, right_run in zip(runs, runs[1:]):
            left, right = left_run[0]["gene"], right_run[0]["gene"]
            a = [s for s in left_run if s.get("start") is not None and s.get("end") is not None]
            b = [s for s in right_run if s.get("start") is not None and s.get("end") is not None]
            if a and b:
                row["locus_joins"].append({
                    "from_gene": left, "to_gene": right,
                    # An overlapping third annotation must not hide a displayed
                    # join between two disjoint genes. The whole-row status above
                    # remains unchanged; this metadata controls reference rows.
                    "locus_status": classify_loci([models[left], models[right]]),
                    "from_coord": (min(s["start"] for s in a) if models[left].get("strand") == "-"
                                   else max(s["end"] for s in a)),
                    "to_coord": (max(s["end"] for s in b) if models[right].get("strand") == "-"
                                 else min(s["start"] for s in b)),
                })
        if len(participating) > 1:
            row["coordinate_source"] = "exon_annotation_tokens"
    response["gene_models"] = list(models.values())
    response["has_multiple_loci"] = len(models) > 1
    return response


def classify_loci(models):
    if len(models) < 2:
        return "single_gene"
    if any(not m.get("chrom") or not m.get("extent") for m in models):
        return "unresolved_loci"
    def chrom_key(value):
        key = str(value).upper().removeprefix("CHR")
        return "M" if key == "MT" else key
    # Closed genomic intervals: sharing even a boundary is genomic overlap.
    for a, b in combinations(models, 2):
        if chrom_key(a["chrom"]) == chrom_key(b["chrom"]):
            if max(a["extent"][0], b["extent"][0]) <= min(a["extent"][1], b["extent"][1]):
                return "overlapping_gene_annotations"
    return "distinct_loci_candidate"
