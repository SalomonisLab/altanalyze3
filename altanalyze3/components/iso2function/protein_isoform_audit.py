"""Audit a complete ISV panel and map protein features without replacing native ORFs.

Run as a module from the repository parent. All sequence identities are exact. Domain
projection uses amino-acid alignment, rather than genomic overlap, to respect reading
frames. A projected feature is called intact only when its entire reference sequence
occurs contiguously and uniquely in the query. Other projections remain descriptive.
This does not discover novel domains or establish translation in vivo.
"""
import argparse
import csv
import hashlib
import gzip
import io
import itertools
import json
import pathlib
import re

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.file_digest(open(path, "rb"), "sha256").hexdigest()


def exact_locations(sequence, fragment):
    return [m.start() for m in re.finditer("(?=" + re.escape(fragment) + ")", sequence)]


def project_feature(reference, query, start, end, mapping):
    """1-based inclusive feature coordinates; never call a partial motif intact."""
    if not (1 <= start <= end <= len(reference)):
        raise ValueError("Feature coordinates must be 1-based inclusive within the reference protein")
    fragment = reference[start - 1:end]
    hits = exact_locations(query, fragment)
    mapped = [(i, mapping[i]) for i in range(start - 1, end) if i in mapping]
    identical = [(i, q) for i, q in mapped if reference[i] == query[q]]
    result = {"reference_start": start, "reference_end": end,
              "reference_length": end - start + 1,
              "aligned_residues": len(mapped), "identical_residues": len(identical),
              "exact_fraction": len(identical) / (end - start + 1),
              "query_start": min((q + 1 for _, q in mapped), default=""),
              "query_end": max((q + 1 for _, q in mapped), default="")}
    if len(hits) == 1:
        result.update(status="intact_sequence", query_start=hits[0] + 1,
                      query_end=hits[0] + len(fragment), exact_fraction=1.0,
                      identical_residues=len(fragment), aligned_residues=len(fragment))
    elif len(hits) > 1:
        result.update(status="ambiguous_repeated_sequence", query_start="", query_end="")
    elif not mapped:
        result["status"] = "not_aligned"
    else:
        result["status"] = "altered_or_partial"
    return result


def align_mapping(reference, query):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, query)[0]
    mapping = {}
    for (a, b), (c, d) in zip(*alignment.aligned):
        mapping.update(zip(range(int(a), int(b)), range(int(c), int(d))))
    return mapping


FEATURE_TYPES = {"Domain", "Zinc finger", "DNA binding", "Motif", "Region", "Coiled coil"}


def matched_activation(rows, measurements):
    """A blank identifier is never identity evidence, even when both sides are blank."""
    result = []
    for r in rows:
        clones = {x for x in r["tfiso_exact_clones"].split(";") if x}
        for m in measurements:
            if m["clone_id"] and m["clone_id"] in clones:
                if m["gene_symbol"] != r["gene_symbol"]:
                    raise ValueError("Exact clone gene identity disagrees with activation source")
                result.append(dict(m, gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"]))
    return result


def validate_native_inputs(rows, rendered, decision):
    """Reject incomplete/replaced identities and altered native sequences before annotation."""
    identities = [(r["gene_symbol"], r["isoform_id"]) for r in rows]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate required identity")
    if set(identities) != rendered:
        raise ValueError("Complete native render roster and source table do not agree")
    if decision["pending_question"] is not None or decision["actual_user_answer"] is None:
        raise ValueError("Unresolved source decision; cannot execute annotation")
    if [(r["gene_symbol"], r["isoform_id"]) for r in decision["required_roster"]] != identities:
        raise ValueError("Expected roster changed")
    for r in rows:
        for name in ("protein_sequence", "cds_sequence", "cdna_sequence"):
            if not r[name]:
                raise ValueError(f"Missing required {name}: {r['isoform_id']}")
        if len(r["cds_sequence"]) % 3:
            raise ValueError(f"Incomplete CDS codon: {r['isoform_id']}")
        if str(Seq(r["cds_sequence"]).translate()).rstrip("*") != r["protein_sequence"]:
            raise ValueError(f"CDS translation differs: {r['isoform_id']}")
        if r["cds_sequence"] not in r["cdna_sequence"]:
            raise ValueError(f"CDS not in native cDNA: {r['isoform_id']}")


def audit(run, out, reference_json, reference_fasta, data):
    out.mkdir(parents=True, exist_ok=True)
    source = run / "isv_pdf/HSC_TF_isoforms/top_isoforms_unique_for_literature.tsv"
    rows = [r for r in read_tsv(source) if "HSC_MPP_pooled:" in r["panels"]]
    rendered = set()
    render_files = sorted((run / "isv_pdf/HSC_TF_isoforms/HSC_MPP_pooled").glob("*.render.json"))
    for path in render_files:
        gene = path.name.split("_")[1]
        for panel in json.loads(path.read_text()).get("panels", []):
            rendered.update((gene, iso) for iso in panel.get("final_ids", []))
    decision = json.loads((out / "DECISION_STATE.json").read_text())
    validate_native_inputs(rows, rendered, decision)

    genes = {r["gene_symbol"] for r in rows}
    entries = {x["genes"][0]["geneName"]["value"]: x
               for x in json.loads(reference_json.read_text())["results"]
               if x.get("genes", [{}])[0].get("geneName", {}).get("value") in genes}
    if set(entries) != genes:
        raise ValueError(f"Required domain reference unresolved: {genes - set(entries)}")
    known = {}
    for rec in SeqIO.parse(reference_fasta, "fasta"):
        match = re.search(r" GN=(\S+)", rec.description)
        if match and match.group(1) in genes:
            known.setdefault((match.group(1), str(rec.seq)), []).append(rec.id.split("|")[1])

    clones = read_tsv(data.parent / "external/SuppTable_CloneList.txt")
    clone_index = {}
    for r in clones:
        clone_index.setdefault((r["gene_symbol"], r["aa_seq"].rstrip("*")), []).append(r["clone_id"])
    yang_index = {}
    from .crosswalk.sequence_map import translate_orf
    for r in read_tsv(data / "paper1_orfs.tsv"):
        sequence = translate_orf(r["Isoform__Open_Reading_Frame_Sequence"])
        yang_index.setdefault((r["Gene_Symbol"], sequence), []).append(r["Isoform_ID"])

    # External accessions are additive evidence; the native transcript remains the key.
    from . import config
    ensembl_index = {}
    needed_proteins = {(r["ensg"], r["protein_sequence"]) for r in rows}
    with gzip.open(config.ENSEMBL_PEP_FASTA, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            gene = re.search(r"gene:(ENSG\d+)", rec.description)
            tid = re.search(r"transcript:(ENST\d+\.\d+)", rec.description)
            if gene and tid and (gene.group(1), str(rec.seq)) in needed_proteins:
                ensembl_index.setdefault((gene.group(1), str(rec.seq)), []).append((tid.group(1), rec.id))
    native = {}
    for line in (out / "native_associations.tsv").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 5:
            if fields[3] in native:
                raise ValueError(f"Duplicate native transcript association: {fields[3]}")
            native[fields[3]] = fields
    if set(native) != {r["isoform_id"] for r in rows}:
        raise ValueError("Native transcript association roster differs from required roster")
    for r in rows:
        if native[r["isoform_id"]][0] != r["ensg"]:
            raise ValueError(f"Native gene identity differs: {r['isoform_id']}")

    feature_rows = []
    extra_path = reference_json.parent / "additional_reference_features.json"
    extra = json.loads(extra_path.read_text()) if extra_path.exists() else []
    for r in rows:
        gene, iso, seq = r["gene_symbol"], r["isoform_id"], r["protein_sequence"]
        entry = entries[gene]
        ref = entry["sequence"]["value"]
        mapping = align_mapping(ref, seq)
        r.update(actual_protein_length=len(seq),
                 protein_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                 cds_translation_exact=True, cds_is_cdna_substring=True,
                 length_metadata_agrees=int(r["protein_length"]) == len(seq),
                 uniprot_exact=";".join(known.get((gene, seq), [])),
                 uniprot_domain_reference=entry["primaryAccession"],
                 canonical_reference_aa=len(ref),
                 tfiso_exact_clones=";".join(clone_index.get((gene, seq), [])),
                 yang2016_exact_clones=";".join(yang_index.get((gene, seq), [])))
        r["native_full_exon_structure"] = native[iso][2]
        r["Ensembl116_exact_protein_equivalent_transcripts"] = ";".join(sorted({x[0] for x in ensembl_index.get((r["ensg"], seq), [])}))
        r["Ensembl116_exact_protein_equivalent_ENSP"] = ";".join(sorted({x[1] for x in ensembl_index.get((r["ensg"], seq), [])}))
        displayed = {i for c in entry.get("comments", []) if c["commentType"] == "ALTERNATIVE PRODUCTS"
                     for isoform in c.get("isoforms", []) if isoform.get("isoformSequenceStatus") == "Displayed"
                     for i in isoform["isoformIds"]}
        protein_keys = set(r["uniprot_exact"].split(";")) - {""}
        if entry["primaryAccession"] in protein_keys:
            protein_keys.update(displayed)
        aliases = [x for x in entry.get("uniProtKBCrossReferences", [])
                   if x["database"] in {"Ensembl", "RefSeq"} and x.get("isoformId") in protein_keys]
        r["protein_equivalent_RefSeq"] = ";".join(sorted({x["id"] for x in aliases if x["database"] == "RefSeq"}))
        r["protein_equivalent_Ensembl_transcripts"] = ";".join(sorted({x["id"] for x in aliases if x["database"] == "Ensembl"}))
        r["protein_equivalent_ENSP"] = ";".join(sorted({p["value"] for x in aliases if x["database"] == "Ensembl"
                                                       for p in x["properties"] if p["key"] == "ProteinId"}))
        annotation_features = [(f, "UniProt reviewed feature") for f in entry.get("features", [])]
        annotation_features += [(f, "Literature-defined construct") for f in extra if f["gene_symbol"] == gene
                                and f["reference_accession"] == entry["primaryAccession"]]
        for f, annotation_source in annotation_features:
            if f["type"] not in FEATURE_TYPES:
                continue
            loc = f["location"]
            if loc["start"].get("modifier") != "EXACT" or loc["end"].get("modifier") != "EXACT":
                continue
            start, end = loc["start"]["value"], loc["end"]["value"]
            if not 1 <= start <= end <= len(ref):
                raise ValueError(f"Invalid feature interval: {gene} {f}")
            feature_rows.append(dict(gene_symbol=gene, isoform_id=iso,
                reference_accession=entry["primaryAccession"], feature_type=f["type"],
                description=f.get("description", ""),
                annotation_source=annotation_source,
                evidence=json.dumps(f.get("evidences", []), separators=(",", ":")),
                **project_feature(ref, seq, start, end, mapping)))

    # Reuse the existing binary, co-tested iso2function implementation unchanged.
    from .network import build
    all_clones = sorted({c for r in rows for c in r["tfiso_exact_clones"].split(";") if c})
    _, edges = build.build_graph(clone_ids=all_clones, data_dir=str(data))
    matched_edges = []
    for r in rows:
        cs = r["tfiso_exact_clones"].split(";")
        for edge in edges.to_dict("records"):
            if edge["source"] in cs:
                matched_edges.append(dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"],
                                          evidence="exact_protein_clone_assay", **edge))
    # Yang 2016 has a separate full tested/negative table; never infer absence from its positives.
    for edge in read_tsv(data / "paper1_ppi.tsv"):
        for r in rows:
            if edge["Isoform_ID"] in r["yang2016_exact_clones"].split(";"):
                matched_edges.append(dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"],
                    evidence="exact_protein_clone_assay", source=edge["Isoform_ID"],
                    target=edge["Interactor_Symbol"] or edge["Interactor_ID"], edge_type="ppi",
                    assay="Yang2016_y2h", detected=edge["detected"], tested=edge["tested"], score=""))
    differences = []
    author_categories = []
    for gene in sorted(genes):
        pair = [r for r in rows if r["gene_symbol"] == gene]
        if len(pair) != 2:
            continue
        a, b = pair
        for ca, cb in itertools.product(a["tfiso_exact_clones"].split(";"), b["tfiso_exact_clones"].split(";")):
            if ca and cb:
                category = build.switch_consequence(ca, cb, data_dir=str(data))
                if category:
                    author_categories.append(dict(category, native_isoform_a=a["isoform_id"], native_isoform_b=b["isoform_id"]))
            for modality in ("ppi", "pdi"):
                base = dict(gene_symbol=gene, isoform_a=a["isoform_id"], isoform_b=b["isoform_id"],
                            clone_a=ca, clone_b=cb, modality=modality,
                            direction="A relative to B; panel order, not a biological switch")
                if ca and cb:
                    # Give the existing function a single-modality table: its tested
                    # denominator otherwise combines PPI and PDI targets even when
                    # edge_type is supplied. No production function or calls are changed.
                    modality_edges = edges[edges["edge_type"] == modality]
                    d = build.differential_interactome(ca, cb, modality_edges, edge_type=modality)
                    differences.append(dict(**base, status="co_tested" if d["co_tested_targets"] else "no_cotested_targets",
                        n_cotested=len(d["co_tested_targets"]), gained=";".join(d["gained"]),
                        lost=";".join(d["lost"]), shared=";".join(d["shared"])))
                else:
                    differences.append(dict(**base, status="no_exact_assayed_clone_for_both",
                                            n_cotested=0, gained="", lost="", shared=""))

    write_tsv(out / "isoform_sequence_accession_audit.tsv", rows)
    write_tsv(out / "reference_mapped_features.tsv", feature_rows)
    write_tsv(out / "measured_interactions_all_tested.tsv", matched_edges)
    write_tsv(out / "measured_pair_differences.tsv", differences)
    baits = read_tsv(data / "pdi_dna_baits.tsv")
    for bait in baits:
        tested_isoforms = sorted({e["isoform_id"] for e in matched_edges
                                  if e["edge_type"] == "pdi" and e["target"] == bait["bait_id"]})
        bait["native_exact_isoforms_tested"] = ";".join(tested_isoforms)
        bait["tested_for_any_native_exact_clone"] = bool(tested_isoforms)
    write_tsv(out / "DNA_baits.tsv", baits)
    m1h = matched_activation(rows, read_tsv(data / "activation_m1h.tsv"))
    write_tsv(out / "measured_activation.tsv", m1h)
    write_tsv(out / "author_functional_categories.tsv", author_categories)
    for field, name in [("protein_sequence", "native_proteins.fasta"), ("cds_sequence", "native_CDS.fasta"),
                        ("cdna_sequence", "native_cDNA.fasta")]:
        with open(out / name, "w") as h:
            for r in rows:
                h.write(f">{r['isoform_id']} gene={r['gene_symbol']} source=native_ISV_run\n{r[field]}\n")
    reference_paths = [source, reference_json, reference_fasta,
        data.parent / "external/SuppTable_CloneList.txt", data / "ppi_y2h.tsv", data / "pdi_ey1h.tsv",
        data / "paper1_orfs.tsv", data / "paper1_ppi.tsv", data / "pdi_dna_baits.tsv",
        data / "activation_m1h.tsv", data / "functional_categories.tsv", out / "native_associations.tsv",
        pathlib.Path(config.ENSEMBL_PEP_FASTA), pathlib.Path(__file__), pathlib.Path(build.__file__)]
    if extra_path.exists():
        reference_paths.append(extra_path)
    manifest = dict(status="ANNOTATION_COMPLETE_ALL_NATIVE_IDENTITIES_RETAINED",
        isoforms=len(rows), genes=len(genes), excluded_isoforms=[], exact_CDS_translation_count=len(rows),
        native_source=str(source),
        source_hashes={str(p): digest(p) for p in reference_paths},
        alignment="Biopython global BLOSUM62, gap open -10, gap extension -0.5; descriptive coordinate projection",
        intact_rule="Entire reference feature sequence contiguous, exact and unique in native protein",
        method_relation="Native ISV ORFs unchanged; existing iso2function binary/co-tested assay method reused unchanged, with each modality table supplied separately so the denominator contains only that assay. New amino-acid feature projection does not use genomic domain_impact coverage as a protein-domain call.",
        cautions=["No evidence of protein translation/stability is created by an ORF prediction.",
          "No novel domains predicted; reference feature nonretention does not prove all activity absent.",
          "Panel order is not an expression switch or young-versus-aged comparison.",
          "Exact protein matches do not establish identical UTR or transcript structure.",
          "HLF metadata reports 144 aa while its CDS translates exactly to the exported 169-aa protein; retain both values."])
    (out / "ANALYSIS_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"isoforms":len(rows), "genes":len(genes), "features":len(feature_rows),
        "TFIso_matched_isoforms":sum(bool(r["tfiso_exact_clones"]) for r in rows),
        "Yang_matched_isoforms":sum(bool(r["yang2016_exact_clones"]) for r in rows),
        "measured_edges":len(matched_edges), "pair_differences":len(differences)}, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "out", "reference-json", "reference-fasta", "data"):
        p.add_argument("--" + name, type=pathlib.Path, required=True)
    a = p.parse_args()
    audit(a.run, a.out, a.reference_json, a.reference_fasta, a.data)


if __name__ == "__main__":
    main()
