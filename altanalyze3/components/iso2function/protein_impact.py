"""Reproducible native-protein impact assessment using AltAnalyze reference tables.

There is no model or network call here. Every feature instance is retained, including
invalid or unresolved annotations. AltAnalyze2 substring and AltAnalyze3 genomic
coverage baselines are reported separately from 1-based AA verification. These are
sequence annotations and putative impacts, never native-isoform functional assays.
"""
import csv
import gzip
import hashlib
import itertools
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from ..annotation import domain_impact as di
from . import config
from .protein_isoform_audit import align_mapping, project_feature
from .query import _axis, _row, _selection, selected_metadata, _table

METHOD = "altanalyze-native-protein-impact-v1"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def op(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def table(path):
    return _table(path)


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def pinned_bundle(path):
    """An explicit release and SHA are mandatory; analysis never downloads references."""
    path = Path(path)
    b = json.loads(path.read_text())
    if b.get("schema") != METHOD or not b.get("release_decision"):
        raise ValueError("Reference bundle requires schema and recorded release decision")
    resources = {}
    for key, r in b["resources"].items():
        if not all(r.get(x) for x in ("path", "sha256", "source", "release")):
            raise ValueError("Incomplete pinned reference: " + key)
        p = Path(config.resolve_shared(str(path.parent / r["path"])))
        if sha(p) != r["sha256"]:
            raise ValueError("Pinned reference hash differs: " + key)
        resources[key] = p
    required = {"protein_relationships", "protein_features", "uniprot_features",
                "protein_coordinates", "protein_sequences", "ensembl_peptides"}
    if required - resources.keys():
        raise ValueError("Missing built-in reference resources: " + repr(required - resources.keys()))
    return b, resources


def partner_genes(row):
    """Other genes a multi-gene (trans-spliced or read-through) native model joins.

    Optional ';'-separated partner_ensg column. Empty for ordinary single-gene models, whose
    feature selection and outputs are unchanged.
    """
    return [g.strip() for g in str(row.get("partner_ensg") or "").split(";") if g.strip()]


def native_models(source, roster_path):
    contract = json.loads(Path(roster_path).read_text())
    if sha(source) != contract["source_table_sha256"]:
        raise ValueError("Native source hash differs from independent roster contract")
    required = [(r["gene_symbol"], r["isoform_id"]) for r in contract["required_roster"]]
    if len(required) != len(set(required)) or not required:
        raise ValueError("Empty or duplicate required native roster")
    full = table(source)
    keys = [(r["gene_symbol"], r["isoform_id"]) for r in full]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate native source identity")
    index = dict(zip(keys, full))
    if set(required) - index.keys():
        raise ValueError("Incomplete required native roster: " + repr(set(required) - index.keys()))
    rows = [index[k] for k in required]
    for r in rows:
        if not all(r.get(k) for k in ("protein_sequence", "cds_sequence", "cdna_sequence", "ensg")):
            raise ValueError("Missing native sequence or gene: " + r["isoform_id"])
        cds, protein = r["cds_sequence"], r["protein_sequence"]
        if len(cds) % 3 or str(Seq(cds).translate()).removesuffix("*") != protein or "*" in protein:
            raise ValueError("Native CDS/protein mismatch: " + r["isoform_id"])
        if cds not in r["cdna_sequence"]:
            raise ValueError("Native CDS not in cDNA: " + r["isoform_id"])
        partners = partner_genes(r)
        if (any(not re.fullmatch(r"ENSG\d+", g) for g in partners) or r["ensg"] in partners
                or len(partners) != len(set(partners))):
            raise ValueError("Invalid partner_ensg for multi-gene native model: " + r["isoform_id"])
        symbols = [s.strip() for s in str(r.get("partner_gene_symbol") or "").split(";") if s.strip()]
        if symbols and len(symbols) != len(partners):
            raise ValueError("partner_gene_symbol must align with partner_ensg: " + r["isoform_id"])
    audit = [dict(gene_symbol=k[0], isoform_id=k[1],
                  scope="required" if k in set(required) else "outside_independently_defined_roster") for k in keys]
    return rows, audit


def read_support(path, metadata, states, conditions, donors, models):
    import h5py
    import numpy as np
    _, selected = selected_metadata(metadata, conditions, donors)
    with h5py.File(path, "r") as h:
        obs, var = _axis(h, "obs"), _axis(h, "var")
        selection, states = _selection(obs, selected, states, "-isoform", "library")
        vi = {v: i for i, v in enumerate(var)}
        wanted = [r["ensg"] + ":" + r["isoform_id"] for r in models]
        if set(wanted) - vi.keys():
            raise ValueError("Required native models missing from read matrix: " + repr(set(wanted) - vi.keys()))
        totals = np.zeros(len(models))
        indices = [vi[k] for k in wanted]
        for i in selection.values():
            values = _row(h, i, len(var))
            if not np.equal(values, np.rint(values)).all():
                raise ValueError("Read support requires assigned integer counts")
            totals += values[indices]
    return dict(zip(wanted, totals.tolist())), dict(
        states=states, conditions=conditions, donors=donors,
        observations=[obs[i] for i in selection.values()],
        source_libraries=[r["library"] for r in selected],
        obs_sha256=hashlib.sha256(("\n".join(obs) + "\n").encode()).hexdigest(),
        var_sha256=hashlib.sha256(("\n".join(var) + "\n").encode()).hexdigest())


def native_genomic_models(directory, coding_file, models):
    """Stream the established GFF/CDS join without loading millions of unused models."""
    from .crosswalk.sequence_map import iter_literal_proteins
    needed = {r["isoform_id"]: r for r in models}
    fa = Path(directory) / "protein_sequences.fasta.gz"
    verified = set()
    for name, gene, seq in iter_literal_proteins(fa, set(needed)):
        if name in verified or gene != needed[name]["ensg"] or seq != needed[name]["protein_sequence"]:
            raise ValueError("Native protein FASTA identity/sequence differs: " + name)
        verified.add(name)
    if verified != needed.keys():
        raise ValueError("Incomplete native protein FASTA roster: " + repr(needed.keys() - verified))
    exons = defaultdict(list)
    genes = defaultdict(set)
    gff = Path(directory) / "combined.gff.gz"
    with op(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            m = re.search(r'transcript_id[ =]"?([^";]+)', p[8])
            if not m:
                continue
            tid = m[1].split(".")[0] if m[1].startswith("ENST") else m[1]
            if tid in needed:
                gm = re.search(r'gene_id[ =]"?([^";]+)', p[8])
                genes[tid].add(gm[1].split(".")[0] if gm else "")
                exons[tid].append((int(p[3]), int(p[4])))
    coding = {}
    with op(coding_file) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tid = r["Transcript ID"]
            if tid in needed:
                if tid in coding:
                    raise ValueError("Duplicate native coding record: " + tid)
                coding[tid] = r
    if exons.keys() != needed.keys() or coding.keys() != needed.keys():
        raise ValueError("Incomplete required native GFF/CDS roster")
    blocks, audit = {}, []
    for tid, r in needed.items():
        c = coding[tid]
        if genes[tid] != {r["ensg"]} or c["Gene ID"].split(".")[0] != r["ensg"]:
            raise ValueError("Native coding gene mismatch: " + tid)
        blocks[tid] = di.coding_blocks(exons[tid], int(c["CDS Genomic Min"]), int(c["CDS Genomic Max"]))
        audit.append(dict(gene_symbol=r["gene_symbol"], isoform_id=tid,
                          source_table_aa=r["protein_length"], protein_fasta_aa=len(r["protein_sequence"]),
                          coding_regions_aa=c["Protein Length"], source_table_nmd=r["nmd_status"],
                          coding_regions_nmd=c["NMD Status"],
                          cds_nt=len(r["cds_sequence"]), genomic_cds_nt=sum(b-a+1 for a,b in blocks[tid]),
                          native_exon_count=len(exons[tid]),
                          length_agrees=int(c["Protein Length"]) == len(r["protein_sequence"]) == int(r["protein_length"])))
    return blocks, audit, [fa, gff, Path(coding_file)]


def builtin_features(resources, models):
    """Retain every selected-gene source row, including malformed feature coordinates."""
    genes = {r["ensg"] for r in models} | {g for r in models for g in partner_genes(r)}
    relations = defaultdict(set)
    with open(resources["protein_relationships"]) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 3 and p[0] in genes and p[2]:
                relations[p[2]].add(p[0])
    if any(len(g) != 1 for g in relations.values()):
        raise ValueError("Ambiguous built-in protein/gene relationships")
    sequences = {}
    sequence_sources = {}
    with open(resources["protein_sequences"]) as f:
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if p[0] in relations:
                if len(p) != 2 or (p[0] in sequences and sequences[p[0]] != p[1]):
                    raise ValueError("Conflicting reference protein sequence: " + p[0])
                sequences[p[0]] = p[1]
                sequence_sources[p[0]] = "AltAnalyze_builtin_protein_sequences"
    # The distributed exon-comparison sequence database is not the complete
    # proteome. Resolve missing IDs against an explicitly pinned peptide FASTA
    # of the SAME feature release. Never replace a conflicting sequence.
    conflicts = set()
    if "feature_peptides" in resources:
        with op(resources["feature_peptides"]) as f:
            for rec in SeqIO.parse(f, "fasta"):
                protein = rec.id.split(".")[0]
                if protein not in relations:
                    continue
                gene = re.search(r"gene:(ENSG\d+)", rec.description)
                if not gene or gene[1] not in relations[protein]:
                    raise ValueError("Feature peptide gene mismatch: " + protein)
                seq = str(rec.seq)
                if protein in sequences and sequences[protein] != seq:
                    conflicts.add(protein)
                elif protein not in sequences:
                    sequences[protein] = seq
                    sequence_sources[protein] = "same_release_feature_peptides"
    ref_blocks = di.load_protein_coding_blocks(resources["protein_coordinates"])
    rows, counts = [], {}
    for source in ("protein_features", "uniprot_features"):
        n = 0
        with open(resources[source]) as f:
            next(f)
            for lineno, line in enumerate(f, 2):
                n += 1
                p = line.rstrip("\r\n").split("\t")
                if p[0] not in relations:
                    continue
                if len(p) != 8:
                    raise ValueError("Malformed selected-gene feature source row: " + str(lineno))
                r = dict(provider=source, feature_id=source + ":" + str(lineno),
                         ensg=next(iter(relations[p[0]])), reference_protein=p[0],
                         start=p[1], end=p[2], genomic_start=p[3], genomic_end=p[4],
                         name=p[5], accession=p[6], description=p[7],
                         reference_sequence="" if p[0] in conflicts else sequences.get(p[0], ""),
                         reference_sequence_source=sequence_sources.get(p[0], "missing"),
                         reference_sequence_conflict=p[0] in conflicts)
                try:
                    start, end = int(p[1]), int(p[2])
                    gs, ge = sorted((int(p[3]), int(p[4])))
                    intervals = di.intersect(ref_blocks.get(p[0], []), gs, ge)
                    r["domain"] = di.Domain(p[0], p[6], p[5], p[7], start, end, gs, ge, tuple(intervals))
                except ValueError:
                    r["domain"] = None
                rows.append(r)
        counts[source] = dict(total_source_rows=n, selected_gene_rows=sum(r["provider"] == source for r in rows))
    return rows, counts


def supplemental_features(resources, models):
    """Pinned reviewed annotations and literature constructs, explicitly separate from built-ins."""
    if "reviewed_uniprot" not in resources:
        return [], {}
    data = json.loads(resources["reviewed_uniprot"].read_text())
    results = data.get("results", data) if isinstance(data, dict) else data
    # Both the official REST results list and the earlier gene-keyed snapshot are accepted.
    entries = list(results.values()) if isinstance(results, dict) else results
    entries = [e for e in entries if isinstance(e, dict) and "sequence" in e]
    seqs = {e["primaryAccession"]: e["sequence"]["value"] for e in entries}
    gene_index = {r["gene_symbol"]: r["ensg"] for r in models}
    rows = []
    for e in entries:
        gs = [g.get("geneName", {}).get("value") for g in e.get("genes", [])]
        for g in gs:
            if g not in gene_index:
                continue
            for i, f in enumerate(e.get("features", [])):
                location = f.get("location", {})
                a, b = location.get("start", {}), location.get("end", {})
                rows.append(dict(provider="reviewed_uniprot_snapshot", feature_id=e["primaryAccession"] + ":" + str(i),
                                 ensg=gene_index[g], reference_protein=e["primaryAccession"],
                                 start=a.get("value", ""), end=b.get("value", ""),
                                 name=f["type"], accession=e["primaryAccession"], description=f.get("description", ""),
                                 coordinate_modifier=a.get("modifier", "") + "/" + b.get("modifier", ""),
                                 reference_sequence=e["sequence"]["value"],
                                 evidence=json.dumps(f.get("evidences", []), sort_keys=True, separators=(",", ":"))))
    if "curated_regions" in resources:
        for i, f in enumerate(json.loads(resources["curated_regions"].read_text())):
            if f["gene_symbol"] in gene_index:
                rows.append(dict(provider="literature_construct_snapshot", feature_id="curated:" + str(i),
                                 ensg=gene_index[f["gene_symbol"]], reference_protein=f["reference_accession"],
                                 start=f["location"]["start"]["value"], end=f["location"]["end"]["value"],
                                 name=f["type"], accession=f["reference_accession"], description=f["description"],
                                 reference_sequence=seqs.get(f["reference_accession"], ""),
                                 evidence=json.dumps(f["evidences"], sort_keys=True, separators=(",", ":"))))
    return rows, seqs


def project_instance(feature, protein, mapping=None):
    """Do not turn a failed reference lookup or invalid interval into a lost domain."""
    ref = feature["reference_sequence"]
    try:
        a, b = int(feature["start"]), int(feature["end"])
    except (ValueError, TypeError):
        return dict(status="unresolved_coordinates", sequence_present="", legacy_sequence_present="")
    if not ref:
        return dict(status="unresolved_reference_sequence", sequence_present="", legacy_sequence_present="")
    if not 1 <= a <= b <= len(ref):
        return dict(status="invalid_reference_interval", sequence_present="", legacy_sequence_present="")
    if feature.get("coordinate_modifier", "EXACT/EXACT") != "EXACT/EXACT":
        return dict(status="uncertain_reference_interval", sequence_present="", legacy_sequence_present="")
    mapping = mapping if mapping is not None else align_mapping(ref, protein)
    p = project_feature(ref, protein, a, b, mapping)
    frag = di.legacy_feature_fragment(ref, a, b)
    p.update(sequence_present=di.feature_sequence_present(ref[a-1:b], protein),
             legacy_sequence_present=di.feature_sequence_present(frag, protein) if frag else "",
             legacy_fragment_length=len(frag))
    return p


def comparison(a, b):
    """B relative to A; absent from both is never reported as a differential loss."""
    if not isinstance(a, bool) or not isinstance(b, bool):
        return "unresolved"
    return "shared_present" if a and b else "absent_from_both" if not a and not b else "lost_in_B" if a else "gained_in_B"


def novelty(resources, models):
    lookup = defaultdict(list)
    wanted = {r["protein_sequence"] for r in models}
    def add(seq, source, accession, gene="", transcript=""):
        if seq in wanted:
            lookup[seq].append(dict(source=source, accession=accession, ensg=gene, transcript=transcript))
    builtin_gene = defaultdict(set)
    with open(resources["protein_relationships"]) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 3:
                for accession in (p[1], p[2]):
                    if accession:
                        builtin_gene[accession].add(p[0])
    with open(resources["protein_sequences"]) as f:
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) == 2 and p[1] in wanted:
                gs = builtin_gene.get(p[0].removesuffix("-PEP"), set())
                add(p[1], "altanalyze_builtin_proteins", p[0], next(iter(gs)) if len(gs) == 1 else "")
    for key in ("ensembl_peptides", "feature_peptides", "uniprot_isoforms", "refseq_proteins"):
        if key not in resources:
            continue
        with op(resources[key]) as f:
            for r in SeqIO.parse(f, "fasta"):
                g, t = re.search(r"gene:(ENSG\d+)", r.description), re.search(r"transcript:([^ ]+)", r.description)
                add(str(r.seq), key, r.id, g[1] if g else "", t[1] if t else "")
    if "all_uniprot" in resources:
        gs = {r["gene_symbol"]: r["ensg"] for r in models}
        for r in json.loads(resources["all_uniprot"].read_text())["results"]:
            names = [g.get("geneName", {}).get("value") for g in r.get("genes", [])]
            genes = {gs[n] for n in names if n in gs}
            add(r["sequence"]["value"], "all_uniprot", r["primaryAccession"], next(iter(genes)) if len(genes) == 1 else "")
    rows, matches = [], []
    for r in models:
        hits = [h for h in lookup[r["protein_sequence"]] if not h["ensg"] or h["ensg"] == r["ensg"]]
        for h in hits:
            matches.append(dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"], **h,
                                gene_scope="same_gene" if h["ensg"] else "sequence_identity_gene_not_encoded_in_snapshot_header"))
        rows.append(dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"], actual_aa=len(r["protein_sequence"]),
                         annotated_transcript_in_native_source=r["annotated"], exact_reference_matches=len(hits),
                         protein_catalog_status="known_exact_sequence" if hits else "no_exact_sequence_in_pinned_search",
                         novelty_interpretation="protein_sequence_already_catalogued" if hits else
                         "annotated_transcript_unmatched_exported_peptide" if r["annotated"] == "yes" else
                         "candidate_unannotated_coding_form_not_proven_globally_novel",
                         matching_sources=";".join(sorted({h["source"] for h in hits}))))
    return rows, matches


def run_protein_impact(args):
    required = ("native_models", "required_roster", "reference_bundle", "artifacts_dir")
    if not all(getattr(args, k, None) for k in required):
        raise ValueError("protein-impact requires --native-models --required-roster --reference-bundle --artifacts-dir")
    out = Path(args.artifacts_dir)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Use an empty output directory; previous audits cannot be overwritten")
    manifest = dict(status="incomplete", method=METHOD, direct_functional_evidence=False,
                    interpretation="Putative sequence impacts; neither translation nor altered function established in HSPC",
                    inputs={}, references={}, outputs={})
    dump(out / "impact_manifest.json", manifest)
    bundle, resources = pinned_bundle(args.reference_bundle)
    models, source_scope = native_models(args.native_models, args.required_roster)
    manifest["required_roster"] = [dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"]) for r in models]
    manifest["references"] = bundle
    inputs = [Path(args.native_models), Path(args.required_roster), Path(args.reference_bundle)]
    counts, selection = {}, None
    if args.pseudobulk:
        if not args.metadata:
            raise ValueError("Native read support requires --metadata")
        counts, selection = read_support(args.pseudobulk, args.metadata, args.state, args.condition, args.donor, models)
        inputs.extend([Path(args.pseudobulk), Path(args.metadata)])
    manifest["read_selection"] = selection
    blocks, genomic_audit = {}, []
    if args.gff or args.coding_regions:
        if not args.gff or not args.coding_regions:
            raise ValueError("Genomic baseline requires both --gff directory and --coding-regions")
        blocks, genomic_audit, files = native_genomic_models(args.gff, args.coding_regions, models)
        inputs.extend(files)
    features, source_counts = builtin_features(resources, models)
    supplemental, reviewed_sequences = supplemental_features(resources, models)
    manifest["feature_source_counts"] = source_counts
    # A multi-gene model's partner-gene features are projected onto every model of the same gene symbol,
    # so that each pair is compared over the same feature rows. Without partner_ensg this is {own gene}.
    group_partners = defaultdict(set)
    for m in models:
        group_partners[m["gene_symbol"]].update(partner_genes(m))

    def feature_genes(m):
        return {m["ensg"]} | group_partners[m["gene_symbol"]]

    if any(group_partners.values()):
        manifest["multi_gene_partner_features"] = dict(
            partners={g: sorted(p) for g, p in sorted(group_partners.items()) if p},
            rule="partner_ensg reference features are projected by the same AA alignment and native genomic coverage "
                 "onto every native model of the gene symbol; the join is not assumed to be in frame")
    cache, projected, index = {}, [], {}
    for m in models:
        for f in features + supplemental:
            if f["ensg"] not in feature_genes(m):
                continue
            key = f["reference_sequence"], m["protein_sequence"]
            if key[0] and key not in cache:
                cache[key] = align_mapping(*key)
            p = project_instance(f, m["protein_sequence"], cache.get(key))
            d = f.get("domain")
            coverage = di.domain_coverage(blocks[m["isoform_id"]], d) if blocks and d and d.intervals else ""
            record = dict(gene_symbol=m["gene_symbol"], isoform_id=m["isoform_id"],
                          **{k: v for k,v in f.items() if k not in {"reference_sequence", "domain"}}, **p,
                          genomic_coverage=coverage, genomic_present=coverage >= 1 if coverage != "" else "",
                          reference_sequence_sha256=hashlib.sha256(key[0].encode()).hexdigest() if key[0] else "")
            projected.append(record)
            index[m["isoform_id"], f["feature_id"], f["provider"]] = record
    comparisons, pairs = [], []
    by_gene = defaultdict(list)
    for m in models:
        by_gene[m["gene_symbol"]].append(m)
    for gene, members in by_gene.items():
        # Explicit deterministic orientation: longest protein A, ties by native ID.
        members = sorted(members, key=lambda r: (-len(r["protein_sequence"]), r["isoform_id"]))
        for a, b in itertools.combinations(members, 2):
            calls = Counter()
            for f in features + supplemental:
                if f["ensg"] not in feature_genes(a):
                    continue
                pa = index[a["isoform_id"], f["feature_id"], f["provider"]]
                pb = index[b["isoform_id"], f["feature_id"], f["provider"]]
                call = comparison(pa["sequence_present"], pb["sequence_present"])
                calls[("builtin_" if f in features else "supplemental_") + call] += 1
                comparisons.append(dict(gene_symbol=gene, isoform_A=a["isoform_id"], isoform_B=b["isoform_id"],
                                        provider=f["provider"], feature_id=f["feature_id"], reference_protein=f["reference_protein"],
                                        name=f["name"], accession=f["accession"], description=f["description"],
                                        reference_start=f["start"], reference_end=f["end"],
                                        A_status=pa["status"], B_status=pb["status"],
                                        A_exact_fraction=pa.get("exact_fraction", ""), B_exact_fraction=pb.get("exact_fraction", ""),
                                        AA_comparison=call,
                                        legacy_substring_comparison=comparison(pa["legacy_sequence_present"], pb["legacy_sequence_present"]),
                                        genomic_comparison=comparison(pa["genomic_present"], pb["genomic_present"]),
                                        evidence="sequence_annotation_not_direct_isoform_function"))
            ca, cb = [counts.get(r["ensg"] + ":" + r["isoform_id"], "") for r in (a,b)]
            pairs.append(dict(gene_symbol=gene, isoform_A=a["isoform_id"], isoform_B=b["isoform_id"],
                              A_aa=len(a["protein_sequence"]), B_aa=len(b["protein_sequence"]),
                              A_reads=ca, B_reads=cb, both_at_least_min_reads=ca >= args.min_isoform_reads and cb >= args.min_isoform_reads if ca != "" else "",
                              A_nmd=a["nmd_status"], B_nmd=b["nmd_status"],
                              same_protein=a["protein_sequence"] == b["protein_sequence"],
                              changed_N5=a["protein_sequence"][:5] != b["protein_sequence"][:5],
                              changed_C5=a["protein_sequence"][-5:] != b["protein_sequence"][-5:],
                              **{k: calls[k] for k in [p+c for p in ("builtin_", "supplemental_") for c in
                                                      ("lost_in_B", "gained_in_B", "shared_present", "absent_from_both", "unresolved")]}))
    nov, matches = novelty(resources, models)
    audit = []
    for r in models:
        audit.append(dict(gene_symbol=r["gene_symbol"], isoform_id=r["isoform_id"], ensg=r["ensg"],
                          reported_aa=r["protein_length"], actual_aa=len(r["protein_sequence"]),
                          protein_sha256=hashlib.sha256(r["protein_sequence"].encode()).hexdigest(),
                          cds_translation_exact=True, cds_in_cdna=True,
                          metadata_length_agrees=int(r["protein_length"]) == len(r["protein_sequence"]),
                          assigned_reads=counts.get(r["ensg"] + ":" + r["isoform_id"], ""), nmd=r["nmd_status"],
                          builtin_feature_rows=sum(f["ensg"] in feature_genes(r) for f in features),
                          unresolved_feature_rows=sum(p["isoform_id"] == r["isoform_id"] and p["status"] in
                                                      {"unresolved_coordinates", "invalid_reference_interval", "unresolved_reference_sequence", "uncertain_reference_interval"} for p in projected)))
        if any(group_partners.values()):
            # Carried to protein-exon-map, which labels partner-gene coding exons. Absent for single-gene runs.
            audit[-1].update(partner_ensg=";".join(partner_genes(r)), partner_gene_symbol=r.get("partner_gene_symbol") or "")
    regression = []
    if args.prior_feature_audit:
        inputs.append(Path(args.prior_feature_audit))
        model_index = {r["isoform_id"]: r for r in models}
        for old in table(args.prior_feature_audit):
            m = model_index.get(old["isoform_id"])
            if not m or old["reference_accession"] not in reviewed_sequences:
                raise ValueError("Unresolved earlier feature audit identity/reference")
            ref = reviewed_sequences[old["reference_accession"]]
            f = dict(reference_sequence=ref, start=old["reference_start"], end=old["reference_end"])
            p = project_instance(f, m["protein_sequence"], cache[(ref, m["protein_sequence"])])
            regression.append(dict(gene_symbol=m["gene_symbol"], isoform_id=m["isoform_id"],
                                   reference_accession=old["reference_accession"], feature_type=old["feature_type"], description=old["description"],
                                   start=f["start"], end=f["end"], old_status=old["status"], new_status=p["status"],
                                   status_agrees=old["status"] == p["status"],
                                   old_exact_fraction=old["exact_fraction"], new_exact_fraction=p.get("exact_fraction", "")))
    write(out / "source_roster_audit.tsv", source_scope)
    write(out / "native_model_audit.tsv", audit)
    write(out / "native_genomic_audit.tsv", genomic_audit, list(genomic_audit[0]) if genomic_audit else ["status"])
    write(out / "feature_projections.tsv", projected)
    write(out / "pair_feature_impacts.tsv", comparisons)
    write(out / "pair_summary.tsv", pairs)
    write(out / "protein_catalog_audit.tsv", nov)
    write(out / "exact_reference_matches.tsv", matches)
    write(out / "prior_feature_regression.tsv", regression, list(regression[0]) if regression else ["status"])
    with open(out / "native_proteins.fasta", "w") as f:
        for r in models:
            f.write(">" + r["isoform_id"] + " gene:" + r["ensg"] + "\n" + r["protein_sequence"] + "\n")
    manifest.update(status="complete_with_explicit_annotation_limits", models=len(models), genes=len(by_gene), pairs=len(pairs),
                    min_isoform_reads=args.min_isoform_reads,
                    feature_projection_status_counts=dict(Counter(p["status"] for p in projected)),
                    prior_regression_changes=sum(not r["status_agrees"] for r in regression),
                    methods=dict(AA="1-based inclusive full-feature exact substring; repeated hits retained with ambiguous placement",
                                 legacy="AltAnalyze2 FeatureAlignment slice/short-fragment expansion + ExonAnalyze exact substring",
                                 genomic="AltAnalyze3 domain_impact.domain_coverage on native GFF/CDS; full coding coverage >=1",
                                 partial="existing global BLOSUM62 alignment; gap open -10, extension -0.5; descriptive, not domain recognition",
                                 pair_orientation="descending native protein length, then literal ID; B relative to A",
                                 novelty="literal exact peptide identity in explicitly pinned references only; transcript annotation status separate"),
                    limits=["Feature-row counts are not independent domains or interactions.",
                            "Gene-level reference annotations are projected to every native protein, not assumed to be isoform assays.",
                            "Absence of an exact feature sequence does not prove loss of domain activity; partial and altered calls require interpretation.",
                            "Short motifs and modification sites require positional and biological validation; sequence presence alone cannot establish activity.",
                            "No new domain is assigned to unmatched tails; this workflow transfers existing annotations and audits coverage gaps.",
                            "Reference descriptions are preserved verbatim, including the known final-character truncation in EnsMart100 UniProt tables.",
                            "Low read support, NMD flags, metadata inconsistencies and unresolved feature instances are retained, not excluded.",
                            "Native ORFs have not been demonstrated to be translated or stable proteins by this analysis."])
    manifest["inputs"] = {str(p): sha(p) for p in inputs}
    code = [Path(__file__), Path(di.__file__), Path(__file__).with_name("protein_isoform_audit.py"),
            Path(__file__).with_name("query.py"), Path(__file__).with_name("cli.py")]
    manifest["code_sha256"] = {str(p): sha(p) for p in code}
    import Bio, platform
    manifest["software"] = dict(python=platform.python_version(), biopython=Bio.__version__)
    manifest["outputs"] = {p.name: sha(p) for p in sorted(out.iterdir()) if p.name != "impact_manifest.json"}
    dump(out / "impact_manifest.json", manifest)
    return 0
