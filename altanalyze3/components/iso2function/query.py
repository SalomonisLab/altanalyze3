"""Reusable measured-interaction queries for the standard iso2function CLI.

All within-gene clone pairs are considered, including those not classified in Data_S8.
Default native RNA support requires literal, same-gene protein identity. An explicit
substitution allowance retains and labels nonidentical same-length matches as putative;
aliases are retained independently. Protein partners can be filtered on native gene CPM. DNA baits
remain DNA baits, and never receive fabricated gene-expression values.
"""

import csv
import gzip
import hashlib
import itertools
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from . import config
from .network.build import assay_call, build_graph, differential_interactome

logger = logging.getLogger(__name__)


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _table(path):
    # Native aliases can legitimately exceed the csv module's 128 KiB default.
    # Retain the complete field rather than truncate required identifiers.
    previous = csv.field_size_limit()
    try:
        csv.field_size_limit(sys.maxsize)
        with _open(path) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    finally:
        csv.field_size_limit(previous)


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _roster_sha(values):
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def _write(path, rows, columns):
    with open(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def symbol_index(annotation, mane, extra=None):
    """Additive, explicit symbol mappings. Ambiguity is retained for the caller to resolve."""
    index = defaultdict(set)
    with _open(annotation) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 1 and fields[0].startswith("ENSG") and fields[1]:
                index[fields[1]].add(fields[0].split(".")[0])
    for row in _table(mane):
        index[row["symbol"]].add(row["Ensembl_Gene"].split(".")[0])
    if extra:
        for row in _table(extra):
            if not all(row.get(k) for k in ("symbol", "ensg", "source")):
                raise ValueError("--symbol-map requires nonempty symbol, ensg and source columns")
            index[row["symbol"]].add(row["ensg"].split(".")[0])
    return index


def selected_metadata(path, conditions=None, donors=None):
    rows = _table(path)
    if not rows or not {"library", "groups", "original_library"} <= rows[0].keys():
        raise ValueError("Metadata must contain library, groups, original_library")
    if len({r["library"] for r in rows}) != len(rows):
        raise ValueError("Duplicate source library in metadata")
    if conditions and set(conditions) - {r["groups"] for r in rows}:
        raise ValueError("Requested condition missing from metadata")
    if donors and set(donors) - {r["original_library"] for r in rows}:
        raise ValueError("Requested donor library missing from metadata")
    selected = [r for r in rows if (not conditions or r["groups"] in conditions)
                and (not donors or r["original_library"] in donors)]
    if not selected:
        raise ValueError("No selected libraries")
    # Source-qualified isoform libraries remain separate for read support. A native gene
    # library is counted once in gene-level means, even if ISV has multiple source records.
    for library in {r["original_library"] for r in selected}:
        signatures = {(r["groups"], r.get("original_uid", library)) for r in selected
                      if r["original_library"] == library}
        if len(signatures) != 1:
            raise ValueError("Inconsistent duplicate library provenance: " + library)
    return rows, selected


def _axis(h, name):
    group = h[name]
    values = group[group.attrs.get("_index", "_index")].asstr()[:].tolist()
    if len(values) != len(set(values)):
        raise ValueError("Duplicate h5ad " + name + " identities")
    return values


def _row(h, i, width):
    x = h["X"]
    if hasattr(x, "keys"):
        if x.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError("Sparse query input must use CSR encoding")
        start, end = map(int, x["indptr"][i:i + 2])
        ix, values = x["indices"][start:end], x["data"][start:end]
        if len(ix) != len(set(ix.tolist())) or (ix < 0).any() or (ix >= width).any():
            raise ValueError("Duplicate or invalid sparse feature indices")
        out = np.zeros(width, dtype=float)
        out[ix] = values
    else:
        out = np.asarray(x[i], dtype=float)
    if len(out) != width or not np.isfinite(out).all() or (out < 0).any():
        raise ValueError("Invalid/nonfinite/negative matrix values")
    return out


def _selection(obs, selected, states, suffix, library_column):
    expected_libraries = {r[library_column] for r in selected}
    wanted = {}
    for i, name in enumerate(obs):
        if not name.endswith(suffix) or "." not in name:
            raise ValueError("Unexpected pseudobulk observation identity: " + name)
        state, library = name[:-len(suffix)].split(".", 1)
        if library in expected_libraries and (not states or state in states):
            key = state, library
            if key in wanted:
                raise ValueError("Duplicate state/library observation")
            wanted[key] = i
    states = list(states or sorted({s for s, _ in wanted}))
    missing = {(s, lib) for s in states for lib in expected_libraries} - set(wanted)
    if not states or missing:
        raise ValueError("Incomplete required state/library roster: " + repr(sorted(missing)))
    return wanted, states


def protein_matches(sequence, gene, exact_lookup, length_lookup, max_substitutions=0):
    """Exact identity first; optional same-gene/length substitutions, never indels/stops.

    Every compatible assay ID is retained, with literal residue differences.
    Substitution matches transfer support for a hypothesis, not a native assay.
    """
    exact = exact_lookup.get((gene, sequence), [])
    if exact:
        return {c: dict(match_type="exact", substitutions=[]) for c in exact}
    matches = {}
    if max_substitutions and "*" not in sequence:
        for c, target in length_lookup.get((gene, len(sequence)), []):
            if "*" in target:
                continue
            changes = []
            for i, (native, assayed) in enumerate(zip(sequence, target), 1):
                if native != assayed:
                    changes.append(dict(position=i, native=native, assayed=assayed))
                    if len(changes) > max_substitutions:
                        break
            if 0 < len(changes) <= max_substitutions:
                matches[c] = dict(match_type="AA_substitutions", substitutions=changes)
    return matches


def native_matches(counts, fasta, selected, states, clones, genes, max_substitutions=0):
    """Audit every native feature in the assay genes; retain all literal protein aliases."""
    import h5py
    from .crosswalk.sequence_map import iter_literal_proteins
    gene_set = set(genes.values())
    with h5py.File(counts, "r") as h:
        obs, var = _axis(h, "obs"), _axis(h, "var")
        roster, states = _selection(obs, selected, states, "-isoform", "library")
        totals = np.zeros(len(var))
        for i in roster.values():
            row = _row(h, i, len(var))
            if not np.equal(row, np.rint(row)).all():
                raise ValueError("Isoform read support requires counts, not CPM")
            totals += row
        records = {}
        for i, feature in enumerate(var):
            if ":" not in feature:
                raise ValueError("Native feature must be gene:isoform: " + feature)
            gene, isoform = feature.split(":", 1)
            if gene in gene_set:
                if isoform in records:
                    raise ValueError("Ambiguous native FASTA identity: " + isoform)
                records[isoform] = dict(ensg=gene, isoform=isoform, reads=float(totals[i]),
                                        protein_status="missing", aa="", protein_sha256="", clones="")
        provenance = dict(obs_sha256=_roster_sha(obs), var_sha256=_roster_sha(var),
                          selected_observations=[obs[i] for i in roster.values()], states=states,
                          missing_assay_genes=sorted(gene_set - {r["ensg"] for r in records.values()}))
    lookup = defaultdict(list)
    length_lookup = defaultdict(list)
    for clone in clones:
        lookup[genes[clone["clone_id"]], clone["aa_seq"]].append(clone["clone_id"])
        length_lookup[genes[clone["clone_id"]], len(clone["aa_seq"])].append((clone["clone_id"], clone["aa_seq"]))
    if max_substitutions:
        for rec in records.values():
            rec["match_details"] = "{}"
    seen = set()
    for name, header_gene, seq in iter_literal_proteins(fasta, set(records)):
        if name in seen:
            raise ValueError("Duplicate native protein record: " + name)
        seen.add(name)
        rec = records[name]
        rec.update(aa=len(seq), protein_sha256=hashlib.sha256(seq.encode()).hexdigest())
        rec["protein_status"] = "verified" if header_gene == rec["ensg"] else "gene_header_mismatch"
        if rec["protein_status"] == "verified":
            matches = protein_matches(seq, rec["ensg"], lookup, length_lookup, max_substitutions)
            rec["clones"] = ";".join(sorted(matches))
            if max_substitutions:
                rec["match_details"] = json.dumps(matches, sort_keys=True, separators=(",", ":"))
    return list(records.values()), provenance


def partner_expression(path, counts, selected, states, targets, index):
    """Native gene CPM; state means include zero-expression donors, never missing inputs."""
    import h5py
    records, mapping = [], {}
    with h5py.File(path, "r") as h:
        obs, var = _axis(h, "obs"), _axis(h, "var")
        if any(":" in g for g in var):
            raise ValueError("--gene-cpm must be gene-level, not isoform expression")
        roster, states = _selection(obs, selected, states, "-gene", "original_library")
        vi = {g: i for i, g in enumerate(var)}
        for symbol in sorted(targets):
            candidates = sorted(index.get(symbol, set()))
            available = sorted(set(candidates) & set(var))
            status = ("unresolved_symbol" if not candidates else "ambiguous" if len(available) > 1
                      else "quantified" if available else "not_in_gene_panel")
            mapping[symbol] = dict(target=symbol, candidate_ensg=";".join(candidates),
                                   ensg=available[0] if len(available) == 1 else "", status=status)
        raw = h5py.File(counts, "r") if counts else None
        try:
            if raw and (_axis(raw, "obs") != obs or _axis(raw, "var") != var):
                raise ValueError("Gene counts/CPM rosters differ")
            by_library = {r["original_library"]: r for r in selected}
            for (state, library), i in sorted(roster.items()):
                row = _row(h, i, len(var))
                if raw:
                    rv = _row(raw, i, len(var))
                    if rv.sum() <= 0 or not np.equal(rv, np.rint(rv)).all():
                        raise ValueError("Invalid gene counts")
                    if not np.allclose(row, rv / rv.sum() * 1e6, rtol=1e-5, atol=1e-4):
                        raise ValueError("Native gene CPM disagrees with counts; no renormalization performed")
                for gene in sorted({r["ensg"] for r in mapping.values() if r["ensg"]}):
                    records.append(dict(state=state, library=library,
                                        donor=by_library[library].get("original_uid", library),
                                        ensg=gene, cpm=float(row[vi[gene]])))
        finally:
            if raw:
                raw.close()
        provenance = dict(obs_sha256=_roster_sha(obs), var_sha256=_roster_sha(var),
                          selected_observations=[obs[i] for i in roster.values()], states=states,
                          donors=len(by_library), cpm_verified_against_counts=bool(counts))
    means = defaultdict(list)
    for r in records:
        means[r["ensg"], r["state"]].append(r["cpm"])
    for rec in mapping.values():
        sm = {state: float(np.mean(means[rec["ensg"], state])) for state in states} if rec["ensg"] else {}
        rec["state_mean_cpm"] = json.dumps(sm, sort_keys=True)
        rec["max_state_mean_cpm"] = max(sm.values()) if sm else None
    return list(mapping.values()), records, provenance


def pair_targets(a, b, edges, modality):
    """Include every tested, negative and unknown target alongside the shared differential API."""
    sub = edges[edges.edge_type.eq(modality) & edges.source.isin([a, b])]
    calls = {}
    for (source, target), block in sub.groupby(["source", "target"], sort=True):
        known = {assay_call(r.detected) for r in block.itertuples()
                 if assay_call(r.tested) is True and assay_call(r.detected) is not None}
        if len(known) > 1:
            raise ValueError("Conflicting measured calls: " + repr((source, target)))
        calls[source, target] = next(iter(known)) if known else None
    diff = differential_interactome(a, b, sub, edge_type=modality)
    rows = []
    for target in sorted(set(sub.target)):
        ca, cb = calls.get((a, target)), calls.get((b, target))
        status = ("A_only" if target in diff["gained"] else "B_only" if target in diff["lost"]
                  else "shared_positive" if target in diff["shared"] else
                  "both_negative" if ca is False and cb is False else "unknown")
        rows.append(dict(clone_A=a, clone_B=b, modality=modality, target=target,
                         call_A=ca, call_B=cb, comparison=status))
    return rows


def run_query(args):
    """CLI query driver. Audits are written before any filtered scientific result."""
    out = Path(args.artifacts_dir or config.ARTIFACTS_DIR) / "query" if not args.artifacts_dir else Path(args.artifacts_dir)
    if any((out / name).exists() for name in ("query_manifest.json", "retained_pairs.tsv", "all_pair_decisions.tsv")):
        raise FileExistsError("Use a fresh query output directory; existing results must remain reviewable: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    manifest = dict(status="incomplete", parameters=vars(args).copy(), sources={}, coverage={})
    manifest["code_sha256"] = {str(p.relative_to(Path(__file__).parent)): _sha(p) for p in
                               [Path(__file__), Path(__file__).with_name("cli.py"),
                                Path(__file__).parent / "network" / "build.py",
                                Path(__file__).parent / "network" / "cytoscape.py",
                                Path(__file__).parent / "network" / "query_network.py",
                                Path(__file__).parent / "crosswalk" / "sequence_map.py"]}
    manifest["parameters"].pop("func", None)
    def source(path):
        resolved = config.require(str(path))
        manifest["sources"][str(Path(resolved).absolute())] = _sha(resolved)
        return resolved
    def write_manifest():
        (out / "query_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    try:
        max_substitutions = getattr(args, "max_aa_substitutions", 0)
        if max_substitutions < 0 or args.min_isoform_reads < 0 or args.min_differential_partners < 1:
            raise ValueError("Read threshold must be nonnegative; minimum differences must be >=1")
        if (args.min_partner_cpm is not None and args.min_partner_cpm < 0):
            raise ValueError("CPM threshold must be nonnegative")
        native = bool(args.pseudobulk or args.protein_fasta)
        if native and not all((args.pseudobulk, args.protein_fasta, args.metadata)):
            raise ValueError("Native query requires --pseudobulk, --protein-fasta, --metadata")
        if args.min_partner_cpm is not None and not (args.gene_cpm and args.metadata):
            raise ValueError("Expression filtering requires --gene-cpm and --metadata")
        if args.gene_cpm and not args.metadata:
            raise ValueError("--gene-cpm requires --metadata for the complete donor roster")
        if (args.condition or args.donor or args.state) and not args.metadata:
            raise ValueError("State/condition/donor selection requires --metadata")
        if args.interaction_type == "pdi" and args.min_partner_cpm is not None:
            logger.warning("Gene CPM does not filter DNA baits; PDI is assayed DNA binding, not target regulation")
        clones = _table(source(args.clonelist or config.SUPP_CLONELIST))
        if not clones or len({r["clone_id"] for r in clones}) != len(clones):
            raise ValueError("Empty or duplicate assayed clone catalogue")
        if any(not r.get("aa_seq") for r in clones):
            raise ValueError("Missing assayed clone protein sequence")
        if args.gene:
            clones = [r for r in clones if r["gene_symbol"] == args.gene]
            if not clones:
                raise ValueError("Gene absent from assayed clone catalogue")
        if args.clone:
            missing = set(args.clone) - {r["clone_id"] for r in clones}
            if missing:
                raise ValueError("Unresolved requested clones: " + repr(sorted(missing)))
            clones = [r for r in clones if r["clone_id"] in args.clone]
        index, mane_path = None, None
        if native or args.gene_cpm:
            mane_path = source(args.mane_summary or config.MANE_SUMMARY)
            index = symbol_index(source(args.gene_annotation or config.ENSEMBL91_ANNOT),
                                 mane_path,
                                 source(args.symbol_map) if args.symbol_map else None)
        selected, states = None, args.state
        if args.metadata:
            full_meta, selected = selected_metadata(source(args.metadata), args.condition, args.donor)
            selected_ids = {r["library"] for r in selected}
            _write(out / "metadata_audit.tsv", [dict(r, selected=r["library"] in selected_ids) for r in full_meta],
                   list(full_meta[0]) + ["selected"])
        native_rows, aliases = [], defaultdict(list)
        if native:
            genes, clone_audit = {}, []
            # MANE supplies the gene locus only, never the identity of an alternative
            # protein. Old annotation builds can retain retired duplicate loci.
            mane_genes = defaultdict(set)
            for row in _table(mane_path):
                mane_genes[row["symbol"]].add(row["Ensembl_Gene"].split(".")[0])
            for r in clones:
                ids = sorted(index.get(r["gene_symbol"], set()))
                preferred = sorted(mane_genes.get(r["gene_symbol"], set()))
                resolved = preferred if len(preferred) == 1 else ids
                if len(resolved) == 1:
                    genes[r["clone_id"]] = resolved[0]
                clone_audit.append(dict(clone=r["clone_id"], gene=r["gene_symbol"], candidate_ensg=";".join(ids),
                                        resolved_ensg=genes.get(r["clone_id"], ""),
                                        source="MANE gene locus" if len(preferred) == 1 else "annotation symbol"))
            _write(out / "clone_gene_audit.tsv", clone_audit, ["clone", "gene", "candidate_ensg", "resolved_ensg", "source"])
            if len(genes) != len(clones):
                raise ValueError("Unresolved/ambiguous assay gene mapping; see clone_gene_audit.tsv")
            native_rows, manifest["coverage"]["native"] = native_matches(
                source(args.pseudobulk), source(args.protein_fasta), selected, states, clones, genes, max_substitutions)
            for r in native_rows:
                for clone in r["clones"].split(";"):
                    if clone:
                        alias = {"isoform": r["isoform"], "reads": r["reads"]}
                        if max_substitutions:
                            alias.update(json.loads(r["match_details"])[clone], protein_sha256=r["protein_sha256"])
                        aliases[clone].append(alias)
            _write(out / "native_isoform_audit.tsv", native_rows,
                   ["ensg", "isoform", "reads", "protein_status", "aa", "protein_sha256", "clones"]
                   + (["match_details"] if max_substitutions else []))
            unresolved = [r for r in native_rows if r["reads"] >= args.min_isoform_reads
                          and r["protein_status"] != "verified"]
            if unresolved:
                raise ValueError("Supported native protein sources unresolved; see native_isoform_audit.tsv")
        by_gene = defaultdict(list)
        catalogue = {r["clone_id"]: r for r in clones}
        for r in clones:
            by_gene[r["gene_symbol"]].append(r["clone_id"])
        pairs = [(g, a, b) for g, cs in sorted(by_gene.items())
                 for a, b in itertools.combinations(sorted(cs), 2)]
        if args.pairs:
            supplied = _table(source(args.pairs))
            pairs = [(r["gene"], r.get("A_clone") or r.get("clone_A"),
                      r.get("B_clone") or r.get("clone_B")) for r in supplied]
            if len({tuple(sorted((a, b))) for _, a, b in pairs}) != len(pairs):
                raise ValueError("Duplicate query pairs")
            if any(a not in catalogue or b not in catalogue or a == b or
                   catalogue[a]["gene_symbol"] != g or catalogue[b]["gene_symbol"] != g for g, a, b in pairs):
                raise ValueError("Pair identity does not match assayed clone catalogue")
        data_dir = args.out_dir or config.DATA_DIR
        modalities = ["ppi", "pdi"] if args.interaction_type == "all" else [args.interaction_type]
        for modality in modalities:
            source(Path(data_dir) / ("ppi_y2h.tsv" if modality == "ppi" else "pdi_ey1h.tsv"))
        _, edges = build_graph(clone_ids=catalogue, data_dir=data_dir, crosswalk_structure=pd.DataFrame())
        if "pdi" in modalities:
            bait_path = Path(data_dir) / "pdi_dna_baits.tsv"
            # A small sequence inventory preserves physical bait identity, including mutants.
            # Optional for older ingests; missing annotation is explicit, not a gene substitute.
            baits = {r["bait_id"]: r.get("seq", "") for r in _table(source(bait_path))} if bait_path.exists() else {}
            _write(out / "DNA_bait_audit.tsv", [dict(bait_id=b, sequence=baits.get(b, ""),
                    annotation_status="sequence_available" if baits.get(b) else "sequence_unavailable")
                    for b in sorted(set(edges.loc[edges.edge_type.eq("pdi"), "target"]))],
                   ["bait_id", "sequence", "annotation_status"])
        calls = []
        calls_by_pair = defaultdict(list)
        for gene, a, b in pairs:
            for modality in modalities:
                rows = [dict(gene=gene, **r) for r in pair_targets(a, b, edges, modality)]
                calls.extend(rows)
                calls_by_pair[a, b].extend(rows)
        def native_support(a, b):
            if not native:
                return True
            if max_substitutions:
                return any(x["reads"] >= args.min_isoform_reads and y["reads"] >= args.min_isoform_reads
                           and x["isoform"] != y["isoform"] and x["protein_sha256"] != y["protein_sha256"]
                           for x in aliases[a] for y in aliases[b])
            return (max((r["reads"] for r in aliases[a]), default=0) >= args.min_isoform_reads and
                    max((r["reads"] for r in aliases[b]), default=0) >= args.min_isoform_reads)
        supported = {(a, b) for _, a, b in pairs if native_support(a,b)}
        ppi_diff_pairs = {(r["clone_A"], r["clone_B"]) for r in calls
                         if r["modality"] == "ppi" and r["comparison"] in ("A_only", "B_only")
                         and (r["clone_A"], r["clone_B"]) in supported}
        expression = {}
        if args.gene_cpm:
            mapping, records, manifest["coverage"]["expression"] = partner_expression(
                source(args.gene_cpm), source(args.gene_counts) if args.gene_counts else None,
                selected, states, {r["target"] for r in calls if r["modality"] == "ppi"
                                   and (r["clone_A"], r["clone_B"]) in ppi_diff_pairs}, index)
            _write(out / "partner_expression_audit.tsv", mapping,
                   ["target", "candidate_ensg", "ensg", "status", "state_mean_cpm", "max_state_mean_cpm"])
            _write(out / "gene_expression_donor_state.tsv", records, ["state", "library", "donor", "ensg", "cpm"])
            expression = {r["target"]: r for r in mapping}
            if any(r["status"] in ("unresolved_symbol", "ambiguous") for r in mapping):
                raise ValueError("Unresolved partner mapping; see partner_expression_audit.tsv")
        for r in calls:
            ex = expression.get(r["target"], {}) if r["modality"] == "ppi" else {}
            r["partner_ensg"] = ex.get("ensg", "")
            r["max_state_mean_cpm"] = ex.get("max_state_mean_cpm")
            r["expression_status"] = ex.get("status", "DNA_bait_not_gene_expression" if r["modality"] == "pdi" else
                                             "not_required_for_unsupported_or_unchanged_pair" if args.gene_cpm else "not_requested")
            r["passes_expression"] = (r["modality"] == "pdi" or args.min_partner_cpm is None or
                                      (r["max_state_mean_cpm"] is not None and r["max_state_mean_cpm"] > args.min_partner_cpm))
            r["qualifying_difference"] = r["passes_expression"] and r["comparison"] in ("A_only", "B_only")
        target_columns = ["gene", "clone_A", "clone_B", "modality", "target", "call_A", "call_B", "comparison",
                          "partner_ensg", "max_state_mean_cpm", "expression_status", "passes_expression", "qualifying_difference"]
        _write(out / "all_pair_target_calls.tsv", calls, target_columns)
        activation_path = Path(data_dir) / "activation_m1h.tsv"
        activation = {r["clone_id"]: r for r in _table(source(activation_path))} if activation_path.exists() else {}
        summaries = []
        for gene, a, b in pairs:
            qualifying = [r for r in calls_by_pair[a, b] if r["qualifying_difference"]]
            ma = max((r["reads"] for r in aliases[a]), default=0) if native else None
            mb = max((r["reads"] for r in aliases[b]), default=0) if native else None
            support = native_support(a,b)
            rec = dict(gene=gene, clone_A=a, clone_B=b, aa_A=len(catalogue[a]["aa_seq"]), aa_B=len(catalogue[b]["aa_seq"]),
                       native_aliases_A=json.dumps(aliases[a]), native_aliases_B=json.dumps(aliases[b]),
                       max_single_native_reads_A=ma, max_single_native_reads_B=mb, native_supported=support,
                       M1H_A=activation.get(a, {}).get("M1H_mean", ""), M1H_B=activation.get(b, {}).get("M1H_mean", ""),
                       qualifying_PPI=sum(r["modality"] == "ppi" for r in qualifying),
                       qualifying_PDI=sum(r["modality"] == "pdi" for r in qualifying),
                       retained=support and len(qualifying) >= args.min_differential_partners)
            if max_substitutions:
                for side, clone in (("A", a), ("B", b)):
                    eligible = [r for r in aliases[clone] if r["reads"] >= args.min_isoform_reads]
                    rec["native_evidence_"+side] = ("not_requested" if not native else "exact" if any(
                        r["match_type"] == "exact" for r in eligible) else "AA_substitutions" if eligible else "unsupported")
                rec["native_pair_evidence"] = ("exact_proteins" if rec["native_evidence_A"] == rec["native_evidence_B"] == "exact"
                                               else "putative_transfer_with_AA_substitutions" if support and native else "unsupported")
            summaries.append(rec)
        columns = ["gene", "clone_A", "clone_B", "aa_A", "aa_B", "native_aliases_A", "native_aliases_B",
                   "max_single_native_reads_A", "max_single_native_reads_B", "native_supported", "M1H_A", "M1H_B",
                   "qualifying_PPI", "qualifying_PDI", "retained"]
        if max_substitutions:
            columns += ["native_evidence_A", "native_evidence_B", "native_pair_evidence"]
        _write(out / "all_pair_decisions.tsv", summaries, columns)
        retained = [r for r in summaries if r["retained"]]
        _write(out / "retained_pairs.tsv", retained, columns)
        retained_keys = {(r["clone_A"], r["clone_B"]) for r in retained}
        divergent = [r for r in calls if r["qualifying_difference"] and (r["clone_A"], r["clone_B"]) in retained_keys]
        _write(out / "divergent_interactions.tsv", divergent, target_columns)
        manifest.update(status="complete", pairs_considered=len(summaries), pairs_retained=len(retained),
                        divergent_pair_targets=len(divergent),
                        interpretation=("Measured atlas assays; nonidentical same-gene, same-length protein matches with <= "
                                        + str(max_substitutions) + " AA substitutions support putative functional transfer only. " if max_substitutions else "")
                                       + "Measured atlas assays mapped to native RNA proteins; not assays performed in selected cells. "
                                       "PDI is DNA-bait binding, not demonstrated regulation; gene CPM applies only to PPI partners.")
        if args.combined_network:
            from .network.query_network import export_query_network
            manifest["network"] = export_query_network(divergent, retained, out, manifest)
        write_manifest()
        logger.info("[iso2function.query] %d/%d pairs retained; %d divergent pair-target rows -> %s",
                    len(retained), len(summaries), len(divergent), out)
        return 0
    except Exception as exc:
        manifest.update(status="incomplete", error=str(exc))
        write_manifest()
        raise
