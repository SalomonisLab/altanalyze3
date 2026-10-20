"""Audit native proteins against specified assayed proteins without transferring mismatched assays.

Uses the standard query action for measurements; every proposed native pair remains in
an audit table regardless of exact identity, RNA support or interaction differences.
"""
import copy
import hashlib
import json
from collections import Counter
from difflib import SequenceMatcher
from pathlib import Path

from . import config
from .crosswalk.sequence_map import iter_literal_proteins
from .protein_impact import dump, read_support, sha, table, write
from .query import run_query


def sequence_differences(native, assayed):
    """Deterministic literal AA edit blocks; native and assay coordinates are 1-based inclusive.

    Empty fragments represent insertions/deletions. Their start/end columns are
    blank, and the zero-based half-open boundaries retain the exact placement.
    """
    return [dict(change=tag, native_start=i+1 if j > i else "", native_end=j if j > i else "",
                 assayed_start=k+1 if l > k else "", assayed_end=l if l > k else "",
                 native_0based_halfopen=f"{i}:{j}", assayed_0based_halfopen=f"{k}:{l}",
                 native_sequence=native[i:j], assayed_sequence=assayed[k:l])
            for tag,i,j,k,l in SequenceMatcher(None,native,assayed,autojunk=False).get_opcodes() if tag != "equal"]


def run_assay_match(args):
    if not args.pairs or not args.protein_fasta or not args.artifacts_dir:
        raise ValueError("assay-match requires --pairs --protein-fasta --artifacts-dir")
    out = Path(args.artifacts_dir); out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Use an empty assay-match output directory")
    manifest = dict(status="incomplete", method="literal-native-assayed-protein-match-v1",
                    interpretation="An assay on a differing protein is supporting evidence only; it does not directly measure the native protein.")
    dump(out/'assay_match_manifest.json', manifest)
    pairs = table(args.pairs)
    fields = {"gene", "ensg", "native_A", "native_B", "clone_A", "clone_B"}
    if not pairs or any(not fields <= r.keys() or not all(r[k] for k in fields) for r in pairs):
        raise ValueError("Pair TSV requires nonempty gene,ensg,native_A,native_B,clone_A,clone_B")
    keys = [(r['gene'],r['clone_A'],r['clone_B']) for r in pairs]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate requested assay pair")
    clpath = config.require(args.clonelist or config.SUPP_CLONELIST)
    clones = table(clpath)
    if len({r['clone_id'] for r in clones}) != len(clones):
        raise ValueError("Duplicate assayed protein ID")
    ci = {r['clone_id']:r for r in clones}
    needed = {}
    for p in pairs:
        for side in ('A','B'):
            clone, native = p['clone_'+side], p['native_'+side]
            if clone not in ci or ci[clone]['gene_symbol'] != p['gene'] or not ci[clone]['aa_seq']:
                raise ValueError("Unresolved assayed protein/gene: " + clone)
            if native in needed and needed[native] != p['ensg']:
                raise ValueError("Conflicting native gene: " + native)
            needed[native] = p['ensg']
    sequences = {}
    for native,gene,seq in iter_literal_proteins(args.protein_fasta,set(needed)):
        if native in sequences or gene != needed[native] or not seq:
            raise ValueError("Native FASTA identity/gene/sequence mismatch: " + native)
        sequences[native] = seq
    if sequences.keys() != needed.keys():
        raise ValueError("Incomplete required native protein roster: " + repr(needed.keys()-sequences.keys()))
    counts, selection = {}, None
    if args.pseudobulk:
        if not args.metadata:
            raise ValueError("Native read support requires --metadata")
        models = [dict(ensg=g,isoform_id=n) for n,g in needed.items()]
        counts,selection = read_support(args.pseudobulk,args.metadata,args.state,args.condition,args.donor,models)
    # The established query computes assays, including unknown calls and optional
    # partner CPM. Native identity/read checks are attached separately, never
    # weakened to make a mismatched protein qualify as a measured native pair.
    q = copy.copy(args)
    q.action = 'query'; q.artifacts_dir = str(out/'assay_calls')
    q.protein_fasta = None; q.pseudobulk = None
    q.clone = sorted({p['clone_'+s] for p in pairs for s in ('A','B')})
    q.combined_network = False
    run_query(q)
    queried = {(r['gene'],r['clone_A'],r['clone_B']):r for r in table(out/'assay_calls/all_pair_decisions.tsv')}
    if set(queried) != set(keys):
        raise ValueError("Assay query failed to retain the complete requested pair roster")
    native_rows, edits, summaries = [], [], []
    for p in pairs:
        row = dict(p)
        for side in ('A','B'):
            native,clone = p['native_'+side],p['clone_'+side]
            ns,cs = sequences[native],ci[clone]['aa_seq']
            differences = sequence_differences(ns,cs)
            for d in differences:
                edits.append(dict(gene=p['gene'],native=native,assayed_protein=clone,**d))
            native_rows.append(dict(gene=p['gene'],native=native,assayed_protein=clone,
                                    native_aa=len(ns),assayed_aa=len(cs),exact_sequence_match=ns==cs,
                                    native_sha256=hashlib.sha256(ns.encode()).hexdigest(),
                                    assayed_sha256=hashlib.sha256(cs.encode()).hexdigest(),
                                    differences=json.dumps(differences,sort_keys=True,separators=(',',':')),
                                    assigned_reads=counts.get(p['ensg']+':'+native,'')))
            row['exact_native_'+side] = ns==cs
            row['native_aa_'+side],row['assayed_aa_'+side] = len(ns),len(cs)
            row['native_reads_'+side] = counts.get(p['ensg']+':'+native,'')
        row['both_native_proteins_exact'] = row['exact_native_A'] and row['exact_native_B']
        row['assay_applies_to_native_pair'] = 'exact_protein_identity' if row['both_native_proteins_exact'] else 'related_protein_sequences_only'
        row['both_native_read_threshold'] = all(row['native_reads_'+s] >= args.min_isoform_reads for s in ('A','B')) if counts else ''
        calls = queried[p['gene'],p['clone_A'],p['clone_B']]
        for k in ('M1H_A','M1H_B','qualifying_PPI','qualifying_PDI'):
            row[k] = calls[k]
        summaries.append(row)
    targets = table(out/'assay_calls/all_pair_target_calls.tsv')
    activation_path = Path(args.out_dir or config.DATA_DIR)/'activation_m1h.tsv'
    activation = {r['clone_id']:r for r in table(activation_path)} if activation_path.exists() else {}
    for row in summaries:
        for side in ('A','B'):
            row['M1H_'+side+'_call'] = activation.get(row['clone_'+side],{}).get('m1h_call','unmeasured')
        for modality in ('ppi','pdi'):
            for call in ('A_only','B_only','shared_positive','both_negative','unknown'):
                values = [r['target'] for r in targets if r['gene']==row['gene'] and r['clone_A']==row['clone_A']
                          and r['clone_B']==row['clone_B'] and r['modality']==modality and r['comparison']==call]
                row[modality+'_'+call] = len(values)
                row[modality+'_'+call+'_targets'] = ';'.join(sorted(values))
            row[modality+'_qualifying_difference_targets'] = ';'.join(sorted(r['target'] for r in targets
                if r['gene']==row['gene'] and r['clone_A']==row['clone_A'] and r['clone_B']==row['clone_B']
                and r['modality']==modality and r['qualifying_difference']=='True'))
    write(out/'native_assayed_protein_comparison.tsv',native_rows)
    write(out/'amino_acid_differences.tsv',edits, list(edits[0]) if edits else ['gene','native','assayed_protein','change'])
    write(out/'all_native_pair_evidence.tsv',summaries)
    files = [Path(args.pairs),Path(args.protein_fasta),Path(clpath)]
    if args.pseudobulk: files.extend([Path(args.pseudobulk),Path(args.metadata)])
    manifest.update(status='complete',pairs=len(pairs),native_proteins=len(needed),read_selection=selection,
                    exact_native_pairs=sum(r['both_native_proteins_exact'] for r in summaries),
                    inputs={str(p):sha(p) for p in files},
                    assay_query_manifest_sha256=sha(out/'assay_calls/query_manifest.json'),
                    code_sha256={str(Path(__file__)):sha(__file__)},
                    outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='assay_match_manifest.json'})
    dump(out/'assay_match_manifest.json',manifest)
    return 0
