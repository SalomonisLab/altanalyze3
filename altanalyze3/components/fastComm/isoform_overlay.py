"""Overlay long-read isoform and splicing changes on differential fastComm communication.

For every significant edge of a fastComm differential (cellHarmony_differential DEG_detailed table, features
'<sender>|<ligand>-><receptor>' per receiver population), the components are
  ligand      the ligand gene(s), checked in the SENDER state
  receptor    each receptor subunit, checked in the RECEIVER state
  downstream  the top response genes of the ligand's NicheNet signature (else the CellChat pathway's), checked in the
              RECEIVER state (fastComm response_signatures.tsv)
A component is altered with age in a state when, in the long-read covariate tables of that state, it has
  event     a splicing event with ebayesPval < rawp, |dPSI| > dpsi and >= min_n donors per group (dPSI-covariate-events)
  isoform   an isoform with ebayesPval < rawp, |log2 fold| > log2(fold) and >= min_n donors per group
            (diff-covariate-isoform), with its protein outcome (protein coding, NMD, same length as the gene's longest
            protein, no protein) from the table's Protein Length and Predicted NMD columns
  switch    a called aged-vs-young isoform switch (isoform_state_switches table), when given
The same evidence is computed for the components of every candidate edge (scored but not necessarily significant),
which is the background for the enrichment of altered components among significant edges.
"""
import glob
import json
import os

import numpy as np
import pandas as pd

STROMAL_PATTERNS = ("MSC", "Fibroblast", "Osteo", "Stromal", "Endotheli", "Adipo", "Chondro", "Pericyte", "Smooth")


def is_stromal(state):
    return any(p.lower() in str(state).lower() for p in STROMAL_PATTERNS)


def _state_tables(folder, case, control):
    out = {}
    for path in sorted(glob.glob(os.path.join(folder, "*_stats-annotated.txt"))):
        stem = os.path.basename(path)[:-len("_stats-annotated.txt")]
        for x, y in ((case, control), (control, case)):
            if stem.startswith(f"{x}-{y}-"):
                out[stem[len(f"{x}-{y}-"):]] = path
                break
    if not out:
        raise FileNotFoundError(f"no {case}-{control}-<state>_stats-annotated.txt tables in {folder}")
    return out


def _oriented(d, case, control, path):
    mc = pd.to_numeric(d[f"Mean_{case}"], errors="coerce")
    mk = pd.to_numeric(d[f"Mean_{control}"], errors="coerce")
    fcol = [c for c in d.columns if c.startswith("DiffMeans_")][0]
    stored = pd.to_numeric(d[fcol], errors="coerce")
    sign = 1.0 if fcol == f"DiffMeans_{case}_m_{control}" else -1.0
    ok = (mc - mk).notna() & stored.notna()
    if not np.allclose((mc - mk)[ok], sign * stored[ok], atol=1e-5):
        raise ValueError(f"{path}: {fcol} does not equal the difference of its Mean columns")
    return mc - mk


def event_evidence(event_dir, case, control, rawp=0.05, dpsi=0.1, min_n=3, log=print):
    """Per state: tested gene symbols and, per altered symbol, events (count, largest |dPSI| with its type)."""
    tested, altered = {}, {}
    for state, path in _state_tables(event_dir, case, control).items():
        d = pd.read_csv(path, sep="\t", usecols=lambda c: c in ("Gene_Symbol", "Feature", "Event_Type", f"N_{case}",
                                                                 f"N_{control}", f"Mean_{case}", f"Mean_{control}",
                                                                 "ebayesPval") or c.startswith("DiffMeans_"),
                        dtype={"Gene_Symbol": str, "Feature": str, "Event_Type": str}, low_memory=False)
        diff = _oriented(d, case, control, path)
        p = pd.to_numeric(d.ebayesPval, errors="coerce")
        ok = (pd.to_numeric(d[f"N_{case}"], errors="coerce") >= min_n) & \
             (pd.to_numeric(d[f"N_{control}"], errors="coerce") >= min_n) & p.notna() & diff.notna()
        tested[state] = set(d.Gene_Symbol[ok].dropna())
        hit = d[ok & (p < rawp) & (diff.abs() > dpsi)].assign(dpsi=diff[ok & (p < rawp) & (diff.abs() > dpsi)])
        per = {}
        for sym, g in hit.groupby("Gene_Symbol"):
            top = g.loc[g.dpsi.abs().idxmax()]
            per[sym] = dict(events=int(len(g)), largest_dpsi=round(float(top.dpsi), 4), largest_type=str(top.Event_Type))
        altered[state] = per
        log(f"[isoform-overlay] events {state}: {len(per):,} altered of {len(tested[state]):,} tested genes")
    return tested, altered


def isoform_evidence(isoform_dir, case, control, rawp=0.05, fold=1.5, min_n=3, log=print):
    """Per state: tested gene symbols and, per altered symbol, isoforms higher in case / control by protein outcome."""
    tested, altered = {}, {}
    cut = float(np.log2(fold))
    for state, path in _state_tables(isoform_dir, case, control).items():
        d = pd.read_csv(path, sep="\t", usecols=lambda c: c in ("Gene_Symbol", "Feature", "Protein Length", "Predicted NMD",
                                                                 f"N_{case}", f"N_{control}", f"Mean_{case}",
                                                                 f"Mean_{control}", "ebayesPval") or c.startswith("DiffMeans_"),
                        dtype={"Gene_Symbol": str, "Feature": str, "Protein Length": str, "Predicted NMD": str},
                        low_memory=False)
        diff = _oriented(d, case, control, path)
        p = pd.to_numeric(d.ebayesPval, errors="coerce")
        ok = (pd.to_numeric(d[f"N_{case}"], errors="coerce") >= min_n) & \
             (pd.to_numeric(d[f"N_{control}"], errors="coerce") >= min_n) & p.notna() & diff.notna()
        tested[state] = set(d.Gene_Symbol[ok].dropna())
        sel = ok & (p < rawp) & (diff.abs() > cut)
        h = d[sel].copy(); h["diff"] = diff[sel]
        L = pd.to_numeric(h["Protein Length"].str.split("/").str[0], errors="coerce")
        M = pd.to_numeric(h["Protein Length"].str.split("/").str[1], errors="coerce")
        nmd = h["Predicted NMD"].fillna("").str.contains("NMD") & ~h["Predicted NMD"].fillna("").str.startswith("Not")
        h["outcome"] = np.where(~(L > 0), "no protein", np.where(L == M, "same length",
                                                                 np.where(nmd, "NMD", "protein coding")))
        h["side"] = np.where(h["diff"] > 0, case, control)
        per = {}
        for sym, g in h.groupby("Gene_Symbol"):
            c = g.groupby(["side", "outcome"]).size()
            per[sym] = {f"{s} {o}": int(v) for (s, o), v in c.items()}
        altered[state] = per
        log(f"[isoform-overlay] isoforms {state}: {len(per):,} altered of {len(tested[state]):,} tested genes")
    return tested, altered


def switch_evidence(switches):
    sw = pd.read_csv(switches, sep="\t")
    out = {}
    for r in sw.itertuples():
        out.setdefault(r.state, set()).add(r.gene_symbol)
    return out


def load_signatures(response_matrix, top_n):
    m = pd.read_csv(response_matrix, sep="\t", index_col=0)
    sig = {}
    for name, row in m.iterrows():
        row = row[row > 0].sort_values(ascending=False)
        sig[str(name)] = list(row.index[:top_n])
    return sig


def edge_components(ligand, receptor, pathway, signatures):
    lig = [g for g in str(ligand).split("+") if g]
    rec = [g for g in str(receptor).split("+") if g]
    down_src = ligand if ligand in signatures else (pathway if pathway in signatures else None)
    down = [g for g in signatures.get(down_src, []) if g not in lig and g not in rec] if down_src else []
    return lig, rec, down, down_src


def _component_rows(edge, role, genes, state, ev_t, ev_a, iso_t, iso_a, sw):
    rows = []
    for g in genes:
        e = ev_a.get(state, {}).get(g)
        i = iso_a.get(state, {}).get(g)
        tested = (g in ev_t.get(state, set())) or (g in iso_t.get(state, set()))
        rows.append(dict(edge, role=role, gene=g, state_checked=state,
                         state_in_long_read=(state in ev_t) or (state in iso_t), tested=tested,
                         event_altered=e is not None, isoform_altered=i is not None,
                         switch=bool(sw and g in sw.get(state, set())),
                         events=e["events"] if e else 0, largest_dpsi=e["largest_dpsi"] if e else np.nan,
                         largest_event_type=e["largest_type"] if e else "",
                         isoform_changes=";".join(f"{k}:{v}" for k, v in sorted(i.items())) if i else ""))
    return rows


def overlay(deg_detailed, candidate_edges, lr_table, response_matrix, event_dir, isoform_dir, case, control,
            outdir, alpha=0.1, use_rawp=False, rawp=0.05, dpsi=0.1, fold=1.5, min_n=3, top_downstream=50,
            switches=None, log=print):
    os.makedirs(outdir, exist_ok=True)
    lr = pd.read_csv(lr_table, sep="\t", dtype=str)
    pathway = {(r.ligand, r.receptor): r.pathway for r in lr.itertuples()}
    signatures = load_signatures(response_matrix, top_downstream)
    ev_t, ev_a = event_evidence(event_dir, case, control, rawp=rawp, dpsi=dpsi, min_n=min_n, log=log)
    iso_t, iso_a = isoform_evidence(isoform_dir, case, control, rawp=rawp, fold=fold, min_n=min_n, log=log)
    sw = switch_evidence(switches) if switches else None

    det = pd.read_csv(deg_detailed, sep="\t")
    metric = det.pval if (use_rawp and "pval" in det.columns) else det.fdr
    det = det[metric < alpha].copy()
    parts = det.gene.str.split("|", n=1)
    det["sender"] = parts.str[0]; det["lr"] = parts.str[1]
    det["ligand"] = det.lr.str.split("->").str[0]; det["receptor"] = det.lr.str.split("->").str[1]
    det["receiver"] = det.population.astype(str)
    det["direction"] = np.where(det.log2fc > 0, f"higher in {case}", f"higher in {control}")
    cand = pd.read_csv(candidate_edges, sep="\t", dtype=str)
    cand = cand.rename(columns={"sender_state": "sender", "receiver_state": "receiver"})
    sig_keys = set(zip(det.sender, det.receiver, det.ligand, det.receptor))

    rows = []
    for src, frame in (("significant", det), ("candidate", cand)):
        for r in frame.drop_duplicates(["sender", "receiver", "ligand", "receptor"]).itertuples():
            pw = pathway.get((r.ligand, r.receptor), "")
            lig, rec, down, down_src = edge_components(r.ligand, r.receptor, pw, signatures)
            edge = dict(edge_set=src, sender=r.sender, receiver=r.receiver, ligand=r.ligand, receptor=r.receptor,
                        pathway=pw, significant=(r.sender, r.receiver, r.ligand, r.receptor) in sig_keys,
                        direction=getattr(r, "direction", ""), log2fc=getattr(r, "log2fc", np.nan),
                        fdr=getattr(r, "fdr", np.nan), sender_stromal=is_stromal(r.sender),
                        receiver_stromal=is_stromal(r.receiver), downstream_source=down_src or "")
            rows += _component_rows(edge, "ligand", lig, r.sender, ev_t, ev_a, iso_t, iso_a, sw)
            rows += _component_rows(edge, "receptor", rec, r.receiver, ev_t, ev_a, iso_t, iso_a, sw)
            rows += _component_rows(edge, "downstream", down, r.receiver, ev_t, ev_a, iso_t, iso_a, sw)
    comp = pd.DataFrame(rows)
    comp["altered"] = comp.event_altered | comp.isoform_altered | comp.switch
    sig_comp = comp[comp.edge_set == "significant"]
    sig_comp.to_csv(os.path.join(outdir, "overlay_components_significant_edges.tsv"), sep="\t", index=False)

    keycols = ["sender", "receiver", "ligand", "receptor"]
    per_edge = []
    for k, g in sig_comp.groupby(keycols, sort=False):
        first = g.iloc[0]
        rec_ = g[g.role == "receptor"]; lig_ = g[g.role == "ligand"]; dn_ = g[g.role == "downstream"]
        per_edge.append(dict(zip(keycols, k), pathway=first.pathway, direction=first.direction, log2fc=first.log2fc,
                             fdr=first.fdr, sender_stromal=first.sender_stromal, receiver_stromal=first.receiver_stromal,
                             ligand_tested=int(lig_.tested.sum()), ligand_altered=";".join(lig_.gene[lig_.altered]),
                             receptor_tested=int(rec_.tested.sum()), receptor_altered=";".join(rec_.gene[rec_.altered]),
                             downstream_tested=int(dn_.tested.sum()), downstream_altered=int(dn_.altered.sum()),
                             downstream_altered_genes=";".join(dn_.gene[dn_.altered][:15])))
    edges = pd.DataFrame(per_edge)
    edges.to_csv(os.path.join(outdir, "overlay_significant_edges.tsv"), sep="\t", index=False)

    # enrichment: altered fraction among tested (state, gene, role) of significant edges vs all candidate edges
    def uniq(frame, role):
        x = frame[(frame.role == role) & frame.tested].drop_duplicates(["state_checked", "gene"])
        return int(x.altered.sum()), int(len(x))
    from scipy import stats
    enr = {}
    for role in ("ligand", "receptor", "downstream"):
        for subset, sel in (("all", slice(None)), ("stromal sender", comp.sender_stromal)):
            c = comp[sel] if subset != "all" else comp
            a_s, n_s = uniq(c[c.edge_set == "significant"], role)
            a_c, n_c = uniq(c[c.edge_set == "candidate"], role)
            orr, p = stats.fisher_exact([[a_s, n_s - a_s], [a_c, n_c - a_c]]) if n_s and n_c else (np.nan, np.nan)
            enr[f"{role}|{subset}"] = dict(significant_altered=a_s, significant_tested=n_s, candidate_altered=a_c,
                                           candidate_tested=n_c, odds_ratio=orr, fisher_p=p)
    summary = dict(significant_edges=int(len(edges)), candidate_edges=int(cand.drop_duplicates(keycols).shape[0]),
                   significant_by_direction=edges.direction.value_counts().to_dict(),
                   stromal_sender_edges=int(edges.sender_stromal.sum()), stromal_receiver_edges=int(edges.receiver_stromal.sum()),
                   edges_receiver_in_long_read=int(sig_comp[sig_comp.role == "receptor"].drop_duplicates(keycols).state_in_long_read.sum()),
                   edges_with_altered_receptor=int((edges.receptor_altered != "").sum()),
                   edges_with_altered_ligand=int((edges.ligand_altered != "").sum()),
                   edges_with_altered_downstream=int((edges.downstream_altered > 0).sum()),
                   enrichment_vs_candidate_edges=enr,
                   parameters=dict(alpha=alpha, use_rawp=use_rawp, rawp=rawp, dpsi=dpsi, fold=fold, min_n=min_n,
                                   top_downstream=top_downstream, case=case, control=control,
                                   deg_detailed=os.path.abspath(deg_detailed), candidate_edges=os.path.abspath(candidate_edges),
                                   lr_table=os.path.abspath(lr_table), response_matrix=os.path.abspath(response_matrix),
                                   event_dir=os.path.abspath(event_dir), isoform_dir=os.path.abspath(isoform_dir),
                                   switches=os.path.abspath(switches) if switches else None))
    with open(os.path.join(outdir, "overlay_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log(json.dumps({k: v for k, v in summary.items() if k != "parameters"}, indent=1, default=str))
    return edges, sig_comp, summary
