"""fastComm isoform-overlay: long-read splicing and isoform evidence for the components of differential edges."""
import json

import pandas as pd

from altanalyze3.components.fastComm import isoform_overlay as io


def _write_events(path, rows):
    cols = ["Gene_Symbol", "Feature", "Event_Type", "N_aged", "N_young", "Mean_aged", "Mean_young",
            "DiffMeans_aged_m_young", "ebayesPval"]
    pd.DataFrame([[g, f, t, na, ny, ma, my, ma - my, p] for g, f, t, na, ny, ma, my, p in rows],
                 columns=cols).to_csv(path, sep="\t", index=False)


def _write_iso(path, rows):
    cols = ["Gene_Symbol", "Feature", "Protein Length", "Predicted NMD", "N_aged", "N_young", "Mean_aged", "Mean_young",
            "DiffMeans_aged_m_young", "ebayesPval"]
    pd.DataFrame([[g, f, L, n, 5, 3, ma, my, ma - my, p] for g, f, L, n, ma, my, p in rows],
                 columns=cols).to_csv(path, sep="\t", index=False)


def test_overlay_flags_components_and_background(tmp_path):
    ev = tmp_path / "events"; ev.mkdir(); iso = tmp_path / "iso"; iso.mkdir()
    # receiver HSC: receptor R1 has an aging event; downstream T1 has an aged-up NMD isoform; R2 tested, unchanged
    _write_events(ev / "young-aged-HSC_stats-annotated.txt",
                  [("R1", "R1:E1|R1:E2", "Cassette Exon", 5, 3, 0.6, 0.3, 0.01),
                   ("R2", "R2:E1|R2:E2", "Cassette Exon", 5, 3, 0.5, 0.48, 0.6)])
    _write_iso(iso / "young-aged-HSC_stats-annotated.txt",
               [("T1", "T1:nov", "80/200", "NMD", 2.0, 0.5, 0.001), ("R2", "R2:ENST", "300/300", "Not-NMD", 1.0, 1.0, 0.9)])
    _write_events(ev / "young-aged-MSC_stats-annotated.txt", [("L1", "L1:E1|L1:E2", "Intron Retention", 5, 3, 0.2, 0.5, 0.02)])
    _write_iso(iso / "young-aged-MSC_stats-annotated.txt", [("L1", "L1:ENST", "100/100", "Not-NMD", 1.0, 1.0, 0.9)])
    det = pd.DataFrame({"gene": ["MSC Fibroblast|L1->R1+R2", "HSC|L2->R2"], "log2fc": [0.3, -0.2], "fdr": [0.01, 0.05],
                        "pval": [0.001, 0.01], "population": ["HSC", "HSC"]})
    dp = tmp_path / "det.tsv"; det.to_csv(dp, sep="\t", index=False)
    cand = pd.DataFrame({"sender_state": ["MSC Fibroblast", "HSC", "HSC"], "receiver_state": ["HSC", "HSC", "HSC"],
                         "ligand": ["L1", "L2", "L3"], "receptor": ["R1+R2", "R2", "R2"]})
    cp = tmp_path / "cand.tsv"; cand.to_csv(cp, sep="\t", index=False)
    lr = pd.DataFrame({"ligand": ["L1", "L2", "L3"], "receptor": ["R1+R2", "R2", "R2"], "pathway": ["PW", "PW", "PW"]})
    lp = tmp_path / "lr.tsv"; lr.to_csv(lp, sep="\t", index=False)
    rm = pd.DataFrame({"T1": [0.9, 0.5], "T2": [0.1, 0.4]}, index=["L1", "PW"])
    rp = tmp_path / "resp.tsv"; rm.to_csv(rp, sep="\t")
    edges, comp, summ = io.overlay(str(dp), str(cp), str(lp), str(rp), str(ev), str(iso), "aged", "young",
                                   str(tmp_path / "out"), alpha=0.1, top_downstream=1, log=lambda *a: None)
    e = edges.set_index("ligand")
    assert e.loc["L1", "receptor_altered"] == "R1" and e.loc["L1", "downstream_altered"] == 1
    assert e.loc["L1", "sender_stromal"] and e.loc["L1", "ligand_tested"] == 0   # MSC Fibroblast is not a long-read state
    assert e.loc["L2", "receptor_altered"] == "" and e.loc["L2", "direction"] == "higher in young"
    t1 = comp[(comp.gene == "T1") & (comp.ligand == "L1")].iloc[0]
    assert t1.isoform_changes == "aged NMD:1" and t1.role == "downstream"
    assert summ["significant_edges"] == 2 and summ["candidate_edges"] == 3
    r = summ["enrichment_vs_candidate_edges"]["receptor|all"]
    assert r["significant_tested"] == 2 and r["significant_altered"] == 1
    assert json.load(open(tmp_path / "out" / "overlay_summary.json"))["stromal_sender_edges"] == 1
