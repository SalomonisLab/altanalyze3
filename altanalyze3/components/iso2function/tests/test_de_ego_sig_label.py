"""plot_regulator_states_de legend names the significance value its threshold gates (sig_label)."""
import os, subprocess, shutil
import pytest
from altanalyze3.components.iso2function.network import coexpr, de_ego


@pytest.mark.skipif(shutil.which("pdftotext") is None, reason="pdftotext not installed")
def test_sig_label_appears_in_legend(tmp_path):
    edges = coexpr.load_edges()
    s2g = coexpr.load_symbol_to_ensg(edges)
    isos, targets, _ = de_ego.build_regulator_structure("MAX", {"ENST00000358402", "ENST00000358664"}, edges, s2g,
                                                        key_by_final_isoform=True)
    assert isos and targets
    ensg = next(v[1] for _, _, t in isos for v in t.values())
    ide = {"S1": {"ENST00000358402": (2.0, 0.01)}, "S2": {"ENST00000358402": (0.1, 0.5)}}
    gde = {"S1": {ensg: (2.0, 0.01)}, "S2": {ensg: (0.0, 0.9)}}
    for label in ("FDR", "rawP"):
        out = tmp_path / f"{label}.pdf"
        kw = {} if label == "FDR" else {"sig_label": label}
        de_ego.plot_regulator_states_de("MAX", ["S1", "S2"], str(out), gde, ide, {"ENST00000358402", "ENST00000358664"},
                                        edges=edges, sym2ensg=s2g, sig_thresh=0.05, min_log2fc=1.0, direction="up",
                                        key_by_final_isoform=True, **kw)
        text = subprocess.run(["pdftotext", str(out), "-"], capture_output=True, text=True).stdout
        assert f"{label}<0.05" in text
        other = "rawP" if label == "FDR" else "FDR"
        assert f"{other}<0.05" not in text
