"""GO-Elite on a flat gene-set collection (WikiPathways via the AltAnalyze gene-mapp layout)."""
import logging

import pytest

from altanalyze3.components.goelite.parser import build_flat_gene_sets
from altanalyze3.components.goelite.runner import EnrichmentSettings, GOEliteRunner
from altanalyze3.components.goelite import main as goelite_main


def _mapp_file(tmp_path):
    p = tmp_path / "Ensembl-MAPP.txt"
    rows = ["UID\tSystemCode\tMAPP"]
    for i in range(12):
        rows.append(f"G{i}\tEn\tStatin Pathway:WP430")
    for i in range(8, 20):
        rows.append(f"G{i}\tEn\tWnt Signaling:WP428")
    for i in range(30, 36):
        rows.append(f"G{i}\tEn\tOther:WP999")
    p.write_text("\n".join(rows) + "\n")
    return p


def test_mapp_layout_parses_into_named_flat_terms(tmp_path):
    parsed = build_flat_gene_sets(_mapp_file(tmp_path))
    assert set(parsed.term_to_genes) == {"WP430", "WP428", "WP999"}
    assert parsed.tree.get("WP430").name == "Statin Pathway" and parsed.tree.get("WP430").namespace == "WikiPathways"
    assert len(parsed.term_to_genes["WP430"]) == 12 and parsed.tree.descendants("WP430") == set()


def test_two_column_layout_and_bad_header(tmp_path):
    p = tmp_path / "sets.tsv"
    p.write_text("gene\tset\tname\nA\tS1\tSet one\nB\tS1\tSet one\nC\tS2\tSet two\n")
    parsed = build_flat_gene_sets(p, label="Custom")
    assert parsed.term_to_genes == {"S1": {"A", "B"}, "S2": {"C"}} and parsed.tree.get("S1").name == "Set one"
    bad = tmp_path / "bad.tsv"
    bad.write_text("x\ty\nA\tB\n")
    with pytest.raises(ValueError):
        build_flat_gene_sets(bad)


def test_enrichment_recovers_the_enriched_pathway(tmp_path):
    parsed = build_flat_gene_sets(_mapp_file(tmp_path))
    runner = GOEliteRunner(parsed, logger=logging.getLogger("t"), settings=EnrichmentSettings(min_term_size=3, max_term_size=2000))
    background = [f"G{i}" for i in range(40)]
    results = runner.run([f"G{i}" for i in range(0, 8)], background)
    assert results and results[0].term_id == "WP430" and results[0].overlap == 8


def test_cli_runs_with_gene_sets_and_writes_results(tmp_path):
    mapp = _mapp_file(tmp_path)
    q = tmp_path / "q.txt"; q.write_text("\n".join(f"G{i}" for i in range(0, 8)) + "\n")
    b = tmp_path / "b.txt"; b.write_text("\n".join(f"G{i}" for i in range(40)) + "\n")
    out = tmp_path / "out"
    rc = goelite_main.main(["--gene-sets", str(mapp), "--query", str(q), "--background", str(b), "--allow-custom-background",
                            "--outdir", str(out), "--name", "demo", "--min-term-size", "3"])
    assert rc == 0
    tsv = out / "demo_goelite_results.tsv"
    assert tsv.exists()
    text = tsv.read_text()
    assert "WP430" in text and "Statin Pathway" in text


def test_custom_background_is_refused_without_the_override(tmp_path):
    mapp = _mapp_file(tmp_path)
    q = tmp_path / "q.txt"; q.write_text("\n".join(f"G{i}" for i in range(0, 8)) + "\n")
    b = tmp_path / "b.txt"; b.write_text("\n".join(f"G{i}" for i in range(40)) + "\n")
    with pytest.raises(SystemExit):
        goelite_main.main(["--gene-sets", str(mapp), "--query", str(q), "--background", str(b),
                           "--outdir", str(tmp_path / "o1"), "--min-term-size", "3"])
    rc = goelite_main.main(["--gene-sets", str(mapp), "--query", str(q), "--background", str(b), "--allow-custom-background",
                            "--outdir", str(tmp_path / "o2"), "--name", "demo", "--min-term-size", "3"])
    assert rc == 0 and (tmp_path / "o2" / "demo_goelite_results.tsv").exists()
