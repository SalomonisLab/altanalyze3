"""Scientific query contracts: measured calls, complete rosters, native identity and CPM."""
import argparse
import json

import h5py
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.iso2function.cli import add_arguments, main
from altanalyze3.components.iso2function.network.build import build_graph, differential_interactome
from altanalyze3.components.iso2function.network.cytoscape import to_cytoscape_elements
from altanalyze3.components.iso2function.query import native_matches, partner_expression, selected_metadata
from altanalyze3.components.iso2function.network.query_network import query_graph


def _h5(path, obs, var, matrix):
    with h5py.File(path, "w") as h:
        for name, values in [("obs", obs), ("var", var)]:
            group = h.create_group(name)
            group.attrs["_index"] = "_index"
            group.create_dataset("_index", data=values, dtype=h5py.string_dtype())
        h.create_dataset("X", data=np.array(matrix, dtype=float))
    return str(path)


def _tsv(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return str(path)


def test_unknown_and_other_modality_do_not_imply_negative(tmp_path):
    _tsv(tmp_path / "ppi_y2h.tsv", [
        dict(ad_clone_id="A", db_gene_symbol="X", detected=True, tested=True),
        dict(ad_clone_id="B", db_gene_symbol="X", detected="", tested=True),
        dict(ad_clone_id="A", db_gene_symbol="Y", detected=True, tested=True),
        dict(ad_clone_id="B", db_gene_symbol="Y", detected=False, tested=False),
        dict(ad_clone_id="A", db_gene_symbol="Z", detected=True, tested=True)])
    _tsv(tmp_path / "pdi_ey1h.tsv", [dict(clone_id="B", bait_id="Z", detected=False, tested=True)])
    nodes, edges = build_graph(data_dir=str(tmp_path))
    result = differential_interactome("A", "B", edges, "ppi")
    assert result["gained"] == []
    assert result["co_tested_targets"] == []
    assert len(edges) == 6
    assert edges[edges.source.eq("B") & edges.edge_type.eq("ppi")].tested.eq(False).all()
    els = to_cytoscape_elements(nodes, edges)
    assert len(els["edges"]) == 3


def test_roster_deduplicates_gene_sources_but_not_isoform_sources(tmp_path):
    p = _tsv(tmp_path / "metadata.tsv", [
        dict(library="source0.D1", original_library="D1", original_uid="one", groups="healthy"),
        dict(library="source1.D1", original_library="D1", original_uid="one", groups="healthy"),
        dict(library="source0.D2", original_library="D2", original_uid="two", groups="healthy")])
    _, selected = selected_metadata(p, ["healthy"])
    cpm = _h5(tmp_path / "genes.h5ad", ["HSC.D1-gene", "HSC.D2-gene"], ["ENSG1"], [[2], [0]])
    mapping, records, provenance = partner_expression(cpm, None, selected, ["HSC"], {"P", "missing"},
                                                     {"P": {"ENSG1"}, "missing": {"ENSG2"}})
    rows = {r["target"]: r for r in mapping}
    assert rows["P"]["max_state_mean_cpm"] == 1  # strict >1 must fail, zero donor retained
    assert rows["missing"]["max_state_mean_cpm"] is None
    assert rows["missing"]["status"] == "not_in_gene_panel"
    assert len(records) == 2 and provenance["donors"] == 2
    with pytest.raises(ValueError, match="Incomplete required"):
        partner_expression(cpm, None, selected, ["HSC", "MPP"], {"P"}, {"P": {"ENSG1"}})


def test_literal_native_aliases_and_no_alias_read_pooling(tmp_path):
    selected = [dict(library="S.D1", original_library="D1", groups="healthy")]
    counts = _h5(tmp_path / "iso.h5ad", ["HSC.S.D1-isoform"],
                 ["ENSG1:novel.1", "ENSG1:ENST2", "ENSG1:with_stop", "ENSG1:wrong_gene"], [[20, 20, 30, 31]])
    fasta = tmp_path / "protein.fa"
    fasta.write_text(">novel.1 ;Protein;gene_id:ENSG1\nMA\n>ENST2 gene:ENSG1\nMA\n"
                     ">with_stop gene:ENSG1\nMA*\n>wrong_gene gene:ENSG2\nMA\n")
    clones = [dict(clone_id="TF-1", aa_seq="MA")]
    rows, _ = native_matches(counts, str(fasta), selected, ["HSC"], clones, {"TF-1": "ENSG1"})
    assert len(rows) == 4
    matched = [r for r in rows if r["clones"]]
    assert {r["isoform"] for r in matched} == {"novel.1", "ENST2"}
    assert max(r["reads"] for r in matched) == 20
    assert rows[-1]["protein_status"] == "gene_header_mismatch"


def test_query_registers_identically_in_shared_parser():
    standalone = argparse.ArgumentParser()
    add_arguments(standalone)
    top = argparse.ArgumentParser()
    top.add_argument("--loglevel", default="INFO")
    add_arguments(top, include_loglevel=False)
    argv = ["query", "--interaction-type", "all", "--clone", "TF-1", "--clone", "TF-2", "--combined-network"]
    assert vars(standalone.parse_args(argv)) == vars(top.parse_args(argv))


def test_end_to_end_ppi_expression_and_pdi_union(tmp_path):
    data, out = tmp_path / "data", tmp_path / "result"
    data.mkdir()
    clones = _tsv(tmp_path / "clones.tsv", [dict(clone_id="TF-1", gene_symbol="TF", aa_seq="MA"),
                                           dict(clone_id="TF-2", gene_symbol="TF", aa_seq="MB")])
    _tsv(data / "ppi_y2h.tsv", [dict(ad_clone_id="TF-1", db_gene_symbol="P", detected=True),
                               dict(ad_clone_id="TF-2", db_gene_symbol="P", detected=False),
                               dict(ad_clone_id="TF-1", db_gene_symbol="Q", detected=True),
                               dict(ad_clone_id="TF-2", db_gene_symbol="Q", detected=False)])
    _tsv(data / "pdi_ey1h.tsv", [dict(clone_id="TF-1", bait_id="MUT_119", detected=False),
                                dict(clone_id="TF-2", bait_id="MUT_119", detected=True)])
    metadata = _tsv(tmp_path / "metadata.tsv", [dict(library="src.D1", original_library="D1", groups="healthy"),
                                               dict(library="src.D2", original_library="D2", groups="healthy")])
    counts = _h5(tmp_path / "iso.h5ad", ["HSC.src.D1-isoform", "HSC.src.D2-isoform"],
                 ["ENSG1:novel.1", "ENSG1:novel.2"], [[15, 20], [15, 10]])
    cpm = _h5(tmp_path / "genes.h5ad", ["HSC.D1-gene", "HSC.D2-gene"], ["ENSG2", "ENSG3"], [[2, 2.2], [0, 0]])
    fasta = tmp_path / "proteins.fa"
    fasta.write_text(">novel.1 gene:ENSG1\nMA\n>novel.2 gene:ENSG1\nMB\n")
    annotation = tmp_path / "genes.tsv"
    annotation.write_text("ENSG1\tTF\nENSG2\tP\nENSG3\tQ\n")
    mane = tmp_path / "mane.tsv"
    mane.write_text("symbol\tEnsembl_Gene\nTF\tENSG1\n")
    code = main(["query", "--out-dir", str(data), "--artifacts-dir", str(out), "--clonelist", clones,
                 "--pseudobulk", counts, "--protein-fasta", str(fasta), "--metadata", metadata,
                 "--state", "HSC", "--condition", "healthy", "--gene-cpm", cpm, "--min-partner-cpm", "1",
                 "--gene-annotation", str(annotation), "--mane-summary", str(mane)])
    assert code == 0
    result = pd.read_csv(out / "retained_pairs.tsv", sep="\t")
    assert len(result) == 1
    assert result.iloc[0].qualifying_PPI == 1
    assert result.iloc[0].qualifying_PDI == 1
    divergent = pd.read_csv(out / "divergent_interactions.tsv", sep="\t")
    assert set(divergent.target) == {"Q", "MUT_119"}
    assert len(pd.read_csv(out / "all_pair_target_calls.tsv", sep="\t")) == 3
    assert json.loads((out / "query_manifest.json").read_text())["status"] == "complete"
    for modality, expected in [("all", {"Q", "MUT_119"}), ("ppi", {"Q"}), ("pdi", {"MUT_119"})]:
        export = tmp_path / (modality + ".txt")
        assert main(["export", "--query-results", str(out), "--out", str(export),
                     "--interaction-type", modality]) == 0
        table = pd.read_csv(export, sep="\t")
        assert set(table.target) == expected
        assert table.native_aliases_A.notna().all()
        assert table.max_single_native_reads_B.eq(30).all()
        assert json.loads((tmp_path / (modality + ".txt.manifest.json")).read_text())["exported_rows"] == len(expected)
    manifest = json.loads((out / "query_manifest.json").read_text())
    manifest["status"] = "incomplete"
    (out / "query_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Only a completed query"):
        main(["export", "--query-results", str(out), "--out", str(tmp_path / "invalid.txt")])


def test_cpm_count_discrepancy_rejects(tmp_path):
    selected = [dict(library="S.D1", original_library="D1", groups="healthy")]
    cpm = _h5(tmp_path / "cpm.h5ad", ["S.D1-gene"], ["ENSG1"], [[10]])
    counts = _h5(tmp_path / "counts.h5ad", ["S.D1-gene"], ["ENSG1"], [[10]])
    with pytest.raises(ValueError, match="disagrees with counts"):
        partner_expression(cpm, counts, selected, ["S"], {"P"}, {"P": {"ENSG1"}})


def test_combined_network_keeps_negative_calls_and_separates_baits_from_genes():
    pair = dict(gene="TF", clone_A="TF-1", clone_B="TF-2", aa_A=100, aa_B=90,
                native_aliases_A='[{"isoform":"novel.1","reads":30}]', native_aliases_B="[]",
                max_single_native_reads_A=30, max_single_native_reads_B=35)
    divergent = [dict(clone_A="TF-1", clone_B="TF-2", modality=m, target="same_label", call_A=True,
                      call_B=False, max_state_mean_cpm=2 if m == "ppi" else None) for m in ["ppi", "pdi"]]
    nodes, edges = query_graph(divergent, [pair])
    assert set(nodes.id) == {"TF-1", "TF-2", "ppi:same_label", "pdi:same_label"}
    assert len(edges) == 4
    assert edges.detected.sum() == 2
    assert edges.tested.all()
    payload = to_cytoscape_elements(nodes, edges, include_undetected=True)
    json.dumps(payload, allow_nan=False)  # missing PDI/source gene CPM must be valid JSON


def test_conflicting_calls_reject():
    rows = [dict(source=s, target="X", edge_type="ppi", tested=True, detected=d)
            for s, d in [("A", True), ("A", False), ("B", False)]]
    with pytest.raises(ValueError, match="Conflicting"):
        differential_interactome("A", "B", pd.DataFrame(rows), "ppi")


def test_substitution_support_is_explicit_same_gene_length_and_exact_first():
    from altanalyze3.components.iso2function.query import protein_matches
    exact={('ENSG1','MACDEF'):['G-1']}
    by_length={('ENSG1',6):[('G-1','MACDEF'),('G-2','MAYDEF'),('G-3','MAYDEY'),('stop','MAYDE*')],
               ('ENSG2',6):[('OTHER','MACDEY')]}
    # Exact-first prevents inflating support for a second assay from the same native protein.
    assert list(protein_matches('MACDEF','ENSG1',exact,by_length,2))==['G-1']
    assert protein_matches('MACDEY','ENSG1',exact,by_length,0)=={}
    hit=protein_matches('MACDEY','ENSG1',exact,by_length,1)
    assert set(hit)=={'G-1','G-3'}
    assert hit['G-1']['substitutions']==[dict(position=6,native='Y',assayed='F')]
    assert 'OTHER' not in hit and 'stop' not in hit
    assert protein_matches('MACDEFY','ENSG1',exact,by_length,2)=={}
    assert protein_matches('MACDE*','ENSG1',exact,by_length,2)=={}


@pytest.mark.parametrize('distinct_native_proteins',[True,False])
def test_cli_substitution_transfer_requires_two_supported_different_proteins(tmp_path,distinct_native_proteins):
    data=tmp_path/'data';data.mkdir()
    clones=_tsv(tmp_path/'clones.tsv',[dict(clone_id='TF-1',gene_symbol='TF',aa_seq='MA'),dict(clone_id='TF-2',gene_symbol='TF',aa_seq='MB')])
    _tsv(data/'ppi_y2h.tsv',[dict(ad_clone_id='TF-1',db_gene_symbol='P',detected=True,tested=True),dict(ad_clone_id='TF-2',db_gene_symbol='P',detected=False,tested=True)])
    meta=_tsv(tmp_path/'meta.tsv',[dict(library='L',original_library='L',groups='healthy')])
    var=['ENSG1:sub.source001']+(['ENSG1:exact'] if distinct_native_proteins else [])
    counts=_h5(tmp_path/'counts.h5ad',['HSC.L-isoform'],var,[[30]*len(var)])
    fasta=tmp_path/'proteins.fa';fasta.write_text('>sub.source001 gene:ENSG1\nMC\n'+('>exact gene:ENSG1\nMA\n' if distinct_native_proteins else ''))
    ann=tmp_path/'genes.tsv';ann.write_text('ENSG1\tTF\n')
    mane=tmp_path/'mane.tsv';mane.write_text('symbol\tEnsembl_Gene\nTF\tENSG1\n')
    out=tmp_path/'result'
    assert main(['query','--out-dir',str(data),'--artifacts-dir',str(out),'--clonelist',clones,
                 '--pseudobulk',counts,'--protein-fasta',str(fasta),'--metadata',meta,'--state','HSC',
                 '--gene-annotation',str(ann),'--mane-summary',str(mane),'--max-aa-substitutions','1','--interaction-type','ppi'])==0
    decisions=pd.read_csv(out/'all_pair_decisions.tsv',sep='\t')
    assert bool(decisions.iloc[0].native_supported)==distinct_native_proteins
    assert bool(decisions.iloc[0].retained)==distinct_native_proteins
    if distinct_native_proteins:
        assert decisions.iloc[0].native_pair_evidence=='putative_transfer_with_AA_substitutions'
        alias=json.loads(decisions.iloc[0].native_aliases_B)[0]
        assert alias['substitutions']==[dict(position=2,native='C',assayed='B')]
        assert alias['reads']==30


def test_large_native_alias_fields_are_never_truncated(tmp_path):
    import csv
    from altanalyze3.components.iso2function.query import _table
    from altanalyze3.components.iso2function.protein_impact import table
    aliases=json.dumps([dict(isoform='novel.'+str(i),reads=30,substitutions=[dict(position=20,native='G',assayed='S')]) for i in range(2500)])
    assert len(aliases)>131072
    path=tmp_path/'large.tsv'
    with path.open('w') as f:
        writer=csv.DictWriter(f,['gene','native_aliases_A'],delimiter='\t');writer.writeheader();writer.writerow(dict(gene='G',native_aliases_A=aliases))
    limit=csv.field_size_limit()
    assert _table(path)[0]['native_aliases_A']==aliases
    assert table(path)[0]['native_aliases_A']==aliases
    assert csv.field_size_limit()==limit
