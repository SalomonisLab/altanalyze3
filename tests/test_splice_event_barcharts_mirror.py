"""splice_event_barcharts --mirror: per-comparison counts, the sample and event filters, and unchanged default outputs."""
import os

import pandas as pd

from altanalyze3.components.visualization import splice_event_barcharts as seb

HEADER = ["Gene_Symbol", "Isoform_1|Length", "Isoform_2|Length", "Coord_1", "Coord_2", "Event_Type", "Event_Direction",
          "Trans_Splicing_Status", "Trans_Splicing_Evidence", "Feature", "N_aged", "N_young", "Median_aged",
          "Median_young", "Mean_aged", "Mean_young", "DiffMeans_aged_m_young", "ebayesStat", "ebayesSign",
          "ebayesPval", "ebayesAdjPval"]


def _write(path, rows, max_n=(5, 4)):
    out = []
    for k, (etype, direction, dpsi, p, n) in enumerate(rows):
        out.append(["G", f"ENST{k}|100", f"ENST{k + 50}|90(NMD)", "c1", "c2", etype, direction, "", "", f"G:E{k}",
                    n, n, 0.5, 0.5, 0.5 + dpsi, 0.5, dpsi, 1, 1, p, 1])
    # a null-direction row carrying the comparison's largest N per condition (never an event)
    out.append(["G", "ENSTX|10", "ENSTY|20", "c1", "c2", "Cassette Exon", "null", "", "", "G:EX",
                max_n[0], max_n[1], 0, 0, 0, 0, 0, 1, 1, 1, 1])
    pd.DataFrame(out, columns=HEADER).to_csv(path, sep="\t", index=False)


def test_mirror_counts_filters_and_defaults(tmp_path):
    d = tmp_path / "events"; d.mkdir()
    hsc = [(t, dr, 0.3, 0.01, 3) for t, dr in [("Cassette Exon", "inclusion"), ("Intron Retention", "exclusion"),
                                               ("AltPromoter", "proximal"), ("Alt5' SS", "distal")]]
    _write(d / "young-aged-HSC-1_stats-annotated.txt", hsc)
    _write(d / "young-aged-MPP-1_stats-annotated.txt", hsc[:1], max_n=(5, 4))       # 1 event
    _write(d / "young-aged-ERP-8_stats-annotated.txt", [r[:4] + (2,) for r in hsc], max_n=(3, 2))  # 2 young samples
    files = sorted(str(x) for x in d.glob("*_stats-annotated.txt"))
    order = ["HSC-1", "MPP-1", "ERP-8"]
    outs = {}
    for name, kw in (("plain", {}), ("mirror", dict(mirror=True, mirror_min_samples=3, mirror_min_events=2))):
        o = tmp_path / name
        seb.analyze_splicing_events(files, output_dir=str(o), pval_threshold=0.05, diff_means_threshold=0.1,
                                    sample_min=2, cluster_order=order, **kw)
        outs[name] = {f: open(o / f).read() for f in sorted(os.listdir(o)) if f.endswith("_significant-results.tsv")}
    assert outs["plain"] == outs["mirror"]
    m = pd.read_csv(tmp_path / "mirror" / "splicing_events_inclusion_exclusion_mirror.tsv", sep="\t",
                    keep_default_na=False).set_index("cell_state")
    assert list(m.index) == order
    assert m.loc["HSC-1", ["samples_aged", "samples_young", "total_events", "inclusion_proximal_events",
                           "exclusion_distal_events"]].tolist() == [5, 4, 4, 2, 2]
    assert m.loc["HSC-1", "kept"] in (True, "True")
    assert m.loc["MPP-1", "excluded_because"] == "events 1 < 2"
    assert m.loc["ERP-8", "excluded_because"].startswith("young samples 2 < 3")
    assert (tmp_path / "mirror" / "splicing_events_inclusion_exclusion_mirror.pdf").stat().st_size > 1000
    assert not (tmp_path / "plain" / "splicing_events_inclusion_exclusion_mirror.pdf").exists()


def test_strict_cutoffs_and_nonfinite_do_not_leak_into_event_counts(tmp_path):
    d=tmp_path/'events';d.mkdir()
    _write(d/'aged-young-HSC-1_stats-annotated.txt',[
        ('Cassette Exon','inclusion',.2,.049,3),
        ('Cassette Exon','inclusion',.2,.05,3),
        ('Cassette Exon','inclusion',.1,.01,3),
        ('Cassette Exon','inclusion',.2,float('nan'),3),
        ('Cassette Exon','inclusion',float('nan'),.01,3)])
    out=tmp_path/'out'
    seb.analyze_splicing_events([str(d/'aged-young-HSC-1_stats-annotated.txt')],
        output_dir=str(out),sample_min=3,cluster_order=['HSC-1'],mirror=True,
        pval_threshold=.05,diff_means_threshold=.1,strict_thresholds=True)
    table=pd.read_csv(out/'splicing_events_inclusion_exclusion_mirror.tsv',sep='\t')
    assert table.total_events.tolist()==[1]
