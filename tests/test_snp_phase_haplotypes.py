"""altanalyze3.components.snp.phase_haplotypes: recover simulated haplotypes from reads with errors."""
import numpy as np

from altanalyze3.components.snp.phase_haplotypes import ParityUnionFind, phase


def test_parity_union_find():
    uf = ParityUnionFind()
    for x in "abcd":
        uf.add(x)
    assert uf.union("a", "b", 1) == "joined"
    assert uf.union("b", "c", 1) == "joined"
    assert uf.union("a", "c", 0) == "consistent"
    assert uf.union("a", "c", 1) == "conflict"
    assert uf.parity_of("a") ^ uf.parity_of("c") == 0
    assert uf.find("d") == "d"


def test_phase_recovers_simulated_haplotypes():
    rng = np.random.default_rng(5)
    sites = ["chr1:{}:A:G".format(p) for p in range(100, 2000, 100)]          # 19 sites, one block
    truth = {s: int(rng.integers(0, 2)) for s in sites}                       # REF on haplotype truth[s]
    reads, true_hap = {}, {}
    for i in range(400):
        h = int(rng.integers(0, 2))
        start = int(rng.integers(0, len(sites) - 3))
        cover = sites[start:start + int(rng.integers(2, 6))]
        calls = []
        for s in cover:
            is_ref = (truth[s] == h)
            if rng.random() < 0.03:                                           # 3% allele errors
                is_ref = not is_ref
            calls.append((s, is_ref))
        reads["r{}".format(i)] = calls
        true_hap["r{}".format(i)] = h
    site_phase, edges, read_haps = phase(reads)
    blocks = {b for b, _ in site_phase.values()}
    assert len(blocks) == 1
    # phase bits equal the truth up to one global flip
    flip = site_phase[sites[0]][1] ^ truth[sites[0]]
    assert all((site_phase[s][1] ^ truth[s]) == flip for s in sites)
    # reads: >= 95% assigned the true haplotype (up to the same flip)
    ok = sum(1 for k, hs in read_haps.items() if hs[0][1] >= 0 and (hs[0][1] ^ flip) == true_hap[k])
    assert ok / len(reads) > 0.95
    assert any(e["status"] in ("consistent", "conflict") for e in edges)


def test_longcallr_check_concordance(tmp_path):
    import gzip, os
    from altanalyze3.components.snp import het_panel
    from altanalyze3.components.snp.phase_haplotypes import longcallr_check
    tmp = str(tmp_path)
    meta = os.path.join(tmp, "meta.txt")
    with open(meta, "w") as fh:
        fh.write("uid\tbam\tlibrary\treverse\tgroups\nU0\t{}/x.bam\tLIBA\tTRUE\tyoung\n".format(tmp))
    dmap = os.path.join(tmp, "dmap.tsv")
    with open(dmap, "w") as fh:
        fh.write("uid\toriginal_uid\tdonor\nLIBA\tU0\tD\n")
    os.makedirs(os.path.join(tmp, "panel", "libraries"))
    gts = {100: "0|1", 200: "1|0", 300: "1|0", 400: "0|1"}
    with gzip.open(os.path.join(tmp, "panel", "libraries", "00__LIBA.sites.tsv.gz"), "wt") as fh:
        fh.write("\t".join(het_panel.LIB_FIELDS) + "\n")
        for pos, gt in gts.items():
            row = {k: "" for k in het_panel.LIB_FIELDS}
            row.update(chrom="chr1", pos=pos, ref="A", alt="G", vtype="SNV", longcallr_gt=gt, longcallr_ps="100")
            fh.write("\t".join(str(row[k]) for k in het_panel.LIB_FIELDS) + "\n")
    # our phase: bits 0,1,1 agree with longcallR; site 400 placed with the wrong bit
    site_phase = {"chr1:100:A:G": ("B", 0), "chr1:200:A:G": ("B", 1), "chr1:300:A:G": ("B", 1),
                  "chr1:400:A:G": ("B", 1)}
    rows = longcallr_check(os.path.join(tmp, "panel"), meta, dmap, "D", "chr1", site_phase)
    assert rows[0]["concordant"] == 2 and rows[0]["discordant"] == 1
