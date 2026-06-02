#!/usr/bin/env python3
"""Tier 1 -- per-sample structure reduction (run per sample, in series or parallel).

For ONE sample's molecule-level ``transcript_associations.txt`` (one row per read:
``gene  strand  structure  molecule_id  source``), reduce the millions of reads to the sample's
DISTINCT structures with read counts, choosing a representative molecule id per structure:

    <sample>.exemplars.tsv     gene, strand, structure, count, rep_molecule_id
    <sample>.mol2struct.tsv     molecule_id, gene, structure

Tier 1 does ONLY exact-structure aggregation -- identical structures unify (hash to the same key);
NO containment folding happens here. The containment fold (the validated ``scored_collapse.collapse_gene``)
runs once in Tier 2 on the cross-sample SUMMED counts, so the Tier1->Tier2 result is identical to the
previously-validated single pooled collapse (same per-structure counts feed collapse_gene), while
memory stays bounded by a sample's distinct structures (not its read count).

``rep_molecule_id`` is the first molecule id seen for that exact structure in this sample. Tier 2
merges per-sample reps; the final representative molecule (sample-tagged) becomes the isoform id, and
the ``mol2struct`` map lets Stage 3 build the per-sample isoform h5ad from molecule-level counts.
"""

from __future__ import annotations

import collections
import os
import time


def collapse_sample(sample_name, ta_path, out_dir, log=print):
    """Tier 1: reduce one sample's per-read TA to distinct structures + counts + rep molecule.
    Returns (exemplars_path, mol2struct_path)."""
    t = time.time()
    os.makedirs(out_dir, exist_ok=True)
    exemplars_path = os.path.join(out_dir, f"{sample_name}.exemplars.tsv")
    mol2_path = os.path.join(out_dir, f"{sample_name}.mol2struct.tsv")

    # streaming aggregation: (gene) -> {structure: read_count}, first molecule per structure, strand.
    gene_struct = collections.defaultdict(collections.Counter)
    gene_struct_mol = collections.defaultdict(dict)
    gene_strand = {}
    n_reads = 0
    with open(ta_path) as f, open(mol2_path, 'w') as mol_out:
        mol_out.write("molecule_id\tgene\tstructure\n")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            gene, strand, struct, mol = p[0], p[1], p[2], p[3]
            if not struct or 'UNK' in gene:
                continue
            gene_struct[gene][struct] += 1
            d = gene_struct_mol[gene]
            if struct not in d:
                d[struct] = mol
            gene_strand[gene] = strand
            mol_out.write(f"{mol}\t{gene}\t{struct}\n")
            n_reads += 1

    n_struct = 0
    with open(exemplars_path, 'w') as o:
        o.write("gene\tstrand\tstructure\tcount\trep_molecule_id\n")
        for gene, sr in gene_struct.items():
            strand = gene_strand[gene]
            for struct, cnt in sr.items():
                o.write(f"{gene}\t{strand}\t{struct}\t{cnt}\t{gene_struct_mol[gene][struct]}\n")
                n_struct += 1

    log(f"[tier1:{sample_name}] {n_reads:,} reads -> {n_struct:,} distinct structures "
        f"({len(gene_struct):,} genes) {time.time()-t:.1f}s")
    return exemplars_path, mol2_path
