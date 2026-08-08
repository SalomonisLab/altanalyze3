'''
scTriangulate, integrated into AltAnalyze3.

Reconciles conflicting single cell cluster annotations by scoring each candidate
cluster for biological stability and letting the annotations compete cell by
cell through a Shapley-style game.

Upstream: https://github.com/frankligy/scTriangulate (version 0.13.0,
commit 8b9598cf, 2026-03-30). See README.md in this directory for what was
changed, what was measured, and how equivalence was checked.

Entry points:
    from altanalyze3.components.sctriangulate import ScTriangulate
    sctri = ScTriangulate(dir=outdir, adata=adata, query=['ann1','ann2'])
    sctri.lazy_run()

  or from a shell:
    python3 -m altanalyze3.components.sctriangulate.cli --h5ad in.h5ad \\
        --query ann1,ann2 --outdir out
'''

from .main_class import ScTriangulate, sctriangulate_setting
from .metrics import (
    tf_idf1_for_cluster, tf_idf5_for_cluster, tf_idf10_for_cluster,
    marker_gene, reassign_score, SCCAF_score, doublet_compute,
)
from .shapley import shapley_batch, which_to_take_batch

__all__ = [
    'ScTriangulate', 'sctriangulate_setting',
    'tf_idf1_for_cluster', 'tf_idf5_for_cluster', 'tf_idf10_for_cluster',
    'marker_gene', 'reassign_score', 'SCCAF_score', 'doublet_compute',
    'shapley_batch', 'which_to_take_batch',
]

__upstream_version__ = '0.13.0'
__upstream_commit__ = '8b9598cfcbf269856f41740514cadcd75f1ee2c6'
