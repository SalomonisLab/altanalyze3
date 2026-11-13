"""Configuration adapter for the validated cross-chromosomal reporting scripts."""
import json
import os
from pathlib import Path

CONFIG = json.loads(Path(os.environ['ALTANALYZE3_NONCANONICAL_CONFIG']).read_text())
def config_path(name):
    return Path(CONFIG[name])
X = config_path('output') / 'trans_splicing'
# Historical symbols retained for imports; all scientific resource paths below are explicit.
REVIO = X
KINNEX = X
GENOME = config_path('genome')
EXON = config_path('exon_model')
REF_GFF = config_path('reference_gff')
SAMPLES = X / 'config/samples.tsv'
BARCODE_CLUSTERS = config_path('barcode_clusters')
def inside_x(path):
    path = Path(path).resolve()
    root = X.resolve()
    if path != root and root not in path.parents:
        raise ValueError('Output outside configured branch: ' + str(path))
    return path
def refuse_existing(*paths):
    for path in paths:
        if Path(path).exists():
            raise FileExistsError(path)
