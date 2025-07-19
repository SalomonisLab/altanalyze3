import types
from pathlib import Path
from altanalyze3.components.aggregate.main import aggregate
import annotate
# Simulated CLI arguments
args = types.SimpleNamespace()
args.juncounts = Path("/Users/saljh8/Dropbox/Revio/BoneMarrow/test_env/altanalyze3/tests/Cal27P5-3_juncounts2.bed")
args.intcounts = Path("/Users/saljh8/Dropbox/Revio/BoneMarrow/test_env/altanalyze3/tests/Cal27P5-3_intcounts2.bed")
args.output = "/Users/saljh8/Dropbox/Revio/BoneMarrow/test_env/altanalyze3/tests/Cal27P5_3_aggregated_counts-new"
args.chr = ["chr6"]  # Limit for testing
args.threads = 2
args.cpus = 4

aggregate(args)

h5ad_path = args.output + '.h5ad'
exon_annotation_file = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'
output_path = args.output+'-dense.txt'
annotate.annotate_junctions_and_export(h5ad_path, exon_annotation_file, output_path)