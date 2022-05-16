import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

CONTROL_MD5_SUMS = {
    "Mm_Ensembl_exon_Mm10.tsv": "6a3a934a6520b5f07e10f5aa8ef7f7cf",
    "Mm_Ensembl_exon_Mm10_filtered.bed.gz": "5a11a3d012655c47b0436bc047193c9f",
    "Mm_Ensembl_exon_Mm10_filtered.bed.gz.tbi": "05500da64f0bbd43c1de19caeb3ba843"
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"