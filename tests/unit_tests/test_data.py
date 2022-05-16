import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

CONTROL_MD5_SUMS = {
    "Mm_Ensembl_exon_Mm10.tsv": "6a3a934a6520b5f07e10f5aa8ef7f7cf"
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"