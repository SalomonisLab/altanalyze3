import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

CONTROL_MD5_SUMS = {
    "gene_model_all.tsv": "2b8531420d7916b9da21102586a96015",
    "Cal27P5-1.bam": "6e418cb604ea82c90d4b3e832afb4997",
    "Cal27P5-1.bam.bai": "6f0d6aceabade33f91ea69e5f036125c",
    "Cal27P5-2.bam": "8ce0bc51544144fd9b9bb1c66c172a61",
    "Cal27P5-2.bam.bai": "28bb5431b5d3bd1d6b8dc7a30a128011",
    "Cal27P5-3.bam": "3d6fdef5e90c16deed6c36b418c597c4",
    "Cal27P5-3.bam.bai": "c92fbf32ffcf874979c41051f0d2e16c",
    "Cal27P5-1-copy.bam": "6e418cb604ea82c90d4b3e832afb4997",
    "Cal27P5_1_aggregated_counts.bed.gz": "8299eed6c0f471850d67939f35c576bc"
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        if location.is_file():
            collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"