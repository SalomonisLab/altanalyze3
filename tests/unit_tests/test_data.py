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
    "Cal27P5_1_aggregated_counts.bed.gz": "faee5946f7d18064773fc67e01e4ce5b",
    "Cal27P5_1_intcounts.bed": "6f5250f1a841438fadcbe1419966f133",
    "Cal27P5_1_juncounts.bed": "24f6a6ac5f831c7709102c235db03ec5",
    "Cal27P5_2_intcounts.bed": "2d9889a6520f96012ba7d4be5cdcfea3",
    "Cal27P5_2_juncounts.bed": "0ab1f2283ba349f9b08dc76d1f29dfe9",
    "Cal27P5_3_intcounts.bed": "be222bf62519f5b6c601377955e4ef76",
    "Cal27P5_3_juncounts.bed": "91863b9f8ce7c7d8fc2177b50ef0175e"
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        if location.is_file():
            collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"