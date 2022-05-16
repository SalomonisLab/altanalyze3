import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

CONTROL_MD5_SUMS = {
    "hg19_ref.tsv": "cb3093cbaeb5a27d95fc09add352641a",
    "hg19_ref_introns.bed.gz": "a9f441a78f87d47100fb4f18c827a976",
    "hg19_ref_introns.bed.gz.tbi": "007a8d670d47f2c54d97f3b26c9c054c",
    "hg19_pe.bam": "016704992a644f8c48d943c09d8bc2ce",
    "hg19_pe.bam.bai": "566547efad2c85af4c2dc465d23b1471",
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"