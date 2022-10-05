import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

CONTROL_MD5_SUMS = {
    "hg19_ref.tsv": "cb3093cbaeb5a27d95fc09add352641a",
    "hg19_ref_introns.bed.gz": "a9f441a78f87d47100fb4f18c827a976",
    "hg19_pe.bam": "b94ee144a2bc4f9764e54aa367a5a42e",
    "hg19_pe.bam.bai": "efc5314e764a8ac3ce7ad060c92fad3f",
    "hg19_pe_copy.bam": "b94ee144a2bc4f9764e54aa367a5a42e",
    "hg19_pe_copy.bam.bai": "efc5314e764a8ac3ce7ad060c92fad3f"
}


def test_all_md5_sums():
    collected_md5_sums = {}
    for location in DATA_FOLDER.glob("**/*"):
        if location.is_file():
            collected_md5_sums[location.name] = get_md5_sum(location)
    for filename, control_md5_sum in CONTROL_MD5_SUMS.items():
        assert filename in collected_md5_sums, "missing file"
        assert collected_md5_sums[filename] == control_md5_sum, "wrong MD5 sum"