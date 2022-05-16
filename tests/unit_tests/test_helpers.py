import pytest
import pathlib
from altanalyze3.utilities.helpers import get_md5_sum


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")


@pytest.mark.parametrize(
    "location, control_md5_sum",
    [
        ("Mm_Ensembl_exon_Mm10.tsv", "6a3a934a6520b5f07e10f5aa8ef7f7cf")
    ]
)
def test_get_md5_sum(location, control_md5_sum):
    calculated_md5_sum = get_md5_sum(DATA_FOLDER.joinpath(location))
    assert calculated_md5_sum == control_md5_sum