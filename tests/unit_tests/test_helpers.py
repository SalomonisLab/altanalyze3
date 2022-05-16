import pytest
import pathlib
from altanalyze3.utilities.helpers import (
    get_md5_sum,
    get_tmp_suffix
)


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")


@pytest.mark.parametrize(
    "location, control_md5_sum",
    [
        ("hg19_ref.tsv", "cb3093cbaeb5a27d95fc09add352641a")
    ]
)
def test_get_md5_sum(location, control_md5_sum):
    calculated_md5_sum = get_md5_sum(DATA_FOLDER.joinpath(location))
    assert calculated_md5_sum == control_md5_sum


@pytest.mark.parametrize(
    "length_without_dot",
    [10, 12, 3]
)
def test_get_tmp_suffix(length_without_dot):
    tmp_suffix = get_tmp_suffix(length_without_dot)
    assert len(tmp_suffix) == length_without_dot + 1