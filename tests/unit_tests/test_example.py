import pytest


@pytest.mark.parametrize(
    "input, control",
    [
        (
            "alpha",
            "alpha"
        ),
        (
            "beta",
            "beta"
        )
    ]
)
def test_dummy(input, control):
    assert input==control, \
        "Dummy test failed"