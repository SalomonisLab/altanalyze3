import pytest
import pathlib
from altanalyze3.utilities.io import (
    get_all_bam_chr,
    get_all_ref_chr,
    is_bam_paired,
    is_bam_indexed,
    get_indexed_references
)


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

@pytest.mark.parametrize(
    "location, control_ref_chr",
    [
        (
            "gene_model_all.tsv",
            [
                "chr1", "chr10", "chr11", "chr12", "chr13",
                "chr14", "chr15", "chr16", "chr17", "chr18",
                "chr19", "chr2", "chr20", "chr21", "chr22",
                "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                "chrMT", "chrX", "chrY", "chrKI270750.1",
                "chrGL000009.2", "chrGL000194.1", "chrGL000195.1",
                "chrGL000205.2", "chrGL000213.1", "chrGL000216.2",
                "chrGL000218.1", "chrGL000219.1", "chrGL000220.1",
                "chrGL000225.1", "chrKI270442.1", "chrKI270711.1",
                "chrKI270713.1", "chrKI270721.1", "chrKI270726.1",
                "chrKI270727.1", "chrKI270728.1", "chrKI270731.1",
                "chrKI270733.1", "chrKI270734.1", "chrKI270744.1",
            ]
        )
    ]
)
def test_get_all_ref_chr(tmp_path, location, control_ref_chr):                           # this will also check guard_chr function
    calculated_ref_chr = get_all_ref_chr(
        get_indexed_references(
            DATA_FOLDER.joinpath(location),
            tmp_path                                                                     # https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html#tmp-path
        ),
        1
    )
    assert sorted(calculated_ref_chr) == sorted(control_ref_chr)


@pytest.mark.parametrize(
    "location, control_bam_chr",
    [
        (
            "Cal27P5-1.bam",
            [
                "chr1", "chr10", "chr11", "chr12", "chr13", "chr14",
                "chr15", "chr16", "chr17", "chr18", "chr19", "chr2",
                "chr20", "chr21", "chr22", "chr3", "chr4", "chr5",
                "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"
            ]
        )
    ]
)
def test_get_all_bam_chr(location, control_bam_chr):                           # this will also check guard_chr function
    calculated_bam_chr = get_all_bam_chr(DATA_FOLDER.joinpath(location), 1)
    assert sorted(calculated_bam_chr) == sorted(control_bam_chr)


@pytest.mark.parametrize(
    "location, control_is_paired",
    [
        ("Cal27P5-1.bam", True)
    ]
)
def test_is_bam_paired(location, control_is_paired):                           # this will also check skip_bam_read function
    calculated_is_paired = is_bam_paired(DATA_FOLDER.joinpath(location), 1)
    assert calculated_is_paired == control_is_paired


@pytest.mark.parametrize(
    "location, control_is_indexed",
    [
        ("Cal27P5-1.bam", True)
    ]
)
def test_is_bam_indexed(location, control_is_indexed):
    calculated_is_indexed = is_bam_indexed(DATA_FOLDER.joinpath(location), 1)
    assert calculated_is_indexed == control_is_indexed