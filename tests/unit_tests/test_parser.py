import pytest
import pathlib
from altanalyze3.utilities.parser import (
    ArgsParser
)


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")


@pytest.mark.parametrize(
    "args, control_chr",
    [
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--ref", "gene_model_all.tsv"
            ],
            [
                "chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                "chrM", "chrMT", "chrX", "chrY"                                   # this will still include chrM and chrMT, because we use MAIN_CRH as default
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--ref", "gene_model_all.tsv",
                "--chr", "chr1"
            ],
            [
                "chr1"
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--ref", "gene_model_all.tsv",
                "--chr", "chr1", "chr2", "chr3"
            ],
            [
                "chr1", "chr2", "chr3"
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--ref", "gene_model_all.tsv",
                "--chr", "1", "2", "3"
            ],
            [
                "chr1", "chr2", "chr3"
            ]
        )
    ]
)
def test_chromosome_selection(monkeypatch, args, control_chr):
    monkeypatch.chdir(DATA_FOLDER)
    args = ArgsParser(args)
    assert sorted(args.chr) == sorted(control_chr)
