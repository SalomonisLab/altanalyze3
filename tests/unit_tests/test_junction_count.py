import pytest
import pathlib
from altanalyze3.utilities.parser import ArgsParser
from altanalyze3.utilities.helpers import get_md5_sum
from altanalyze3.components.junction_count.main import get_jobs


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")


@pytest.mark.parametrize(
    "args, control_chr",
    [
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
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
                "--chr", "1", "2", "3"
            ],
            [
                "chr1", "chr2", "chr3"
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "chr1", "chr2", "chr3", "chr99"
            ],
            [
                "chr1", "chr2", "chr3"
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "chr99"
            ],
            []
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam"
            ],
            [
                "chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                "chrX", "chrY"
            ]
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1-copy.bam"
            ],
            [
                "chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                "chrX", "chrY"
            ]
        )
    ]
)
def test_get_jobs(monkeypatch, args, control_chr):
    monkeypatch.chdir(DATA_FOLDER)
    args = ArgsParser(args)
    calculated_chr = [job.contig for job in get_jobs(args)]
    assert sorted(calculated_chr) == sorted(control_chr)


@pytest.mark.parametrize(
    "args, control_md5_sum",
    [
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "chr1", "chr2", "chr3"
            ],
            "849c5a8479892e291e17af2aaa09dedd"
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "1", "2", "3"
            ],
            "849c5a8479892e291e17af2aaa09dedd"
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "chr1", "chr2", "chr3", "chr99"
            ],
            "849c5a8479892e291e17af2aaa09dedd"
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam",
                "--chr", "chr99"
            ],
            "d41d8cd98f00b204e9800998ecf8427e"
        ),
        (
            [
                "juncount",
                "--bam", "Cal27P5-1.bam"
            ],
            "24f6a6ac5f831c7709102c235db03ec5"
        )
    ]
)
def test_count_junctions(monkeypatch, args, control_md5_sum):             # tests all function from junction_count.main file
    monkeypatch.chdir(DATA_FOLDER)
    args = ArgsParser(args)
    args.func(args)
    calculated_md5_sum = get_md5_sum(args.output.with_suffix(".bed"))
    args.output.with_suffix(".bed").unlink()
    assert calculated_md5_sum == control_md5_sum