import enum
from collections import namedtuple


MAIN_CRH = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM"
]


class IntRetCat (enum.IntEnum):
    """
    A class to include any categorical infromation to be used intead of
    string or other data types comparison for intron retention component.
    """
    # to define categories of overlap
    PRIME_5 = enum.auto()
    PRIME_3 = enum.auto()
    INTRON = enum.auto()
    DISCARD = enum.auto()
    # to define strand specificity of RNA library
    AUTO = enum.auto()
    FORWARD = enum.auto()
    REVERSE = enum.auto()
    UNSTRANDED = enum.auto()

    def __str__(self):
        return self.name


Job = namedtuple(
    "Job",
    "contig location"
)


IntRetRawData = namedtuple(
    "IntRetRawData",
    "contig intron_start intron_end intron_name intron_strand read_start read_end read_strand xs_strand read_name read_1 read_2"
)


Annotation = namedtuple(
    "Annotation",
    "gene exon strand position name"
)