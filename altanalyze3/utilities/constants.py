import enum
from collections import namedtuple


ChrConverter = lambda c: c if c.startswith("chr") else f"chr{c}"


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
    "chrM",
    "chrMT"    # somehow in the gene model file the chromosome is named MT
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


class AnnMatchCat (enum.IntEnum):
    EXON_START = enum.auto()
    EXON_END = enum.auto()
    EXON_MID = enum.auto()
    INTRON_MID = enum.auto()
    CLOSEST = enum.auto()
    def __str__(self):
        return self.name


Annotation = namedtuple(
    "Annotation",
    "gene exon strand position match"
)


IntronsParams = {
    "usecols": [0, 1, 2, 3, 4, 5],
    "names": ["chr", "start", "end", "name", "score", "strand"],
    "index_col": ["chr", "start", "end", "name", "strand"],
    "converters": {"chr": ChrConverter},
    "dtype": {
        "start": "int32",
        "end": "int32",
        "name": "string",
        "score": "int32",
        "strand": "category"
    },
    "sep": "\t"
}


JunctionsParams = {
    "usecols": [0, 1, 2, 3, 4],
    "names": ["chr", "start", "end", "name", "score"],
    "index_col": ["chr", "start", "end", "name"],
    "converters": {"chr": ChrConverter},
    "dtype": {
        "start": "int32",
        "end": "int32",
        "name": "string",
        "score": "int32"
    },
    "sep": "\t"
}


ReferencesParams = {
    "usecols": [0, 1, 2, 3, 4, 5],
    "names": ["gene", "chr", "strand", "exon", "start", "end"],
    "index_col": ["chr", "start", "end"],
    "converters": {
        "chr": ChrConverter,
        "start": lambda c: int(c)-1                      # need to make [start, end) intervals
    },
    "dtype": {
        "gene": "string",
        "strand": "category",
        "exon": "category",
        "end": "int32"
    },
    "sep": "\t"
}

AnnotationsParams = {
    "usecols": [0, 1, 2, 3, 4, 5],
    "names": ["chr", "start", "end", "name", "annotation", "strand"],
    "index_col": ["chr", "start", "end", "name"],
    "converters": {"chr": ChrConverter},
    "dtype": {
        "start": "int32",
        "end": "int32",
        "name": "string",
        "annotation": "string",
        "strand": "category"
    },
    "sep": "\t"
}
