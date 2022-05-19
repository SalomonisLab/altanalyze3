import enum
from collections import namedtuple


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