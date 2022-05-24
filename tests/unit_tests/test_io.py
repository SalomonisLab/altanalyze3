import pytest
import pathlib
from altanalyze3.utilities.io import (
    get_all_bam_chr,
    get_all_ref_chr,
    is_bam_paired
)


DATA_FOLDER = pathlib.Path(__file__).resolve().parents[1].joinpath("data")

@pytest.mark.parametrize(
    "location, control_ref_chr",
    [
        (
            "hg19_ref_introns.bed.gz",
            [
                "chr1", "chr10", "chr11", "chr12", "chr13",
                "chr14", "chr15", "chr16", "chr17", "chr18",
                "chr19", "chr2", "chr20", "chr21", "chr22",
                "chr3", "chr4", "chr5", "chr6", "chr7",
                "chr8", "chr9", "chrGL000192.1", "chrGL000223.1",
                "chrHG104_HG975_PATCH", "chrHG1063_PATCH",
                "chrHG1079_PATCH", "chrHG115_PATCH",
                "chrHG1211_PATCH", "chrHG122_PATCH",
                "chrHG1257_PATCH", "chrHG1287_PATCH",
                "chrHG1292_PATCH", "chrHG1304_PATCH",
                "chrHG1308_PATCH", "chrHG1426_PATCH",
                "chrHG1433_PATCH", "chrHG1436_HG1432_PATCH",
                "chrHG1439_PATCH", "chrHG1440_PATCH", "chrHG1442_PATCH",
                "chrHG144_PATCH", "chrHG1462_PATCH", "chrHG1463_PATCH",
                "chrHG1486_PATCH", "chrHG1497_PATCH", "chrHG1501_PATCH",
                "chrHG1502_PATCH", "chrHG1592_PATCH", "chrHG183_PATCH",
                "chrHG185_PATCH", "chrHG186_PATCH", "chrHG19_PATCH",
                "chrHG243_PATCH", "chrHG280_PATCH", "chrHG306_PATCH",
                "chrHG311_PATCH", "chrHG344_PATCH", "chrHG385_PATCH",
                "chrHG388_HG400_PATCH", "chrHG418_PATCH",
                "chrHG506_HG507_HG1000_PATCH", "chrHG544_PATCH",
                "chrHG686_PATCH", "chrHG729_PATCH", "chrHG79_PATCH",
                "chrHG7_PATCH", "chrHG865_PATCH", "chrHG873_PATCH",
                "chrHG905_PATCH", "chrHG946_PATCH", "chrHG962_PATCH",
                "chrHG971_PATCH", "chrHG989_PATCH", "chrHSCHR15_1_CTG4",
                "chrHSCHR17_4_CTG4", "chrHSCHR19LRC_COX1_CTG1",
                "chrHSCHR19LRC_COX2_CTG1", "chrHSCHR19LRC_LRC_I_CTG1",
                "chrHSCHR19LRC_LRC_J_CTG1", "chrHSCHR19LRC_LRC_S_CTG1",
                "chrHSCHR19LRC_LRC_T_CTG1", "chrHSCHR19LRC_PGF1_CTG1",
                "chrHSCHR19LRC_PGF2_CTG1", "chrHSCHR19_1_CTG3", "chrHSCHR1_1_CTG31",
                "chrHSCHR21_4_CTG1_1", "chrHSCHR3_1_CTG1", "chrHSCHR6_MHC_APD",
                "chrHSCHR6_MHC_COX", "chrHSCHR6_MHC_DBB", "chrHSCHR6_MHC_MANN",
                "chrHSCHR6_MHC_MCF", "chrHSCHR6_MHC_QBL", "chrHSCHR6_MHC_SSTO",
                "chrHSCHR7_1_CTG6", "chrHSCHR9_3_CTG35", "chrLRG_1", "chrLRG_10",
                "chrLRG_135", "chrLRG_137", "chrLRG_156", "chrLRG_162", "chrLRG_179",
                "chrLRG_189", "chrLRG_190", "chrLRG_196", "chrLRG_198", "chrLRG_201",
                "chrLRG_230", "chrLRG_231", "chrLRG_249", "chrLRG_251", "chrLRG_267",
                "chrLRG_272", "chrLRG_283", "chrLRG_286", "chrLRG_287", "chrLRG_293",
                "chrLRG_331", "chrLRG_334", "chrLRG_336", "chrLRG_363", "chrLRG_377",
                "chrLRG_383", "chrLRG_392", "chrLRG_41", "chrLRG_411", "chrLRG_426",
                "chrLRG_53", "chrLRG_55", "chrLRG_56", "chrLRG_6", "chrLRG_7",
                "chrLRG_76", "chrLRG_94", "chrX", "chrY"
            ]
        )
    ]
)
def test_get_all_ref_chr(location, control_ref_chr):                           # this will also check guard_chr function
    calculated_ref_chr = get_all_ref_chr(DATA_FOLDER.joinpath(location), 1)
    assert sorted(calculated_ref_chr) == sorted(control_ref_chr)


@pytest.mark.parametrize(
    "location, control_bam_chr",
    [
        (
            "hg19_pe.bam",
            [
                "chr1", "chr10", "chr11", "chr12", "chr13", "chr14",
                "chr15", "chr16", "chr17", "chr18", "chr19", "chr2",
                "chr20", "chr21", "chr22", "chr3", "chr4", "chr5",
                "chr6", "chr7", "chr8", "chr9", "chrM", "chrX", "chrY"
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
        ("hg19_pe.bam", True)
    ]
)
def test_is_bam_paired(location, control_is_paired):                           # this will also check skip_bam_read function
    calculated_is_paired = is_bam_paired(DATA_FOLDER.joinpath(location), 1)
    assert calculated_is_paired == control_is_paired