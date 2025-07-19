import os
import pandas as pd
import logging
from pathlib import Path
from altanalyze3.utilities.constants import AnnotationsParams
from altanalyze3.components.aggregate.main import get_annotation, Job

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Test junction coordinates
TEST_JUNCTIONS = [
    "chr16:176811-176928",
    "chr8:7016862-7018221",
    "chr16:173007-173124",
    "chr8:6997763-6999122",
    "chr16:173329-173471",
    "chr8:6978646-6980012",
    "chr16:177133-177282",
    "chr12:69348544-69350107",
    "chr1:153390243-153390394",
    "chr14:105742222-105742340",
    "chr1:153358433-153360643",
    "chr8:6977879-6978459",
    "chr14:105741795-105741892",
    "chr14:105742385-105742776",
    "chr11:5226799-5226929",
    "chr11:5225726-5226576",
    "chr10:69097231-69103870",
    "chr8:6996996-6997576",
    "chr11:57389317-57389886",
    "chr12:69352298-69353152"
]

# Mock reference data (example format - we'll need to populate this with real data)
MOCK_REFERENCE = """
chr16	176800	176900	E1.1	ENSG00000123456	+
chr8	7016800	7018300	E2.1	ENSG00000789012	+
chr16	173000	173200	E3.1	ENSG00000123456	+
chr8	6997700	6999200	E4.1	ENSG00000789012	+
chr16	173300	173500	E5.1	ENSG00000123456	+
chr8	6978600	6980100	E6.1	ENSG00000789012	+
chr16	177100	177300	E7.1	ENSG00000123456	+
chr12	69348500	69350100	E8.1	ENSG00000987654	+
chr1	153390200	153390400	E9.1	ENSG00000111111	+
chr14	105742200	105742400	E10.1	ENSG00000222222	+
chr1	153358400	153360700	E11.1	ENSG00000111111	+
chr8	6977800	6978500	E12.1	ENSG00000789012	+
chr14	105741700	105741900	E13.1	ENSG00000222222	+
chr14	105742300	105742800	E14.1	ENSG00000222222	+
chr11	5226700	5227000	E15.1	ENSG00000333333	+
chr11	5225700	5226600	E16.1	ENSG00000333333	+
chr10	69097200	69103900	E17.1	ENSG00000444444	+
chr8	6996900	6997600	E18.1	ENSG00000789012	+
chr11	57389300	57389900	E19.1	ENSG00000333333	+
chr12	69352200	69353200	E20.1	ENSG00000987654	+
"""

def create_test_files(tmp_dir):
    """Create temporary test files for junctions and references"""
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(exist_ok=True)
    
    # Create reference file
    ref_file = tmp_dir / "reference.bed"
    with open(ref_file, "w") as f:
        f.write(MOCK_REFERENCE.strip())
    
    # Create junction file
    jun_file = tmp_dir / "junctions.bed"
    with open(jun_file, "w") as f:
        for i, junction in enumerate(TEST_JUNCTIONS):
            chrom, coords = junction.split(":")
            start, end = coords.split("-")
            f.write(f"{chrom}\t{start}\t{end}\tJUNC_{i}\t1\t.\n")
    
    return ref_file, jun_file

def test_junction_annotation():
    """Test the junction annotation process"""
    # Create test files
    tmp_dir = Path("test_junction_annotation")
    ref_file, jun_file = create_test_files(tmp_dir)
    
    # Create a job for testing
    job = Job(contig="chr16", location=tmp_dir / "output.bed")
    
    # Process each junction
    logger.info("Starting junction annotation test")
    for junction in TEST_JUNCTIONS:
        chrom, coords = junction.split(":")
        start, end = map(int, coords.split("-"))
        
        logger.info(f"\nProcessing junction: {junction}")
        logger.info(f"Chromosome: {chrom}")
        logger.info(f"Start: {start}")
        logger.info(f"End: {end}")
        
        # Get annotation
        annotations = get_annotation(job, jun_file, ref_file)
        
        # Log results
        logger.info("Annotation results:")
        for annotation in annotations:
            if annotation is None:
                logger.info("No annotation found")
            else:
                logger.info(f"Gene: {annotation.gene}")
                logger.info(f"Exon: {annotation.exon}")
                logger.info(f"Strand: {annotation.strand}")
                logger.info(f"Position: {annotation.position}")
                logger.info(f"Match type: {annotation.match}")
    
    # Cleanup
    import shutil
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    test_junction_annotation() 