"""Execute core AltAnalyze3 stages without loading unrelated CLI components."""
import argparse
import logging
from pathlib import Path
from types import SimpleNamespace


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['index', 'junctions', 'introns', 'aggregate'])
    p.add_argument('--gtf'); p.add_argument('--bam'); p.add_argument('--reference')
    p.add_argument('--junctions', nargs='+'); p.add_argument('--introns', nargs='+')
    p.add_argument('--output', required=True)
    p.add_argument('--cpus', type=int, default=1)
    p.add_argument('--strandness', default='auto', choices=['auto', 'forward', 'reverse', 'unstranded'])
    p.add_argument('--novel-gene-mode', default='corrected', choices=['corrected', 'legacy'])
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    if a.stage == 'index':
        from altanalyze3.components.gene_model.gene_model_index import build_index
        build_index(SimpleNamespace(gtf=a.gtf, output=a.output))
    elif a.stage in ('junctions', 'introns'):
        import pysam
        from altanalyze3.utilities.io import get_all_bam_chr
        from altanalyze3.utilities.constants import IntRetCat
        bam = Path(a.bam)
        with pysam.AlignmentFile(str(bam), 'rb') as handle:
            if not handle.has_index():
                raise ValueError('BAM must be coordinate-sorted and indexed')
        args = SimpleNamespace(bam=bam, output=Path(a.output), ref=Path(a.reference) if a.reference else None,
            cpus=a.cpus, threads=1, tmp=Path('tmp_'+a.stage), chr=get_all_bam_chr(bam, 1),
            loglevel=logging.INFO, span=0, strandness=IntRetCat[a.strandness.upper()], savereads=False)
        args.tmp.mkdir(exist_ok=True)
        if a.stage == 'junctions':
            from altanalyze3.components.junction_count.main import count_junctions
            count_junctions(args)
        else:
            from altanalyze3.components.intron_count.main import count_introns
            count_introns(args)
    else:
        from altanalyze3.components.aggregate.main import aggregate
        reference = Path(a.reference)
        aggregate(SimpleNamespace(juncounts=[Path(x) for x in a.junctions or []],
                  intcounts=[Path(x) for x in a.introns or []], ref=reference,
                  output=a.output, chr=[], novel_gene_mode=a.novel_gene_mode))

if __name__ == '__main__':
    main()
