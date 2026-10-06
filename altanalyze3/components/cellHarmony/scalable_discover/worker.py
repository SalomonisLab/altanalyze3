"""Disposable worker entry point for one scALABLE-discover pipeline run."""
import argparse
from pathlib import Path

from altanalyze3.components.cellHarmony.flask.job_manager import JobStore

from .tasks import DiscoverJobRunner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('storage', type=Path)
    parser.add_argument('registry', type=Path)
    parser.add_argument('job_id')
    parser.add_argument('task', choices=['pipeline'])
    parser.add_argument('compression', choices=['lzf', 'gzip', 'none'])
    parser.add_argument('export_pdfs', type=int, choices=[0, 1])
    args = parser.parse_args()
    runner = DiscoverJobRunner(JobStore(args.storage), args.registry, h5ad_compression=args.compression,
                               export_approx_pdfs=bool(args.export_pdfs))
    try:
        runner._run_pipeline(args.job_id)
    finally:
        runner.executor.shutdown(wait=True)


if __name__ == '__main__':
    main()
