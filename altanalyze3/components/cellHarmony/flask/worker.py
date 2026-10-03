"""Disposable worker entry point for upload and differential analyses."""
import argparse
from pathlib import Path

from .job_manager import JobStore
from .tasks import JobRunner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('storage', type=Path)
    parser.add_argument('registry', type=Path)
    parser.add_argument('job_id')
    parser.add_argument('task', choices=['pipeline', 'differential'])
    parser.add_argument('compression', choices=['lzf', 'gzip', 'none'])
    parser.add_argument('export_pdfs', type=int, choices=[0, 1])
    args = parser.parse_args()
    runner = JobRunner(JobStore(args.storage), args.registry,
                       h5ad_compression=args.compression, export_approx_pdfs=bool(args.export_pdfs))
    try:
        if args.task == 'pipeline':
            runner._run_pipeline(args.job_id)
        else:
            runner._run_differential(args.job_id)
    finally:
        runner.executor.shutdown(wait=True)


if __name__ == '__main__':
    main()
