"""Compare legacy status-log reads with bounded snapshots on a large saved log."""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','baseline','bounded'])
    parser.add_argument('path',type=Path)
    parser.add_argument('--rounds',type=int,default=20)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.mode=='prepare':
        with args.path.open('w') as fh:
            fh.write('adata shape: (128388, 32738)\nCells remaining after min_genes 500 filtering: 125000\n[INFO] Applied min_alignment_score=0.4. Excluded 5000 cells, kept 120000.\n')
            chunk='[INFO] Processing a completed gene feature block for saved results.\n'*1000
            for _ in range(1500):
                fh.write(chunk)
        return
    from altanalyze3.components.cellHarmony.webapp.log_snapshot import read_pipeline_log
    def baseline(path):
        lines=path.read_text(encoding='utf-8').splitlines(True)
        markers=('adata shape:','Cells remaining after','Applied min_alignment_score=','Auto-selected rho for library')
        return lines[:80],lines[-200:],[line for line in lines if any(marker in line for marker in markers)]
    reader=baseline if args.mode=='baseline' else read_pipeline_log
    times=[]
    for _ in range(args.rounds):
        start=time.perf_counter()
        result=reader(args.path)
        times.append(time.perf_counter()-start)
    divisor=1024**3 if sys.platform=='darwin' else 1024**2
    report=dict(mode=args.mode,log_bytes=args.path.stat().st_size,rounds=args.rounds,first_read_seconds=times[0],
                total_seconds=sum(times),repeat_mean_seconds=sum(times[1:])/max(1,len(times)-1),
                peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/divisor,
                snapshot_sha256=hashlib.sha256(json.dumps(result).encode()).hexdigest())
    if args.output:
        args.output.write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
