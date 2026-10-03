"""Compare repeated H5AD group-menu construction without retaining AnnData."""
import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['baseline','cached'])
    parser.add_argument('h5ad',type=Path)
    parser.add_argument('--rounds',type=int,default=20)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from altanalyze3.components.cellHarmony.flask.pipeline import _candidate_group_fields
    from altanalyze3.components.cellHarmony.webapp.metadata_options import cached_group_fields
    reader=_candidate_group_fields if args.mode=='baseline' else cached_group_fields
    times=[]
    for _ in range(args.rounds):
        start=time.perf_counter()
        result=reader(args.h5ad,preferred=['Library','group','sample'],max_categories=None)
        times.append(time.perf_counter()-start)
    divisor=1024**3 if sys.platform=='darwin' else 1024**2
    report=dict(mode=args.mode,rounds=args.rounds,total_seconds=sum(times),first_read_seconds=times[0],
                repeat_mean_seconds=sum(times[1:])/max(1,len(times)-1),
                peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/divisor,
                menu_sha256=hashlib.sha256(json.dumps(result).encode()).hexdigest())
    args.output.write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
