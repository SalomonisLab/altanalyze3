"""Read-only simulated Chat request benchmark; never launches differential analyses.

python tests/benchmark_chat_queries.py --base http://127.0.0.1:8062 --job Hs-Lung-COPD-metacells --out /tmp/chat-benchmark.json
Use --questions to supply a JSON list of dataset-specific natural-language questions.
"""
import argparse
import json
from pathlib import Path
from statistics import median
from time import perf_counter
import requests

DEFAULT_QUESTIONS=[
 'Which genes in AT2 track FEV1?',
 'Which genes in AT2 track age?',
 'Which genes coexpress with SFTPC in AT2?',
 'Which cell states differ in abundance by copd_status?',
 'Show SFTPC across GOLD stages in AT2',
 'Which cell type is most affected in COPD vs non-COPD?',
 'Is the COPD vs non-COPD signature in AT2 present in all donors or a subset?',
 'Compare annotation agreement in AT2',
 'Show GO processes up in AT2',
]

def benchmark(base,job,questions,repeats=3):
    records=[]
    for question in questions:
        for repeat in range(repeats):
            started=perf_counter()
            response=requests.post(f'{base.rstrip("/")}/api/jobs/{job}/chat',json={'question':question},timeout=60)
            elapsed=perf_counter()-started;response.raise_for_status();payload=response.json()
            records.append({'question':question,'repeat':repeat,'seconds':round(elapsed,6),
                'intent':payload.get('intent'),'status':payload.get('status'),
                'reading':payload.get('reading'),'performance':payload.get('performance'),
                'plot':(payload.get('plot') or {}).get('kind')})
    return {'requests':records,'median_seconds':median(r['seconds'] for r in records),
            'first_request_median_seconds':median(r['seconds'] for r in records if r['repeat']==0),
            'repeat_median_seconds':median(r['seconds'] for r in records if r['repeat']>0)}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',required=True);parser.add_argument('--job',required=True)
    parser.add_argument('--questions',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();questions=json.loads(args.questions.read_text()) if args.questions else DEFAULT_QUESTIONS
    result=benchmark(args.base,args.job,questions);args.out.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='requests'},indent=2))
