"""General source preparation for the opt-in terminal-evidence collapse path."""
from pathlib import Path
import csv,json,multiprocessing,os
from .. import io_utils as IO
from .source_evidence import load_reference_models
from .source_stream import build_event_index,scan_sample
from .stream_context import build_contexts


_REGIONS = _INDEXES = None


def _scan_task(args):
    sample,regions,indexes,folder,policy=args
    return scan_sample(sample,regions,indexes,folder,policy)


def _fork_scan_task(args):
    sample,folder,policy=args
    return scan_sample(sample,_REGIONS,_INDEXES,folder,policy)


def prepare_streaming(samples,exon_reference,enst_cache,ref_gff,workdir,policy,*,nproc=1,log=print):
    """Measure complete source inputs, then form classes without altering WTA."""
    names=[s[0] for s in samples]
    if not names or len(set(names))!=len(names):
        raise ValueError('Distinct source identifiers required; preserve and disambiguate every original metadata row')
    work=Path(workdir);work.mkdir(parents=True,exist_ok=True)
    regions={}
    with IO.smart_open(str(exon_reference)) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            a,b=int(r['exon-region-start(s)']),int(r['exon-region-stop(s)'])
            value=min(a,b),max(a,b);g=r['gene'];k=r['exon-id']
            if k in regions.setdefault(g,{}) and regions[g][k]!=value:
                raise ValueError('Conflicting required native annotation region')
            regions[g][k]=value
    refs=load_reference_models(str(enst_cache),str(ref_gff))
    indexes={g:build_event_index(m,{s:r['exons'] for s,r in refs.get(g,{}).items()}) for g,m in regions.items()}
    tasks=[(s,regions,indexes,str(work/f'sample_{i:04d}'),policy) for i,s in enumerate(samples)]
    if nproc>1 and os.name=='posix':
        # On Linux fork shares the immutable annotation/ref index. Existing
        # platform callers can select serial mode when parallelism is unavailable.
        global _REGIONS, _INDEXES
        _REGIONS, _INDEXES = regions, indexes
        try:
            with multiprocessing.get_context('fork').Pool(min(nproc,len(samples))) as pool:
                scans=pool.map(_fork_scan_task,[(s,str(work/f'sample_{i:04d}'),policy)
                                              for i,s in enumerate(samples)],chunksize=1)
        finally:
            _REGIONS = _INDEXES = None
    else:scans=[_scan_task(a) for a in tasks]
    if [s['library'] for s in scans]!=names or not all(s['complete'] for s in scans):
        raise ValueError('Complete ordered source identity coverage failed')
    context=build_contexts(scans,regions,refs,work/'pooled',policy,log=log)
    context['reference_models']=refs
    registry={s['library']:s['database'] for s in scans}
    return context,registry
