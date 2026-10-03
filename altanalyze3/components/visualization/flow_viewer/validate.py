"""Check every cell space, array shape, label code and gate mask before serving."""
import argparse,json,sys
from pathlib import Path
import numpy as np


def validate(root):
    root=Path(root);m=json.loads((root/'manifest.json').read_text());fail=[];checks=0
    spaces=m.get('spaces') or {'flow':dict(n=m['n_events'],features=m['channels'],features_file=m['channels_file'],embeddings=m['embeddings'],labels=m['labels'])}
    def check(name,ok):
        nonlocal checks
        checks+=1
        if not ok:fail.append(name);print('FAIL',name)
    def array(spec,dtype,shape):
        p=root/'arrays'/spec['file'];expected=int(np.prod(shape))*np.dtype(dtype).itemsize
        if not p.exists() or p.stat().st_size!=expected:return None
        if not expected:return np.empty(shape,dtype)
        return np.memmap(p,dtype,'r',shape=shape)
    for name,s in spaces.items():
        n=int(s['n']);features=s['features'];check(name+' unique features',len(features)==len(set(features)))
        x=array({'file':s['features_file']},np.float32,(n,len(features)))
        check(name+' feature matrix',x is not None and np.isfinite(x).all())
        for en,e in s['embeddings'].items():
            x=array(e,np.float32,(n,2));check(name+' embedding '+en,x is not None and np.isfinite(x).all())
        for ln,l in s['labels'].items():
            x=array(l,np.int16,(n,));check(name+' labels '+ln,x is not None and ((x>=-1)&(x<len(l['levels']))).all())
        for key in ['default_embedding','default_label']:
            check(name+' '+key,key not in s or s[key] in s['embeddings' if key=='default_embedding' else 'labels'])
        if 'cell_ids_file' in s:
            ids=(root/'arrays'/s['cell_ids_file']).read_text().splitlines();check(name+' cell IDs',len(ids)==n and len(set(ids))==n)
        for ln,prefix in s.get('transfer_links',{}).items():
            check(name+' transfer link '+ln,ln in s['labels'] and any(k.startswith(prefix) for k in spaces['flow']['labels']))
    for gn,g in m.get('gatesets',{}).items():
        n=spaces[g.get('space','flow')]['n']
        for nd in g['nodes']:
            x=array({'file':nd['mask_file']},np.uint8,(n,));check(gn+' mask '+nd['path'],x is not None and ((x==0)|(x==1)).all())
    print(f'{checks-len(fail)} / {checks} bundle checks passed')
    return fail


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True);a=p.parse_args();sys.exit(bool(validate(a.bundle)))
if __name__=='__main__':main()
