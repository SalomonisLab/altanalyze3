"""Compare synthetic serving payloads with a specified committed implementation.

Prepare ROOT using benchmark_plot_builds.py first. No biological analysis.
"""
import argparse, ast, importlib, json, subprocess, sys, time
from pathlib import Path
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from altanalyze3.components.cellHarmony.webapp.job_bundle import _StoreMatrix
w=importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('root',type=Path)
parser.add_argument('--baseline-ref',default='HEAD')
args=parser.parse_args()
root=args.root; dims=json.loads((root/'dimensions.json').read_text())
n=dims['cells']; genes=dims['genes']
idx,data,ptr=[np.load(root/f'{name}.npy',mmap_mode='r') for name in ('indices','data','indptr')]
owner=SimpleNamespace(n_obs=n,n_vars=genes,_indices=idx,_data=data,_indptr=ptr,_to_h5ad=None,_sparse=True)
owner._dense=lambda j:np.bincount(idx[ptr[j]:ptr[j+1]],weights=data[ptr[j]:ptr[j+1]],minlength=n).astype(np.float32)
a=ad.AnnData(sp.csr_matrix((n,genes),dtype=np.float32),obs=pd.read_pickle(root/'obs.pkl'),var=pd.DataFrame(index=[f'G{j}' for j in range(genes)]))
a.__dict__['_X']=_StoreMatrix(owner)
cache=dict(adata=a,var_names=a.var_names.to_numpy(),obs_names=a.obs_names.to_numpy(),populations=a.obs.state.astype(str).to_numpy(),umap_x=np.zeros(n),umap_y=np.zeros(n),obs_filter_values={},cluster_key='state',sample_field='Library')
w._get_expression_cache=lambda *args,**kwargs:cache
# Use the committed, pre-change implementation as the comparison baseline.
code=subprocess.check_output(['git','show',args.baseline_ref+':altanalyze3/components/cellHarmony/webapp/app.py'],text=True)
module=ast.parse(code)
fn=next(node for node in module.body if isinstance(node,ast.FunctionDef) and node.name=='_build_expression_payload')
namespace=dict(vars(w));exec(compile(ast.Module(body=[fn],type_ignores=[]),'committed baseline','exec'),namespace)
old=namespace['_build_expression_payload']
report={'scope':'Synthetic expression serving; all input observations retained','cells':n,'features':genes,'baseline_ref':args.baseline_ref}
for label,build in [('baseline',old),('optimized',w._build_expression_payload)]:
 timings=[]
 for repeat in range(3):
  start=time.perf_counter(); payload=build(None,{},'G0',view='violin');timings.append(time.perf_counter()-start)
 report[label+'_violin_seconds']=timings
 if label=='baseline':previous=payload
 else:assert previous==payload
report['plotted_states']=len(payload['violin']);report['plotted_observations']=sum(len(p['values']) for p in payload['violin'])
(root/'details.json').write_text(json.dumps(report,indent=2)+'\n')
fn=next(node for node in module.body if isinstance(node,ast.FunctionDef) and node.name=='_gene_state_stats')
namespace=dict(vars(w));exec(compile(ast.Module(body=[fn],type_ignores=[]),'committed baseline','exec'),namespace)
old=namespace['_gene_state_stats']; wanted=[f'G{i}' for i in range(12)]
for label,build in [('baseline',old),('optimized',w._gene_state_stats)]:
 timings=[]
 for repeat in range(3):
  start=time.perf_counter(); payload=build(cache,wanted);timings.append(time.perf_counter()-start)
 report[label+'_dotplot_seconds']=timings
 if label=='baseline': previous=payload
 else: assert previous==payload
(root/'details.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
