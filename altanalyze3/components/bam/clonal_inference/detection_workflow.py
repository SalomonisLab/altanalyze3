"""Reviewed detection-aware tree workflow; unchanged 2026-11-16 inference core."""
from pathlib import Path
import json,time,platform,os
import numpy as np,pandas as pd
from .integrity import load_inputs,require,numerical_gate,sha256
from .detection_tree import detector_evidence,TreeSearch
from .phylogeny_validation import validate_gain_tree
from .detection_reporting import POLICY, DESCRIPTION, tables
from .evidence_review import evidence_consistency_audit
from .lineage_support import supported_lineage


def specification(seed=222,sex=None,pooling_strength=2.,concentration=2.):
    return dict(method='detection-tree',seed=seed,concentration=float(concentration),FP=.01,FN=.30,
                pooling_strength=float(pooling_strength),sex=sex,max_tree_rounds=30,
                tree='gain-only, full-state search then BIC sampled-state pruning',
                frequency=DESCRIPTION, reporting_policy=POLICY,
                core_sha256=sha256(Path(__file__).with_name('detection_tree.py')))


def run(args):
    spec=specification(args.seed,args.sex,args.pooling_strength,args.rna_concentration)
    require(args.detector is not None,'--detector is required: no categorical substitute for missing detector provenance')
    require(args.sex in ['M','F','U'],'--sex M/F/U required; U permitted only with no chrX event')
    X,c,e,contract=load_inputs(args.matrix,args.cells,args.events,args.contract,'detection-tree')
    require(contract.get('algorithm')==spec,'Detection-tree algorithm differs from explicit contract')
    require(sha256(args.detector)==contract.get('detector_sha256'),'Detector hash differs from contract')
    require(bool(contract.get('authorization')),'Record authorization for the detection-tree workflow')
    z=np.load(args.matrix);dz=np.load(args.detector)
    for k in ['cell_ids','features']:require(np.array_equal(z[k],dz[k]),'Detector identity/order mismatch: '+k)
    require({'A','T'}<=set(z.files),'Direct ALT/total counts required; no synthetic reads')
    A=z['A'];T=z['T'];Z=dz['Z'];available=dz['available']
    require(A.shape==T.shape==Z.shape==available.shape==X.shape,'Count/detector shape mismatch')
    require(np.isfinite(A).all() and np.isfinite(T).all() and ((A>=0)&(T>=A)).all(),'Invalid read counts')
    require(np.isin(Z,[0,1]).all() and available.dtype==bool,'Invalid detector matrix')
    require(not((Z==1)&~available).any(),'Detection outside source availability')
    require(not((Z==1)&(T==0)&(X!=1)).any(),'Detector positives do not reconcile to accepted calls')
    require('chrom' in e and 'event_label' in e,'Chromosome and event labels required')
    require('patient' not in c or c.patient.nunique()==1,'Fit one complete patient series at a time')
    require(not (args.sex=='U' and e.chrom.isin(['X','chrX']).any()),'Unknown sex with chrX event; resolve metadata')
    numerical_gate();out=Path(args.output);require(not out.exists(),'Refuse to overwrite an existing analysis directory')
    samples=c['sample'].drop_duplicates().tolist();sample=c['sample'].map({s:i for i,s in enumerate(samples)}).to_numpy();mu=np.full(len(e),.5)
    if args.sex=='M':mu[e.chrom.isin(['X','chrX']).to_numpy()]=.99
    out.mkdir(parents=True);start=time.time()
    c.to_csv(out/'cell_roster.tsv',sep='\t',index=False);e.to_csv(out/'event_roster.tsv',sep='\t',index=False)
    manifest=dict(algorithm=spec,contract=contract,cells=len(c),events=len(e),samples=samples,
                  python=platform.python_version(),numpy=np.__version__,OPENBLAS_CORETYPE=os.environ.get('OPENBLAS_CORETYPE'),
                  code_sha256={p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')},
                  interpretation='User-reviewed default. Inferred genotypes and point-estimated errors remain model-dependent; no external biological validation claimed.')
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2))
    l0,l1,calibration=detector_evidence(X,A,T,mu,sample,Z,available,spec['concentration'],spec['FP'],spec['FN'],spec['pooling_strength'])
    cal=pd.DataFrame(calibration)
    if len(cal):
        cal['sample']=cal.sample_index.map(dict(enumerate(samples)));cal['event']=cal.event_index.map(dict(enumerate(e.event_label)))
        cal['sparse_mutant_calibration']=cal.effective_mutant_validation_cells<5
        cal['sparse_WT_calibration']=cal.effective_WT_validation_cells<5
        cal['detector_direction_reversed']=cal.TPR<=cal.FPR
    cal.to_csv(out/'detector_calibration.tsv',sep='\t',index=False)
    t=TreeSearch(l0,l1,sample,prune=False);fit=t.run(args.seed,max_rounds=30)
    require(fit['converged'],'Tree search hit round limit; no completed result')
    collapsed=TreeSearch(l0,l1,sample,prune=True);final=collapsed.evaluate(fit['parents'])
    require(final['max_gap']<2e-5,'Mixture optimality check failed')
    validate_gain_tree(final['G'],final['parents'],list(range(len(final['G']))),list(e.feature))
    P,ll=collapsed.posterior(final);require(np.allclose(P.sum(1),1),'Posterior normalization failed')
    G=final['G'][final['active']]
    np.savez_compressed(out/'model.npz',**{k:final[k] for k in ['G','parents','active','weights']},posterior=P,log_evidence=ll,cell_ids=z['cell_ids'],features=z['features'])
    pd.DataFrame(final['G'],columns=e.feature).assign(node=np.arange(len(final['G'])),parent=final['parents']).to_csv(out/'tree_genotypes.tsv',sep='\t',index=False)
    model=dict(final,posterior=P,cell_ids=z['cell_ids'],features=z['features'])
    Q, assignments, F, audit = tables(X,c,e,model)
    F.to_csv(out/'clone_frequencies.tsv',sep='\t',index=False)
    assignments.to_csv(out/'cell_assignments.tsv.gz',sep='\t',index=False)
    audit.to_csv(out/'assignment_evidence_audit.tsv',sep='\t',index=False)
    review_model = dict(model, log_evidence=ll)
    event_review, pair_review, review_status = evidence_consistency_audit(X,c,e,review_model,assignments)
    event_review.to_csv(out/'clone_event_review.tsv',sep='\t',index=False)
    pair_review.to_csv(out/'clone_pair_review.tsv',sep='\t',index=False)
    # Non-mutating cross-check: this does not replace assignments or frequencies.
    support = supported_lineage(X, final['G'], final['active'])
    lineage_review = assignments[['cell_id', 'sample', 'node']].copy()
    lineage_review['minimum_supported_node'] = support['node']
    lineage_review['incompatible_positive_combination'] = support['incompatible']
    lineage_review['ambiguous_lineage'] = support['ambiguous_minimum']
    lineage_review['reference_RNA_conflicts'] = support['negative_evidence_conflicts']
    lineage_review.to_csv(out/'minimum_supported_lineage_audit.tsv.gz',sep='\t',index=False)
    review_status['minimum_supported_lineage_disagreements'] = int((assignments.node.to_numpy() != support['node']).sum())
    review_status['incompatible_positive_combinations'] = int(support['incompatible'].sum())
    review_status['minimum_lineage_is_not_exclusive_genotype'] = True
    (out/'biological_review.json').write_text(json.dumps(review_status,indent=2))
    np.savez_compressed(out/'reporting_memberships.npz',memberships=Q,nodes=final['active'],
                        cell_ids=z['cell_ids'],features=z['features'])
    (out/'reporting_policy.json').write_text(json.dumps(dict(policy=POLICY,description=DESCRIPTION),indent=2))
    summary=dict(score=final['score'],loglik=final['loglik'],max_gap=final['max_gap'],search=fit['trace'],tree_evaluations=fit['evaluations'],active_nodes=final['active'].tolist(),elapsed_seconds=time.time()-start,
                 cells=len(c),events=len(e),samples=samples,method='detection-tree',seed=args.seed,globally_uncalled_events=e.feature[(X>=0).sum(0)==0].tolist(),
                 calibration_flags='Descriptive flags, not exclusions or confidence thresholds: effective class size <5; TPR<=FPR',
                 reporting_policy=POLICY, interpretation=DESCRIPTION+'; source calls and native model posterior unchanged; missing genotypes not filled.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2));(out/'COMPLETE.json').write_text(json.dumps(summary,indent=2));(out/'COMPLETE').write_text('Input, numerical, ancestry and denominator checks passed; external biological validation not implied.')
    return summary
