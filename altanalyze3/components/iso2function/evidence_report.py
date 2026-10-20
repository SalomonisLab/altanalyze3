"""One reproducible report for native proteins, assayed proteins and predicted features."""
import html
import json
from collections import defaultdict
from pathlib import Path

from .protein_impact import dump, sha, table, write


def run_report(args):
    if not args.artifacts_dir or not any((args.impact_results,args.query_results,args.assay_results)):
        raise ValueError('report requires --artifacts-dir and at least one completed result directory')
    out=Path(args.artifacts_dir);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise ValueError('Use an empty report output directory')
    sources={};data={};manifests={}
    for label,path,manifest_name,files in [
        ('domains',args.impact_results,'impact_manifest.json',['native_model_audit.tsv','protein_catalog_audit.tsv','pair_summary.tsv','pair_feature_impacts.tsv','prior_feature_regression.tsv']),
        ('interactions',args.query_results,'query_manifest.json',['all_pair_decisions.tsv','divergent_interactions.tsv']),
        ('assay_comparisons',args.assay_results,'assay_match_manifest.json',['all_native_pair_evidence.tsv','native_assayed_protein_comparison.tsv','amino_acid_differences.tsv'])]:
        if not path: continue
        base=Path(path);mfile=base/manifest_name;m=json.loads(mfile.read_text())
        if not str(m.get('status','')).startswith('complete'):
            raise ValueError('Incomplete report prerequisite: '+str(mfile))
        manifests[label]=m;sources[str(mfile)]=sha(mfile);data[label]={}
        for name in files:
            p=base/name;sources[str(p)]=sha(p)
            if name in m.get('outputs',{}) and sources[str(p)]!=m['outputs'][name]:
                raise ValueError('Changed report source: '+str(p))
            data[label][name]=table(p)
    genes=set()
    native=defaultdict(list);catalog={};features=defaultdict(list);interactions=defaultdict(list);assays=defaultdict(list)
    supported_pairs=defaultdict(list);query_proteins=defaultdict(list);seen_proteins=set()
    for r in data.get('domains',{}).get('native_model_audit.tsv',[]):
        genes.add(r['gene_symbol']);native[r['gene_symbol']].append(r)
    for r in data.get('domains',{}).get('protein_catalog_audit.tsv',[]):catalog[r['isoform_id']]=r
    for r in data.get('domains',{}).get('pair_feature_impacts.tsv',[]):features[r['gene_symbol']].append(r)
    for r in data.get('interactions',{}).get('all_pair_decisions.tsv',[]):
        if r['native_supported']=='True':
            genes.add(r['gene']);supported_pairs[r['gene']].append(r)
            threshold=manifests['interactions']['parameters']['min_isoform_reads']
            for side in ('A','B'):
                for alias in json.loads(r['native_aliases_'+side]):
                    key=(r['gene'],r['clone_'+side],alias['isoform'])
                    if alias['reads'] < threshold or key in seen_proteins:continue
                    seen_proteins.add(key)
                    changes='; '.join(str(x['position'])+': '+x['native']+' to '+x['assayed'] for x in alias.get('substitutions',[]))
                    query_proteins[r['gene']].append(dict(native=alias['isoform'],assayed_protein=r['clone_'+side],reads=alias['reads'],
                                                         sequence_match=alias.get('match_type','exact'),changes=changes))
    for r in data.get('interactions',{}).get('divergent_interactions.tsv',[]):interactions[r['gene']].append(r)
    for r in data.get('assay_comparisons',{}).get('all_native_pair_evidence.tsv',[]):
        genes.add(r['gene']);assays[r['gene']].append(r)
    summary=[]
    for gene in sorted(genes):
        edges=interactions[gene]
        ppi=sorted({r['target'] for r in edges if r['modality']=='ppi'})
        pdi=sorted({r['target'] for r in edges if r['modality']=='pdi'})
        summary.append(dict(gene=gene,original_panel_proteins=len(native[gene]),
                            divergent_PPI_partners=len(ppi),PPI_partners=';'.join(ppi),
                            divergent_DNA_baits=len(pdi),DNA_baits=';'.join(pdi),
                            requested_assay_comparisons=len(assays[gene]),
                            assay_pairs_with_AA_differences=sum(r['both_native_proteins_exact']=='False' for r in assays[gene])))
    write(out/'gene_evidence_summary.tsv',summary)
    esc=lambda x:html.escape(str(x))
    def tbl(rows,cols):
        if not rows:return '<p>No records in the supplied results.</p>'
        return '<table><thead><tr>'+''.join('<th>'+esc(label)+'</th>' for key,label in cols)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(r.get(key,''))+'</td>' for key,label in cols)+'</tr>' for r in rows)+'</tbody></table>'
    parts=['<!doctype html><html><meta charset="utf-8"><title>HSPC TF isoform evidence</title><style>body{font:15px/1.5 system-ui;margin:28px;color:#17212b}table{border-collapse:collapse;width:100%;font-size:12px;margin:16px 0}th,td{padding:7px;border:1px solid #cad2da;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#e9f0f6}tr:nth-child(even){background:#f7f9fa}h2{margin-top:32px}p{max-width:1200px}code{overflow-wrap:anywhere}</style><h1>HSPC transcription-factor isoforms: assays and predicted protein changes</h1>',
           '<p>Experimental results describe the proteins used in the atlas assays. A matching native protein supports applying that evidence to the native sequence; small amino-acid differences support a prediction, not a direct measurement of that native protein. None of these assays was performed in the selected HSPCs.</p>',
           '<p>The interaction query uses the recorded read-support and partner-expression thresholds. Protein partner CPM is a donor mean in the selected HSC/MPP states. DNA baits are assayed sequences, have no gene CPM, and do not establish target activation or repression. The complete query comparison universe remains in all_pair_decisions.tsv; this report displays genes with native RNA support, all original protein-impact models, and all explicitly requested native/assayed protein comparisons. No requested comparison is removed for an amino-acid mismatch.</p>',
           '<p>Protein-feature calls are predictions from pinned annotation databases. An altered or partial domain is not the complete absence of that domain; a missing feature sequence does not prove loss of activity. Features absent from both proteins are not differential losses. Feature counts refer to annotation records, not independent domains.</p>',
           tbl(summary,[('gene','Gene'),('original_panel_proteins','Original protein models'),('divergent_PPI_partners','Different PPI partners passing expression filter'),('divergent_DNA_baits','Different DNA baits'),('requested_assay_comparisons','Specified native/assayed comparisons'),('assay_pairs_with_AA_differences','Comparisons containing AA mismatches')])]
    for gene in sorted(genes):
        parts.append('<h2>'+esc(gene)+'</h2>')
        nr=[]
        for r in native[gene]:
            c=catalog.get(r['isoform_id'],{})
            nr.append(dict(r,novelty=c.get('novelty_interpretation','')))
        parts.append(tbl(nr,[('isoform_id','Native isoform'),('actual_aa','AA'),('assigned_reads','Assigned reads'),('nmd','NMD annotation'),('metadata_length_agrees','Reported length agrees'),('novelty','Reference database comparison')]))
        parts.append('<h3>Measured interaction differences and native RNA support</h3>')
        parts.append(tbl(supported_pairs[gene],[('clone_A','Assayed protein A'),('clone_B','Assayed protein B'),
            ('max_single_native_reads_A','Most reads for one native A isoform'),('max_single_native_reads_B','Most reads for one native B isoform'),
            ('native_evidence_A','A sequence match'),('native_evidence_B','B sequence match'),('native_pair_evidence','Native protein evidence'),
            ('qualifying_PPI','PPI differences passing expression filter'),('qualifying_PDI','DNA bait differences'),('retained','Passes interaction filter')]))
        parts.append(tbl(query_proteins[gene],[('native','Native isoform meeting read threshold'),('assayed_protein','Assayed protein'),
            ('reads','Assigned reads'),('sequence_match','Protein sequence match'),('changes','AA position: native to assayed residue')]))
        parts.append(tbl(interactions[gene],[('clone_A','Assayed protein A'),('clone_B','Assayed protein B'),('modality','Assay type'),('target','Protein partner or DNA bait'),('call_A','A assay call'),('call_B','B assay call'),('max_state_mean_cpm','Maximum state mean partner CPM'),('expression_status','Expression status')]))
        parts.append('<h3>Native proteins compared with the proteins used in published assays</h3>')
        parts.append(tbl(assays[gene],[('native_A','Native A'),('native_B','Native B'),('clone_A','Assayed A'),('clone_B','Assayed B'),('exact_native_A','A exact AA match'),('exact_native_B','B exact AA match'),('native_reads_A','A reads'),('native_reads_B','B reads'),('qualifying_PPI','PPI differences passing expression filter'),('qualifying_PDI','DNA bait differences'),('M1H_A','A reporter mean'),('M1H_B','B reporter mean'),('M1H_A_call','A reporter category'),('M1H_B_call','B reporter category')]))
        differences=[r for r in data.get('assay_comparisons',{}).get('amino_acid_differences.tsv',[]) if r['gene']==gene]
        parts.append(tbl(differences,[('native','Native isoform'),('assayed_protein','Assayed protein'),('native_start','Native AA start'),('native_end','Native AA end'),('native_sequence','Native residues'),('assayed_start','Assayed AA start'),('assayed_end','Assayed AA end'),('assayed_sequence','Assayed residues')]))
        # Retain all non-null feature changes; duplicate source rows are folded for
        # display only, never in the exported instance-level audit.
        fs=[];seen=set()
        for r in features[gene]:
            if r['AA_comparison'] not in ('lost_in_B','gained_in_B','unresolved'):continue
            key=tuple(r.get(k,'') for k in ['isoform_A','isoform_B','provider','name','accession','description','reference_start','reference_end','A_status','B_status','AA_comparison'])
            if key not in seen:fs.append(r);seen.add(key)
        parts.append('<h3>Predicted changes to annotated protein features (B relative to A)</h3>')
        parts.append(tbl(fs,[('isoform_A','Native A'),('isoform_B','Native B'),('provider','Annotation source'),('name','Feature'),('accession','Accession'),('description','Description'),('reference_start','Reference AA start'),('reference_end','Reference AA end'),('A_status','A sequence status'),('B_status','B sequence status'),('AA_comparison','Full feature sequence comparison'),('genomic_comparison','Genomic coding-coverage comparison')]))
    parts.append('<h2>Sources and reproducibility</h2><pre>'+esc(json.dumps(sources,indent=2,sort_keys=True))+'</pre></html>')
    (out/'iso2function_evidence.html').write_text(''.join(parts))
    dump(out/'report_manifest.json',dict(status='complete',genes=len(genes),sources=sources,
                                        code_sha256=sha(__file__),outputs={p.name:sha(p) for p in out.iterdir()}))
    return 0
