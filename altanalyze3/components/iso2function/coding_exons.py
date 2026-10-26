"""Audited coding-exon maps using the existing registered-structure ISV index.

Raw GFF/CDS files are read once to attest the required indexed records. Plotting
subsequently reads only the small, hashed map. No sequence-change color overrides.
"""
import colorsys
import csv
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

from ..annotation.domain_impact import coding_blocks
from ..visualization.isv_web.precompute_gff_structures import StructureIndex,parse_attrs,strip_version,_norm_chrom
from . import config
from .protein_impact import dump,sha,table,write,op

SCHEMA='iso2function-coding-exon-map-v1'


def assign_exon_labels(model, reference):
    """Annotate existing CDS blocks with AltAnalyze reference exon-group IDs.

    E4.1/E4.2 are regions of exon group E4, not different physical exons.
    A block spanning retained intronic sequence keeps both E and I identifiers.
    Never invent a sequential number for a novel or unresolved interval.
    A reference row carrying label_prefix belongs to a partner gene of a multi-gene
    model; its IDs stay gene-qualified (CLK1:E4, ENSG...:E4.1) so they are never
    confused with the registered gene's own exon numbers.
    """
    reference=[r for r in reference if _norm_chrom(r['chrom'])==_norm_chrom(model['chrom'])
               and r['strand']==model['strand']]
    for block in model['segments']:
        lo,hi=block['visible_genomic_start'],block['visible_genomic_end']
        hits=sorted([r for r in reference if r['start']<=hi and r['end']>=lo],
                    key=lambda r:(r['start'],r['end'],r['exon_id']),reverse=model['strand']=='-')
        covered=lo
        for start,end in sorted((max(lo,r['start']),min(hi,r['end'])) for r in hits):
            if start>covered:break
            covered=max(covered,end+1)
        parents=[]
        for r in hits:
            match=re.fullmatch(r'([EIU]\d+)\.\d+',r['exon_id'])
            if not match:raise ValueError('Unrecognized AltAnalyze exon-region ID: '+r['exon_id'])
            parent=r['label_prefix']+':'+match[1] if r.get('label_prefix') else match[1]
            if parent not in parents:parents.append(parent)
        block.update(altanalyze_region_ids=[r['gene']+':'+r['exon_id'] if r.get('label_prefix') else r['exon_id'] for r in hits],
                     altanalyze_exon_ids=parents,
                     exon_label='/'.join(parents) if covered>hi else '',
                     exon_label_status='assigned' if covered>hi else 'unresolved_reference_coverage')
    return model


def _partners(model):
    """(partner ENSG, display symbol) pairs of a multi-gene model; empty for single-gene models."""
    genes=[g for g in str(model.get('partner_ensg') or '').split(';') if g]
    symbols=[s for s in str(model.get('partner_gene_symbol') or '').split(';') if s]
    if symbols and len(symbols)!=len(genes):raise ValueError('partner_gene_symbol must align with partner_ensg: '+model['isoform_id'])
    return list(zip(genes,symbols or genes))


def annotate_exon_map(mapped,index_path):
    """Read required genes from ISV's gene_model schema without AnnData loading."""
    index_path=Path(config.require(index_path,'ISV AltAnalyze gene-model index'))
    references={}
    connection=sqlite3.connect(index_path.resolve().as_uri()+'?mode=ro',uri=True)
    try:
        for gene in sorted({r['ensg'] for r in mapped}|{g for r in mapped for g,_ in _partners(r)}):
            # Same per-gene query and normalized intervals as ISV's
            # read_gene_model_from_sqlite; no full database or matrix read.
            rows=connection.execute('SELECT exon_id,chrom,strand,start,end FROM gene_model WHERE gene=?',(gene,)).fetchall()
            segments=[dict(gene=gene,exon_id=e,chrom=c,strand=s,start=min(int(a),int(b)),
                           end=max(int(a),int(b)),type=e[0].upper()) for e,c,s,a,b in rows]
            references[gene]=sorted(segments,key=lambda r:(r['exon_id'],r['start'],r['end']))
    finally:connection.close()
    for model in mapped:
        assign_exon_labels(model,references[model['ensg']]+[dict(r,label_prefix=symbol)
                                                            for gene,symbol in _partners(model) for r in references[gene]])
    digest=__import__('hashlib').sha256(json.dumps(references,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return dict(index_path=str(index_path),queried_gene_rows_sha256=digest,
                reference_rows=sum(map(len,references.values())),
                unresolved_blocks=sum(e['exon_label_status']!='assigned' for m in mapped for e in m['segments']),
                label_rule='AltAnalyze exon-group IDs (E4 from E4.1/E4.2); I marks reference intronic sequence; combined blocks retain every group; never renumbered per isoform'),references


def peptide_segments(exons,cds_min,cds_max,strand,protein_aa):
    """Clip CDS with the existing baseline, then remove only the terminal stop.

    Protein boundaries are continuous 0-based AA edges. Fractional edges retain
    codons split across exons instead of rounding them to different residues.
    """
    if strand not in ('+','-'):raise ValueError('Unresolved coding strand')
    blocks=sorted(coding_blocks(exons,cds_min,cds_max),reverse=strand=='-')
    expected=3*int(protein_aa);total=sum(hi-lo+1 for lo,hi in blocks)
    if total not in (expected,expected+3):raise ValueError('Genomic coding length differs from verified peptide')
    offset=0;out=[]
    for lo,hi in blocks:
        nt=min(hi-lo+1,max(0,expected-offset))
        if nt:
            visible_lo,visible_hi=(lo,lo+nt-1) if strand=='+' else (hi-nt+1,hi)
            out.append(dict(genomic_start=lo,genomic_end=hi,visible_genomic_start=visible_lo,
                visible_genomic_end=visible_hi,protein_nt_start0=offset,protein_nt_end0=offset+nt,
                aa_start0=offset/3,aa_end0=(offset+nt)/3))
            offset+=nt
    if offset!=expected:raise ValueError('Coding-exon map does not cover complete peptide')
    return out,total-expected


def color_exons(models):
    """Shared overlapping genomic coding-exon groups, ordered 5' to 3'."""
    groups=defaultdict(list)
    for m in models:groups[(m['gene_symbol'],m['chrom'],m['strand'])].append(m)
    for (_,_,strand),members in groups.items():
        intervals=sorted({(e['genomic_start'],e['genomic_end']) for m in members for e in m['segments']})
        merged=[]
        for lo,hi in intervals:
            if merged and lo<=merged[-1][1]+1:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
            else:merged.append((lo,hi))
        if strand=='-':merged.reverse()
        for m in members:
            for e in m['segments']:
                hits=[i for i,(lo,hi) in enumerate(merged) if lo<=e['genomic_start']<=e['genomic_end']<=hi]
                if len(hits)!=1:raise ValueError('Ambiguous coding-exon color group')
                i=hits[0];rgb=colorsys.hsv_to_rgb((.98+.88*i/max(len(merged)-1,1))%1,.65,.88)
                e.update(exon_group=i+1,color='#'+''.join(f'{round(c*255):02x}' for c in rgb))
    return models


def _source(impact,name):
    hits=[(p,h) for p,h in impact['inputs'].items() if Path(p).name==name]
    if len(hits)!=1:raise ValueError('Exactly one pinned source required: '+name)
    raw,digest=hits[0];path=Path(config.require(raw,name))
    return path,digest


def build_coding_exon_map(base,impact,models,sequences,index_path):
    """Use ISV geometry; independently attest it once against the pinned GFF."""
    audit_path=base/'native_genomic_audit.tsv'
    if sha(audit_path)!=impact['outputs'].get(audit_path.name):raise ValueError('Changed native genomic audit')
    audited={r['isoform_id']:r for r in table(audit_path)}
    wanted={r['isoform_id']:r for r in models}
    if audited.keys()!=wanted.keys():raise ValueError('Incomplete genomic audit roster')
    index=StructureIndex(str(config.require(index_path,'registered ISV structure index')))
    indexed={}
    try:
        for tid,m in wanted.items():
            # Do not permit an arbitrary first row or a bare-molecule collision.
            rows=index._con.execute('SELECT chrom,strand,exons FROM structures WHERE gene=? AND tid_base=?',
                                    (m['ensg'],strip_version(tid))).fetchall()
            if len(set(rows))!=1:raise ValueError('Exact indexed native identity unresolved or ambiguous: '+tid)
            geometry=index.geometry(m['ensg'],tid)
            if not geometry:raise ValueError('ISV registered geometry unresolved: '+tid)
            indexed[tid]=geometry
    finally:index._con.close()
    gff,gff_hash=_source(impact,'combined.gff.gz');coding,coding_hash=_source(impact,'coding_regions.txt')
    for p,d in [(gff,gff_hash),(coding,coding_hash)]:
        if sha(p)!=d:raise ValueError('Changed source for indexed map: '+str(p))
    # One-time attestation, never performed by the plotting endpoint.
    registered=defaultdict(list);metadata={}
    with op(gff) as fh:
        for line in fh:
            if line.startswith('#'):continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9 or p[2]!='exon':continue
            attr=parse_attrs(p[8]);tid=strip_version(attr.get('transcript_id',''))
            if tid not in wanted:continue
            if strip_version(attr.get('gene_id',''))!=wanted[tid]['ensg']:raise ValueError('Registered gene differs: '+tid)
            info=(p[0],p[6])
            if tid in metadata and metadata[tid]!=info:raise ValueError('Mixed genomic loci for native protein: '+tid)
            metadata[tid]=info;registered[tid].append((int(p[3]),int(p[4])))
    cds={}
    with op(coding) as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            tid=r['Transcript ID']
            if tid in wanted:
                if tid in cds:raise ValueError('Duplicate native CDS record: '+tid)
                cds[tid]=r
    if registered.keys()!=wanted.keys() or cds.keys()!=wanted.keys():raise ValueError('Incomplete original native GFF/CDS roster')
    out=[]
    for tid,m in wanted.items():
        g=indexed[tid];exons=sorted(registered[tid]);c=cds[tid];old=audited[tid]
        if sorted(g['exons'])!=exons or (g['chrom'],g['strand'])!=metadata[tid]:
            raise ValueError('ISV indexed geometry differs from pinned native GFF: '+tid)
        if (strip_version(c['Gene ID'])!=m['ensg'] or _norm_chrom(c['Chromosome'])!=_norm_chrom(metadata[tid][0])
                or c['Strand']!=metadata[tid][1]):
            raise ValueError('Native CDS locus differs from registered model: '+tid)
        aa=len(sequences[tid])
        if (int(c['Protein Length'])!=int(old['coding_regions_aa']) or c['NMD Status']!=old['coding_regions_nmd']
                or aa!=int(old['protein_fasta_aa']) or m['nmd']!=old['source_table_nmd']):
            raise ValueError('Native CDS/protein metadata differs from the existing audit: '+tid)
        segments,stop=peptide_segments(g['exons'],int(c['CDS Genomic Min']),int(c['CDS Genomic Max']),g['strand'],aa)
        if len(exons)!=int(old['native_exon_count']) or aa*3+stop!=int(old['genomic_cds_nt']):
            raise ValueError('Existing genomic baseline differs: '+tid)
        out.append(dict(gene_symbol=m['gene_symbol'],isoform_id=tid,ensg=m['ensg'],chrom=g['chrom'],strand=g['strand'],
                        protein_aa=aa,reported_aa=m['reported_aa'],coding_regions_aa=c['Protein Length'],
                        source_nmd=m['nmd'],coding_regions_nmd=c['NMD Status'],
                        protein_sha256=m['protein_sha256'],terminal_stop_nt=stop,segments=segments))
        if m.get('partner_ensg'):
            out[-1].update(partner_ensg=m['partner_ensg'],partner_gene_symbol=m.get('partner_gene_symbol',''))
    return color_exons(out),dict(gff_sha256=gff_hash,coding_sha256=coding_hash,index_path=str(index_path),
        indexed_rows_attested=len(out),index_row_sha256=__import__('hashlib').sha256(
            json.dumps(indexed,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        source_scans='One-time pinned GFF/CDS attestation; subsequent plots read only this map')


def run_protein_exon_map(args):
    from .plotting.protein_features import load_audit
    cached=getattr(args,'coding_exon_map',None)
    if not args.impact_results or not args.artifacts_dir or not (cached or args.isv_structure_index):
        raise ValueError('protein-exon-map requires --impact-results --artifacts-dir and --isv-structure-index or --coding-exon-map')
    if cached and not getattr(args,'gene_model_index',None):
        raise ValueError('Cached exon-map enrichment requires --gene-model-index')
    out=Path(args.artifacts_dir);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use an empty exon-map output directory')
    dump(out/'exon_map_manifest.json',dict(schema=SCHEMA,status='incomplete'))
    base,impact,models,sequences,_,_=load_audit(args.impact_results)
    if cached:
        previous,parent,parent_manifest=load_coding_exon_map(cached,base,models)
        mapped=list(previous.values());sources=dict(parent_manifest['sources'])
        sources['parent_map']=dict(path=str(parent),manifest_sha256=sha(parent/'exon_map_manifest.json'))
        sources['source_scans']='None: verified complete parent coding-exon map reused'
    else:
        mapped,sources=build_coding_exon_map(base,impact,models,sequences,args.isv_structure_index)
    if getattr(args,'gene_model_index',None):
        sources['altanalyze_exon_labels'],reference_rows=annotate_exon_map(mapped,args.gene_model_index)
        dump(out/'altanalyze_gene_model_rows.json',reference_rows)
    dump(out/'coding_exon_map.json',mapped)
    rows=[{k:json.dumps(v,separators=(',',':')) if isinstance(v,list) else v
           for k,v in dict(**{k:v for k,v in m.items() if k!='segments'},**e).items()}
          for m in mapped for e in m['segments']]
    write(out/'coding_exon_coordinates.tsv',rows)
    dump(out/'exon_map_manifest.json',dict(schema=SCHEMA,status='complete',required_roster=impact['required_roster'],
        source_manifest_sha256=sha(base/'impact_manifest.json'),sources=sources,code_sha256=sha(__file__),
        color_rule='Overlapping coding-exon groups on the same gene/chromosome/strand; no sequence-change palette overrides',
        coordinates='1-based inclusive genomic intervals; 0-based continuous AA edges, retaining split codons',
        outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.name!='exon_map_manifest.json'}))
    return 0


def load_coding_exon_map(path,base,models):
    root=Path(config.require(path,'completed coding-exon map'));m=json.loads((root/'exon_map_manifest.json').read_text())
    if m.get('status')!='complete' or m['source_manifest_sha256']!=sha(base/'impact_manifest.json'):
        raise ValueError('Exon map is incomplete or belongs to another impact audit')
    for name,digest in m['outputs'].items():
        if sha(root/name)!=digest:raise ValueError('Changed coding-exon map')
    rows=json.loads((root/'coding_exon_map.json').read_text())
    if [(r['gene_symbol'],r['isoform_id']) for r in rows]!=[(r['gene_symbol'],r['isoform_id']) for r in models]:
        raise ValueError('Complete native coding-exon roster not preserved')
    for old,row in zip(models,rows):
        if old['protein_sha256']!=row['protein_sha256'] or int(old['actual_aa'])!=row['protein_aa']:
            raise ValueError('Coding-exon map protein differs from native model')
    return {r['isoform_id']:r for r in rows},root,m
