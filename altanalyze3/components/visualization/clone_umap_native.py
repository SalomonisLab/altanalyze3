"""Run the installed AltAnalyze2 UMAP interface on saved coordinates, with verified group joins.

AltAnalyze2 is an external, separately installed dependency. No UMAP fitting, expression
replacement, or mutation inference occurs here. The Python-2-compatible driver only adds
SVG export and verifies the exact plotted coordinate multiset and group counts.
"""
from pathlib import Path
import argparse
import collections
import hashlib
import json
import os
import shutil
import subprocess
import pandas as pd


def render(inputs, expression, coordinates, altanalyze2, python2, output, groups=None):
    from ..bam.clonal_inference.integrity import require, sha256
    source, out = Path(inputs), Path(output)
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source / 'mapping_validation.json').read_text())
    require(all(sha256(source / f) == h for f, h in manifest['render_inputs_sha256'].items()), 'Mapping inputs changed')
    d = pd.read_csv(source / 'mapped_umap_cells.tsv.gz', sep='\t').set_index('umap_id')
    with open(expression) as handle:
        header = [x.split(':')[-1] for x in handle.readline().rstrip('\r\n').split('\t')[1:]]
    require(len(header) == len(set(header)) == len(d) and set(header) == set(d.index), 'Expression header and complete coordinate roster differ')
    coords = pd.read_csv(coordinates, sep='\t', header=None, names=['umap_id', 'UMAP1', 'UMAP2'])
    require(coords.umap_id.tolist() == d.index.tolist(), 'Coordinate order changed')
    require((coords[['UMAP1','UMAP2']].to_numpy() == d[['UMAP1','UMAP2']].to_numpy()).all(), 'Coordinate values changed')
    labels = d.loc[header, 'plot_clone'].replace({'Aged healthy':'Healthy'})
    if groups is not None:
        annotation=pd.read_csv(groups,sep='\t',header=None,names=['cell','group','label'],dtype=str).set_index('cell')
        require(annotation.index.is_unique and set(annotation.index)==set(header),'Annotation must cover every coordinate cell exactly once')
        labels=annotation.loc[header,'label']
    require(not labels.str.contains(':').any(), 'Invalid native group label')
    annotation_source=Path(groups) if groups is not None else None
    groups = list(dict.fromkeys(labels))
    group_path = out / 'groups.KINNEX-clones.txt'
    with group_path.open('w') as handle:
        for uid, label in labels.items():
            handle.write('%s\t%s\t%s\n' % (uid, groups.index(label)+1, label))
    # Native AltAnalyze2 locates its groups file using the expression filename.
    exp = out / 'exp.KINNEX-clones.txt'
    if exp.exists() or exp.is_symlink():
        require(exp.is_symlink() and exp.resolve() == Path(expression).resolve(), 'Refuse to replace an existing expression file')
    else:
        exp.symlink_to(Path(expression).resolve())
    shutil.copyfile(coordinates, out / 'coordinates.tsv')
    config = dict(base=str(out.resolve()), app=str(Path(altanalyze2).resolve()), expected_groups={str(k):int(v) for k,v in labels.value_counts().items()})
    (out / 'native_config.json').write_text(json.dumps(config,indent=2))
    for f in ('native_plot_validation.json',):
        if (out/f).exists(): (out/f).unlink()  # stale validation cannot satisfy this run
    env = dict(os.environ, MPLBACKEND='Agg', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', LC_ALL='C', LANG='C')
    driver = Path(__file__).with_name('clone_umap_native_py2.py')
    result = subprocess.run([str(python2),str(driver),str(out/'native_config.json')],env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (out/'native.log').write_text(result.stdout);(out/'native.err').write_text(result.stderr)
    require(result.returncode == 0, 'AltAnalyze2 failed; see native.err and native.log')
    require('Performing UMAP' not in result.stdout, 'AltAnalyze2 unexpectedly fitted UMAP')
    report = json.loads((out/'native_plot_validation.json').read_text())
    require(report['exact_coordinate_multiset_equal'] and report['native_point_count']==len(d), 'Native plotted roster differs')
    require(report['native_group_sizes']==config['expected_groups'], 'Native group counts differ')
    app = Path(altanalyze2)
    report['source_sha256'] = {str(p):sha256(p) for p in [app,app.parent/'visualization_scripts/clustering.py',app.parent/'UI.py',driver,Path(expression),Path(coordinates),group_path]}
    if annotation_source is not None:report['source_sha256'][str(annotation_source)]=sha256(annotation_source)
    (out/'native_plot_validation.json').write_text(json.dumps(report,indent=2))
    return report


def add_arguments(p):
    for name in ['inputs','expression','coordinates','altanalyze2','python2','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--groups',help='Optional complete original annotation table; default is clone assignments')
def run(a):return render(a.inputs,a.expression,a.coordinates,a.altanalyze2,a.python2,a.output,a.groups)
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);add_arguments(p);return run(p.parse_args(argv))
if __name__=='__main__':main()
