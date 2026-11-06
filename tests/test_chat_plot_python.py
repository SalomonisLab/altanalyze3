"""Generated plotting programs get evidence and plotting capabilities, not Python execution."""
import copy
import xml.etree.ElementTree as ET
import pytest
from altanalyze3.components.cellHarmony.webapp.plot_python import Program, render


CODE = '''import matplotlib.pyplot as plt
from matplotlib_venn import venn2
fig, ax = plt.subplots(figsize=(7, 5))
counts = tables["q1"]["count"]
venn2(subsets=(counts[0], counts[1], counts[2]), set_labels=("A", "B"), ax=ax)
ax.set_title("Exact RNA feature overlap")
'''


def test_native_plot_uses_all_exact_counts_and_preserves_sources():
    tables={'q1':{'count':[10,20,3]}};before=copy.deepcopy(tables)
    svg=Program(tables).run(CODE,'Exact feature identities; A only 10, B only 20, shared 3. No statistical test.')
    root=ET.fromstring(svg);texts=[t.text for t in root.iter() if t.tag.endswith('text')]
    assert all(label in texts for label in ['10','20','3','A','B'])
    assert tables==before


def test_plotting_runs_in_bounded_child_and_returns_vector_figure():
    svg=render(CODE,{'q1':{'count':[1,1,1]}},'All requested identities retained; descriptive overlap.')
    assert '<svg' in svg and 'All requested identities retained' in svg


@pytest.mark.parametrize('code',[
    'import os', 'import subprocess', 'import pathlib',
    'open("/tmp/forbidden", "w")', '__import__("os")',
    'import matplotlib.pyplot as plt\nplt.savefig("/tmp/forbidden")',
    'import matplotlib.pyplot as plt\nplt.imread("/etc/passwd")',
    'import seaborn as sns\nsns.load_dataset("penguins")',
    'import matplotlib.pyplot as plt\nplt.__dict__',
    'import matplotlib.pyplot as plt\nfig, ax = plt.subplots()\nfig.savefig("/tmp/forbidden")',
    'from matplotlib_venn import __dict__',
    'while True: pass', 'def f(): pass',
    'import matplotlib.pyplot as plt\nplt.figure(figsize=(10000,10000))',
    'import matplotlib.pyplot as plt\nplt.figure(dpi=2400)',
    'x = [1] * 1000000000', 'x = 10 ** 1000000000',
    'tables["q1"]["count"] = [999,999,999]',
    'from matplotlib_venn import venn2\nvenn2(subsets=(999,999,999))',
])
def test_files_network_reflection_unbounded_allocations_and_fabricated_plots_are_unavailable(code):
    with pytest.raises((ValueError,KeyError)):
        Program({'q1':{'count':[1,1,1]}}).run(code,'Test')


def test_histogram_uses_full_evidence_vector():
    values=list(range(1000))
    code='import matplotlib.pyplot as plt\nplt.figure(figsize=(6,4))\nplt.hist(tables["q1"]["values"],bins=20)\nplt.xlabel("Stored values")'
    svg=Program({'q1':{'values':values}}).run(code,'All 1,000 supplied values; no sampling or test.')
    assert 'Stored values' in svg and 'All 1,000 supplied values' in svg


def test_literal_vectors_require_a_unique_exact_ordered_complete_source():
    code='from matplotlib_venn import venn2\nvenn2(subsets=(10,20,3))'
    program=Program({'q1':{'count':[10,20,3]}})
    program.run(code,'Verified complete vector.')
    assert program.literal_bindings==[dict(table='q1',column='count',rows=3,verification='exact ordered complete-column match')]
    for tables in ({'q1':{'count':[10,3,20]}}, {'q1':{'count':[10,20,3,7]}},
                   {'q1':{'count':[10,20,3]},'q2':{'other':[10,20,3]}}):
        with pytest.raises(ValueError,match='not invented'):
            Program(tables).run(code,'Unreconciled vectors must not be plotted.')


def test_full_verified_groups_support_standard_plotting_and_aliases():
    code='from matplotlib_venn import venn2\na = {"MRC1","SFTPC"}\nb = {"MRC1","SFTPA1"}\nvenn2([a,b],set_labels=("A","B"))'
    groups={'A':['id1','id2'],'B':['id1','id3']}
    aliases={'A':['MRC1','SFTPC'],'B':['MRC1','SFTPA1']}
    program=Program({},groups,aliases)
    svg=program.run(code,'All requested source identities verified.')
    assert len(program.literal_bindings)==2 and 'A' in svg and 'B' in svg
    with pytest.raises(ValueError,match='not invented'):
        Program({},groups,{}).run(code,'Unresolved alias lists must not be plotted.')


def test_child_seeds_only_native_font_cache_and_keeps_private_configuration(tmp_path, monkeypatch):
    import json
    import matplotlib
    from pathlib import Path
    from altanalyze3.components.cellHarmony.webapp import plot_python as module
    (tmp_path/'fontlist-test.json').write_text('{}')
    (tmp_path/'matplotlibrc').write_text('text.usetex: True')
    monkeypatch.setattr(matplotlib,'get_cachedir',lambda:str(tmp_path))
    class Child:
        returncode=0
        def __init__(self,args,**kwargs):
            assert kwargs['start_new_session']
            directory=Path(kwargs['env']['MPLCONFIGDIR'])
            assert directory==Path(kwargs['env']['HOME'])
            assert (directory/'fontlist-test.json').read_text()=='{}'
            assert not (directory/'matplotlibrc').exists()
        def communicate(self,*args,**kwargs):return json.dumps({'svg':'<svg/>'}),''
    monkeypatch.setattr(module.subprocess,'Popen',Child)
    assert module.render('test',{},'test')=='<svg/>'


def test_render_timeout_stops_its_private_process_group(monkeypatch):
    from altanalyze3.components.cellHarmony.webapp import plot_python as module
    calls=[]
    class Child:
        pid=12345
        def __init__(self,*args,**kwargs):self.first=True
        def communicate(self,*args,**kwargs):
            if self.first:
                self.first=False
                raise module.subprocess.TimeoutExpired('plot worker',30)
            return '',''
    monkeypatch.setattr(module.subprocess,'Popen',Child)
    monkeypatch.setattr(module.os,'killpg',lambda pid,sig:calls.append((pid,sig)))
    with pytest.raises(ValueError,match='30-second'):
        module.render('test',{},'test')
    assert calls==[(12345,module.signal.SIGKILL)]


def test_literal_mappings_require_complete_ordered_data_and_empty_native_figures_are_rejected():
    code="venn2(subsets={'10':10,'01':20,'11':3})"
    program=Program({'q1':{'count':[10,20,3]}})
    assert '<svg' in program.run(code,'Complete source counts.')
    assert program.literal_bindings[0]['keys']==['10','01','11']
    with pytest.raises(ValueError,match='not invented'):
        Program({'q1':{'count':[10,20,3,9]}}).run(code,'Missing source counts cannot be hidden.')
    with pytest.raises(ValueError,match='zero area'):
        Program({'q1':{'count':[10,20,3]}}).run("venn2(subsets={'101':10,'011':20,'111':3})",'Invalid native keys must not become an empty figure.')
