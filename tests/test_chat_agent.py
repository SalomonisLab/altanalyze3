"""The general fallback retrieves real data; it does not need plot-specific intents."""
import json
import numpy as np
import pandas as pd
import pytest
from test_grn_web_release import uploaded, seed_comparisons, web
from altanalyze3.components.cellHarmony.webapp import chat_agent as agent


@pytest.fixture
def tools(uploaded):
    seed_comparisons(uploaded)
    return agent.JobTools(uploaded.app, uploaded.store.get_job(uploaded.job))


def test_schema_and_all_cell_descriptive_fractions(tools):
    schema = tools.schema()
    assert schema['cells'] == 192 and schema['features'] == 5
    assert 'condition' in {f['name'] for f in schema['obs']}
    assert len(schema['comparisons']) == 3 and 'X_umap' in schema['embeddings']
    result = tools.query(dict(source='obs', group_by=['Library', 'cell_type'], aggregate='count', fraction_within=['Library']))
    f = tools.frames[result['id']]
    assert result['matched_rows'] == 192 and len(f) == 16
    assert (f.value == 12).all() and (f.denominator == 24).all() and (f.fraction == .5).all()
    # Descriptive plots need no extra biological-replicate count or new test.
    one = tools.query(dict(source='obs', filters=[dict(column='Library', op='eq', value='donor0')],
                           group_by=['cell_type'], aggregate='count', fraction_within=[]))
    assert one['matched_rows'] == 24
    assert tools.frames[one['id']].fraction.tolist() == [.5, .5]


def test_sets_use_complete_tables_not_previews_and_keep_missing(tools):
    tools.frames['large_a'] = pd.DataFrame({'gene': [f'G{i}' for i in range(100)] + ['G0', None]})
    tools.frames['large_b'] = pd.DataFrame({'gene': [f'G{i}' for i in range(50, 150)]})
    a = tools.query({'source': 'large_a'})
    b = tools.query({'source': 'large_b'})
    assert a['preview_truncated'] and len(a['preview']) == 30
    result = tools.sets({'sets': [dict(query=a['id'], column='gene', label='A'), dict(query=b['id'], column='gene', label='B')]})
    assert result['data']['count'] == [50, 50, 50]
    assert result['unique_union'] == 150 and sum(result['data']['count']) == 150
    assert result['missing_ids'] == {'A': 1, 'B': 0} and result['duplicates'] == {'A': 1, 'B': 0}
    # References also resolve complete columns, not the preview.
    assert len(tools.resolve({'$ref': a['id'] + '/data/gene'})) == 102


def test_retrieves_saved_differentials_and_unmodified_expression(tools):
    q = tools.query(dict(source='differential:run-rna', filters=[dict(column='population', op='eq', value='AT1')]))
    assert set(tools.frames[q['id']].gene) == {'TF1', 'TARGET1'}
    original = np.asarray(tools.cache['adata'].X).copy()
    q = tools.query(dict(source='expression:rna', features=['TF1', 'TARGET3'],
                        group_by=['cell_type'], aggregate='mean', value='TF1'))
    expected = pd.DataFrame({'state': tools.cache['adata'].obs.cell_type, 'TF1': original[:, 0]}).groupby('state', sort=False).TF1.mean()
    np.testing.assert_allclose(tools.frames[q['id']].value, expected)
    np.testing.assert_array_equal(tools.cache['adata'].X, original)
    with pytest.raises(ValueError, match='Unresolved features'):
        tools.query(dict(source='expression:rna', features=['invented']))
    with pytest.raises(ValueError, match='Unknown source'):
        tools.query(dict(source='result:/etc/passwd'))


def venn_actions():
    # Model-authored layout: the server has no Venn renderer or Venn intent.
    return [
        dict(tool='query', args=dict(source='var:rna', filters=[dict(column='_feature', op='in', value=['TF1', 'TARGET1'])])),
        dict(tool='query', args=dict(source='var:rna', filters=[dict(column='_feature', op='in', value=['TF1', 'TARGET2'])])),
        dict(tool='sets', args={'sets': [dict(query='q1', column='_feature', label='A'), dict(query='q2', column='_feature', label='B')]}),
        dict(tool='finish', args=dict(answer='The two requested feature sets share TF1.',
            caption='Exact feature identities from the uploaded RNA feature panel. All selected features retained; schematic circles are not area-proportional. Descriptive counts, no significance test.',
            source_ids=['q3'], table_query='q3', figure=dict(data=[], layout=dict(title={'text': 'Feature-set overlap'},
                shapes=[dict(type='circle', x0=0, x1=2, y0=0, y1=2), dict(type='circle', x0=1, x1=3, y0=0, y1=2)],
                annotations=[dict(x=x, y=1, text={'$ref': 'q3/data/count', 'index': i}, showarrow=False)
                             for i, x in enumerate([.5, 2.5, 1.5])]))))]


def test_unanticipated_diagram_returns_grounded_figure_with_progress(uploaded, monkeypatch):
    actions = iter(venn_actions()); calls = []
    def model(messages, timeout):
        calls.append(messages)
        assert 'standard_result' not in json.loads(messages[1]['content']), 'A failed intent must not bias the fresh fallback question.'
        return next(actions)
    monkeypatch.setattr(agent, '_complete', model)
    monkeypatch.setattr(web, '_chat_read_question', lambda *args: {'intent': 'clarify'})
    r = uploaded.client.post(f'/api/jobs/{uploaded.job}/chat', json={'question': 'Make a Venn diagram of A (TF1,TARGET1) and B (TF1,TARGET2)', 'stream': True})
    assert r.status_code == 200 and 'application/x-ndjson' in r.headers['content-type']
    events = [json.loads(line) for line in r.text.splitlines()]
    assert sum(e['type'] == 'progress' for e in events) >= 4
    result = events[-1]['result']
    assert result['intent'] == 'data_assistant' and result['plot']['kind'] == 'agent_figure'
    assert [a['text'] for a in result['plot']['figure']['layout']['annotations'][:3]] == [1, 1, 1]
    assert result['table']['rows'] == [['A only', 1], ['B only', 1], ['A & B', 1]]
    assert len(result['provenance']['tools']) == 3 and len(calls) == 4


def test_existing_answers_do_not_call_fallback(uploaded, monkeypatch):
    monkeypatch.setattr(agent, '_complete', lambda *a, **k: pytest.fail('A working standard query must stay unchanged'))
    monkeypatch.setattr(web, '_chat_read_question', lambda *args: {'intent': 'expression', 'genes': ['TF1']})
    r = uploaded.client.post(f'/api/jobs/{uploaded.job}/chat', json={'question': 'Where is TF1 expressed?', 'stream': True})
    assert r.status_code == 200 and 'application/json' in r.headers['content-type']
    assert r.json()['intent'] == 'expression' and r.json()['table']['rows']


def test_failures_are_explicit_and_model_can_repair_tool_errors(uploaded, monkeypatch):
    actions = iter([dict(tool='query', args={'source': 'obs', 'columns': ['invented']}),
                    dict(tool='query', args={'source': 'obs', 'group_by': ['cell_type'], 'aggregate': 'count'}),
                    dict(tool='finish', args={'answer': 'Counts for every stored cell type.', 'table_query': 'q1'})])
    calls = []
    def model(messages, timeout):
        calls.append(json.loads(json.dumps(messages)))
        return next(actions)
    monkeypatch.setattr(agent, '_complete', model)
    events = list(agent.answer_events(uploaded.app, uploaded.store.get_job(uploaded.job), 'Count cells', {'status': 'not_covered'}))
    assert 'invented' in calls[1][-1]['content'] and 'error' in calls[1][-1]['content']
    assert events[-1]['result']['table']['rows'] == [['AT1', 96], ['AT2', 96]]
    monkeypatch.setattr(agent, '_complete', lambda *a: (_ for _ in ()).throw(TimeoutError('model unavailable')))
    result = list(agent.answer_events(uploaded.app, uploaded.store.get_job(uploaded.job), 'question', {'answer': 'No standard result.'}))[-1]['result']
    assert result['fallback']['status'] == 'unavailable' and 'model unavailable' in result['answer']


def test_custom_plot_rejects_invented_values_and_remote_data(tools):
    q = tools.query(dict(source='obs', group_by=['cell_type'], aggregate='count'))
    args = dict(answer='Counts.', caption='All retained cells, descriptive only.', source_ids=[q['id']],
                figure={'data': [{'type': 'bar', 'x': [1, 2], 'y': [500, 999]}]})
    with pytest.raises(ValueError, match='must reference'):
        tools.finish(args)
    args['figure'] = dict(data=[dict(type='bar', x={'$ref': q['id'] + '/data/cell_type'}, y={'$ref': q['id'] + '/data/value'})])
    result = tools.finish(args)
    assert result['plot']['figure']['data'][0]['y'] == [96, 96]
    args['figure']['layout'] = {'images': [{'source': 'https://example.com/track'}]}
    with pytest.raises(ValueError, match='local retrieved'):
        tools.finish(args)


def test_requested_graphic_not_silently_replaced_by_marker_lookup():
    result = {'intent': 'markers', 'table': {'rows': [['TF1']]}, 'plot': {'kind': 'dotplot'}}
    assert agent.needs_fallback(result, 'Make a Venn diagram of marker genes')
    assert agent.needs_fallback(result, 'Make a histogram of gene values')
    assert not agent.needs_fallback(result, 'What are the marker genes?')
    assert not agent.needs_fallback(result, 'Show a dotplot of marker genes')


def test_established_view_uses_native_plot_payload(tools):
    view = tools.view(dict(kind='dotplot', genes=['TF1'], modality='rna'))
    result = tools.finish(dict(answer='Stored TF1 expression.', view_id=view['id']))
    assert result['plot']['view'] == 'dotplot'
    assert result['plot']['payload']['genes'] == ['TF1']
    freq=tools.view(dict(kind='frequency',modality='rna',group_by='Library'))
    assert freq['sample_field']=='Library' and freq['counted_cells']==192
    result=tools.finish(dict(answer='Descriptive frequency only; no significance test.',view_id=freq['id']))
    assert result['plot']['view']=='frequency' and result['plot']['payload']['query']==[]


def test_verified_alias_mapping_preserves_native_features(tools):
    from altanalyze3.components.cellHarmony.webapp.feature_lookup import FeatureLookup
    cache = tools.cache
    original = cache['var_names'].copy()
    cache['adata'].var['gene_symbols'] = ['FACTOR_A', 'FACTOR_B', 'GENE_A', 'GENE_B', 'GENE_C']
    cache['feature_lookup'] = FeatureLookup(cache['var_names'], cache['adata'].var)
    q = tools.query({'source': 'var:rna', 'features': ['FACTOR_A', 'GENE_A']})
    assert tools.frames[q['id']]._feature.tolist() == ['TF1', 'TARGET1']
    assert tools.schema('FACTOR_A')['mentioned_features'] == {'FACTOR_A': 'TF1'}
    np.testing.assert_array_equal(cache['var_names'], original)
    with pytest.raises(ValueError, match='Unresolved features'):
        tools.query({'source': 'var:rna', 'features': ['FACTOR_A', 'MISSING']})


def test_model_json_mode_receives_actual_schema_bindings_and_tolerates_text_newlines(tools, monkeypatch):
    import io
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost:9/v1/chat/completions')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','test')
    captured=[]
    def respond(request,timeout):
        captured.append(json.loads(request.data))
        return io.BytesIO(json.dumps({'choices':[{'message':{'content':'{"tool":"finish","args":{"answer":"Test"}}'}}]}).encode())
    monkeypatch.setattr(agent.urllib.request,'urlopen',respond)
    messages=[dict(role='system',content=agent.INSTRUCTIONS),dict(role='user',content=json.dumps({'dataset':tools.schema(),'question':'Make a plot'}))]
    agent._complete(messages,2)
    assert captured[0]['response_format']=={'type':'json_object'}
    schema=json.loads(captured[0]['messages'][1]['content'])['dataset']
    assert 'var:rna' in schema['available_sources'] and 'var:MODALITY' not in schema['available_sources']
    assert captured[0]['scalable_timeout_seconds']==1
    messages.append(dict(role='user',content='Tool result: {"id":"q1","columns":["gene"]}'))
    agent._complete(messages,2)
    assert 'q1_gene' in captured[-1]['messages'][1]['content']
    messages.extend([dict(role='assistant',content=json.dumps({'tool':'finish','args':{'python':'broken_plot()'}})),
        dict(role='user',content='Tool result: {"error":"Unknown plotting name"}')])
    agent._complete(messages,2)
    assert 'Previous draft to repair:' in captured[-1]['messages'][1]['content']
    assert 'broken_plot()' in captured[-1]['messages'][1]['content']
    native=messages[:2]+[dict(role='user',content='Tool result: {"id":"v1","ready":true,"kind":"frequency"}')]
    agent._complete(native,2)
    assert '"view_id":"v1"' in captured[-1]['messages'][0]['content']
    assert 'Do not generate Python' in captured[-1]['messages'][0]['content']
    rejected=messages[:2]+[dict(role='assistant',content='{"tool":"view","args":{"kind":"wrong"}}'),
        dict(role='user',content='Tool result: {"error":"Use marker_heatmap with only kind and modality"}')]
    agent._complete(rejected,2)
    repair=captured[-1]['messages']
    assert len(repair)==3 and 'Replace it using the explicit correction' in repair[0]['content']
    assert 'only kind and modality' in repair[-1]['content']
    assert not any(m['role']=='assistant' for m in repair)
    def multiline(request,timeout):
        content='{"tool":"finish","args":{"answer":"Two lines\nfrom local model"}}'
        return io.BytesIO(json.dumps({'choices':[{'message':{'content':content}}]}).encode())
    monkeypatch.setattr(agent.urllib.request,'urlopen',multiline)
    assert agent._complete(messages,2)['args']['answer']=='Two lines\nfrom local model'


def test_direct_explicit_sets_resolve_every_alias_and_plan_completes_in_one_call(tools, uploaded, monkeypatch):
    direct = {'sets': [dict(source='var:rna', members=['TF1', 'TF2'], label='A'),
                       dict(source='var:rna', members=['TF2', 'TARGET1'], label='B')]}
    result = tools.sets(direct)
    assert result['data']['count'] == [1, 1, 1] and result['unique_union'] == 3
    assert result['relationships'][0]['relationship']=='partially overlapping; neither contains the other'
    with pytest.raises(ValueError, match='Unresolved features'):
        tools.sets({'sets': [dict(source='var:rna', members=['MISSING'], label='A'), direct['sets'][1]]})
    calls = []
    def model(*args):
        calls.append(1)
        return dict(tool='plan', args={'actions': [dict(tool='sets', args=direct), dict(tool='finish', args={
            'answer': 'Exact overlap of the requested sets.', 'table_query': 'q1'})]})
    monkeypatch.setattr(agent, '_complete', model)
    result = list(agent.answer_events(uploaded.app, uploaded.store.get_job(uploaded.job), 'Compare lists', {}))[-1]['result']
    assert result['table']['rows'] == [['A only', 1], ['B only', 1], ['A & B', 1]] and len(calls) == 1


def test_query_cannot_filter_out_explicitly_requested_features(tools):
    with pytest.raises(ValueError, match='Filters removed explicitly requested features'):
        tools.query(dict(source='var:rna', features=['TF1', 'TF2'], filters=[dict(column='_feature', op='eq', value='TF1')]))
    empty=tools.sets({'sets':[dict(source='var:rna',members=[],label='empty'),dict(source='var:rna',members=['TF1'],label='one')]})
    assert empty['sizes']=={'empty':0,'one':1} and empty['unique_union']==1


def test_generic_plot_program_uses_verified_sets_and_library_geometry(tools):
    q = tools.sets({'sets': [dict(source='var:rna', members=['TF1', 'TF2'], label='A'),
                           dict(source='var:rna', members=['TF2', 'TARGET1'], label='B')]})
    tools.question='Draw a diagram of the sets.'
    args=dict(answer='Exact set overlap.',caption='Unique features from the complete requested lists; descriptive counts.',
        source_ids=[q['id']],python="from matplotlib_venn import venn2\nvenn2((groups['A'], groups['B']), set_labels=('A','B'))")
    result=tools.finish(args)
    assert result['plot']['kind']=='agent_svg'
    assert result['provenance']['plot_versions']['matplotlib-venn']
    assert result['table']['rows']==[['A only',1],['B only',1],['A & B',1]]
    with pytest.raises(ValueError,match='actual figure'):
        tools.finish(dict(answer='Here is a diagram',caption='descriptive'))
    with pytest.raises(ValueError,match='freehand SVG'):
        tools.finish(dict(args,python=None,svg='<svg/>'))
    with pytest.raises(ValueError,match='custom Plotly'):
        tools.view(dict(kind='dotplot',genes=['TF1']))


def test_explicit_lists_infer_only_unambiguous_feature_source(tools, monkeypatch):
    with pytest.raises(ValueError,match='Choose the feature source'):
        tools.sets({'sets':[dict(label='A',members=['TF1']),dict(label='B',members=['TF2'])]})
    monkeypatch.setattr(tools,'schema',lambda:{'available_sources':['var:rna']})
    result=tools.sets({'sets':[dict(label='A',members=['TF1','TF2']),dict(label='B',members=['TF2','TARGET1'])]})
    assert result['sizes']=={'A':2,'B':2}
    assert tools.audit[-1]['arguments']['sets'][0]['source']=='var:rna'


def test_prepared_native_view_returns_without_an_extra_model_call(uploaded,monkeypatch):
    calls=[]
    def model(*args):
        calls.append(1)
        assert len(calls)==1
        return dict(tool='view',args=dict(kind='frequency',modality='rna',group_by='Library'))
    monkeypatch.setattr(agent,'_complete',model)
    result=list(agent.answer_events(uploaded.app,uploaded.store.get_job(uploaded.job),'Show cell frequencies by sample',{}))[-1]['result']
    assert result['status']=='answered' and result['plot']['view']=='frequency'
    assert '192 retained selected cells' in result['caption']
    assert len(calls)==1


@pytest.mark.parametrize('question,kind', [
    ('Show a UMAP', 'umap_states'), ('Show a violin plot of TF1', 'violin'),
    ('Show a comb plot', 'combplot'), ('Show a heatmap', 'marker_heatmap'),
    ('Show a volcano plot', 'volcano'), ('Show a pathway', 'pathway'),
    ('Show a network', 'network'), ('Show a dot plot', 'dotplot')])
def test_native_plot_requests_cannot_generate_replacement_code(tools, question, kind):
    tools.question=question
    with pytest.raises(ValueError,match='established renderer'):
        tools.finish(dict(answer='Generated.',python='plt.plot([1,2])'))
    with pytest.raises(ValueError,match='established renderer'):
        tools.finish(dict(answer='Generated.',figure={'data':[]}))
    assert kind in tools.schema()['native_views']
    assert not agent.needs_fallback({'plot':{'kind':'existing_view','view':kind}},question)
    assert agent.needs_fallback({'plot':{'kind':'barchart'}},question)


def test_saved_comparison_views_reuse_builders_and_exact_run(tools):
    # The active run is GRN; choosing RNA must not inherit its modality/statistics.
    assert tools.meta['differential']['run_id']=='run-grn'
    from altanalyze3.components.cellHarmony.webapp.grn_data import comparison_metadata
    snapshot=comparison_metadata(tools.meta,'run-rna')
    for kind,builder in [('volcano',web._build_differential_volcano_payload),
                         ('differential_heatmap',web._build_differential_heatmap_payload)]:
        view=tools.view(dict(kind=kind,run_id='run-rna',cell_state='AT1'))
        result=tools.finish(dict(answer='The selected saved RNA comparison.',view_id=view['id']))
        plot=result['plot']
        assert plot['payload']==builder(tools.app,snapshot,'AT1')
        assert plot['modality']=='rna' and plot['contrast']=='run-rna'
        assert plot['state_layer']=='cell_type'
    assert tools.meta['differential']['run_id']=='run-grn'
    with pytest.raises(ValueError,match='exact completed run_id'):
        tools.view(dict(kind='volcano',cell_state='AT1'))
    with pytest.raises(ValueError,match='modality does not match'):
        tools.view(dict(kind='volcano',run_id='run-rna',cell_state='AT1',modality='grn'))
    with pytest.raises(ValueError,match='exact cell_state'):
        tools.view(dict(kind='volcano',run_id='run-rna',cell_state='invented'))


def test_integrated_and_marker_heatmap_delegate_to_existing_components(tools,monkeypatch):
    checked=[]
    monkeypatch.setattr(web,'_get_marker_heatmap_cache_entry',lambda app,meta,**kw:checked.append(kw) or {})
    view=tools.view(dict(kind='marker_heatmap'))
    assert tools.views[view['id']]['view']=='marker_heatmap'
    assert checked==[{'modality':'rna','compact':True}]
    with pytest.raises(ValueError,match='Remove run_id, cell_state and genes'):
        tools.view(dict(kind='marker_heatmap',run_id='analysis/assignments'))
    for kind in ('network','pathway'):
        args=dict(kind=kind,run_id='run-rna',cell_state='AT1',modality='lipid' if kind=='pathway' else 'rna')
        view=tools.view(args)
        plot=tools.finish(dict(answer='Existing integrated component.',view_id=view['id']))['plot']
        assert plot['contrast']=='run-rna' and plot['source']=='differential' and plot['cell_state']=='AT1'
        assert plot['payload'] is None  # Its cached native endpoint is fetched on access.
    with pytest.raises(ValueError,match='lipid or metabolite'):
        tools.view(dict(kind='pathway',run_id='run-rna',cell_state='AT1',modality='rna'))


def test_regeneration_attempt_repairs_to_native_volcano(uploaded,monkeypatch):
    seed_comparisons(uploaded)
    actions=iter([
        dict(tool='finish',args=dict(answer='Generated volcano.',python='plt.scatter([1],[2])')),
        dict(tool='view',args=dict(kind='volcano',run_id='run-rna',cell_state='AT1'))])
    monkeypatch.setattr(agent,'_complete',lambda *args:next(actions))
    from altanalyze3.components.cellHarmony.webapp import plot_python
    monkeypatch.setattr(plot_python,'render',lambda *a,**k:pytest.fail('An existing plot must never execute generated code.'))
    result=list(agent.answer_events(uploaded.app,uploaded.store.get_job(uploaded.job),'Show a volcano plot for AT1',{}))[-1]['result']
    assert result['plot']['view']=='volcano' and len(result['plot']['payload']['points'])==2
    assert 'established renderer' in result['provenance']['assistant_actions'][0]['errors'][0]
