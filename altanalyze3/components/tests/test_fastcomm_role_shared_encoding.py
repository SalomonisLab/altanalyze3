"""Visible topology colours and shared context never alter native differential edges."""
from copy import deepcopy
import importlib
from types import SimpleNamespace
import pandas as pd
import pytest

W = importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')


def test_receiver_sender_balanced_and_unconnected_colours_preserve_source():
    elements = [{'data': {'id': n, 'color': '#old'}} for n in ['S', 'B', 'R', 'U']]
    elements += [{'data': {'id': 'e1', 'source': 'S', 'target': 'B', 'delta_score': -2.0}},
                 {'data': {'id': 'e2', 'source': 'B', 'target': 'R', 'delta_score': 3.0}}]
    before = deepcopy(elements); actual = W._fastcomm_role_encoding(elements)
    assert elements == before
    nodes = {e['data']['id']: e['data'] for e in actual if 'source' not in e['data']}
    assert nodes['S']['sender_fraction'] == 1 and nodes['S']['color'] == '#018571'
    assert nodes['R']['sender_fraction'] == 0 and nodes['R']['color'] == '#5e3c99'
    assert nodes['B']['sender_fraction'] == .5 and nodes['B']['color'] == '#f7f7f7'
    assert nodes['U']['sender_fraction'] is None and nodes['U']['color'] == '#bdbdbd'
    assert actual[-2:] == before[-2:]


def test_self_loop_counts_once_in_each_direction():
    e = [{'data': {'id': 'A'}}, {'data': {'id': 'e', 'source': 'A', 'target': 'A'}}]
    a = W._fastcomm_role_encoding(e)[0]['data']
    assert a['incoming_arrows'] == a['outgoing_arrows'] == 1 and a['sender_fraction'] == .5


def test_shared_context_never_recolours_or_rescales_highlighted_native_changes(monkeypatch):
    groups = [{'label': 'MDS', 'splits': ['M1','M2','M3']}, {'label': 'aged', 'splits': ['A1','A2','A3']}]
    layout = dict(comparison={'groups': groups}, score_column='lr_expression_score', states=['S','R','T'],
                  diff_pairs=[('S','R')], positions={n: {'x': i, 'y': 0} for i,n in enumerate(['S','R','T'])},
                  diff_pairs_per_direction=(1,0,1,0), comparison_modes=['native unchanged'], n_significant=1, n_tested=12,
                  shared_context_pairs=[dict(sender='S',receiver='R',case_total_score=4.,control_total_score=2.),
                                        dict(sender='S',receiver='T',case_total_score=2.,control_total_score=3.)])
    detail = pd.DataFrame([dict(sender_state='S',receiver_state='R',delta_score=.4)])
    original_edges = [{'data': dict(id='ccd0',source='state::S',target='state::R',delta_score=.4,
                                    edge_color='#dc2626',weight=7.,edge_opacity=.6)}]
    nodes = [{'data': dict(id='state::'+n,label=n), 'position': layout['positions'][n]} for n in layout['states']]
    monkeypatch.setattr(W, '_fastcomm_shared_layout', lambda *args: deepcopy(layout))
    monkeypatch.setattr(W, '_get_differential_detail_table', lambda *args: detail.copy())
    monkeypatch.setattr(W.pipeline_mod, 'cell_communication_significant', lambda d: pd.Series([True],index=d.index))
    monkeypatch.setattr(W, '_cell_communication_diff_elements', lambda *args, **kwargs: (deepcopy(nodes+original_edges),detail.copy()))
    p = W._fastcomm_shared_plot_payload(SimpleNamespace(),{},'differential_network',pd.DataFrame(),None,None,5)
    edges = [e for e in p['elements'] if 'source' in e['data']]
    assert edges[0] == original_edges[0]
    assert len(edges) == 2 and edges[1]['data']['context_only']
    assert edges[1]['data']['edge_color'] == '#d8ac00' and edges[1]['data']['weight'] == 1.5
    assert 'delta_score' not in edges[1]['data']  # no fabricated zero or "unchanged" statistic
    assert p['summary']['n_tested'] == 12 and p['summary']['n_significant'] == 1
    assert p['summary']['n_displayed_shared_context_pairs'] == 1
    assert 'does not establish equal activity' in p['caption']
    assert p['communication_legend']['node_stops'][0]['label'] == 'Receiver'
