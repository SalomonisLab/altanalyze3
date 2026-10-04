"""Flask routes for the flow_viewer. All data access goes through data_api.FlowBundle.

Every data route takes an optional ?space=; omitted, it serves the default space. A route
that cannot resolve its space returns 404 rather than quietly serving the wrong array.
"""
from __future__ import annotations

import numpy as np
from functools import lru_cache
import json
import threading
from flask import Flask, Response, jsonify, render_template, request

from .data_api import FlowBundle

BIN = "application/octet-stream"


def create_app(bundle_root: str) -> Flask:
    app = Flask(__name__)
    bundle = FlowBundle(bundle_root)
    app.config['FLOW_BUNDLE'] = bundle
    gate_lock = threading.Lock()

    @app.errorhandler(ValueError)
    def invalid(e):
        return jsonify(error=str(e)), 400

    @app.errorhandler(KeyError)
    def unknown(e):
        return jsonify(error=str(e)), 404

    @lru_cache(maxsize=64)
    def fit_gate(key):
        from ...rna2flow.gate_search import optimize_population, optimize_signed_gate
        q = json.loads(key); space = q.get('space', 'flow')
        names = bundle.features(space)
        from ...rna2flow.crosswalk import normalize_marker
        def resolve_channel(channel):
            if channel in names:return channel
            hits=[c for c in names if not c.startswith('RNA:') and normalize_marker(c)==normalize_marker(channel)]
            if len(hits)==1:return hits[0]
            if len(hits)>1:raise ValueError('Ambiguous antibody alias; select a measured channel: '+str(channel))
            raise ValueError('Parent / required marker unavailable in this panel: '+str(channel))
        if q.get('parent_conditions'):
            q['parent_conditions']=[dict(c,channel=resolve_channel(c['channel'])) for c in q['parent_conditions']]
        if q.get('signs'):q['signs']={resolve_channel(c):op for c,op in q['signs'].items()}
        # Scatter, viability and reporter channels are not surface phenotype markers.
        allowed = q.get('channels') or [c for c in names if not
                   (c.startswith(('FSC', 'SSC', 'RNA:')) or c in ('DEAD', 'tdTomato-A'))]
        if len(allowed) != len(set(allowed)) or not set(allowed).issubset(names):
            raise ValueError('Unknown or duplicate channels')
        search_allowed=list(allowed)
        for condition in q.get('parent_conditions') or []:
            channel=condition.get('channel')
            if channel not in names:raise ValueError('Parent marker unavailable in this panel: '+str(channel))
            if channel not in allowed:allowed.append(channel)
        X = bundle.feature_matrix(space)[:, [names.index(c) for c in allowed]]
        levels = bundle.levels(q['label_set'], space)
        codes = bundle.labels(q['label_set'], space)
        labs = np.array([levels[c] if 0 <= c < len(levels) else 'unassigned' for c in codes])
        parent=q.get('parent');parent_mask=None
        if parent:
            if bundle.gateset(parent['gateset']).get('space','flow')!=space:
                raise ValueError('Parent gate must use the same flow events')
            parent_mask=bundle.gate_mask(parent['gateset'],parent['path'])
        extra=dict(parent_conditions=q.get('parent_conditions'),parent_mask=parent_mask)
        if q.get('signs'):
            return dict(optimize_signed_gate(X, labs, allowed, q['populations'], q['signs'],float(q.get('beta',1)),min_purity=float(q.get('min_purity',0)),min_recovery=float(q.get('min_recovery',0)),**extra),parent=parent)
        return dict(optimize_population(X, labs, allowed, q['populations'],
                                   beta=float(q.get('beta', 1)), seed=42,max_markers=int(q.get('max_markers',6)),search_features=search_allowed,include_hypergate=bool(q.get('include_hypergate',False)),min_purity=float(q.get('min_purity',0)),min_recovery=float(q.get('min_recovery',0)),**extra),parent=parent)

    def sp():
        return request.args.get("space") or None

    app.config["TEMPLATES_AUTO_RELOAD"] = True   # an edited template serves without a restart

    @app.route("/")
    def index():
        return render_template("interactive.html")

    @app.route("/isolation")
    def isolation():
        return render_template("isolation.html")

    @app.route("/favicon.ico")
    def favicon():
        return Response(status=204)

    @app.route("/api/manifest")
    def manifest():
        return jsonify(bundle.catalog())

    @app.route("/api/embedding/<path:name>")
    def embedding(name):
        try:
            return Response(np.asarray(bundle.embedding(name, sp())).tobytes(), mimetype=BIN)
        except KeyError as e:
            return jsonify({"error": str(e)}), 404

    @app.route("/api/labels/<path:name>")
    def labels(name):
        try:
            return Response(np.asarray(bundle.labels(name, sp())).tobytes(), mimetype=BIN)
        except KeyError as e:
            return jsonify({"error": str(e)}), 404

    @app.route("/api/channel/<path:name>")
    def channel(name):
        try:
            return Response(bundle.feature(name, sp()).tobytes(), mimetype=BIN)
        except KeyError as e:
            return jsonify({"error": str(e)}), 404

    @app.route("/api/gatemask/<gateset>/<path:path>")
    def gatemask(gateset, path):
        """The population's membership as one byte per event, for instant overlay."""
        try:
            return Response(bundle.gate_mask(gateset, path).astype(np.uint8).tobytes(),
                            mimetype=BIN)
        except KeyError as e:
            return jsonify({"error": str(e)}), 404

    @app.route("/api/gate", methods=["POST"])
    def gate():
        q = request.get_json(force=True)
        return jsonify(bundle.gate(q["polygon"], q.get("mode", q.get("space_mode", "channel")),
                                   q.get("label_sets", []), q.get("x"), q.get("y"),
                                   q.get("embedding"), space=q.get("space"),
                                   parent=q.get("parent")))

    @app.route("/api/concordance", methods=["POST"])
    def concordance():
        q = request.get_json(force=True)
        return jsonify(bundle.crosstab(q["a"], q["b"], q.get("polygon"), q.get("x"),
                                       q.get("y"), q.get("space"),
                                       q.get("mode", "channel"), q.get("embedding")))

    @app.route('/api/optimize', methods=['POST'])
    def optimize():
        q = request.get_json(force=True)
        if not isinstance(q.get('populations'), list) or not q['populations']:
            raise ValueError('Select at least one population')
        key = json.dumps({k: q[k] for k in ('space','label_set','populations','beta','channels','signs','parent_conditions','parent','max_markers','include_hypergate','min_purity','min_recovery') if k in q}, sort_keys=True)
        with gate_lock:
            result = fit_gate(key)
        result = dict(result, display_transform=bundle.space(q.get('space','flow')).get('display_transform',{}),
                      flow_provenance=bundle.space(q.get('space','flow')).get('provenance',{}),
                      measurement_provenance=bundle.space(q.get('space','flow')).get('provenance',{}),
                      space=q.get('space','flow'),label_set=q['label_set'])
        return jsonify(result)

    @app.route('/api/capture_templates')
    def capture_templates():
        from ...rna2flow.crosswalk import normalize_marker
        features=bundle.features(sp());templates=[]
        definitions=[('initial','Paper initial progenitor capture',[('Lin','<='),('CD117','>'),('CD34','>'),('CD115','<='),('Ly6C','<=')]),
                     ('multilin','Paper refined MultiLin capture',[('Lin','<='),('Sca-1','<='),('CD117','>'),('CD27','>')])]
        for key,name,signs in definitions:
            conditions=[];missing=[]
            for marker,op in signs:
                hits=[f for f in features if not f.startswith('RNA:') and normalize_marker(f)==normalize_marker(marker)]
                conditions.append(dict(channel=hits[0] if len(hits)==1 else marker,operator=op,threshold=None))
                if len(hits)!=1:missing.append(marker)
            templates.append(dict(id=key,name=name,conditions=conditions,missing=missing))
        return jsonify(templates=templates,source='nihms-2119165.pdf',thresholds='Require user control/workspace thresholds; capture strategies are optional and target-specific')

    @app.route('/api/populations/<path:name>')
    def populations(name):
        codes = bundle.labels(name, sp())
        counts = np.bincount(codes[codes >= 0], minlength=len(bundle.levels(name, sp())))
        return jsonify([dict(label=p, n=int(counts[i])) for i,p in enumerate(bundle.levels(name, sp()))])

    @app.route('/api/strategy', methods=['POST'])
    def strategy():
        from ...rna2flow.gate_search import membership
        q = request.get_json(force=True); space = q.get('space', 'flow')
        X = bundle.feature_matrix(space); names = bundle.features(space)
        sel = membership(X, names, q['conditions'])
        if q.get('parent'):
            if bundle.gateset(q['parent']['gateset']).get('space','flow')!=space:
                raise ValueError('Parent gate must use the same flow events')
            sel &= bundle.gate_mask(q['parent']['gateset'],q['parent']['path'])
        out = dict(n_gated=int(sel.sum()), n_total=len(sel), percent=100*float(sel.mean()), composition={})
        for ls in q.get('label_sets', []):
            levels=bundle.levels(ls, space); vals,counts=np.unique(bundle.labels(ls, space)[sel],return_counts=True)
            out['composition'][ls]=[dict(label=levels[int(v)] if 0<=v<len(levels) else 'unassigned',n=int(n),percent=100*int(n)/max(int(sel.sum()),1)) for v,n in zip(vals,counts)]
        return jsonify(out)

    @app.route('/api/export', methods=['POST'])
    def export():
        from ...rna2flow.gate_search import membership
        q = request.get_json(force=True); space=q.get('space','flow')
        sel = membership(bundle.feature_matrix(space), bundle.features(space), q['conditions'])
        if q.get('parent'):
            if bundle.gateset(q['parent']['gateset']).get('space','flow')!=space:
                raise ValueError('Parent gate must use the same flow events')
            sel &= bundle.gate_mask(q['parent']['gateset'],q['parent']['path'])
        import io, csv
        buf=io.StringIO();w=csv.writer(buf,delimiter='\t');w.writerow(['event_index_1based','label'])
        levels=bundle.levels(q['label_set'],space);codes=bundle.labels(q['label_set'],space)
        for i in np.flatnonzero(sel):w.writerow([int(i)+1,levels[int(codes[i])] if 0<=codes[i]<len(levels) else 'unassigned'])
        return Response(buf.getvalue(),mimetype='text/tab-separated-values',headers={'Content-Disposition':'attachment; filename=gated_events.tsv'})

    @app.route('/api/isolation/studies')
    def isolation_studies():
        from pathlib import Path
        root=Path(bundle.root)/'validation_20261003'/'ML_CLP_prospective'
        rows=[]
        for path in sorted(root.glob('selected__balanced_F1__*.json')):
            result=json.loads(path.read_text())
            rows.append(dict(file=path.name,space=result['space'],population=result['populations'][0],
                             source=result.get('source_space',result['space']),panel=result.get('panel','flow_panel'),
                             hypothesis=result['hypothesis'],test=result['best']['test']))
        return jsonify(studies=rows)

    @app.route('/api/isolation/studies/<filename>')
    def isolation_study(filename):
        from pathlib import Path
        if not filename.startswith('selected__balanced_F1__') or not filename.endswith('.json') or Path(filename).name!=filename:
            raise ValueError('Select an available study strategy')
        path=Path(bundle.root)/'validation_20261003'/'ML_CLP_prospective'/filename
        if not path.is_file():raise KeyError('Study strategy not found')
        return jsonify(json.loads(path.read_text()))

    @app.route('/api/isolation/audit', methods=['POST'])
    def isolation_audit():
        """Descriptive virtual experiment; edited thresholds never inherit held-out scores."""
        from ...rna2flow.gate_search import membership, metrics
        from matplotlib.path import Path as PolygonPath
        q=request.get_json(force=True);space=q.get('space','flow')
        X=bundle.feature_matrix(space);features=bundle.features(space)
        levels=bundle.levels(q['label_set'],space);codes=bundle.labels(q['label_set'],space)
        missing=set(q['populations'])-set(levels)
        if missing:raise ValueError('Unknown target populations: '+', '.join(sorted(missing)))
        target=np.isin(codes,[levels.index(p) for p in q['populations']]);current=np.ones(len(X),bool)
        baseline=float(target.mean());rows=[];rng=np.random.default_rng(42)
        def append_stage(name, before, conditions=None, polygon=None, x=None, y=None):
            nonlocal current
            cs=conditions or []
            if cs:current &= membership(X,features,cs)
            if polygon is not None:
                vertices=np.asarray(polygon,float)
                if vertices.ndim!=2 or vertices.shape[1]!=2 or len(vertices)<3 or not np.isfinite(vertices).all():raise ValueError('Polygon requires at least three finite x/y vertices')
                if x not in features or y not in features:raise ValueError('Select measured polygon axes')
                current &= PolygonPath(vertices).contains_points(X[:,[features.index(x),features.index(y)]])
            stats=metrics(current,target);nt=int((before&target).sum())
            stats.update(stage=name,n_entering=int(before.sum()),target_entering=nt,
                         target_lost=int((before&target&~current).sum()),
                         stage_recovery=float((current&target).sum()/nt) if nt else 0.,
                         enrichment=stats['precision']/baseline if baseline else 0.)
            x=x or q.get('plot_x') or (cs[0]['channel'] if cs else q.get('x',features[0]))
            y=y or q.get('plot_y') or next((c['channel'] for c in cs if c['channel']!=x),q.get('y',features[min(1,len(features)-1)]))
            if x==y:y=next((c for c in [q.get('x'),q.get('y')]+features if c in features and c!=x),x)
            if x not in features or y not in features:raise ValueError('Unknown plot channel')
            # A reproducible stratified display sample; counts always use every cell/event.
            display=[]
            for is_target,is_pass in [(True,True),(True,False),(False,True),(False,False)]:
                idx=np.flatnonzero(before&(target==is_target)&(current==is_pass))
                if len(idx)>800:idx=rng.choice(idx,800,replace=False)
                display.extend([[float(X[i,features.index(x)]),float(X[i,features.index(y)]),int(is_target),int(is_pass)] for i in idx])
            stats.update(x=x,y=y,points=display,conditions=cs,polygon=polygon,
                         composition=[dict(label=levels[int(v)] if 0<=v<len(levels) else 'unassigned',n=int(n)) for v,n in zip(*np.unique(codes[current],return_counts=True))])
            rows.append(stats)
        append_stage('Starting population',current.copy())
        if q.get('parent'):
            pa=q['parent']
            if bundle.gateset(pa['gateset']).get('space','flow')!=space:raise ValueError('Capture parent and analysis must use the same events')
            before=current.copy();current &= bundle.gate_mask(pa['gateset'],pa['path'])
            append_stage('Capture: '+pa['path'],before)
        for i,step in enumerate(q.get('steps',[])):
            append_stage(step.get('name','Gate '+str(i+1)),current.copy(),step.get('conditions'),step.get('polygon'),step.get('x'),step.get('y'))
        return jsonify(stages=rows,space=space,scale=bundle.space(space).get('feature_scale',bundle.space(space).get('display_transform',{})),
                       limitation='Descriptive agreement with the selected annotation on all events/cells. Manual edits and stage reordering do not constitute new held-out validation. CITE-only ADT thresholds require calibration on a flow panel; RNA features are not sort channels.',
                       sampling='Up to 800 points per target/pass category per plot. Counts use all events; points are cells/events, not biological replicates.')

    @app.route('/api/profiles', methods=['POST'])
    def profiles():
        q=request.get_json(force=True);populations=q['populations']
        rows=[]
        for space,annotation in [('marrow_reference_RNA','Population'),('cite_grimes','Mm-MarrowAtlas-L4'),('cite_chinese','Mm-MarrowAtlas-L4')]:
            if space not in bundle.man['spaces']:continue
            levels=bundle.levels(annotation,space);codes=bundle.labels(annotation,space)
            selected=np.isin(codes,[levels.index(p) for p in populations if p in levels])
            features=bundle.features(space);X=bundle.feature_matrix(space)
            genes=[f for f in features if f.startswith('RNA:')]
            for gene in genes:
                v=X[selected,features.index(gene)]
                rows.append(dict(space=space,gene=gene,n_cells=int(selected.sum()),
                                 median=float(np.median(v)) if len(v) else None,
                                 fraction_detected=float((v>0).mean()) if len(v) else None))
        return jsonify(rows=rows,limitation='Measured RNA in cells assigned these marrow states. Cell-level descriptive summaries; no donor-level inference. Reference and query preprocessing may differ. State assignment is not evidence of stem-cell function.')

    @app.route('/api/marker_profiles', methods=['POST'])
    def marker_profiles():
        from ...rna2flow.crosswalk import normalize_marker
        q=request.get_json(force=True);populations=q['populations'];rows=[]
        for space in ['cite_marrow_ADT195','cite_marrow_ADT112','cite_grimes','cite_chinese']:
            if space not in bundle.man['spaces']:continue
            annotation='Population' if space.startswith('cite_marrow') else 'Mm-MarrowAtlas-L4'
            levels=bundle.levels(annotation,space);codes=bundle.labels(annotation,space);features=bundle.features(space)
            X=bundle.feature_matrix(space)
            for population in populations:
                mask=codes==levels.index(population) if population in levels else np.zeros(len(codes),bool)
                for marker in q.get('markers',['Sca-1','CD11c','CD11b','CD27','CD117','CD25','CD4','CD8','RNA:Spi1','RNA:Irf4','RNA:Irf8']):
                    # The flow channel is called CD8; CITE panels explicitly measure CD8a.
                    # Keep CD8b distinct, as in the transfer crosswalk's explicit alias.
                    marker_key=normalize_marker('CD8a' if marker=='CD8' else marker)
                    hits=[f for f in features if f==marker or (not f.startswith('RNA:') and not marker.startswith('RNA:') and normalize_marker(f)==marker_key)]
                    if not hits:
                        rows.append(dict(space=space,population=population,marker=marker,n_cells=int(mask.sum()),available=False));continue
                    feature=hits[0];v=X[mask,features.index(feature)];allv=X[:,features.index(feature)]
                    quartiles=np.quantile(v,[.25,.5,.75]) if len(v) else [None]*3
                    rows.append(dict(space=space,population=population,marker=marker,feature=feature,n_cells=int(mask.sum()),available=bool(len(v)),
                       q25=float(quartiles[0]) if len(v) else None,median=float(quartiles[1]) if len(v) else None,q75=float(quartiles[2]) if len(v) else None,
                       median_percentile_within_panel=float(100*(allv<=quartiles[1]).mean()) if len(v) else None))
        return jsonify(rows=rows,limitation='Measured CITE-seq ADT and RNA, with panels kept separate. ADT units are TotalVI for marrow/Grimes and DSB for Chinese; RNA is log-normalized. Percentiles are within each complete panel and depend on capture composition. Low signal does not establish experimental negativity without controls; cells are not independent biological replicates.')

    @app.route('/api/validation')
    def validation():
        from pathlib import Path
        import pandas as pd
        root=Path(bundle.root)/'validation_20261002'
        p=root/'held_out_marker_summary.tsv'
        rows=pd.read_csv(p,sep='\t').replace({np.nan:None}).to_dict('records') if p.exists() else []
        p=root/'marrow_holdout'/'held_out_marker_summary.tsv'
        if p.exists():rows+=pd.read_csv(p,sep='\t').replace({np.nan:None}).to_dict('records')
        p=root/'cellharmony_holdout'/'held_out_marker_summary.tsv'
        if p.exists():rows+=pd.read_csv(p,sep='\t').replace({np.nan:None}).to_dict('records')
        p=root/'DN_holdout'/'held_out_marker_summary.tsv'
        if p.exists():rows+=pd.read_csv(p,sep='\t').replace({np.nan:None}).to_dict('records')
        return jsonify(rows=rows,limitation='Half of shared antibodies withheld from assignment in each fold. This tests cross-marker consistency, not functional identity or donor replication. Missing populations are retained in the accompanying population table.')

    return app
