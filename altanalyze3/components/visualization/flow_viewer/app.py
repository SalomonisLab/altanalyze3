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
        # Scatter, viability and reporter channels are not surface phenotype markers.
        allowed = q.get('channels') or [c for c in names if not
                   (c.startswith(('FSC', 'SSC')) or c in ('DEAD', 'tdTomato-A'))]
        if len(allowed) != len(set(allowed)) or not set(allowed).issubset(names):
            raise ValueError('Unknown or duplicate channels')
        X = bundle.feature_matrix(space)[:, [names.index(c) for c in allowed]]
        levels = bundle.levels(q['label_set'], space)
        codes = bundle.labels(q['label_set'], space)
        labs = np.array([levels[c] if 0 <= c < len(levels) else 'unassigned' for c in codes])
        if q.get('signs'):
            return optimize_signed_gate(X, labs, allowed, q['populations'], q['signs'],float(q.get('beta',1)))
        return optimize_population(X, labs, allowed, q['populations'],
                                   beta=float(q.get('beta', 1)), seed=42)

    def sp():
        return request.args.get("space") or None

    app.config["TEMPLATES_AUTO_RELOAD"] = True   # an edited template serves without a restart

    @app.route("/")
    def index():
        return render_template("interactive.html")

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
        key = json.dumps({k: q[k] for k in ('space','label_set','populations','beta','channels','signs') if k in q}, sort_keys=True)
        with gate_lock:
            result = fit_gate(key)
        result = dict(result, display_transform=bundle.space(q.get('space','flow')).get('display_transform',{}),
                      flow_provenance=bundle.space(q.get('space','flow')).get('provenance',{}))
        return jsonify(result)

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
        import io, csv
        buf=io.StringIO();w=csv.writer(buf,delimiter='\t');w.writerow(['event_index_1based','label'])
        levels=bundle.levels(q['label_set'],space);codes=bundle.labels(q['label_set'],space)
        for i in np.flatnonzero(sel):w.writerow([int(i)+1,levels[int(codes[i])] if 0<=codes[i]<len(levels) else 'unassigned'])
        return Response(buf.getvalue(),mimetype='text/tab-separated-values',headers={'Content-Disposition':'attachment; filename=gated_events.tsv'})

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
