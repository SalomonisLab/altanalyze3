"""Read side of the flow_viewer, in the shape of scalable_viewer.data_api.

Every method reads a precomputed bundle. Nothing opens an FCS, and nothing recomputes a
statistic that precompute wrote.

A bundle holds SPACES. A space is one set of objects with its own count: `flow` holds FCS
events, a `cite_*` space holds CITE-seq cells. Embeddings, label sets and features resolve
inside a space, so a flow embedding is never drawn against CITE labels.

Loading is lazy and memory-mapped: the manifest is read once, and an array is mapped on the
first request that needs it.
"""
from __future__ import annotations

import json
import os
import numpy as np

__all__ = ["FlowBundle"]


class FlowBundle:
    def __init__(self, root: str):
        self.root = os.path.abspath(root)
        with open(os.path.join(self.root, "manifest.json")) as fh:
            self.man = json.load(fh)
        self._cache = {}
        if "spaces" not in self.man:      # pre-spaces bundle: present it as one space
            self.man["spaces"] = {"flow": {
                **{k:self.man[k] for k in ('display_transform','provenance') if k in self.man},
                "n": self.man["n_events"], "features": self.man["channels"],
                "features_file": self.man["channels_file"],
                "embeddings": self.man["embeddings"], "labels": self.man["labels"]}}
        self.default_space = "flow" if "flow" in self.man["spaces"] else \
            next(iter(self.man["spaces"]))

    # ---- spaces ---------------------------------------------------------------------
    def space(self, name=None):
        s = self.man["spaces"].get(name or self.default_space)
        if s is None:
            raise KeyError("unknown space %r" % name)
        return s

    def n(self, space=None):
        return int(self.space(space)["n"])

    def features(self, space=None):
        return list(self.space(space)["features"])

    # ---- arrays, memory-mapped on first touch ---------------------------------------
    def _map(self, fname, dtype, shape=None):
        key = (fname, str(dtype))
        if key not in self._cache:
            a = np.memmap(os.path.join(self.root, "arrays", fname), dtype=dtype, mode="r")
            self._cache[key] = a.reshape(shape) if shape else a
        return self._cache[key]

    def feature_matrix(self, space=None):
        s = self.space(space)
        return self._map(s["features_file"], np.float32, (int(s["n"]), len(s["features"])))

    def feature(self, name, space=None):
        s = self.space(space)
        if name not in s["features"]:
            raise KeyError("unknown feature %r in space %r" % (name, space))
        return np.ascontiguousarray(self.feature_matrix(space)[:, s["features"].index(name)])

    channel = feature          # flow vocabulary

    def embedding(self, name, space=None):
        s = self.space(space)
        return self._map(s["embeddings"][name]["file"], np.float32, (int(s["n"]), 2))

    def labels(self, name, space=None):
        return self._map(self.space(space)["labels"][name]["file"], np.int16)

    def levels(self, name, space=None):
        return self.space(space)["labels"][name]["levels"]

    def mask(self, fname):
        """A boolean population mask written by precompute (uint8 on disk)."""
        return np.asarray(self._map(fname, np.uint8)).astype(bool)

    # ---- gate sets -------------------------------------------------------------------
    def gatesets(self):
        return self.man.get("gatesets", {})

    def gateset(self, name):
        g = self.gatesets().get(name)
        if g is None:
            raise KeyError("unknown gate set %r" % name)
        return g

    def gate_mask(self, gateset, path):
        for nd in self.gateset(gateset)["nodes"]:
            if nd["path"] == path:
                return self.mask(nd["mask_file"])
        raise KeyError("unknown population %r in gate set %r" % (path, gateset))

    # ---- catalog ---------------------------------------------------------------------
    def catalog(self):
        strip = lambda d: {k: v for k, v in d.items() if k != "file"}
        spaces = {}
        for name, s in self.man["spaces"].items():
            spaces[name] = {
                **{k: v for k, v in s.items() if k in ('description', 'modality', 'provenance', 'transfer_links', 'default_embedding', 'default_label', 'feature_scale', 'cell_ids_file', 'display_transform')},
                "n": int(s["n"]), "features": list(s["features"]),
                "embeddings": {k: strip(v) for k, v in s["embeddings"].items()},
                "labels": {k: strip(v) for k, v in s["labels"].items()}}
        gs = {}
        for name, g in self.gatesets().items():
            gs[name] = {"space": g.get("space", "flow"), "n_events": g.get("n_events"),
                        "kind": g.get("kind", "manual"),
                        "workspace": os.path.basename(g.get("workspace", "")
                                                      or g.get("source", "")),
                        "nodes": [{k: v for k, v in nd.items() if k != "mask_file"}
                                  for nd in g["nodes"]]}
        return {"spaces": spaces, "default_space": self.default_space, "gatesets": gs,
                "audit": self.man.get('audit', {}),
                "n_events": self.n(), "channels": self.features()}

    # ---- gating ----------------------------------------------------------------------
    @staticmethod
    def inside(poly, X, Y):
        """Even-odd ray casting, vectorized over every event."""
        p = np.asarray(poly, dtype=np.float64)
        inside = np.zeros(len(X), dtype=bool)
        j = len(p) - 1
        for i in range(len(p)):
            xi, yi = p[i]
            xj, yj = p[j]
            inside ^= ((yi > Y) != (yj > Y)) & (X < (xj - xi) * (Y - yi) / (yj - yi + 1e-30) + xi)
            j = i
        return inside

    def _xy(self, space, mode, x=None, y=None, embedding=None):
        if mode == "embedding":
            E = self.embedding(embedding, space)
            return np.asarray(E[:, 0], np.float64), np.asarray(E[:, 1], np.float64)
        return (self.feature(x, space).astype(np.float64),
                self.feature(y, space).astype(np.float64))

    def gate(self, poly, mode, label_sets, x=None, y=None, embedding=None, top=15,
             space=None, parent=None):
        """Gate inside an optional parent population, as FlowJo does."""
        space = space or self.default_space
        n = self.n(space)
        X, Y = self._xy(space, mode, x, y, embedding)
        sel = self.inside(poly, X, Y)
        base = np.ones(n, bool)
        if parent:
            base = self.gate_mask(parent["gateset"], parent["path"])
            sel &= base
        k = int(sel.sum())
        out = {"n_gated": k, "n_total": n, "n_parent": int(base.sum()),
               "percent": round(100.0 * k / max(n, 1), 3),
               "percent_of_parent": round(100.0 * k / max(int(base.sum()), 1), 3),
               "composition": {}}
        if not k:
            return out
        labs = self.space(space)["labels"]
        for ls in label_sets:
            if ls not in labs:
                continue
            lev = self.levels(ls, space)
            vals, cnt = np.unique(np.asarray(self.labels(ls, space))[sel], return_counts=True)
            order = np.argsort(-cnt)[:top]
            out["composition"][ls] = [
                {"label": lev[vals[i]] if 0 <= vals[i] < len(lev) else "NA",
                 "n": int(cnt[i]), "percent": round(100.0 * cnt[i] / k, 2)} for i in order]
        return out

    def crosstab(self, a, b, poly=None, x=None, y=None, space=None, mode="channel",
                 embedding=None):
        space = space or self.default_space
        la, lb = self.levels(a, space), self.levels(b, space)
        aa, bb = np.asarray(self.labels(a, space)), np.asarray(self.labels(b, space))
        sel = (aa >= 0) & (aa < len(la)) & (bb >= 0) & (bb < len(lb))
        if poly:
            X, Y = self._xy(space, mode, x, y, embedding)
            sel &= self.inside(poly, X, Y)
        M = np.zeros((len(la), len(lb)), dtype=np.int64)
        np.add.at(M, (np.asarray(self.labels(a, space))[sel].astype(int),
                      np.asarray(self.labels(b, space))[sel].astype(int)), 1)
        from sklearn.metrics import adjusted_rand_score
        return {"rows": la, "cols": lb, "counts": M.tolist(), 'n_compared': int(sel.sum()),
                'n_excluded': int((~sel).sum()),
                'ARI': float(adjusted_rand_score(aa[sel], bb[sel])) if sel.any() else None}
