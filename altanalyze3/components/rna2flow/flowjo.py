#!/usr/bin/env python3
"""Read a FlowJo workspace: the authors' real gating strategy, as they drew it.

A FlowJo workspace (.wsp) or template (.wspt) stores a gate tree in Gating-ML. This module
reads that tree and applies it to an event matrix, so the manual strategy and a derived
virtual strategy are scored on the SAME events with the SAME definitions.

Two coordinate spaces matter and the module never conflates them:

  raw      what the gate coordinates are written in, and what the gate is applied in.
           A FlowJo polygon on FSC-A x SSC-A carries values like 163414.4.
  display  the logicle values the viewer draws. Gate vertices are mapped through the SAME
           per-channel logicle the viewer used, so an outline lands on its own events.

Gate types read: PolygonGate, RectangleGate (2-D box and 1-D interval), EllipsoidGate.
A gate naming a channel the data lacks is reported unevaluable, never silently skipped.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np

NS = {"gating": "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
      "data-type": "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes",
      "transforms": "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"}
GT = "{%s}" % NS["gating"]
DT = "{%s}" % NS["data-type"]
TR = "{%s}" % NS["transforms"]

__all__ = ["Gate", "Population", "parse_workspace", "apply_tree", "population_table"]


@dataclass
class Gate:
    kind: str                      # polygon | rectangle | ellipsoid
    dims: list                     # channel names, in gate order
    vertices: list = field(default_factory=list)      # polygon
    bounds: list = field(default_factory=list)        # rectangle: (min, max) per dim
    mean: list = field(default_factory=list)          # ellipsoid
    cov: list = field(default_factory=list)
    dist2: float = 1.0

    def mask(self, cols):
        """cols: list of 1-D arrays, one per self.dims, in RAW space."""
        if self.kind == "polygon":
            return point_in_polygon(np.asarray(self.vertices, float), cols[0], cols[1])
        if self.kind == "rectangle":
            m = np.ones(len(cols[0]), bool)
            for c, (lo, hi) in zip(cols, self.bounds):
                m &= (c >= lo) & (c <= hi)
            return m
        if self.kind == "ellipsoid":
            d = np.column_stack(cols) - np.asarray(self.mean, float)
            inv = np.linalg.pinv(np.asarray(self.cov, float))
            return np.einsum("ij,jk,ik->i", d, inv, d) <= self.dist2
        raise ValueError("unsupported gate kind %r" % self.kind)


@dataclass
class Population:
    name: str
    path: str
    gate: Gate | None
    count_in_workspace: int | None
    children: list = field(default_factory=list)


def point_in_polygon(poly, X, Y):
    """Even-odd ray casting, vectorized over all events."""
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    inside = np.zeros(len(X), bool)
    j = len(poly) - 1
    for i in range(len(poly)):
        xi, yi = poly[i]; xj, yj = poly[j]
        inside ^= ((yi > Y) != (yj > Y)) & (X < (xj - xi) * (Y - yi) / (yj - yi + 1e-30) + xi)
        j = i
    return inside


def _dim_name(d):
    f = d.find(DT + "fcs-dimension")
    return f.get(DT + "name") if f is not None else d.get(DT + "name")


def _read_gate(gate_el):
    for el in gate_el:
        tag = el.tag.replace(GT, "")
        dims = [_dim_name(d) for d in el.findall(GT + "dimension")]
        if tag == "PolygonGate":
            verts = [[float(c.get(DT + "value"))
                      for c in v.findall(GT + "coordinate")]
                     for v in el.findall(GT + "vertex")]
            return Gate("polygon", dims, vertices=verts)
        if tag == "RectangleGate":
            b = []
            for d in el.findall(GT + "dimension"):
                lo = d.get(GT + "min"); hi = d.get(GT + "max")
                b.append((float(lo) if lo is not None else -np.inf,
                          float(hi) if hi is not None else np.inf))
            return Gate("rectangle", dims, bounds=b)
        if tag == "EllipsoidGate":
            mean = [float(c.get(DT + "value"))
                    for c in el.find(GT + "mean").findall(GT + "coordinate")]
            rows = [[float(e.get(DT + "value")) for e in r.findall(GT + "entry")]
                    for r in el.find(GT + "covarianceMatrix").findall(GT + "row")]
            dq = el.find(GT + "distanceSquare")
            return Gate("ellipsoid", dims, mean=mean, cov=rows,
                        dist2=float(dq.get(DT + "value")) if dq is not None else 1.0)
    return None


def _walk(node, parent_path, out):
    for pop in node.findall("Population"):
        name = pop.get("name")
        path = parent_path + "/" + name if parent_path else name
        g_el = pop.find("Gate")
        cnt = pop.get("count")
        p = Population(name, path, _read_gate(g_el) if g_el is not None else None,
                       int(cnt) if cnt not in (None, "-1") else None)
        out.append(p)
        sub = pop.find("Subpopulations")
        if sub is not None:
            _walk(sub, path, out)


def parse_transforms(path):
    """Per-channel FlowJo transform parameters, for reference and for display mapping."""
    root = ET.parse(path).getroot()
    out = {}
    for el in root.iter():
        tag = el.tag
        if not tag.startswith(TR):
            continue
        kind = tag.replace(TR, "")
        if kind not in ("biex", "linear", "log", "fasinh"):
            continue
        par = el.find(DT + "parameter")
        if par is None:
            continue
        out[par.get(DT + "name")] = dict(
            kind=kind, **{k.replace(TR, ""): v for k, v in el.attrib.items()})
    return out


def parse_workspace(path):
    """Return (populations, transforms). Populations are depth-first, parents before children."""
    root = ET.parse(path).getroot()
    pops = []
    for sub in root.iter("Subpopulations"):
        parent = [p for p in root.iter() if sub in list(p)]
        # only walk top-level Subpopulations (those whose parent is a GroupNode/SampleNode)
        if parent and parent[0].tag in ("GroupNode", "SampleNode", "Sample"):
            _walk(sub, "", pops)
    if not pops:                       # fall back: any Subpopulations block
        for sub in root.iter("Subpopulations"):
            _walk(sub, "", pops)
            break
    seen, uniq = set(), []
    for p in pops:
        if p.path not in seen:
            seen.add(p.path); uniq.append(p)
    return uniq, parse_transforms(path)


def apply_tree(pops, data, channel_lookup=None):
    """Apply the gate tree to raw event data.

    data: dict channel name -> 1-D raw array, or (matrix, names).
    Returns {path: {mask, n, parent, unevaluable}} with each child intersected with its parent.
    """
    if isinstance(data, tuple):
        M, names = data
        data = {n: np.asarray(M[:, i]) for i, n in enumerate(names)}
    get = channel_lookup or (lambda c: data.get(c))
    n = len(next(iter(data.values())))
    res = {}
    for p in pops:
        parent = p.path.rsplit("/", 1)[0] if "/" in p.path else None
        pm = res[parent]["mask"] if parent in res else np.ones(n, bool)
        if p.gate is None:
            res[p.path] = dict(mask=pm, n=int(pm.sum()), parent=parent, unevaluable=None)
            continue
        cols, missing = [], []
        for c in p.gate.dims:
            v = get(c)
            (cols.append(np.asarray(v, float)) if v is not None else missing.append(c))
        if missing:
            res[p.path] = dict(mask=np.zeros(n, bool), n=0, parent=parent,
                               unevaluable="missing channel(s): " + ", ".join(missing))
            continue
        m = pm & p.gate.mask(cols)
        res[p.path] = dict(mask=m, n=int(m.sum()), parent=parent, unevaluable=None)
    return res


def population_table(pops, applied, n_total):
    """One row per population: workspace count, recomputed count, and % of parent."""
    rows = []
    for p in pops:
        a = applied[p.path]
        par = a["parent"]
        denom = applied[par]["n"] if par in applied else n_total
        rows.append(dict(
            path=p.path, name=p.name, depth=p.path.count("/"),
            gate=p.gate.kind if p.gate else "none",
            dims="|".join(p.gate.dims) if p.gate else "",
            n=a["n"], pct_of_parent=round(100.0 * a["n"] / denom, 3) if denom else 0.0,
            pct_of_total=round(100.0 * a["n"] / n_total, 3),
            workspace_count=p.count_in_workspace, unevaluable=a["unevaluable"]))
    return rows


def build_lookup(data_channels, antibodies=None):
    """Resolve a gate's channel name against an FCS parameter list.

    A FlowJo gate names the bare detector ("BV421-A"). This FCS writes $PnN as
    "BV421-A :: CD135". Three keys are registered per column, each exact: the full $PnN,
    the detector prefix before " :: ", and the $PnS antibody. A key that two columns would
    claim is dropped rather than resolved by preference, so no gate is ever applied to a
    channel picked by guesswork.

    Returns (lookup dict name -> column index, ambiguous list).
    """
    reg, dup = {}, {}
    for i, name in enumerate(map(str, data_channels)):
        keys = {name, name.split(" :: ")[0].strip()}
        if antibodies is not None:
            keys.add(str(antibodies[i]))
        for k in keys:
            if k in reg and reg[k] != i:
                dup.setdefault(k, {reg[k]}).add(i)
            else:
                reg[k] = i
    for k in dup:
        reg.pop(k, None)
    return reg, sorted(dup)
