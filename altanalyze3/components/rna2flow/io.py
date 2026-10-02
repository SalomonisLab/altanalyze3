"""Readers for flow-cytometry inputs. Formats vary between facilities, so each reader is
independent and every one returns the same thing: (events x channels float32 matrix,
channel labels, per-channel metadata).

read_fcs            FCS 3.0/3.1 binary, pure numpy, no external package
read_flowjo_rds     the inputCSV slot of a FlowJo FJObj, via Rscript
read_matrix_csv     a plain events x channels CSV

A FlowJo channel is named "<detector> :: <antibody>", e.g. "BUV395-A :: CD11b". `antibody_name`
returns the antibody, which is what a CITE-seq panel can be matched on.
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

RSCRIPT = "/Library/Frameworks/R.framework/Resources/bin/Rscript"

__all__ = ["read_fcs", "read_flowjo_rds", "read_matrix_csv", "antibody_name", "FlowData"]


class FlowData:
    """events x channels, with labels. Deliberately thin."""

    def __init__(self, X: np.ndarray, channels, source: str, meta=None):
        self.X = np.asarray(X, dtype=np.float32)
        self.channels = list(channels)
        self.source = source
        self.meta = dict(meta or {})
        if self.X.shape[1] != len(self.channels):
            raise ValueError("channel count %d does not match matrix width %d"
                             % (len(self.channels), self.X.shape[1]))

    @property
    def antibodies(self):
        return [antibody_name(c) for c in self.channels]

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.X, columns=self.channels)

    def __repr__(self):
        return "FlowData(%d events x %d channels, %s)" % (self.X.shape[0], self.X.shape[1],
                                                          os.path.basename(self.source))


def antibody_name(channel: str) -> str:
    """'BUV395-A :: CD11b' -> 'CD11b'. A channel with no '::' returns itself, trimmed."""
    s = str(channel)
    if "::" in s:
        s = s.split("::", 1)[1]
    return s.strip()


def _parse_fcs_text(raw: bytes):
    """TEXT segment: <delim>KEY<delim>VALUE<delim>... The delimiter is the first byte.
    A doubled delimiter is a literal delimiter inside a value (FCS 3.1 section 3.2)."""
    delim = raw[:1]
    body = raw[1:]
    parts = body.split(delim)
    # Re-join the pairs produced by an escaped (doubled) delimiter.
    merged, i = [], 0
    while i < len(parts):
        cur = parts[i]
        while i + 1 < len(parts) and parts[i + 1] == b"":
            cur = cur + delim + (parts[i + 2] if i + 2 < len(parts) else b"")
            i += 2
        merged.append(cur)
        i += 1
    merged = [p for p in merged if p != b""]
    return {merged[i].decode("latin-1").strip(): merged[i + 1].decode("latin-1").strip()
            for i in range(0, len(merged) - 1, 2)}


def read_fcs(path: str) -> FlowData:
    """Minimal FCS 3.0/3.1 reader. Raises on anything it cannot honour exactly."""
    with open(path, "rb") as fh:
        header = fh.read(58)
        if not header.startswith(b"FCS"):
            raise ValueError("%s is not an FCS file (magic %r)" % (path, header[:6]))
        t0, t1, d0, d1 = (int(header[i:i + 8]) for i in (10, 18, 26, 34))
        fh.seek(t0)
        text = _parse_fcs_text(fh.read(t1 - t0 + 1))

        # DATA offsets of 0 in the HEADER mean "see $BEGINDATA/$ENDDATA" (large files).
        if d0 == 0 or d1 == 0:
            d0, d1 = int(text["$BEGINDATA"]), int(text["$ENDDATA"])

        n_events = int(text["$TOT"])
        n_par = int(text["$PAR"])
        mode = text.get("$MODE", "L")
        if mode != "L":
            raise ValueError("only $MODE L (list) is supported, got %r" % mode)
        dtype_code = text["$DATATYPE"].upper()
        byteord = text["$BYTEORD"]
        little = byteord.startswith("1,2")

        widths = [int(text["$P%dB" % (i + 1)]) for i in range(n_par)]
        if len(set(widths)) != 1:
            raise ValueError("mixed channel widths %s are not supported" % sorted(set(widths)))
        bits = widths[0]
        if dtype_code == "F":
            np_dt = np.dtype(np.float32 if bits == 32 else np.float64)
        elif dtype_code == "D":
            np_dt = np.dtype(np.float64)
        elif dtype_code == "I":
            np_dt = np.dtype({8: np.uint8, 16: np.uint16, 32: np.uint32, 64: np.uint64}[bits])
        else:
            raise ValueError("$DATATYPE %r is not supported" % dtype_code)
        np_dt = np_dt.newbyteorder("<" if little else ">")

        fh.seek(d0)
        buf = fh.read(d1 - d0 + 1)
        need = n_events * n_par * np_dt.itemsize
        if len(buf) < need:
            raise ValueError("DATA segment holds %d bytes, need %d" % (len(buf), need))
        X = np.frombuffer(buf[:need], dtype=np_dt).reshape(n_events, n_par)

        names, meta = [], {}
        for i in range(n_par):
            n = text.get("$P%dN" % (i + 1), "P%d" % (i + 1))
            s = text.get("$P%dS" % (i + 1), "")
            names.append("%s :: %s" % (n, s) if s else n)
            meta["$P%dR" % (i + 1)] = text.get("$P%dR" % (i + 1))
        meta["text"] = text
    return FlowData(X.astype(np.float32), names, path, meta)


def read_flowjo_rds(path: str, rscript: str = RSCRIPT) -> FlowData:
    """Read the inputCSV slot of a FlowJo FJObj. R is used only to open an R-native file."""
    if not os.path.exists(rscript):
        raise FileNotFoundError("Rscript not found at %s" % rscript)
    out = tempfile.mkdtemp(prefix="rna2flow_rds_")
    npy, cols = os.path.join(out, "m.csv"), os.path.join(out, "c.txt")
    script = (
        'x <- readRDS("%s"); m <- attr(x,"inputCSV")[[1]];'
        'write.table(m, "%s", sep=",", row.names=FALSE, col.names=FALSE);'
        'writeLines(colnames(m), "%s")' % (path, npy, cols)
    )
    r = subprocess.run([rscript, "-e", script], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("Rscript failed reading %s:\n%s" % (path, r.stderr[-2000:]))
    X = np.loadtxt(npy, delimiter=",", dtype=np.float32)
    channels = [l.rstrip("\n") for l in open(cols)]
    return FlowData(X, channels, path, {"reader": "flowjo_rds"})


def read_matrix_csv(path: str) -> FlowData:
    df = pd.read_csv(path)
    num = df.select_dtypes(include=[np.number])
    return FlowData(num.to_numpy(np.float32), list(num.columns), path, {"reader": "csv"})
