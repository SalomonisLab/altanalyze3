"""ISV web data layer.

Wraps the validated, matplotlib-free isoform_structure_view engine + its precompute caches to serve
fast, JSON-able isoform-structure queries. No drawing here -- the frontend renders SVG from this data.

Reused engine functions (do NOT reimplement):
  * load_isoform_counts_indexed   -> per-(gene x sample x cell_states) records from gene_indexes_v2 caches
  * build_plotted_isoforms        -> attach token structures, trim, normalize, rank to max_isoforms
  * group_structures              -> cluster isoforms (exposes thresholds + filter_junctions)
  * build_isoform_segments        -> per-isoform exon/intron segments with genomic coords (mouseover)
  * read_gene_model / *_from_sqlite -> exon_lookup
  * load_transcript_associations_auto -> transcript token structures (uses transcripts.db if present)
"""
from __future__ import annotations

import os
import threading
from collections import defaultdict
from functools import lru_cache

import numpy as np

from .. import isoform_structure_view as isv
from ...annotation.junction_isoform import load_gene_symbols


# --------------------------------------------------------------------------- run context

class RunContext:
    """Holds everything loaded once at startup for fast queries: per-sample h5ad paths + their
    gene_indexes_v2 dirs, the gene model exon_lookup, transcript-structure source, protein indexes,
    and the gene-symbol map. All lookups below read precompute caches, never the full h5ad."""

    def __init__(self, run_dir, metadata_file, gene_model_path, gene_symbol_path,
                 gff_output_dir=None):
        self.run_dir = os.path.abspath(run_dir)
        self.metadata_file = metadata_file
        self.gene_model_path = gene_model_path
        self.gene_symbol_path = gene_symbol_path
        self.gff_output_dir = gff_output_dir or os.path.join(self.run_dir, "gff-output")

        # sample -> {h5ad, library, uid, group, reverse}
        self.samples = {}              # keyed by library name (the column id in pseudobulks)
        self.sample_order = []
        self.groups = {}               # group -> [library, ...]
        self.cell_types = []           # union of obs['cluster'] across barcode cluster files
        self.barcode_series = {}       # library -> pd.Series(barcode -> cluster)

        # gene model / transcripts
        self.exon_lookup = {}          # (gene, exon_id) -> {chrom, strand, start, end, type}
        self.transcripts_db = None     # path to transcripts.db (fast per-gene token lookup) or None
        self.transcript_assoc_path = None

        # gene symbols
        self.gene_symbol_dict = {}     # ENSG -> symbol
        self.symbol_to_gene = {}       # SYMBOL(upper) -> ENSG
        self.all_genes = []            # ENSG list seen in var_names

        # protein layer
        self._protein_meta = {}        # transcript_id -> {protein_length, nmd_status, intron_retention, longest}
        self._fasta_index = {}         # transcript_id -> (path, offset, length) for protein_sequences.fasta
        self._mrna_index = {}          # transcript_id -> (path, offset, length) for transcript_sequences.fasta
        self._fasta_lock = threading.Lock()
        self._gene_index_cache = {}    # h5ad path -> (gene_to_pos, indptr, indices), loaded once
        self._struct_cache = {}        # gene -> {transcript_id: [tokens]} (known + novel)

        # combined pseudobulk matrix (THE fast counts source): one h5ad, obs = '<cellType>.<library>-isoform'
        # columns, var = '<ENSG>:<isoform>' isoforms. Composition = slice the gene's vars + sum the
        # selected columns. Loaded once at startup.
        self._adata = None             # the combined AnnData (X = obs x var counts)
        self._X = None                 # CSC for fast column(=isoform) slicing
        self._gene_to_vars = {}        # ENSG -> [var indices] (per-gene isoform columns)
        self._col_meta = []            # per obs row: {"col": name, "sample": library, "cell_type": ct}

    # -- loaders -----------------------------------------------------------

    def load_all(self):
        self._load_metadata_and_clusters()
        self._load_combined_matrix()
        self._load_gene_model()
        self._load_gene_symbols()
        self._index_genes_from_var_names()
        self._index_protein_summary()
        self._index_protein_fasta()
        return self

    def _load_metadata_and_clusters(self):
        from ...long_read.isoform_automate import import_metadata
        from ...long_read.isoform_matrix import import_barcode_clusters
        from ...long_read import cli as _cli

        sample_dict = import_metadata(self.metadata_file)
        for uid, libs in sample_dict.items():
            for lib in libs:
                library = lib.get("library", uid)
                group = lib.get("groups") or lib.get("group") or "ungrouped"
                h5 = lib.get("matrix")
                if h5 is None:
                    # per-sample molecule h5ad sits beside the BAM as <library>-isoform.h5ad
                    cand = _find_sample_h5ad(self.run_dir, library)
                    h5 = cand
                self.samples[library] = {
                    "h5ad": h5, "library": library, "uid": uid,
                    "group": group, "reverse": bool(lib.get("reverse", False)),
                }
                self.sample_order.append(library)
                self.groups.setdefault(group, []).append(library)

        # cell types from the per-sample barcode->cluster TSVs
        try:
            bcd = _cli._discover_barcode_cluster_dirs(self.metadata_file, None)
            bsd = import_barcode_clusters(bcd)
            cell_types = set()
            for sample_key, series in bsd.items():
                self.barcode_series[sample_key] = series
                cell_types.update(str(v) for v in series.unique())
            self.cell_types = sorted(cell_types)
        except Exception as exc:  # pragma: no cover - menu degrades gracefully
            print(f"[isv_web] cell-type discovery failed: {exc}")
            self.cell_types = []

    def _load_combined_matrix(self):
        """Load the combined isoform pseudobulk h5ad once. obs = '<cellType>.<library>-isoform'
        columns; var = '<ENSG>:<isoform>'. Build a CSC view for fast column slicing and a per-gene
        var-index map. This is the fast counts source -- every query slices + sums columns in-memory."""
        import anndata as ad
        import numpy as np
        from scipy import sparse
        path = _find_combined_h5ad(self.run_dir)
        if not path:
            print("[isv_web] combined pseudobulk h5ad not found; counts queries will be unavailable")
            return
        a = ad.read_h5ad(path)
        self._adata = a
        X = a.X
        self._X = X.tocsc() if sparse.issparse(X) else sparse.csc_matrix(np.asarray(X))
        # per-gene var indices
        g2v = {}
        for vi, v in enumerate(a.var_names):
            v = str(v)
            gid = v.split(":", 1)[0] if ":" in v else None
            if gid:
                g2v.setdefault(gid, []).append(vi)
        self._gene_to_vars = g2v
        # parse each obs column into (sample/library, cell_type)
        self._col_meta = [self._parse_column(str(c)) for c in a.obs_names]

    def _parse_column(self, col):
        """'<cellType>.<library>-isoform' -> {'col', 'cell_type', 'sample'(library)}.
        The library may itself contain '.'/'-', so split cell_type at the FIRST '.', then strip a
        trailing dataType suffix from the remainder."""
        name = col
        for suf in ("-isoform", "-junction"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        if "." in name:
            cell_type, sample = name.split(".", 1)
        else:
            cell_type, sample = "", name
        return {"col": col, "cell_type": cell_type, "sample": sample}

    def _load_gene_model(self):
        # Prefer the transcripts.db (fast per-gene); keep gene_model_path for exon_lookup fallback.
        db = _find_transcripts_db(self.run_dir, self.gff_output_dir)
        self.transcripts_db = db
        self.transcript_assoc_path = os.path.join(self.gff_output_dir, "transcript_associations.txt")
        # exon_lookup: full gene model loaded once (coords for every exon region). It's a flat TSV
        # read; do it at startup so per-query mouseover coords are in-memory.
        if self.gene_model_path and os.path.exists(self.gene_model_path):
            _, self.exon_lookup = isv.read_gene_model(self.gene_model_path)

    def _load_gene_symbols(self):
        if self.gene_symbol_path and os.path.exists(self.gene_symbol_path):
            self.gene_symbol_dict = load_gene_symbols(self.gene_symbol_path) or {}
            for ensg, sym in self.gene_symbol_dict.items():
                if sym:
                    self.symbol_to_gene[str(sym).upper()] = ensg

    def _index_genes_from_var_names(self):
        # Genes present = keys of the combined matrix's per-gene var map (O(#genes), already built).
        self.all_genes = sorted(self._gene_to_vars.keys())

    def _index_protein_summary(self):
        path = os.path.join(self.gff_output_dir, "protein_summary.txt")
        if not os.path.exists(path):
            return
        with open(path) as fh:
            header = fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                tid = p[1].strip()
                self._protein_meta[tid] = {
                    "protein_length": _to_int(p[2]),
                    "nmd_status": p[3].strip() if len(p) > 3 else "",
                    "intron_retention": p[4].strip() if len(p) > 4 else "",
                    "longest_isoform_length": _to_int(p[5]) if len(p) > 5 else None,
                }

    def _index_protein_fasta(self):
        self._fasta_index = _build_fasta_offset_index(
            os.path.join(self.gff_output_dir, "protein_sequences.fasta"))
        # mRNA / full transcript sequences (header '>TRANSCRIPT_ID ;transcript;gene_id:ENSG...')
        self._mrna_index = _build_fasta_offset_index(
            os.path.join(self.gff_output_dir, "transcript_sequences.fasta"))

    # -- protein / mRNA lookups (tolerant to ENST vs novel molecule.sample, version, (NMD)) ---------

    def protein_meta(self, isoform_id):
        for key in _id_lookup_keys(isoform_id):
            m = self._protein_meta.get(key)
            if m:
                return m
        return None

    def protein_fasta(self, isoform_id):
        return self._seek_fasta(self._fasta_index, isoform_id)

    def mrna_fasta(self, isoform_id):
        return self._seek_fasta(getattr(self, "_mrna_index", {}), isoform_id)

    def _seek_fasta(self, index, isoform_id):
        """Return the FASTA record text for an isoform id (trying tolerant key variants), or None."""
        entry = None
        for key in _id_lookup_keys(isoform_id):
            entry = index.get(key)
            if entry:
                break
        if not entry:
            return None
        path, start, length = entry
        with self._fasta_lock, open(path, "rb") as fh:
            fh.seek(start)
            return fh.read(length).decode("utf-8", "replace")

    # -- gene search -------------------------------------------------------

    def resolve_gene(self, query):
        """A user query (symbol or ENSG) -> ENSG present in the data (or None)."""
        q = (query or "").strip()
        if not q:
            return None
        if q in self.all_genes:
            return q
        ensg = self.symbol_to_gene.get(q.upper())
        if ensg and ensg in self.all_genes:
            return ensg
        # bare ENSG without version
        if q.upper().startswith("ENSG"):
            base = q.split(".")[0]
            if base in self.all_genes:
                return base
        return None

    def search_genes(self, prefix, limit=25):
        p = (prefix or "").strip().upper()
        if not p:
            return []
        out = []
        for ensg in self.all_genes:
            sym = self.gene_symbol_dict.get(ensg, "")
            if sym and sym.upper().startswith(p):
                out.append({"gene": ensg, "symbol": sym})
            elif ensg.upper().startswith(p):
                out.append({"gene": ensg, "symbol": sym})
            if len(out) >= limit:
                break
        return out


# --------------------------------------------------------------------------- query

def query_isoforms(ctx: RunContext, gene, samples, groups, cell_types, combine_by="group",
                   cluster_similarity_threshold=0.85, min_split_fraction=0.05, min_count=1,
                   max_isoforms=1500, cluster_strategy="subsequence", cluster_mode="block",
                   filter_junctions=None, include_introns=True):
    """Core query. Composes counts across the selected (sample x cell_type) atoms, clusters, and
    returns JSON-able isoforms with exon_segments + per-column expression + protein meta.

    `combine_by`: 'group' -> columns are the metadata groups; 'cell_type' -> columns are cell types.
    Returns dict {gene, symbol, columns[], isoforms[...], cluster_count}.
    """
    gene = ctx.resolve_gene(gene) or gene
    sel_samples = _resolve_selected_samples(ctx, samples, groups)
    cell_states = list(cell_types) if cell_types else None

    # transcript token structures for this gene (fast: from transcripts.db if present).
    # load_transcript_associations_auto returns a 2-tuple (structures_dict, supplementary_dict); we
    # need the first element {transcript_id: [tokens]}.
    transcript_structures = _load_structures(ctx, gene)

    # define output COLUMNS as structured objects {key, label, group}. `group` lets the frontend draw
    # covariate block labels + separators. Three modes:
    #   group               -> one column per metadata group (covariate)
    #   cell_type           -> one column per selected cell type
    #   cell_type_x_covariate -> one column per (covariate, cell_type), grouped/ordered by covariate
    sel_groups = _sorted_groups(ctx, sel_samples)
    if combine_by == "cell_type_x_covariate" and cell_states:
        col_defs = []
        for grp in sel_groups:
            for ct in cell_states:
                col_defs.append({"key": f"{grp}|{ct}", "label": ct, "group": grp})
    elif combine_by == "cell_type" and cell_states:
        col_defs = [{"key": ct, "label": ct, "group": None} for ct in cell_states]
    else:
        combine_by = "group"
        col_defs = [{"key": g, "label": g, "group": g} for g in sel_groups]

    col_counts, total_counts = _compose_counts(ctx, gene, sel_samples, cell_states, combine_by, col_defs)
    columns = col_defs

    # build the isoforms list the engine expects (one per isoform id, total count for ranking)
    isoform_records = [{"isoform_id": iid, "count": cnt} for iid, cnt in total_counts.items()]
    plotted, _, _ = isv.build_plotted_isoforms(
        isoform_records, transcript_structures, gene,
        cluster_mode=cluster_mode, cluster_features="tokens",
        min_count=min_count, max_isoforms=max_isoforms,
    )
    grouped = isv.group_structures(
        plotted, gene, cluster_features="tokens", cluster_strategy=cluster_strategy,
        cluster_similarity_threshold=cluster_similarity_threshold,
        min_split_fraction=min_split_fraction, include_introns=include_introns,
        cluster_mode=cluster_mode, exon_lookup=ctx.exon_lookup,
        filter_junctions=list(filter_junctions) if filter_junctions else None,
    )

    # serialize. group_structures returns: List[cluster] where cluster = List[structure_group]; each
    # structure_group has 'structure_key' (shared tokens) and 'items' (the isoforms that share it).
    out_isoforms = []
    for cluster_id, cluster in enumerate(grouped):
        structure_groups = cluster if isinstance(cluster, list) else [cluster]
        for sg in structure_groups:
            sg_tokens = list(sg.get("structure_key") or sg.get("structure_key_trimmed") or [])
            for it in sg.get("items", []):
                iid = it.get("resolved_id") or it.get("isoform_id")
                tokens = it.get("tokens") or it.get("tokens_trimmed") or sg_tokens
                segs = isv.build_isoform_segments(tokens, ctx.exon_lookup, gene)
                pmeta = ctx.protein_meta(iid) or {}
                out_isoforms.append({
                    "isoform_id": iid,
                    "known": bool(str(iid).startswith("ENST")),
                    "cluster_id": cluster_id,
                    "exon_segments": [_seg_json(seg) for seg in segs],
                    "expression": {c["key"]: round(col_counts[c["key"]].get(iid, 0.0), 3) for c in columns},
                    "total_count": round(total_counts.get(iid, 0.0), 3),
                    "protein_length": pmeta.get("protein_length"),
                    "nmd_status": pmeta.get("nmd_status"),
                    "intron_retention": pmeta.get("intron_retention"),
                })

    return {
        "gene": gene,
        "symbol": ctx.gene_symbol_dict.get(gene, gene),
        "combine_by": combine_by,
        "columns": columns,
        "isoforms": out_isoforms,
        "cluster_count": len(grouped),
    }


def _load_structures(ctx: RunContext, gene):
    """{transcript_id: [tokens]} for a gene, INCLUDING novel isoforms.

    NOTE: the prebuilt transcripts.db indexes only the ENST reference transcripts (novel
    'molecule.sample' ids are dropped), but the combined matrix references those novel ids. So read
    from transcript_associations.txt (which contains both known + novel) -- it's a per-gene scan,
    cached on ctx. Returns element [0] of the 2-tuple."""
    cached = ctx._struct_cache.get(gene)
    if cached is not None:
        return cached
    s = {}
    path = ctx.transcript_assoc_path
    if path and os.path.exists(path):
        # Read the TSV DIRECTLY (load_transcript_associations, NOT *_auto) -- *_auto silently redirects
        # to transcripts.db, which indexes only ENST references and drops the novel 'molecule.sample'
        # isoforms the combined matrix references. The TSV has both known + novel.
        res = isv.load_transcript_associations(path, gene)
        s = res[0] if isinstance(res, tuple) else (res or {})
    ctx._struct_cache[gene] = s
    return s


def gene_junctions(ctx: RunContext, gene):
    """Return the gene's junction tokens (for the 'restrict by junctions' control)."""
    gene = ctx.resolve_gene(gene) or gene
    ts = _load_structures(ctx, gene)
    junctions = set()
    for tid, tokens in ts.items():
        trimmed = isv.strip_terminal_coords(tokens, gene)
        toks = isv.filter_exon_intron_tokens(trimmed, gene, normalize_mode="block")
        # adjacent exon pairs are the junctions; keep individual exon/intron region tokens too
        for t in toks:
            junctions.add(t)
    return sorted(junctions)


# --------------------------------------------------------------------------- helpers

def _compose_counts(ctx, gene, sel_samples, cell_states, combine_by, col_defs):
    """Slice the gene's isoform columns from the combined matrix and sum the selected obs rows per
    output column. `col_defs` is the list of {key,label,group}. Returns
    (col_counts {col_key -> {isoform_id: count}}, total_counts {isoform_id: count}). All in-memory."""
    import numpy as np
    keys = [c["key"] for c in col_defs]
    col_counts = {k: {} for k in keys}
    total_counts = {}
    if ctx._X is None:
        return col_counts, total_counts
    var_idx = ctx._gene_to_vars.get(gene, [])
    if not var_idx:
        return col_counts, total_counts

    sel_set = set(sel_samples)
    cs_set = set(cell_states) if cell_states else None
    keyset = set(keys)

    # which obs rows (the pseudobulk's celltype.sample columns) feed each output column key
    rows_for_output = {k: [] for k in keys}
    for ri, m in enumerate(ctx._col_meta):
        if m["sample"] not in sel_set:
            continue
        if cs_set is not None and m["cell_type"] not in cs_set:
            continue
        grp = ctx.samples.get(m["sample"], {}).get("group")
        if combine_by == "cell_type_x_covariate":
            key = f"{grp}|{m['cell_type']}"
        elif combine_by == "cell_type":
            key = m["cell_type"]
        else:
            key = grp
        if key in keyset:
            rows_for_output[key].append(ri)

    var_names = [str(ctx._adata.var_names[vi]) for vi in var_idx]
    Xg = ctx._X[:, var_idx].tocsr()          # obs x (gene isoforms); CSR for row(obs) slicing
    for k in keys:
        rows = rows_for_output[k]
        if not rows:
            continue
        colsum = np.asarray(Xg[rows, :].sum(axis=0)).ravel()
        for j, v in enumerate(colsum):
            if v:
                iid = _isoform_id_from_var(var_names[j], gene)
                col_counts[k][iid] = col_counts[k].get(iid, 0.0) + float(v)
                total_counts[iid] = total_counts.get(iid, 0.0) + float(v)
    return col_counts, total_counts


def _sorted_groups(ctx, sel_samples):
    return sorted({ctx.samples[s]["group"] for s in sel_samples})


def _isoform_id_from_var(var, gene):
    """'<ENSG>:<isoform>' -> '<isoform>' (the id used for structures / protein lookup)."""
    return var.split(":", 1)[1] if ":" in var else var


def _resolve_selected_samples(ctx, samples, groups):
    sel = set()
    if samples:
        for s in samples:
            if s in ctx.samples:
                sel.add(s)
    if groups:
        for g in groups:
            for lib in ctx.groups.get(g, []):
                sel.add(lib)
    if not sel:
        sel = set(ctx.samples.keys())
    return sorted(sel)


def _seg_json(seg):
    return {
        "exon_id": seg.get("exon_id") or seg.get("label") or seg.get("base"),
        "label": seg.get("label"),
        "type": seg.get("type"),
        "start": int(seg["start"]) if seg.get("start") is not None else None,
        "end": int(seg["end"]) if seg.get("end") is not None else None,
        "strand": seg.get("strand"),
    }


def _id_lookup_keys(isoform_id):
    """Yield the candidate keys an isoform id may be stored under in protein_summary / the FASTAs,
    most-specific first. Handles:
      * 'ENSG...:<id>' matrix prefix  -> strip it
      * ENST version  'ENST...x.5'    -> also try the unversioned 'ENST...x'
      * novel 'molecule.sample' (e.g. '3947043.WF40_ND21-251_HSC_3k') -> also try the BARE molecule
        number '3947043' (that's how protein_summary.txt and the protein/transcript FASTAs key novel
        isoforms -- the .sample suffix is NOT in those files).
    """
    s = str(isoform_id)
    keys = [s]
    if ":" in s:
        s = s.split(":", 1)[1]
        keys.append(s)
    if s.startswith("ENST") and "." in s:
        keys.append(s.split(".")[0])
    elif "." in s:
        # novel molecule.sample -> bare molecule id (left of the FIRST dot)
        keys.append(s.split(".", 1)[0])
    # de-dup, keep order
    seen = set()
    out = []
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def _build_fasta_offset_index(path):
    """Byte-offset index {record_id: (path, start, length)} for a FASTA. Seek-read on demand; never
    holds sequences in RAM. record_id = the token after '>' up to the first space."""
    index = {}
    if not path or not os.path.exists(path):
        return index
    offset = 0
    cur_id = None
    cur_start = 0
    with open(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                if cur_id is not None:
                    index[cur_id] = (path, cur_start, offset - cur_start)
                cur_id = raw[1:].split(b" ", 1)[0].decode("utf-8", "replace").strip()
                cur_start = offset
            offset += len(raw)
        if cur_id is not None:
            index[cur_id] = (path, cur_start, offset - cur_start)
    return index


def _to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def _find_sample_h5ad(run_dir, library):
    import glob
    for pat in (f"**/{library}-isoform.h5ad", f"**/{library}.h5ad"):
        hits = glob.glob(os.path.join(run_dir, pat), recursive=True)
        if hits:
            return hits[0]
    return None


def _find_combined_h5ad(run_dir):
    """Prefer the full combined counts h5ad; fall back to the filtered one."""
    for name in ("isoform_combined_pseudo_cluster_counts.h5ad",
                 "isoform_combined_pseudo_cluster_counts-filtered.h5ad"):
        cand = os.path.join(run_dir, name)
        if os.path.exists(cand):
            return cand
    import glob
    hits = glob.glob(os.path.join(run_dir, "**", "isoform_combined_pseudo_cluster_counts*.h5ad"), recursive=True)
    return hits[0] if hits else None


def _find_transcripts_db(run_dir, gff_output_dir):
    import glob
    for cand in (
        os.path.join(gff_output_dir, "gene_indexes_v2", "transcripts.db"),
        os.path.join(gff_output_dir, "transcripts.db"),
    ):
        if os.path.exists(cand):
            return cand
    hits = glob.glob(os.path.join(run_dir, "**", "transcripts.db"), recursive=True)
    return hits[0] if hits else None
