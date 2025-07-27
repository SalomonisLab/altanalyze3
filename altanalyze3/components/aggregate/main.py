import sys
import pandas as pd
import numpy as np
import anndata
import scipy.sparse as sp
from pathlib import Path
from typing import List, Literal
import logging
import glob
from collections import defaultdict

# Step 1: Collect global UID metadata
def collect_all_uids(files, selected_chr):
    all_uids = set()
    meta_data = {}
    for file in files:
        df = pd.read_csv(file, sep='\t', header=None,
                         names=["chr", "start", "end", "name", "score", "strand"],
                         dtype={"chr": str})
        df = df[df["chr"].isin(selected_chr)]
        df["uid"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        for row in df.itertuples(index=False):
            uid = row.uid
            all_uids.add(uid)
            if uid not in meta_data:
                meta_data[uid] = {"chr": row.chr, "start": int(row.start), "end": int(row.end), "strand": row.strand}

    meta_df = pd.DataFrame.from_dict(meta_data, orient="index")
    meta_df.index.name = "uid"
    return meta_df.sort_index()

def collect_sample_counts(files, selected_chr):
    sample_uid_counts = defaultdict(lambda: defaultdict(int))
    for file in files:
        df = pd.read_csv(file, sep='\t', header=None,
                         names=["chr", "start", "end", "name", "score", "strand"],
                         dtype={"chr": str})
        df = df[df["chr"].isin(selected_chr)]
        df["uid"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        sample = file.stem.replace("_juncounts", "").replace("_intcounts", "")
        for row in df.itertuples(index=False):
            sample_uid_counts[sample][row.uid] += int(row.score)
    return sample_uid_counts

def build_sparse_matrix(uid_index, sample_uid_counts):
    sample_names = sorted(sample_uid_counts.keys())
    row, col, data = [], [], []
    for col_idx, sample in enumerate(sample_names):
        for uid, count in sample_uid_counts[sample].items():
            if uid in uid_index:
                row.append(col_idx)
                col.append(uid_index[uid])
                data.append(count)
    X = sp.csr_matrix((data, (row, col)), shape=(len(sample_names), len(uid_index)), dtype=np.int32)
    return X, sample_names

def build_adata(junc_df: pd.DataFrame, sample_names: List[str]) -> anndata.AnnData:
    """Convert combined junction DataFrame to AnnData."""
    var = junc_df[["chr", "start", "end", "strand", "type"]].copy()
    X = sp.csr_matrix(junc_df[sample_names].values.T)  # TRANSPOSE for correct shape
    obs = pd.DataFrame(index=sample_names)
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    return adata

def sort_adata_by_genomic_position(adata: anndata.AnnData) -> anndata.AnnData:
    sorter = adata.var.sort_values(by=["chr", "start", "end"]).index
    return adata[:, sorter].copy()  # Force copy to avoid ImplicitModificationWarning

def write_h5ad_compressed(adata: anndata.AnnData, out_path: Path):
    for col in ["chr", "strand", "type"]:
        if col in adata.var.columns:
            adata.var[col] = adata.var[col].astype(str)
    adata.write(out_path, compression="gzip")

def aggregate(args):
    def resolve_files(x):
        if isinstance(x, list):
            return x
        p = Path(x)
        if p.is_dir():
            return sorted(p.glob("*.bed"))
        elif p.exists():
            return [p]
        raise FileNotFoundError(f"Invalid path: {x}")

    junc_files = resolve_files(args.juncounts) if args.juncounts else []
    intron_files = resolve_files(args.intcounts) if args.intcounts else []
    all_files = junc_files + intron_files
    chromosomes = args.chr if hasattr(args, "chr") and args.chr else [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # Collect and align UID metadata
    uid_meta_df = collect_all_uids(junc_files + intron_files, chromosomes)
    all_uids = uid_meta_df.index.tolist()
    uid_index = {uid: i for i, uid in enumerate(all_uids)}

    # Collect counts across files
    sample_uid_counts = collect_sample_counts(junc_files + intron_files, chromosomes)

    # Build sparse matrix
    X, sample_names = build_sparse_matrix(uid_index, sample_uid_counts)

    # Final AnnData
    uid_meta_df["type"] = "mixed"
    adata = anndata.AnnData(X=X, obs=pd.DataFrame(index=sample_names), var=uid_meta_df)
    adata = sort_adata_by_genomic_position(adata)

    output_path = Path(args.output).with_suffix(".h5ad")
    write_h5ad_compressed(adata, output_path)

    logging.info(f"Exported aggregated counts to {output_path}")

    # === Annotate junctions and export dense matrix ===
    from altanalyze3.components.aggregate.annotate import annotate_junctions, export_dense_matrix

    # Assume exon annotation file is named like "gene_model_all.tsv"
    exon_file = Path(args.ref).with_name(Path(args.ref).stem.replace(".bed", "_all.tsv"))
    if not exon_file.exists():
        raise FileNotFoundError(f"Exon annotation file not found: {exon_file}")

    adata = anndata.read_h5ad(output_path)
    annotate_junctions(adata, exon_file)

    annotated_path = output_path.with_name(output_path.stem + "_annotated.h5ad")
    adata.write(annotated_path)
    export_dense_matrix(adata, annotated_path.with_suffix(".tsv"))
    logging.info(f"Annotated file written to {annotated_path}")
