import pandas as pd
import numpy as np
import anndata
import scipy.sparse as sp
from pathlib import Path
from typing import List, Literal
import logging
import glob


def load_counts_bed(files: List[Path], sample_names: List[str], splice_type: Literal["exon", "intron"], selected_chr: List[str]) -> pd.DataFrame:
    """
    Load junction/intron BED files, create sparse matrix, and attach UID index.
    Each file must contain: chr, start, end, name, score, strand
    """
    df_all = None
    for file, sample in zip(files, sample_names):
        df = pd.read_csv(file, sep='\t', header=None, names=["chr", "start", "end", "name", "score", "strand"], dtype={"chr": str})
        df = df[df["chr"].isin(selected_chr)]

        # UID: chr:start-end
        df["uid"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        df.set_index("uid", inplace=True)

        # Set up sample column with score
        df = df[["chr", "start", "end", "strand"]].copy().assign(**{sample: df["score"].astype(np.int32)})

        if df_all is None:
            df_all = df
        else:
            df_all = df_all.join(df, how="outer")

    df_all.fillna(0, inplace=True)
    df_all[sample_names] = df_all[sample_names].astype(pd.SparseDtype("int32", 0))
    df_all["type"] = splice_type
    return df_all


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
    adata.write(out_path, compression="gzip")


def aggregate(args):

    def resolve_files(path_str):
        p = Path(path_str)
        return sorted(p.glob("*.bed")) if p.is_dir() else [p]

    junc_files = resolve_files(args.juncounts) if args.juncounts else []
    intron_files = resolve_files(args.intcounts) if args.intcounts else []

    # Derive sample names
    sample_files = junc_files or intron_files
    sample_names = [f.stem.replace("_juncounts", "").replace("_intcounts", "") for f in sample_files]

    chromosomes = args.chr if hasattr(args, "chr") and args.chr else [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    logging.info("Loading exon-exon junctions")
    junc_df = load_counts_bed(junc_files, sample_names, splice_type="exon", selected_chr=chromosomes) if junc_files else None

    logging.info("Loading exon-intron junctions")
    intron_df = load_counts_bed(intron_files, sample_names, splice_type="intron", selected_chr=chromosomes) if intron_files else None

    if junc_df is not None and intron_df is not None:
        combined_df = pd.concat([junc_df, intron_df])
    else:
        combined_df = junc_df or intron_df

    adata = build_adata(combined_df, sample_names)
    adata = sort_adata_by_genomic_position(adata)

    output_path = Path(args.output).with_suffix(".h5ad")
    write_h5ad_compressed(adata, output_path)

    logging.info(f"Exported aggregated counts to {output_path}")
