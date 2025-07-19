import anndata
import pandas as pd
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2]))  # Up to altanalyze3
import altanalyze3.components.long_read.gff_process as gff_process


def annotate_junctions(adata, exon_file):
    exonCoordinates, geneData, strandData = gff_process.importEnsemblGenes(exon_file,include_introns=True)
    annotations = []

    for idx, row in adata.var.iterrows():
        chr = row["chr"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        # Assign donor and acceptor based on strand
        if strand == "+":
            uid1 = (chr, start, strand, 2)  # Donor
            uid2 = (chr, end, strand, 1)    # Acceptor
            pos1, pos2 = start, end
        else:
            uid1 = (chr, end, strand, 2)    # Donor on negative strand
            uid2 = (chr, start, strand, 1)  # Acceptor on negative strand
            pos1, pos2 = end, start

        def get_annot(uid, pos):
            if uid in exonCoordinates:
                gene, exon, _ = exonCoordinates[uid]
                return f"{gene}:{exon}"
            else:
                for gene in geneData:
                    if geneData[gene][0][0] <= pos <= geneData[gene][-1][1] and strandData[gene] == strand:
                        return f"{gene}:{gff_process.findNovelSpliceSite(gene, pos, strand)}"
                return f"NA:{pos}"

        annot1 = get_annot(uid1, pos1)
        annot2 = get_annot(uid2, pos2)

        # Always report both gene:exon components
        annotations.append(f"{annot1}-{annot2}")

    adata.var["annotation"] = annotations


def export_dense_matrix(adata, out_path):
    uids = (
        adata.var["chr"].astype(str) + ":"
        + adata.var["start"].astype(str) + "-"
        + adata.var["end"].astype(str) + ":"
        + adata.var["strand"].astype(str)
    )
    uids.name = "uid"

    df = pd.DataFrame(
        adata.X.toarray().T,
        index=uids,
        columns=adata.obs_names
    )

    df.insert(0, "annotation", adata.var["annotation"].values)
    df.index.name = "uid"
    df.to_csv(out_path, sep="\t")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--exon_file", required=True, help="Exon annotation file")
    parser.add_argument("--export_dense", action="store_true", help="Export dense matrix")
    args = parser.parse_args()

    adata = anndata.read_h5ad(args.h5ad)
    annotate_junctions(adata, args.exon_file)

    annotated_path = Path(args.h5ad).with_name(Path(args.h5ad).stem + "_annotated.h5ad")
    adata.write(annotated_path)

    if args.export_dense:
        export_path = annotated_path.with_suffix(".tsv")
        export_dense_matrix(adata, export_path)
        print(f"Dense matrix exported to: {export_path}")
