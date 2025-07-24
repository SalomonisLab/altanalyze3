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
        if strand == '-':
            start,end = end,start

        def get_annot(chr,pos,strand,site):
            uid = (chr,pos,strand,site)
            if uid in exonCoordinates:
                gene, exon, _ = exonCoordinates[uid]
                return gene, exon #f"{gene}:{exon}"
            else:
                if site == 1: 
                    site = 2 # splice donor
                else: 
                    site = 1 # splice acceptor
                uid = (chr,pos,strand, site)
                if uid in exonCoordinates:
                    gene, exon, _ = exonCoordinates[uid]
                    return gene, exon
                else:
                    return None, None

        def find_completely_novel_annot(chr,pos,strand):
            for gene in geneData:
                if geneData[gene][0][0] <= pos <= geneData[gene][-1][1] and strandData[gene] == strand:
                    return gene, gff_process.findNovelSpliceSite(gene, pos, strand)
            return None, None #f"NA:{pos}"

        gene1,exon1 = get_annot(chr, start, strand, 2)
        gene2,exon2 = get_annot(chr, end, strand, 1)

        if strand == '-' and (end>start):
            # Correct naming and order 
            gene1,exon1,gene2,exon2 = gene2,exon2,gene1,exon1
            start,end = end,start
        if gene1 != gene2: # novel splice site or gene
            if gene1 == None:
                exon1 = gff_process.findNovelSpliceSite(gene2, start, strand)
                annotations.append(f"{gene2}:{exon1}-{exon2}={chr}:{start}-{end}")
            elif gene2 == None:
                exon2 = gff_process.findNovelSpliceSite(gene1, end, strand)
                annotations.append(f"{gene1}:{exon1}-{exon2}={chr}:{start}-{end}")
            else:
                # trans-splicing
                annotations.append(f"{gene1}:{exon1}-{gene2}:{exon2}={chr}:{start}-{end}")
        elif gene1 == None: #both == None
            gene1, exon1 = find_completely_novel_annot(chr,start,strand)
            gene2, exon2 = find_completely_novel_annot(chr,start,strand)
            if gene1 == gene2:
                annotations.append(f"{gene1}:{exon1}-{exon2}={chr}:{start}-{end}")
            else:
                annotations.append(f"{gene1}:{exon1}-{gene2}:{exon2}={chr}:{start}-{end}")
        else:
            annotations.append(f"{gene1}:{exon1}-{exon2}={chr}:{start}-{end}")

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
