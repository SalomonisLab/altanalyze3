import os
import pandas as pd
import glob
from concurrent.futures import ProcessPoolExecutor

# File paths
transcript_assoc_file = 'gff-output/transcript_associations.txt'
gene_symbol_file = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl-annotations.txt'
junction_coords_file = 'junction_combined_pseudo_cluster_counts-filtered.txt'
protein_summary_file = 'protein_summary.txt'
stats_folder = 'dPSI-events'

# 1. Load transcript associations into a nested dictionary {gene: {junction: [isoforms]}}
def load_transcript_associations(file_path, dataType='junction'):
    transcript_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene, _, exon_sequence, isoform, _ = line.strip().split('\t')
            exons = exon_sequence.split('|')
            junctions = [f"{exons[i]}-{exons[i+1]}" for i in range(len(exons) - 1)]

            if gene not in transcript_dict:
                transcript_dict[gene] = {}
            if dataType == 'junction':
                for junction in junctions:
                    transcript_dict[gene].setdefault(junction, []).append(isoform)
            else:
                transcript_dict[gene][isoform] = exon_sequence
    return transcript_dict

# 2. Load gene symbols
def load_gene_symbols(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1]).set_index(0).to_dict()[1]

# 3. Load junction coordinates
def load_junction_coordinates(file_path):
    return dict(line.strip().split('=') for line in open(file_path) if '=' in line)

# 4. Load protein summary
def load_protein_summary(file_path):
    protein_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene_id, isoform, length, nmd_status, _, longest_length = line.strip().split('\t')
            isoform = isoform.split('_')[-1]
            try:
                length = int(length)
            except ValueError:
                continue
            protein_dict[(gene_id, isoform)] = (length, nmd_status, longest_length)
    return protein_dict

# 5. Select the longest isoform
def get_longest_isoform(gene, isoforms, protein_dict):
    max_length = -1
    best_isoform = 'N/A'
    for isoform in isoforms:
        length, _, _ = protein_dict.get((gene, isoform), (0, None, None))
        if length > max_length:
            max_length = length
            best_isoform = f"{isoform}|{length}"
    return best_isoform

# 6. Annotate a stats file
def annotate_stats_file(stats_file, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, output_folder):
    stats_df = pd.read_csv(stats_file, sep='\t')

    annotations = [
        [
            gene_symbol_dict.get(j1.split(':')[0], 'N/A'),  # Gene Symbol
            get_longest_isoform(j1.split(':')[0], transcript_dict.get(j1.split(':')[0], {}).get(j1.split(':')[1], []), protein_dict),
            get_longest_isoform(j2.split(':')[0], transcript_dict.get(j2.split(':')[0], {}).get(j2.split(':')[1], []), protein_dict),
            coord_dict.get(j1, 'N/A'),  # Coord 1
            coord_dict.get(j2, 'N/A')   # Coord 2
        ]
        for j1, j2 in stats_df['Feature'].str.split('|')
    ]

    annotated_df = pd.DataFrame(
        annotations, 
        columns=['Gene_Symbol', 'Isoform_1|Length', 'Isoform_2|Length', 'Coord_1', 'Coord_2']
    )
    final_df = pd.concat([annotated_df, stats_df], axis=1)

    output_file = os.path.join(output_folder, os.path.basename(stats_file).replace('.txt', '-annotated.txt'))
    final_df.to_csv(output_file, sep='\t', index=False)

# 7. Annotate all stats files in the folder
def annotate_all_stats_files(stats_folder, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, dataType='junction'):
    output_folder = stats_folder

    # Get stats files excluding those that already have '-annotated.txt' versions
    stats_files = [
        f for f in glob.glob(f"{stats_folder}/*.txt") 
        if not os.path.exists(f.replace('.txt', '-annotated.txt'))
    ]

    with ProcessPoolExecutor() as executor:
        executor.map(
            annotate_stats_file, stats_files,
            [transcript_dict] * len(stats_files),
            [gene_symbol_dict] * len(stats_files),
            [coord_dict] * len(stats_files),
            [protein_dict] * len(stats_files),
            [stats_folder] * len(stats_files)
        )

def annotate(gene_symbol_file, transcript_assoc_file, junction_coords_file, protein_summary_file, stats_folder, dataType='junction'):
    coord_dict = load_junction_coordinates(junction_coords_file) if junction_coords_file else {}
    transcript_dict = load_transcript_associations(transcript_assoc_file, dataType)
    gene_symbol_dict = load_gene_symbols(gene_symbol_file)
    protein_dict = load_protein_summary(protein_summary_file)

    annotate_all_stats_files(stats_folder, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, dataType)
