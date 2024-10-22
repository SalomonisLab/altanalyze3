import os,sys
import pandas as pd
import glob
from tqdm import tqdm

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
            # Generate junctions by pairing adjacent exons
            junctions = [f"{exons[i]}-{exons[i+1]}" for i in range(len(exons) - 1)]
            # Create an entry for the gene if it doesn't exist
            if gene not in transcript_dict:
                transcript_dict[gene] = {}

            # Associate the isoform with each junction
            if dataType == 'junction':
                for junction in junctions:
                    if junction not in transcript_dict[gene]:
                        transcript_dict[gene][junction] = []
                    transcript_dict[gene][junction].append(isoform)
            else:
                transcript_dict[gene][isoform] = exon_sequence
    return transcript_dict

# 2. Load gene symbols {Ensembl ID: Gene Symbol}
def load_gene_symbols(file_path):
    gene_symbol_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:  # Ensure there are at least two columns
                ensembl_id, gene_symbol = parts[:2]
                gene_symbol_dict[ensembl_id] = gene_symbol
    return gene_symbol_dict


# 3. Load junction coordinates {junction_id: coordinates}
def load_junction_coordinates(file_path):
    coord_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            feature, _ = line.strip().split('\t', 1)
            try: junction_id, coords = feature.split('=')
            except:
                continue # header
            coord_dict[junction_id] = coords
    return coord_dict

# 4. Load protein summary into a dictionary {isoform: (length, NMD status)}
def load_protein_summary(file_path):
    protein_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene_id, isoform, length, nmd_status, _, longest_length = line.strip().split('\t')
            if '_' in isoform:
                isoform = isoform.split('_')[1]  # Remove dataset prefix
            try: length = int(length)
            except:
                continue
            protein_dict[gene_id,isoform] = (length, nmd_status, longest_length)
    return protein_dict

# 5. Select the longest isoform for each junction
def get_longest_isoform(gene,isoforms, protein_dict):
    longest_isoform = None
    max_length = -1
    for isoform in isoforms:
        length, _, _ = protein_dict.get((gene,isoform), (0, None, None))
        #print(isoform,[length]);sys.exit()
        if length > max_length:
            max_length = length
            longest_isoform = f"{isoform}|{length}"
    return longest_isoform or 'N/A'

# 6. Annotate a stats file
def annotate_junction_stats_file(file_path, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, output_folder):
    stats_df = pd.read_csv(file_path, sep='\t')

    annotations = []  # Store annotated rows

    for index, row in stats_df.iterrows():
        feature = row['Feature']
        junction1, junction2 = feature.split('|')

        gene_id = junction1.split(':')[0]  # Extract gene ID
        gene_symbol = gene_symbol_dict.get(gene_id, 'N/A')

        # Get isoforms and select the longest one for each junction
        isoforms1 = transcript_dict.get(gene_id, {}).get(junction1.split(':')[1], [])
        isoforms2 = transcript_dict.get(gene_id, {}).get(junction2.split(':')[1], [])
        longest_isoform1 = get_longest_isoform(gene_id,isoforms1, protein_dict)
        longest_isoform2 = get_longest_isoform(gene_id,isoforms2, protein_dict)
        # Extract coordinates
        coord1 = coord_dict.get(junction1, 'N/A')
        coord2 = coord_dict.get(junction2, 'N/A')

        # Add new annotated row
        annotations.append([
            gene_symbol, longest_isoform1, longest_isoform2, coord1, coord2
        ])

    # Create DataFrame with annotations
    annotated_df = pd.DataFrame(
        annotations, 
        columns=['Gene_Symbol', 'Isoform_1|Length', 'Isoform_2|Length', 'Coord_1', 'Coord_2']
    )

    # Merge with original stats DataFrame
    final_df = pd.concat([annotated_df, stats_df], axis=1)

    # Save annotated file
    output_file = os.path.join(output_folder, os.path.basename(file_path).replace('.txt', '-annotated.txt'))
    final_df.to_csv(output_file, sep='\t', index=False)

# 6. Annotate a stats file
def annotate_iso_stats_file(file_path, transcript_dict, gene_symbol_dict, protein_dict, output_folder):
    stats_df = pd.read_csv(file_path, sep='\t')

    annotations = []  # Store annotated rows

    for index, row in stats_df.iterrows():
        feature = row['Feature']

        gene_id, isoform = feature.split(':')  # Extract gene ID
        gene_symbol = gene_symbol_dict.get(gene_id, 'N/A')

        # Get isoforms and select the longest one for each junction
        junctions = transcript_dict[gene_id][isoform]
        try:
            prot_len, nmd_status, longest = protein_dict[gene_id,isoform]
        except:
            prot_len, nmd_status, longest = 0, 'NMD', 'UNK'

        # Add new annotated row
        annotations.append([gene_symbol, junctions, f'{prot_len}/{longest}', nmd_status])

    # Create DataFrame with annotations
    annotated_df = pd.DataFrame(
        annotations, 
        columns=['Gene_Symbol', 'Junctions', 'Protein Length', 'Predicted NMD']
    )

    # Merge with original stats DataFrame
    final_df = pd.concat([annotated_df, stats_df], axis=1)

    # Save annotated file
    output_file = os.path.join(output_folder, os.path.basename(file_path).replace('.txt', '-annotated.txt'))
    final_df.to_csv(output_file, sep='\t', index=False)

# 7. Annotate all stats files in the folder
def annotate_all_stats_files(stats_folder, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, dataType='junction'):
    output_folder = stats_folder

    # Get a list of all .txt files excluding '-annotated.txt' files
    stats_files = [f for f in glob.glob(f"{stats_folder}/*.txt") if '-annotated' not in f]

    # Use tqdm to show progress over the stats files
    for stats_file in tqdm(stats_files, desc="Annotating Stats Files"):
        # Check if the annotated version already exists
        annotated_file = os.path.join(output_folder, os.path.basename(stats_file).replace('.txt', '-annotated.txt'))
        if os.path.exists(annotated_file):
            #print(f"Skipping {stats_file}, already annotated.")
            continue  # Skip this file

        try:
            if dataType == 'junction':
                annotate_junction_stats_file(stats_file, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, output_folder)
            else:
                annotate_iso_stats_file(stats_file, transcript_dict, gene_symbol_dict, protein_dict, output_folder)
        except ZeroDivisionError as e:
            print(f"Error processing {stats_file}: {e}")
            continue


def annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType='junction'):
    if junction_coords_file == None:
        #isoform expression/ratio
        coord_dict = None
    else:
        coord_dict = load_junction_coordinates(junction_coords_file)
    transcript_dict = load_transcript_associations(transcript_assoc_file,dataType=dataType)
    gene_symbol_dict = load_gene_symbols(gene_symbol_file)
    protein_dict = load_protein_summary(protein_summary_file)

    # Annotate all stats files
    annotate_all_stats_files(stats_folder, transcript_dict, gene_symbol_dict, coord_dict, protein_dict,dataType=dataType)

