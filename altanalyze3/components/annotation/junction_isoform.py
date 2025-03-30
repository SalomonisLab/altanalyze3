import os,sys
import pandas as pd
import glob
from tqdm import tqdm

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import splice_event as se

""" This module assigns isoform annotations to the statistical output from comparisons.py """

# 1. Load transcript associations into a nested dictionary {gene: {junction: [isoforms]}}
def load_transcript_associations(file_path, dataType='junction'):
    transcript_dict = {}

    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f)  # Count total lines for tqdm
        f.seek(0)  # Reset file pointer

        for line in tqdm(f, total=total_lines, desc="Processing transcript associations", unit=" lines"):
            gene, _, exon_sequence, isoform, _ = line.strip().split('\t')
            exons = exon_sequence.split('|')
            junctions = [f"{exons[i]}-{exons[i+1]}" for i in range(len(exons) - 1)]

            if gene not in transcript_dict:
                transcript_dict[gene] = {}

            if dataType == 'junction':
                for junction in junctions:
                    if junction not in transcript_dict[gene]:
                        transcript_dict[gene][junction] = []
                    transcript_dict[gene][junction].append((isoform, exons))
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
def load_junction_coordinates(file_input):
    coord_dict = {}

    if isinstance(file_input, str):
        with open(file_input, 'r') as f:
            total_lines = sum(1 for _ in f)  # Count total lines for tqdm
            f.seek(0)  # Reset file pointer

            for line in tqdm(f, total=total_lines, desc="Processing junction coordinates", unit=" lines"):
                feature, _ = line.split('\t', 1)
                try:
                    junction_id, coords = feature.split('=')
                    coord_dict[junction_id] = coords
                except ValueError:
                    continue  

    else:
        for line in tqdm(file_input, desc="Processing junction coordinates", unit=" lines"):
            feature, _ = line.split('\t', 1)
            try:
                junction_id, coords = feature.split('=')
                coord_dict[junction_id] = coords
            except ValueError:
                continue  

    return coord_dict

# 4. Load protein summary into a dictionary {isoform: (length, NMD status)}
def load_protein_summary(file_path):
    protein_dict = {}

    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f)  # Count total lines for tqdm
        f.seek(0)  # Reset file pointer

        for line in tqdm(f, total=total_lines, desc="Processing protein summary", unit=" lines"):
            gene_id, isoform, length, nmd_status, _, longest_length = line.strip().split('\t')
            if '_' in isoform:
                isoform = isoform.split('_')[1]  
            if ';' in gene_id:
                gene_id = gene_id.split(';')[0]
            try:
                length = int(length)
            except:
                continue
            protein_dict[(gene_id, isoform)] = (length, nmd_status, longest_length)

    return protein_dict

# 5. Select the longest isoform for each junction
def get_longest_isoform(gene,isoforms,protein_dict):
    longest_isoform = 'N/A'
    longest_exons = 'N/A'
    max_length = -1
    for (isoform,exons) in isoforms:
        length, _, _ = protein_dict.get((gene,isoform), (0, None, None))
        #print(isoform,[length]);sys.exit()
        if length > max_length:
            max_length = length
            longest_isoform = f"{isoform}|{length}"
            longest_exons = exons
    return (longest_isoform,longest_exons)

# 6. Annotate a stats file
def annotate_junction_stats_file(file_path, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, output_folder):
    stats_df = pd.read_csv(file_path, sep='\t')

    annotations = []  # Store annotated rows

    for index, row in stats_df.iterrows():
        feature = row['Feature']
        mwuSign = row['mwuSign']

        if feature in event_isoform_db:
            # Do not re-generate prior produced associations from other stat files
            gene_symbol, longest_isoform1, longest_isoform2, coord1, coord2, event_type, proximity = event_isoform_db[feature]
        else:
            junction1, junction2 = feature.split('|')
            gene_id = junction1.split(':')[0]  # Extract gene ID
            gene_symbol = gene_symbol_dict.get(gene_id, 'N/A')

            # Get isoforms and select the longest one for each junction (junction1.split(":", 1)[1] considers only the first ":" for trans-splicing)
            isoforms1 = transcript_dict.get(gene_id, {}).get(junction1.split(":", 1)[1], [])
            isoforms2 = transcript_dict.get(gene_id, {}).get(junction2.split(":", 1)[1], [])

            if len(isoforms1)==0:
                original_gene_id = str(gene_id)
                original_junction1 = str(junction1)
                try:
                    junction1 = f"{gene_id}:{junction1.split(':')[1].split('-')[0]}-{junction1.split(':', 2)[2]}"
                    gene_id = original_junction1.split("-")[1].split(":")[0]
                    isoforms1 = transcript_dict.get(gene_id, {}).get(junction1, [])
                    if len(isoforms1) == 0:
                        isoforms1 = transcript_dict.get(gene_id, {}).get(junction1.split(":", 1)[1], [])
                    longest_isoform1,longest_exons1 = get_longest_isoform(gene_id, isoforms1, protein_dict)
                    gene_id = original_gene_id
                    junction1 = original_junction1

                except: # Errors can occur in re-assigning junction1 or other configurations of trans-splicing events (rare)
                    # In the future, we need to catch these and assign gene IDs as a prefix to the junction and re-assign the gene (may not be indicated in the feature)
                    pass
            else:
                longest_isoform1,longest_exons1 = get_longest_isoform(gene_id, isoforms1, protein_dict)

            if len(isoforms2)==0:
                #gene_id = "ENSG00000002016" > target_gene_id = "ENSG00000250132"
                #source = "ENSG00000002016:E1.1-ENSG00000250132:I3.1_988668" > target = "ENSG00000002016:E1.1-I3.1_988668"
                original_gene_id = str(gene_id)
                original_junction2 = str(junction2)
                try:
                    junction2 = f"{gene_id}:{junction2.split(':')[1].split('-')[0]}-{junction2.split(':', 2)[2]}"
                    gene_id = original_junction2.split("-")[1].split(":")[0]
                    isoforms2 = transcript_dict.get(gene_id, {}).get(junction2, [])
                    if len(isoforms2) == 0:
                        isoforms2 = transcript_dict.get(gene_id, {}).get(junction2.split(":", 1)[1], [])
                    longest_isoform2,longest_exons2 = get_longest_isoform(gene_id, isoforms2, protein_dict)
                    gene_id = original_gene_id
                    junction2 = original_junction2
                except:
                    # In the future, we need to catch these and assign gene IDs as a prefix to the junction and re-assign the gene (may not be indicated in the feature)
                    pass # Errors can occur in re-assigning junction2 or other configurations of trans-splicing events (rare)
            else:
                longest_isoform2,longest_exons2 = get_longest_isoform(gene_id, isoforms2, protein_dict)

            # Extract coordinates
            coord1 = coord_dict.get(junction1, 'N/A')
            coord2 = coord_dict.get(junction2, 'N/A')

            try: 
                strand = "+" if int(coord1.split(":")[1].split("-")[1]) > int(coord1.split(":")[1].split("-")[0]) else "-"
                event_type, proximity = se.classify_splicing_event(feature,strand,longest_exons1,longest_exons2,coord1,coord2)
            except Exception as e:
                event_type, proximity = 'Trans-Splicing', 'null'
            
            event_isoform_db[feature] = gene_symbol, longest_isoform1, longest_isoform2, coord1, coord2, event_type, proximity

        # If dPSI is negative, invert the value of proximity 
        prox_dict = {"distal": "proximal", "proximal": "distal", "inclusion": "exclusion", "exclusion": "inclusion"}
        proximity = prox_dict.get(proximity, proximity) if mwuSign == -1 else proximity

        # Add new annotated row
        annotations.append([
            gene_symbol, longest_isoform1, longest_isoform2, coord1, coord2, event_type, proximity
        ])

    # Create DataFrame with annotations
    annotated_df = pd.DataFrame(
        annotations, 
        columns=['Gene_Symbol', 'Isoform_1|Length', 'Isoform_2|Length', 'Coord_1', 'Coord_2', 'Event_Type', 'Event_Direction']
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
    global event_isoform_db
    event_isoform_db={}
    output_folder = stats_folder

    # Get a list of all .txt files excluding '-annotated.txt' files
    stats_files = [f for f in glob.glob(f"{stats_folder}/*.txt") if '-annotated' not in f]

    # Use tqdm to show progress over the stats files
    for stats_file in tqdm(stats_files, desc="Annotating Stats Files"):
        # Check if the annotated version already exists
        annotated_file = os.path.join(output_folder, os.path.basename(stats_file).replace('.txt', '-annotated.txt'))
        if os.path.exists(annotated_file):
            print(f"Skipping {stats_file}, already annotated.")
            continue  # Skip this file

        if dataType == 'junction':
            annotate_junction_stats_file(stats_file, transcript_dict, gene_symbol_dict, coord_dict, protein_dict, output_folder)
        else:
            annotate_iso_stats_file(stats_file, transcript_dict, gene_symbol_dict, protein_dict, output_folder)

def annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType='junction'):
    if junction_coords_file == None:
        #isoform expression/ratio
        coord_dict = None
    else:
        import zipfile
        # Handle both text and zipped versions of the file
        if junction_coords_file.endswith(".zip"):
            with zipfile.ZipFile(junction_coords_file, "r") as z:
                file_names = z.namelist()  # List files in the ZIP
                assert len(file_names) == 1, "Expected exactly one file in the ZIP archive"
                with z.open(file_names[0], "r") as f:
                    lines = [line.decode("utf-8").strip() for line in f.readlines()]  # Decode bytes to str
            coord_dict = load_junction_coordinates(lines)  # Pass lines directly
        else:
            coord_dict = load_junction_coordinates(junction_coords_file)

    transcript_dict = load_transcript_associations(transcript_assoc_file,dataType=dataType)
    gene_symbol_dict = load_gene_symbols(gene_symbol_file)
    protein_dict = load_protein_summary(protein_summary_file)

    # Annotate all stats files
    annotate_all_stats_files(stats_folder, transcript_dict, gene_symbol_dict, coord_dict, protein_dict,dataType=dataType)

if __name__ == '__main__':
    # File paths
    transcript_assoc_file = 'gff-output/transcript_associations.txt'
    gene_symbol_file = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl-annotations.txt'
    junction_coords_file = 'junction_combined_pseudo_cluster_counts-filtered.txt'
    protein_summary_file = 'protein_summary.txt'
    stats_folder = 'dPSI-covariate-events'

    annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType='junction')
