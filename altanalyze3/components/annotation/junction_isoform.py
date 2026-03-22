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
    best_rank = None
    for (isoform,exons) in isoforms:
        length, nmd_status, _ = protein_dict.get((gene,isoform), (0, None, None))
        # Prefer non-NMD transcripts when available, then maximize protein length.
        is_nmd = _is_potential_nmd_status(nmd_status)
        rank = (is_nmd, -int(length), str(isoform))
        if best_rank is None or rank < best_rank:
            best_rank = rank
            nmd_suffix = "(NMD)" if is_nmd else ""
            longest_isoform = f"{isoform}|{length}{nmd_suffix}"
            longest_exons = exons
    return (longest_isoform,longest_exons)


def _is_ensembl_isoform(isoform_id):
    return str(isoform_id).startswith("ENST")


def _get_protein_length(gene, isoform_id, protein_dict):
    length, _, _ = protein_dict.get((gene, isoform_id), (0, None, None))
    try:
        return int(length)
    except Exception:
        return 0


def _is_potential_nmd_status(nmd_status):
    status = str(nmd_status or "").strip().lower()
    if not status:
        return False
    # Treat explicit non-NMD labels as non-NMD.
    if "not" in status and "nmd" in status:
        return False
    return "nmd" in status


def _is_potential_nmd(gene, isoform_id, protein_dict):
    _, nmd_status, _ = protein_dict.get((gene, isoform_id), (0, None, None))
    return _is_potential_nmd_status(nmd_status)


def _prefer_non_nmd(gene, candidates, protein_dict):
    if not candidates:
        return []
    non_nmd = [(iso, exons) for iso, exons in candidates if not _is_potential_nmd(gene, iso, protein_dict)]
    return non_nmd if non_nmd else candidates


def _format_isoform_with_length(gene, isoform_id, exons, protein_dict):
    length = _get_protein_length(gene, isoform_id, protein_dict)
    nmd_suffix = "(NMD)" if _is_potential_nmd(gene, isoform_id, protein_dict) else ""
    return f"{isoform_id}|{length}{nmd_suffix}", exons


def _candidate_pool_prioritized(isoforms):
    """Always prioritize Ensembl isoforms when present for a junction."""
    if not isoforms:
        return []
    ensembl = [(iso, exons) for iso, exons in isoforms if _is_ensembl_isoform(iso)]
    return ensembl if ensembl else isoforms


def _exon_jaccard_distance(exons_a, exons_b):
    set_a = set(exons_a or [])
    set_b = set(exons_b or [])
    if not set_a and not set_b:
        return 0.0
    union = set_a | set_b
    if not union:
        return 0.0
    inter = set_a & set_b
    return 1.0 - (len(inter) / float(len(union)))


def _top_k_by_protein_length(gene, isoforms, protein_dict, top_k=12):
    if not isoforms:
        return []
    return sorted(
        isoforms,
        key=lambda pair: (
            _is_potential_nmd(gene, pair[0], protein_dict),
            -_get_protein_length(gene, pair[0], protein_dict),
            str(pair[0]),
        ),
    )[:top_k]


def _pick_best_single(gene, candidates, protein_dict, anchor_exons=None):
    if not candidates:
        return "N/A", "N/A"
    if len(candidates) == 1:
        isoform_id, exons = candidates[0]
        return _format_isoform_with_length(gene, isoform_id, exons, protein_dict)

    candidates = _prefer_non_nmd(gene, candidates, protein_dict)

    # Length-first selection by design. Only use exon composition as a tie-breaker.
    lengths = {
        isoform_id: _get_protein_length(gene, isoform_id, protein_dict)
        for isoform_id, _ in candidates
    }
    max_len = max(lengths.values()) if lengths else 0
    longest = [(iso, exons) for iso, exons in candidates if lengths.get(iso, 0) == max_len]

    if anchor_exons is None:
        best_iso, best_exons = sorted(longest, key=lambda x: str(x[0]))[0]
        return _format_isoform_with_length(gene, best_iso, best_exons, protein_dict)

    if len(longest) == 1:
        best_iso, best_exons = longest[0]
        return _format_isoform_with_length(gene, best_iso, best_exons, protein_dict)

    best = None
    best_rank = None
    for isoform_id, exons in longest:
        length = lengths.get(isoform_id, 0)
        dist = _exon_jaccard_distance(exons, anchor_exons)
        rank = (dist, -length, str(isoform_id))
        if best_rank is None or rank < best_rank:
            best_rank = rank
            best = (isoform_id, exons)

    if best is None:
        return "N/A", "N/A"
    return _format_isoform_with_length(gene, best[0], best[1], protein_dict)


def get_prioritized_isoform_pair(gene, isoforms1, isoforms2, protein_dict, top_k=12):
    """Select representative isoforms for competing junctions efficiently.

    Priority:
    1) Ensembl isoforms when present for each junction.
    2) If only one Ensembl isoform exists for a junction, use it directly.
    3) Otherwise, length is primary; exon composition distance is only a tie-breaker
       among equally longest candidates.
    """
    pool1 = _candidate_pool_prioritized(isoforms1)
    pool2 = _candidate_pool_prioritized(isoforms2)

    if not pool1 and not pool2:
        return "N/A", "N/A", "N/A", "N/A"
    if not pool1:
        iso2, ex2 = _pick_best_single(gene, _top_k_by_protein_length(gene, pool2, protein_dict, top_k), protein_dict)
        return "N/A", "N/A", iso2, ex2
    if not pool2:
        iso1, ex1 = _pick_best_single(gene, _top_k_by_protein_length(gene, pool1, protein_dict, top_k), protein_dict)
        return iso1, ex1, "N/A", "N/A"

    # Fast path requested: if one Ensembl candidate exists for a junction, don't over-scan that side.
    if len(pool1) == 1:
        iso1, ex1 = _format_isoform_with_length(gene, pool1[0][0], pool1[0][1], protein_dict)
        iso2, ex2 = _pick_best_single(
            gene,
            _top_k_by_protein_length(gene, pool2, protein_dict, top_k),
            protein_dict,
            anchor_exons=pool1[0][1],
        )
        return iso1, ex1, iso2, ex2
    if len(pool2) == 1:
        iso2, ex2 = _format_isoform_with_length(gene, pool2[0][0], pool2[0][1], protein_dict)
        iso1, ex1 = _pick_best_single(
            gene,
            _top_k_by_protein_length(gene, pool1, protein_dict, top_k),
            protein_dict,
            anchor_exons=pool2[0][1],
        )
        return iso1, ex1, iso2, ex2

    cand1 = _top_k_by_protein_length(gene, pool1, protein_dict, top_k=top_k)
    cand2 = _top_k_by_protein_length(gene, pool2, protein_dict, top_k=top_k)

    if not cand1 and not cand2:
        return "N/A", "N/A", "N/A", "N/A"
    if not cand1:
        iso2, ex2 = _pick_best_single(gene, cand2, protein_dict)
        return "N/A", "N/A", iso2, ex2
    if not cand2:
        iso1, ex1 = _pick_best_single(gene, cand1, protein_dict)
        return iso1, ex1, "N/A", "N/A"

    cand1 = _prefer_non_nmd(gene, cand1, protein_dict)
    cand2 = _prefer_non_nmd(gene, cand2, protein_dict)

    lengths1 = {iso: _get_protein_length(gene, iso, protein_dict) for iso, _ in cand1}
    lengths2 = {iso: _get_protein_length(gene, iso, protein_dict) for iso, _ in cand2}
    max_len1 = max(lengths1.values()) if lengths1 else 0
    max_len2 = max(lengths2.values()) if lengths2 else 0
    top1 = [(iso, exons) for iso, exons in cand1 if lengths1.get(iso, 0) == max_len1]
    top2 = [(iso, exons) for iso, exons in cand2 if lengths2.get(iso, 0) == max_len2]

    best_pair = None
    best_rank = None
    for iso1, ex1 in top1:
        len1 = lengths1.get(iso1, 0)
        for iso2, ex2 in top2:
            len2 = lengths2.get(iso2, 0)
            dist = _exon_jaccard_distance(ex1, ex2)
            rank = (dist, -(len1 + len2), -len1, -len2, str(iso1), str(iso2))
            if best_rank is None or rank < best_rank:
                best_rank = rank
                best_pair = ((iso1, ex1), (iso2, ex2))

    if best_pair is None:
        iso1, ex1 = _pick_best_single(gene, cand1, protein_dict)
        iso2, ex2 = _pick_best_single(gene, cand2, protein_dict)
        return iso1, ex1, iso2, ex2

    (iso1_id, ex1), (iso2_id, ex2) = best_pair
    iso1, exons1 = _format_isoform_with_length(gene, iso1_id, ex1, protein_dict)
    iso2, exons2 = _format_isoform_with_length(gene, iso2_id, ex2, protein_dict)
    return iso1, exons1, iso2, exons2

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
                    gene_id = original_gene_id
                    junction1 = original_junction1

                except: # Errors can occur in re-assigning junction1 or other configurations of trans-splicing events (rare)
                    # In the future, we need to catch these and assign gene IDs as a prefix to the junction and re-assign the gene (may not be indicated in the feature)
                    pass
            else:
                pass

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
                    gene_id = original_gene_id
                    junction2 = original_junction2
                except:
                    # In the future, we need to catch these and assign gene IDs as a prefix to the junction and re-assign the gene (may not be indicated in the feature)
                    pass # Errors can occur in re-assigning junction2 or other configurations of trans-splicing events (rare)
            else:
                pass

            longest_isoform1,longest_exons1,longest_isoform2,longest_exons2 = get_prioritized_isoform_pair(
                gene_id, isoforms1, isoforms2, protein_dict, top_k=12
            )

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
    stats_files = [f for f in glob.glob(f"{stats_folder}/*.txt") if f.endswith("stats.txt")]

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
