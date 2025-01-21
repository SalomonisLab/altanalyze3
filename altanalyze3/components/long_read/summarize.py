import pandas as pd
import os
import re

def determine_splicing_type(exon1, exon2, coord1, coord2):
    """
    Identify splicing event type based on exon patterns and coordinates.
    """
    # Alternative promoter if one of the exons starts with 'U1' or 'E1'
    if exon1.startswith("U1") or exon1.startswith("E1") or exon2.startswith("U1") or exon2.startswith("E1"):
        return "Alternative Promoter"
    
    # Alternative polyadenylation: last exon with different sites marked by '_'
    if "_" in exon1 and "_" in exon2 and exon1.split("_")[0] == exon2.split("_")[0]:
        return "Alternative Polyadenylation"
    
    # Trans-splicing if exons belong to different genes
    gene1, gene2 = exon1.split(":")[0], exon2.split(":")[0]
    if gene1 != gene2:
        return "Trans-Splicing"
    
    # Intron retention: must have 'I' in one of the exons without any '_' in either junction,
    # and coordinates should differ by exactly 1 nt
    if ("I" in exon1 or "I" in exon2) and "_" not in exon1 and "_" not in exon2 and abs(int(coord1.split("-")[1]) - int(coord2.split("-")[0])) == 1:
        return "Intron Retention"
    
    # Novel cassette exon: requires 'I' and '_', indicating a novel exon within an intron
    if ("I" in exon1 or "I" in exon2) and ("_" in exon1 or "_" in exon2):
        return "Novel Cassette Exon"
    
    # Alternative 5' splice site if start exon blocks differ by region
    exon1_blocks, exon2_blocks = exon1.split(":")[-1].split("-"), exon2.split(":")[-1].split("-")
    if exon1_blocks[0][:-2] == exon2_blocks[0][:-2] and exon1_blocks[0] != exon2_blocks[0]:
        return "Alternative 5' Splice Site"
    
    # Alternative 3' splice site if end exon blocks differ by region
    if exon1_blocks[-1][:-2] == exon2_blocks[-1][:-2] and exon1_blocks[-1] != exon2_blocks[-1]:
        return "Alternative 3' Splice Site"
    
    # Cassette exon if a middle exon block appears in one isoform but not the other
    if len(exon1_blocks) != len(exon2_blocks):
        return "Cassette Exon"
    
    return None

def determine_inclusion_exclusion(exon1_blocks, exon2_blocks, diff_means):
    """
    Determine inclusion or exclusion based on exon block proximity.
    """
    # Check for cassette exon events and analyze based on proximity
    if len(exon1_blocks) > len(exon2_blocks):  # Exon in exon1 but not in exon2
        return "Exclusion" if diff_means > 0 else "Inclusion"
    elif len(exon2_blocks) > len(exon1_blocks):  # Exon in exon2 but not in exon1
        return "Inclusion" if diff_means > 0 else "Exclusion"
    
    return "Undetermined"

def filter_and_summarize_differentials(directory_path):
    results = []
    significant_results_dir = os.path.join(directory_path, "significant-results")

    # Create the directory for significant results if it doesn't exist
    if not os.path.exists(significant_results_dir):
        os.makedirs(significant_results_dir)

    for filename in os.listdir(directory_path):
        if filename.endswith("stats-annotated.txt"):
            filepath = os.path.join(directory_path, filename)
            data = pd.read_csv(filepath, sep="\t")

            # Dynamically identify N_* columns
            count_columns = [col for col in data.columns if col.startswith("N_")]
            if len(count_columns) < 2:
                print(f"Skipping {filename}, insufficient N_* columns found.")
                continue
            
            n_col1, n_col2 = count_columns[0], count_columns[1]

            # Filter significant splicing events based on conditions
            filtered_data = data[
                (data['mwuPval'] < 0.05) &
                (data['DiffMeans'].abs() > 0.1) &
                (data[n_col1] > 2) &
                (data[n_col2] > 2)
            ]

            # Count total significant events
            total_significant_events = filtered_data.shape[0]

            # Skip files with no significant events
            if total_significant_events == 0:
                continue

            # Count NMD changes
            nmd_increase = 0
            nmd_decrease = 0

            # Analyze splicing event types and inclusion/exclusion
            event_types = {
                "Alternative Promoter": 0,
                "Alternative Polyadenylation": 0,
                "Trans-Splicing": 0,
                "Intron Retention": 0,
                "Novel Cassette Exon": 0,
                "Alternative 5' Splice Site": 0,
                "Alternative 3' Splice Site": 0,
                "Cassette Exon Exclusion": 0,
                "Cassette Exon Inclusion": 0
            }

            # Annotate each row in filtered_data with the splicing type and inclusion/exclusion
            filtered_data['Splicing_Type'] = ""
            filtered_data['Inclusion_Exclusion'] = ""

            for idx, row in filtered_data.iterrows():
                exon1, exon2 = row['Isoform_1|Length'].split("|")[0], row['Isoform_2|Length'].split("|")[0]
                coord1, coord2 = row['Coord_1'], row['Coord_2']
                
                # Determine splicing event type
                event_type = determine_splicing_type(exon1, exon2, coord1, coord2)
                filtered_data.at[idx, 'Splicing_Type'] = event_type

                if event_type == "Cassette Exon":
                    # Determine inclusion or exclusion for cassette exon events
                    exon1_blocks, exon2_blocks = exon1.split(":")[-1].split("-"), exon2.split(":")[-1].split("-")
                    inclusion_exclusion = determine_inclusion_exclusion(exon1_blocks, exon2_blocks, row['DiffMeans'])
                    filtered_data.at[idx, 'Inclusion_Exclusion'] = inclusion_exclusion

                    if inclusion_exclusion == "Inclusion":
                        event_types["Cassette Exon Inclusion"] += 1
                    elif inclusion_exclusion == "Exclusion":
                        event_types["Cassette Exon Exclusion"] += 1
                elif event_type:
                    event_types[event_type] += 1
                
                # Count NMD changes based on NMD potential in isoforms
                nmd_in_iso1 = "NMD-Potential" in row['Isoform_1|Length'] and row['DiffMeans'] > 0
                nmd_in_iso2 = "NMD-Potential" in row['Isoform_2|Length'] and row['DiffMeans'] < 0

                if nmd_in_iso1 or nmd_in_iso2:
                    nmd_increase += nmd_in_iso1
                    nmd_decrease += nmd_in_iso2

            # Record summary for this file
            result = {
                "Filename": filename,
                "Total Significant Events": total_significant_events,
                "NMD Increase": nmd_increase,
                "NMD Decrease": nmd_decrease,
                "Event Types": event_types
            }
            results.append(result)

            # Export the significant results with annotations to the significant results directory
            output_file = os.path.join(significant_results_dir, f"{filename.replace('.txt', '_significant_results.txt')}")
            filtered_data.to_csv(output_file, sep="\t", index=False)

    # Export the summary results for all files to a single file
    summary_df = pd.DataFrame(results)
    summary_output_file = os.path.join(directory_path, "all_files_summary_results.txt")
    summary_df.to_csv(summary_output_file, sep="\t", index=False)
    print(f"Summary results saved to: {summary_output_file}")
    print(f"Individual significant results saved in directory: {significant_results_dir}")

directory_path = "/path/to/your/directory"
results = filter_and_summarize_differentials(directory_path)
print(results)
