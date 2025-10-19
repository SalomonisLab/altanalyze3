import pysam
from collections import defaultdict, Counter
import os
import numpy as np
from scipy import stats
import pandas as pd
import argparse

"""
SNV or Indel mode to specify mutation or indels calling logic
SNV mode: only examine precise variant genomic location. R (reference) or V (variant) calls. 
Indel mode: check indel within +/- 50bp of given location. Only 0/+1 indels to the location get called by I (indel) otherwise R
Usage: python3 variant_extraction.py --bam scisoseq.mapped.bam --sample C-05-1893 --mutations C-05-1893-variants.txt --reference genome.fa --output variants
Mutation file format: 
chr17	31159086	31159086	NF1c.281T>G	SNV
chr17	31163273	31163273	NF1c.377delA	Indel
"""

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Supervised variant/indel to cell barcode reporting')
    parser.add_argument('--sample', '--s', required=True, help='Sample name (e.g., Y1493)')
    parser.add_argument('--bam', '--b', required=True, help='Path to BAM file')
    parser.add_argument('--mutations', '--m', required=True, help='Path to mutations file')
    parser.add_argument('--reference', '--r', required=True, help='Path to reference genome FASTA')
    parser.add_argument('--output-dir', '--o', required=True, help='Output directory')
    return parser.parse_args()

def parse_mutation_file(mutation_file):
    """Parse mutation file"""
    mutations = []
    with open(mutation_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                chrom = parts[0] if parts[0].startswith('chr') else f'chr{parts[0]}'
                pos = int(parts[1])
                ref = parts[2] if len(parts) > 2 else 'N'
                alt = parts[3] if len(parts) > 3 else 'N'
                mode = parts[4] if len(parts) > 4 else 'SNV'
                gene = parts[5] if len(parts) > 5 else 'NA'
                
                is_indel = mode.lower() == 'indel'
                mutations.append((chrom, pos, ref, alt, gene, mode, is_indel))
    return mutations

def detect_indels_near_position(read, target_pos, window=50):
    """Detect indels within window of target position"""
    indels = []
    current_ref_pos = read.reference_start
    
    for op, length in read.cigartuples:
        if op == 0 or op == 7 or op == 8:  # Match/Mismatch (M/=/X)
            current_ref_pos += length
        elif op == 1:  # Insertion (I) - doesn't advance reference position
            if abs(current_ref_pos - target_pos) <= window:
                indels.append({
                    'type': 'insertion',
                    'pos': current_ref_pos,
                    'length': length,
                    'distance': abs(current_ref_pos - target_pos)
                })
        elif op == 2:  # Deletion (D)
            del_start = current_ref_pos
            del_end = current_ref_pos + length
            distance_to_target = min(abs(del_start - target_pos), abs(del_end - target_pos))
            
            if (distance_to_target <= window or (del_start <= target_pos <= del_end)):
                indels.append({
                    'type': 'deletion',
                    'pos': del_start,
                    'length': length,
                    'distance': distance_to_target
                })
            
            current_ref_pos += length
        elif op == 3:  # Skip/Intron (N)
            current_ref_pos += length
        elif op == 4:  # Soft clip (S) - doesn't advance reference position
            pass
        elif op == 5:  # Hard clip (H) - doesn't advance reference position
            pass
    
    return indels

def get_homopolymer_length(sequence, position):
    """Get homopolymer length around position"""
    if position >= len(sequence) or position < 0:
        return 0
    
    base = sequence[position]
    left = sum(1 for i in range(position-1, -1, -1) if sequence[i] == base)
    right = sum(1 for i in range(position+1, len(sequence)) if sequence[i] == base)
    return left + right + 1

def check_tandem_repeat(sequence):
    """Check for tandem repeats"""
    for unit_size in range(2, min(7, len(sequence)//2 + 1)):
        for start in range(len(sequence) - unit_size + 1):
            unit = sequence[start:start + unit_size]
            repeats = 1
            pos = start + unit_size
            while pos + unit_size <= len(sequence) and sequence[pos:pos + unit_size] == unit:
                repeats += 1
                pos += unit_size
            if repeats >= 3:
                return True
    return False

def calculate_pcr_flags(read_details, target_base):
    """Calculate PCR artifact flags"""
    if not read_details:
        return {'strand_bias_pvalue': 1.0, 'position_bias_score': 0.0, 
                'mean_base_quality': 0.0, 'quality_std': 0.0}
    
    # Strand bias
    forward = sum(1 for r in read_details if not r.get('is_reverse', False))
    reverse = len(read_details) - forward
    strand_bias_pvalue = 1.0
    
    if forward > 0 and reverse > 0:
        try:
            from scipy.stats import fisher_exact
            _, strand_bias_pvalue = fisher_exact([[forward, reverse], [forward, reverse]])
        except:
            pass
    
    # Position bias
    positions = [r.get('rel_position', 0) for r in read_details]
    end_positions = [p for p in positions if p < 0.2 or p > 0.8]
    position_bias_score = len(end_positions) / len(positions) if positions else 0
    
    # Quality stats
    qualities = [r.get('base_quality', 0) for r in read_details if r.get('base_quality', 0) > 0]
    mean_quality = np.mean(qualities) if qualities else 0
    quality_std = np.std(qualities) if len(qualities) > 1 else 0
    
    return {
        'strand_bias_pvalue': strand_bias_pvalue,
        'position_bias_score': position_bias_score,
        'mean_base_quality': mean_quality,
        'quality_std': quality_std
    }

def analyze_comprehensive(bam_file, mutations_list, ref_genome, context_bases=5, min_qual=20):
    """Comprehensive mutation analysis"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    ref = pysam.FastaFile(ref_genome)
    results = []
    
    for chrom, pos, ref_base, alt, gene, variant_id, is_indel in mutations_list:
        pos_0based = pos - 1
        
        # Get reference context
        ref_start = max(0, pos_0based - context_bases)
        ref_end = pos_0based + context_bases + 1
        ref_context = ref.fetch(chrom, ref_start, ref_end).upper()
        ref_single_base = ref.fetch(chrom, pos_0based, pos_0based + 1).upper()
        
        print(f"Processing {chrom}:{pos} ({'indel' if is_indel else 'SNV'} mode)")
        
        # Set fetch region based on mode
        if is_indel:
            fetch_start = pos_0based - 100
            fetch_end = pos_0based + 100
        else:
            fetch_start = pos_0based
            fetch_end = pos_0based + 1
        
        cell_observations = defaultdict(list)
        all_read_details = []
        
        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or read.mapping_quality < min_qual:
                continue
            
            # Get cell barcode
            barcode = None
            for tag in ["CB", "BC", "XC", "UB"]: # possible cell barcode tags
                try:
                    barcode = read.get_tag(tag)
                    break
                except KeyError:
                    continue
            if not barcode:
                continue
            
            # For indel mode: get sequence at reference position if possible
            if is_indel:
                read_positions = read.get_reference_positions(full_length=True)
                if pos_0based in read_positions:
                    query_idx = read_positions.index(pos_0based)
                    observed_base = read.query_sequence[query_idx]
                else:
                    continue
            else:
                # For SNV mode: need exact position overlap
                read_positions = read.get_reference_positions(full_length=True)
                if pos_0based not in read_positions:
                    continue
                query_idx = read_positions.index(pos_0based)
                observed_base = read.query_sequence[query_idx]
            
            # Detect indels within 100bp of target position
            indels = detect_indels_near_position(read, pos_0based, window=100)
            
            # Get base quality
            base_qual = 0
            if read.query_qualities and query_idx < len(read.query_qualities):
                base_qual = read.query_qualities[query_idx]
            
            # Context sequence
            start = max(0, query_idx - context_bases)
            end = min(len(read.query_sequence), query_idx + context_bases + 1)
            obs_seq = read.query_sequence[start:end]
            var_pos_in_context = query_idx - start
            marked_seq = (obs_seq[:var_pos_in_context] + f"[{obs_seq[var_pos_in_context]}]" + 
                         obs_seq[var_pos_in_context+1:])
            
            # Pad sequence
            expected_length = (2 * context_bases) + 1
            if len(obs_seq) < expected_length:
                if start == 0:
                    marked_seq = 'N'*(expected_length - len(obs_seq)) + marked_seq
                else:
                    marked_seq = marked_seq + 'N'*(expected_length - len(obs_seq))
            
            # Additional metrics
            homopolymer_length = get_homopolymer_length(read.query_sequence, query_idx)
            in_tandem_repeat = check_tandem_repeat(obs_seq)
            
            # Frameshift analysis
            has_frameshift = any(indel['length'] % 3 != 0 for indel in indels)
            
            # Store read details
            read_info = {
                'observed_base': observed_base,
                'base_quality': base_qual,
                'is_reverse': read.is_reverse,
                'rel_position': query_idx / len(read.query_sequence) if len(read.query_sequence) > 0 else 0
            }
            all_read_details.append(read_info)
            
            # Store cell observation
            cell_observations[barcode].append({
                'observed_base': observed_base,
                'marked_seq': marked_seq,
                'base_quality': base_qual,
                'mapping_quality': read.mapping_quality,
                'indels': indels,
                'homopolymer_length': homopolymer_length,
                'in_tandem_repeat': in_tandem_repeat,
                'is_reverse': read.is_reverse,
                'template_length': abs(read.template_length) if read.template_length else 0,
                'has_frameshift': has_frameshift
            })
        
        # Calculate PCR flags
        pcr_flags = calculate_pcr_flags(all_read_details, alt)
        
        print(f"Found {len(cell_observations)} cells, {len(all_read_details)} total reads")
        
        # Process each cell
        for barcode, observations in cell_observations.items():
            for obs in observations:
                # Create summaries
                indel_summary = ';'.join([
                    f"{i['type']}@{i['pos']}({i['length']}bp,dist={i['distance']})" 
                    for i in obs['indels']
                ]) if obs['indels'] else 'None'
                
                # Determine allele status based on mutation type
                if is_indel:
                    # For indel positions, only consider indels very close to target (0-1bp)
                    close_indels = [i for i in obs['indels'] if i['distance'] <= 1]
                    
                    if close_indels:
                        allele_status = "I"  # Indel detected at target position
                    else:
                        allele_status = "R"  # Reference (no relevant indel)
                else:
                    # For SNV positions, compare bases
                    allele_status = "V" if obs['observed_base'] != ref_single_base else "R"
                
                # Context changes (simplified)
                context_changes = 0
                
                results.append([
                    barcode, chrom, str(pos), ref_base, alt, gene, variant_id,
                    ref_context, ref_single_base, obs['observed_base'], obs['marked_seq'],
                    allele_status, str(obs['base_quality']), str(obs['mapping_quality']),
                    str(len(obs['indels'])), indel_summary,
                    "Y" if obs['has_frameshift'] else "N", indel_summary,
                    str(context_changes), "None",
                    str(obs['homopolymer_length']), "Y" if obs['in_tandem_repeat'] else "N",
                    "Y" if obs['is_reverse'] else "N", str(obs['template_length']),
                    f"{pcr_flags['strand_bias_pvalue']:.4f}",
                    f"{pcr_flags['position_bias_score']:.3f}",
                    f"{pcr_flags['mean_base_quality']:.1f}",
                    f"{pcr_flags['quality_std']:.1f}",
                    "0.0", "0.0", "1.0000",  # Placeholder PCR flags
                    str(len(all_read_details)), str(len([r for r in all_read_details if r['observed_base'] == alt]))
                ])
    
    bam.close()
    ref.close()
    return results

def reverse_complement(bc):
    """Reverse complement cell barcode"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement.get(base, base) for base in bc[::-1]])

def determine_mutation_status(row):
    """Determine if observation is WT or MUT based on allele_status"""
    if row['allele_status'] == 'V':  # Variant detected (SNV mode)
        return 'MUT'
    elif row['allele_status'] == 'I':  # Indel detected (indel mode)
        return 'MUT'
    elif row['allele_status'] == 'R':  # Reference/wildtype
        return 'WT'
    else:
        return 'UNK'

def genotype_cell(wt_count, mut_count):
    """Genotype cell based on WT and MUT counts"""
    if mut_count > 0:
        return 'MUT'
    elif wt_count > 0:
        return 'WT'
    else:
        return 'UNK'

def create_mutation_matrix(input_file, output_file, sample_name):
    """Create mutation matrix from TSV file"""
    
    print("Loading data...")
    df = pd.read_csv(input_file, sep="\t")
    print(f"Loaded {len(df)} observations from {len(df['cell_barcode'].unique())} cells")
    
    # Reverse complement cell barcodes
    print("Reverse complementing cell barcodes...")
    df['cell_barcode'] = df['cell_barcode'].apply(reverse_complement)
    
    # Determine mutation status
    print("Determining mutation status...")
    df['mut_status'] = df.apply(determine_mutation_status, axis=1)
    
    # Clean up alt column and handle NaN values
    df['alt'] = df['alt'].fillna('Unknown')
    df['alt'] = df['alt'].astype(str)
    
    # Get unique alterations for matrix (using alt column)
    unique_alterations = sorted(df['alt'].unique())
    
    print(f"Found {len(unique_alterations)} unique genetic alterations:")
    for alt in unique_alterations:
        print(f"  {alt}")
    
    # Get all unique cell barcodes
    cell_barcodes = sorted(df['cell_barcode'].unique())
    
    # Initialize results dataframe with placeholders
    results = []
    
    print("\nProcessing each cell...")
    for barcode in cell_barcodes:
        cell_data = {'cell_barcode': barcode, 'genotype': '', 'mutations_present': ''}
        
        for alteration in unique_alterations:
            # Filter for this cell and alteration
            cell_alt_data = df[(df['cell_barcode'] == barcode) & 
                              (df['alt'] == alteration)]
            
            # Count WT and MUT observations
            wt_count = len(cell_alt_data[cell_alt_data['mut_status'] == 'WT'])
            mut_count = len(cell_alt_data[cell_alt_data['mut_status'] == 'MUT'])
            
            # Add counts to cell data
            cell_data[f"{alteration}_WT"] = wt_count
            cell_data[f"{alteration}_MUT"] = mut_count
            
            # Genotype this cell for this alteration
            genotype = genotype_cell(wt_count, mut_count)
            cell_data[f"{alteration}_GENOTYPE"] = genotype
        
        results.append(cell_data)
    
    # Convert to dataframe
    matrix_df = pd.DataFrame(results)
    
    # Now fill in the placeholder columns
    print("Assigning final genotype for each cell...")
    
    def assign_final_genotype(row):
        """Assign final genotype based on all individual genotypes"""
        genotype_cols = [col for col in row.index if col.endswith('_GENOTYPE')]
        individual_genotypes = [row[col] for col in genotype_cols]
        
        # If any alteration is MUT, cell is MUT
        if 'MUT' in individual_genotypes:
            return 'MUT'
        # If any alteration is WT (and none are MUT), cell is WT
        elif 'WT' in individual_genotypes:
            return 'WT'
        # If all are UNK, cell is UNK
        else:
            return 'UNK'
    
    def get_mutations_present(row):
        """Get list of mutations present in this cell"""
        mutations = []
        for alteration in unique_alterations:
            if row[f"{alteration}_GENOTYPE"] == 'MUT':
                mutations.append(alteration)
        return ';'.join(mutations) if mutations else 'None'
    
    # Fill the placeholder columns
    matrix_df['genotype'] = matrix_df.apply(assign_final_genotype, axis=1)
    matrix_df['mutations_present'] = matrix_df.apply(get_mutations_present, axis=1)
    
    # Summary of final genotypes
    final_genotype_counts = matrix_df['genotype'].value_counts()
    print(f"\nFinal genotype summary:")
    for genotype, count in final_genotype_counts.items():
        print(f"  {genotype}: {count} cells ({count/len(matrix_df)*100:.1f}%)")
    
    print(f"\nMatrix summary:")
    print(f"Cells: {len(matrix_df)}")
    print(f"Alterations: {len(unique_alterations)}")
    print(f"Total columns: {len(matrix_df.columns)}")
    
    # Save matrix
    print(f"\nSaving matrix to {output_file}...")
    matrix_df.to_csv(output_file, index=False)
    
    # Create summary file
    summary_file = output_file.replace('.csv', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"Mutation Matrix Summary - {sample_name}\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total cells: {len(matrix_df)}\n")
        f.write(f"Total alterations: {len(unique_alterations)}\n\n")
        
        f.write("Alterations included:\n")
        for alt in unique_alterations:
            f.write(f"  {alt}\n")
        
        f.write(f"\nFinal genotype summary:\n")
        for genotype, count in final_genotype_counts.items():
            f.write(f"  {genotype}: {count} cells ({count/len(matrix_df)*100:.1f}%)\n")
        
        f.write("\nGenotype counts per alteration:\n")
        genotype_cols = [col for col in matrix_df.columns if col.endswith('_GENOTYPE')]
        for col in genotype_cols:
            alteration = col.replace('_GENOTYPE', '')
            genotype_counts = matrix_df[col].value_counts()
            f.write(f"\n{alteration}:\n")
            for genotype, count in genotype_counts.items():
                f.write(f"  {genotype}: {count} cells\n")
    
    print(f"Summary saved to {summary_file}")
    print("Done!")
    
    return matrix_df

def main():
    """Main pipeline function"""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output files based on sample name
    tsv_output = os.path.join(args.output_dir, f"{args.sample}_complete_analysis.tsv")
    matrix_output = os.path.join(args.output_dir, f"{args.sample}_mutation_matrix.csv")
    
    print(f"=== Mutation Analysis Pipeline for {args.sample} ===\n")
    
    # Step 1: Comprehensive mutation analysis
    print("STEP 1: Comprehensive mutation analysis")
    print("Loading mutations...")
    mutations = parse_mutation_file(args.mutations)
    print(f"Loaded {len(mutations)} mutation positions")
    
    print("Processing BAM file...")
    header = [
        "cell_barcode", "chr", "position", "ref", "alt", "gene", "variant_id",
        "reference_sequence", "reference_base", "observed_base", "observed_sequence_marked",
        "allele_status", "base_quality", "mapping_quality",
        "num_indels", "indel_details", "frameshift_detected", "frameshift_details",
        "num_context_changes", "context_changes",
        "homopolymer_length", "in_tandem_repeat", "is_reverse_strand", "template_length",
        "pcr_strand_bias_pvalue", "pcr_position_bias_score", "pcr_mean_base_quality", 
        "pcr_quality_std", "pcr_homopolymer_fraction", "pcr_tandem_repeat_fraction",
        "pcr_template_length_bias_pvalue", "position_total_reads", "position_alt_reads"
    ]
    
    results = analyze_comprehensive(args.bam, mutations, args.reference)
    
    print(f"Found {len(results)} sequence observations")
    
    # Write TSV output
    with open(tsv_output, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in results:
            f.write("\t".join(row) + "\n")
    
    print(f"TSV results written to: {tsv_output}")
    
    # Step 2: Create mutation matrix
    print(f"\nSTEP 2: Creating mutation matrix")
    matrix = create_mutation_matrix(tsv_output, matrix_output, args.sample)
    
    print(f"\n=== Pipeline Complete for {args.sample} ===")
    print(f"TSV output: {tsv_output}")
    print(f"Matrix output: {matrix_output}")
    print(f"Summary: {matrix_output.replace('.csv', '_summary.txt')}")

if __name__ == "__main__":
    main()