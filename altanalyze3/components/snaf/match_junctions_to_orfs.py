import os,sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

# File paths - replace with a main function call to pass these variables when connected with upstream functions
junctions_file = "/data/salomonis-archive/LabFiles/Nathan/Revio/BoneMarrow/SNAF-B_MDS-AML-neojunctions.txt"
genome_fasta = "/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/genome.fa"
transcript_associations_file = "/data/salomonis-archive/LabFiles/Nathan/Revio/Pancancer1/gff-output/transcript_associations.txt"
protein_summary_file = "/data/salomonis-archive/LabFiles/Nathan/Revio/Pancancer1/gff-output/protein_summary.txt"
orf_sequences_file = "/data/salomonis-archive/LabFiles/Nathan/Revio/Pancancer1/gff-output/orf_sequences.fasta"
output_file = "/data/salomonis-archive/LabFiles/Nathan/Revio/Pancancer1/gff-output/SNAF-B_functional_annotations-broad.txt"
uniprot_coordinates = '/data/salomonis2/software/AltAnalyze-100/AltDatabase/EnsMart100/uniprot/Hs/Hs_FeatureCoordinate.txt'
protein_translation_dir = '/data/salomonis2/software/AltAnalyze-100/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_Protein__100_38.tab'
gene_symbol_file = '/data/salomonis2/software/AltAnalyze-100/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl-annotations.txt'
uniprot_seq_file = '/data/salomonis2/software/AltAnalyze-100/AltDatabase/EnsMart100/uniprot/Hs/uniprot_sequence.txt'

# Helper functions
def load_gene_symbols(file_path):
    gene_symbol_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:  # Ensure there are at least two columns
                ensembl_id, gene_symbol = parts[:2]
                gene_symbol_dict[ensembl_id] = gene_symbol
    return gene_symbol_dict

def load_reference_isoform_seq(file_path):
    uniprot_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            protein_seq = parts[2]
            uids = parts[4].split(',')
            gene = next((item for item in uids if item.startswith('ENSG')), None)
            uniprot_dict[gene] = protein_seq
    return uniprot_dict

def load_domain_locations(translation_path,uniprot_coordinates_path,uniprot_dict):
    protein_translation_db={}
    with open(translation_path, 'r') as f:
        for line in f:
            values = line.strip().split('\t')
            if len(values)>2:
                gene,transcript,protein = values
                protein_translation_db[protein] = gene

    extracellular_domain_db={}
    transmembrane_domain_db={}
    transmembrane_domain_seq_db={}
    with open(uniprot_coordinates_path, 'r') as f:
        for line in f:
            ensembl_protein,aa_start,aa_stop,start,stop,name,interpro_id,domain = line.strip().split('\t')
            if ensembl_protein in protein_translation_db:
                gene = protein_translation_db[ensembl_protein]
                if 'TOPO_DOM-Extracell' in domain:
                    try: extracellular_domain_db[gene].append([int(start),int(stop)])
                    except: extracellular_domain_db[gene] = [[int(start),int(stop)]]
                if 'TRANSMEM-Helica' in domain:
                    try: transmembrane_domain_db[gene].append([int(start),int(stop)])
                    except: transmembrane_domain_db[gene] = [[int(start),int(stop)]]
                    if gene in uniprot_dict:
                        prot_sequence = uniprot_dict[gene]
                        tm_seq = prot_sequence[int(aa_start):int(aa_stop)]
                    try: transmembrane_domain_seq_db[gene].append(tm_seq)
                    except: transmembrane_domain_seq_db[gene] = [tm_seq]

    #### Need to add code to import uniprot protein sequences for transmembranes for a TM sequence check in the predicted isoforms
    return extracellular_domain_db, transmembrane_domain_db, transmembrane_domain_seq_db

def parse_junction(junction_line):
    """Parse a junction line into components."""
    junction_name, coordinates = junction_line.split("=")
    gene = junction_name.split(":")[0]
    chrom, coord_range = coordinates.split(":")
    start, end = map(int, coord_range.split("-"))
    return gene, junction_name, chrom, start, end

def reverse_complement(seq):
    """Return the reverse complement of a sequence."""
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def load_protein_summary(file_path):
    protein_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene_id, isoform, length, nmd_status, _, longest_length = line.strip().split('\t')
            if '_' in isoform:
                isoform = isoform.split('_')[1]  # Remove dataset prefix
            if ';' in gene_id:
                gene_id = gene_id.split(';')[0]
            try: length = int(length)
            except:
                continue
            protein_dict[gene_id,isoform] = f"{length}/{longest_length}({nmd_status})"
    return protein_dict

def generate_forward_frames(mrna_sequence,protein_sequence):
    """Generate the 3 forward translation frames for a given mRNA sequence."""
    # Ensure the sequence is uppercase and a valid Seq object
    mrna_seq = Seq(mrna_sequence.upper())
    hit = []
    for offset in range(3):  # Frame offsets: 0, 1, 2
        translated_seq = mrna_seq[offset:].translate(to_stop=False)  # Translate without stopping at stop codon
        if translated_seq in protein_sequence:
            hit.append(translated_seq)
    if len(hit)>0:
        return hit[0]
    else:
        return ""

def domain_architecture_evaluation(ensembl,junc_start,junct_stop):
    extracellular_evidence = 'FALSE'
    transmembrane_evidence = 'FALSE'
    if ensembl in extracellular_domain_db:
        for (ex_start1, ex_stop1) in extracellular_domain_db[ensembl]:
            temp = [ex_start1, ex_stop1]
            temp.sort()
            ex_start, ex_stop = temp
            ex_start1, ex_stop1 = str(ex_start1), str(ex_stop1)
            if junc_start >= ex_start and junc_start <= ex_stop:
                extracellular_evidence = 'TRUE'
            elif junct_stop >= ex_start and junct_stop <= ex_stop:
                extracellular_evidence = 'TRUE'
    if ensembl in transmembrane_domain_db:
        for (ex_start1, ex_stop1) in transmembrane_domain_db[ensembl]:
            temp = [ex_start1, ex_stop1]
            temp.sort()
            ex_start, ex_stop = temp
            ex_start1, ex_stop1 = str(ex_start1), str(ex_stop1)
            if junc_start >= ex_start and junc_start <= ex_stop:
                transmembrane_evidence = 'TRUE'
            elif junct_stop >= ex_start and junct_stop <= ex_stop:
                transmembrane_evidence = 'TRUE'
            #print extracellular_evidence
    return extracellular_evidence, transmembrane_evidence

def sliding_window(sequence, window_size):
    """Generate all subsequences of a given length (window_size) from the sequence."""
    # Properly format the function and ensure parentheses are closed
    return {sequence[i:i + window_size]: i for i in range(len(sequence) - window_size + 1)}

def refine_unique_sequences(seq1, seq2, unique_matches):
    """Merge overlapping unique sequences and refine boundaries."""
    refined_sequences = []
    sorted_matches = sorted(unique_matches.items(), key=lambda x: x[1])  # Sort by position

    merged_sequence = None
    for peptide, position in sorted_matches:
        if merged_sequence is None:
            merged_sequence = [peptide, position, position + len(peptide)]
        else:
            # Check for overlap
            if position <= merged_sequence[2]:  # Overlap detected
                overlap_start = max(merged_sequence[2] - position, 0)
                merged_sequence[0] += peptide[overlap_start:]  # Extend
                merged_sequence[2] = position + len(peptide)  # Update end position
            else:
                # Finalize the previous merged sequence
                refined_sequences.append(merged_sequence)
                merged_sequence = [peptide, position, position + len(peptide)]
    
    # Append the last sequence
    if merged_sequence:
        refined_sequences.append(merged_sequence)

    # Refine the boundaries by removing shared prefix/suffix with seq2
    final_sequences = []
    for merged_peptide, start, end in refined_sequences:
        # Adjust start
        while start < len(seq1) and start < len(seq2) and seq1[start] == seq2[start]:
            start += 1
        # Adjust end
        while end > start and end <= len(seq1) and end <= len(seq2) and seq1[end - 1] == seq2[end - 1]:
            end -= 1
        
        # Extract refined sequence and include flanking
        refined_peptide = seq1[start:end]
        flanking_region = seq1[max(0, start - 10):min(len(seq1), end + 10)]
        final_sequences.append((refined_peptide, flanking_region))

    return final_sequences


def find_unique_sequences_with_refinement(seq1, seq2, window_size=5):
    """Find and refine unique subsequences in seq1 compared to seq2."""
    peptides1 = sliding_window(seq1, window_size)
    peptides2 = sliding_window(seq2, window_size)

    # Ensure peptides1 and peptides2 are dictionaries
    unique_to_seq1 = {peptide: position for peptide, position in peptides1.items() if peptide not in peptides2}

    # Refine unique peptides by merging and boundary adjustment
    refined_sequences = refine_unique_sequences(seq1, seq2, unique_to_seq1)

    return refined_sequences

# Step 1: Import genome fasta and annotations
genome = {rec.id: str(rec.seq) for rec in SeqIO.parse(genome_fasta, "fasta")}
has_chr_prefix = any(key.startswith("chr") for key in genome.keys())
gene_symbol_dict = load_gene_symbols(gene_symbol_file)

# Step 2: Import junctions
junctions = {}
with open(junctions_file, "r") as f:
    for line in f:
        line = line.strip()
        values = line.split("\t")
        if line:
            try: 
                gene, junction_name, chrom, start, end = parse_junction(values[0])
            except:
                continue
            source = junction_name
            if has_chr_prefix:
                pass
            else:
                chrom = chrom.replace('chr','')
            junction_name = junction_name.split(':')[1].replace('-', '|')
            try: 
                junctions[gene].append([junction_name, source, chrom, start, end])
            except: 
                junctions[gene] = [[junction_name, source, chrom, start, end]]

# Step 3: Import transcript associations
transcript_associations = {}
with open(transcript_associations_file, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            gene_id, strand, transcript_structure, transcript_id, _ = parts
            if gene_id not in transcript_associations:
                transcript_associations[gene_id] = []
            transcript_associations[gene_id].append((transcript_structure, transcript_id, strand))

protein_dict = load_protein_summary(protein_summary_file)

# Step 4: Import transmembrane domain genomic positions
uniprot_dict = load_reference_isoform_seq(uniprot_seq_file)
extracellular_domain_db, transmembrane_domain_db, transmembrane_domain_seq_db = load_domain_locations(protein_translation_dir,uniprot_coordinates,uniprot_dict)

# Step 5: Match junctions to transcripts and check sequences
results = []
with open(orf_sequences_file, "r") as orf_handle, open(output_file, "w") as out_handle:
    for record in SeqIO.parse(orf_handle, "fasta"):
        header = record.description  # Includes the full header line
        transcript_id = header.split(";")[0].split("_")[-1]
        transcript_id = transcript_id.split(" ")[0]
        #transcript_id = record.id.split(" ")[0]
        gene_id = None
        if "gene_id:" in header:
            gene_id = header.split("gene_id:")[1].split()[0]
        sequence = str(record.seq).upper()
        if gene_id in junctions:
            for junction_name, source_name, chrom, start, end in junctions[gene_id]:
                if gene_id in transcript_associations:
                    for transcript_structure, assoc_transcript_id, strand in transcript_associations[gene_id]:
                        #print ([[transcript_id,assoc_transcript_id,junction_name,strand]]);sys.exit()
                        if transcript_id == assoc_transcript_id and junction_name in transcript_structure:
                            # Extract junction sequence from genome
                            if strand == "+":
                                left_flank = genome[chrom][start - 10:start]  # 10 nt before the junction
                                right_flank = genome[chrom][end-1:end-1 + 10]  # 10 nt after the junction
                                junction_seq = left_flank + right_flank
                            else:
                                # For the negative strand, reverse-complement and swap left/right logic
                                left_flank = genome[chrom][start-1:start-1 + 10]  # 10 nt after the junction
                                right_flank = genome[chrom][end - 10:end]  # 10 nt before the junction
                                junction_seq = reverse_complement(right_flank+left_flank)
                            if source_name == 'ENSG00000137571:E11.1-E12.1':
                                if 'ENSX' in transcript_id:
                                    print (left_flank,right_flank,junction_seq)
                                    #print ([transcript_id,assoc_transcript_id,junction_name,strand])
                                    print (transcript_id,junction_seq,sequence)
                            if junction_seq in sequence:
                                try: 
                                    protein_info = protein_dict[gene_id,transcript_id]
                                except:
                                    protein_info = ''
                                found = 'Found'
                            else:
                                try: 
                                    protein_info = protein_dict[gene_id,transcript_id]
                                except:
                                    protein_info = ''
                                found = 'Not Found'
                            protein_sequence = Seq(sequence).translate(to_stop=True)
                            peptide = generate_forward_frames(junction_seq,protein_sequence)
                            EM_domain_presence, TM_overlap = domain_architecture_evaluation(gene_id,start,end)

                            ### Identify unique cell surface epitiope minimal sequence in the alternative isoform versus UniProt reference
                            unique_peptide_seq = ''
                            ref_length = ''
                            ref_isoform_seq = ''
                            tm_count = ''
                            tm_ref_count = ''
                            if found and gene_id in uniprot_dict: # Compare the isoform protein sequence to the UniProt reference
                                ref_length = len(uniprot_dict[gene_id])
                                ref_isoform_seq = uniprot_dict[gene_id]
                                try: unique = find_unique_sequences_with_refinement(protein_sequence, uniprot_dict[gene_id])
                                except:
                                    print (protein_sequence)
                                    print (uniprot_dict[gene_id]);sys.exit()
                                for (unique_peptide, flanking) in unique:
                                    if len(peptide)>0 and peptide in flanking:
                                        unique_peptide_seq = unique_peptide
                            
                                ### Report the number of transmembrane domains in the alternative and UniProt reference isoform
                                if gene_id in transmembrane_domain_seq_db:
                                    tm_sequences = transmembrane_domain_seq_db[gene_id]
                                    tm_count = len([peptide for peptide in tm_sequences if peptide in protein_sequence])
                                    tm_ref_count = len([peptide for peptide in tm_sequences if peptide in uniprot_dict[gene_id]])                               
                            results.append((gene_id,source_name,transcript_id,found,protein_info,protein_sequence,junction_seq,
                                    peptide,EM_domain_presence,TM_overlap,ref_isoform_seq,unique_peptide_seq,ref_length,tm_count,tm_ref_count,sequence))

                else:
                    print (gene_id,'failure2');sys.exit()
        else:
             #print (gene_id,'failure1');sys.exit()
             pass

    # Write results to output file
    out_handle.write(f"symbol\tgene\tjunction_name\ttranscript_id\tPeptide in ORF\tprotein_info\tNovel-isoform-sequence\tjunc-seq\tjunction-peptide\tEM_domain_presence\tTM_overlap\tunique_peptide_seq\tref_isoform_seq\tprotein_length)\tref_length\ttm_count\tref_tm_count\ORF\n")

    for gene,junction_name, transcript_id, status, protein_info, protein_sequence, junc_seq, peptide, EM_domain_presence, TM_overlap, ref_isoform_seq, unique_peptide_seq, ref_length, tm_count, tm_ref_count, orf in results:
        if gene in gene_symbol_dict:
            symbol = gene_symbol_dict[gene]
        else:
            symbol = ''
        out_handle.write(f"{symbol}\t{gene}\t{junction_name}\t{transcript_id}\t{status}\t{protein_info}\t{protein_sequence}\t{junc_seq}\t{peptide}\t{EM_domain_presence}\t{TM_overlap}\t{unique_peptide_seq}\t{ref_isoform_seq}\t{len(protein_sequence)}\t{ref_length}\t{tm_count}\t{tm_ref_count}\t{orf}\n")
