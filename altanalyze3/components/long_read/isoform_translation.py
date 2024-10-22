import re, sys, os
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
from collections import defaultdict

def normalize_chromosome_name(chromosome):
    return chromosome.replace('chr', '')

def parse_gff(gff_file):
    transcripts = {}
    with open(gff_file) as file:
        lines = file.readlines()
        for line in tqdm(lines, desc="Parsing GFF"):
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            if feature_type == "exon":
                match = re.search(r'transcript_id[ =]"?([^";]+)"?', attributes)
                gene_match = re.search(r'gene_id[ =]"?([^";]+)"?', attributes)
                if match:
                    transcript_id = match.group(1)
                    gene_id = gene_match.group(1) if gene_match else ""
                    if transcript_id not in transcripts:
                        transcripts[transcript_id] = {"exons": [], "gene_id": gene_id}
                    normalized_chr = normalize_chromosome_name(parts[0])
                    transcripts[transcript_id]["exons"].append((normalized_chr, int(parts[3]), int(parts[4]), parts[6]))
    print(f"{len(transcripts)} transcripts imported from the GFF")
    return transcripts

def get_longest_orf(seq):
    max_orf = ''
    max_orf_start = -1
    strand, nuc = (+1, seq)  # Only consider the 5' to 3' sequence
    for frame in range(3):
        trans = str(nuc[frame:].translate(to_stop=False))
        trans_len = len(trans)
        aa_start = 0
        while aa_start < trans_len:
            aa_start = trans.find('M', aa_start)  # Find the next start codon
            if aa_start == -1:
                break
            aa_end = trans.find('*', aa_start)
            if aa_end == -1:
                aa_end = trans_len
            orf = trans[aa_start:aa_end]
            if len(orf) > len(max_orf):
                max_orf = orf
                max_orf_start = aa_start * 3 + frame
            aa_start = aa_end + 1
    return max_orf, max_orf_start

def is_nmd(transcript_seq, stop_codon_pos, last_exon_start, start_codon_pos, penultimate_exon_start):
    if stop_codon_pos < last_exon_start - 50:
        #"start codon condition", decreased accuracy when tested
        """if stop_codon_pos < start_codon_pos + 100:
           #return False
        else:
        """
        return True
    #predicts NMD when PTC is more than 2 exons away from end of transcript
    elif stop_codon_pos < penultimate_exon_start:
        return True
    return False

def extract_cds_and_protein(transcripts, genome_fasta, ref_first_exons=None, query_transcript_to_gene={}):
    genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    cds_records = []
    transcript_records = []
    protein_records = []
    for transcript_id, data in tqdm(transcripts.items(), desc="Processing Transcripts"):
        #if transcript_id=='ENST00000356073.8':
        exons = data["exons"]
        gene_id = data["gene_id"]
        if transcript_id in query_transcript_to_gene:
            gene_id = query_transcript_to_gene[transcript_id]
        exons.sort(key=lambda x: x[1])
        try:
            if exons[0][3] == '+':
                transcript_seq_list = [str(genome_seq[exon[0]].seq[exon[1]-1:exon[2]]) for exon in exons]
            else:
                transcript_seq_list = [str(genome_seq[exon[0]].seq[exon[1]-1:exon[2]].reverse_complement()) for exon in exons]
                transcript_seq_list = transcript_seq_list[::-1]  # Reverse for negative strand
            transcript_seq = ''.join(transcript_seq_list)

            ref_cds_seq = None
            ref_cds_start_pos = None
            if gene_id in ref_first_exons:
                ref_cds_seq = []
                for (ref_cds_start, ref_cds_end) in ref_first_exons[gene_id]:
                    ref_cds_seq_instance = str(genome_seq[exons[0][0]].seq[ref_cds_start-1:ref_cds_end])
                    if exons[0][3] == '-':
                        ref_seq = Seq(ref_cds_seq_instance)
                        ref_cds_seq_instance = str(ref_seq.reverse_complement())
                    ref_cds_seq.append(ref_cds_seq_instance)

            # Restrict to ORFs with defined start codons if transcript has matching starts
            orf_seq = None
            if ref_cds_seq:
                candidates=[]
                for ref_cds_seq_instance in ref_cds_seq:
                    if len(ref_cds_seq_instance) >= 10:
                        if ref_cds_seq_instance in transcript_seq:
                            orf_start = transcript_seq.find(ref_cds_seq_instance)
                            transcript_seq_alt = Seq(transcript_seq[orf_start:])
                            orf_seq = str(transcript_seq_alt.translate(to_stop=True))
                            candidates.append((len(orf_seq),orf_start,orf_seq))
                if len(candidates)>0:
                    candidates.sort()
                    candidates.reverse()
                    ol,orf_start,orf_seq = candidates[0] # select the most 5' start site

            if orf_seq == None:
                orf_seq, orf_start = get_longest_orf(Seq(transcript_seq))

            if orf_seq:
                cds_seq = transcript_seq[orf_start:orf_start + len(orf_seq) * 3]
                description = f";CDS;gene_id:{gene_id}" if gene_id else "CDS"
                cds_records.append(SeqRecord(Seq(cds_seq), id=transcript_id, description=description))
                description = f";transcript;gene_id:{gene_id}" if gene_id else "CDS"
                transcript_records.append(SeqRecord(Seq(transcript_seq), id=transcript_id, description=description))
                stop_codon_pos = orf_start + len(orf_seq) * 3
                nmd_annotation = ""
                if exons[0][3] == '-':
                    last_exon_length = exons[0][2] - exons[0][1] + 1
                    try: penultimate_exon_length = exons[1][2] - exons[1][1] + 11
                    except: penultimate_exon_legth = 0
                else:
                    last_exon_length = exons[-1][2] - exons[-1][1] + 1
                    try: penultimate_exon_length = exons[-2][2] - exons[-2][1] + 1
                    except: penultimate_exon_length = 0
                last_exon_start_relative = len(transcript_seq) - last_exon_length
                penultimate_exon_start_relative = last_exon_start_relative - penultimate_exon_length #check if last intron length needs to be included?
                nmd_annotation = ";(NMD)" if is_nmd(transcript_seq, stop_codon_pos, last_exon_start_relative, orf_start, penultimate_exon_start_relative) else ""
                protein_description = f";Protein';'gene_id:{gene_id}{nmd_annotation}" if gene_id else f"Protein{nmd_annotation}"
                protein_records.append(SeqRecord(Seq(orf_seq), id=transcript_id, description=protein_description))
                #print (orf_start,(len(orf_seq) * 3),last_exon_length,last_exon_start_relative,stop_codon_pos,transcript_seq,[nmd_annotation])
        except KeyError as e:
            print(f"Error: {e}. Check if chromosome identifiers match between GFF and genome FASTA.")
    return cds_records, transcript_records, protein_records

def parse_transcript_associations(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ['gene_id', 'strand', 'exon_structure', 'transcript_id', 'source']
    transcript_to_gene = {}
    intron_retention_dict = {}

    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Parsing Transcript Associations"):
        new_transcript_id = f"{row['source']}_{row['transcript_id']}"
        transcript_to_gene[new_transcript_id] = row['gene_id']
        transcript_to_gene[row['transcript_id']] = row['gene_id']

        exon_structure = row['exon_structure'].split('|')
        intron_retention = any(part.startswith('I') and '_' not in part for part in exon_structure)
        intron_retention_dict[new_transcript_id] = intron_retention
        intron_retention_dict[row['transcript_id']] = intron_retention

    return transcript_to_gene, intron_retention_dict

def get_reference_first_exons(ref_gff_file):
    ref_first_exons = {}
    transcript_cds = {}
    with open(ref_gff_file) as file:
        for line in tqdm(file, desc="Parsing Reference GFF"):
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attributes = parts[8]
            gene_match = re.search(r'gene_id[ =]"?([^";]+)"?', attributes)
            transcript_match = re.search(r'transcript_id[ =]"?([^";]+)"?', attributes)
            if gene_match:
                gene_id = gene_match.group(1).split('.')[0]
                transcript_id = transcript_match.group(1).split('.')[0]
                cds_start = int(parts[3])
                cds_end = int(parts[4])
                if transcript_id not in transcript_cds: # Ensure it is the first exon for that transcript
                    try: 
                        if (cds_start, cds_end) not in ref_first_exons[gene_id]:
                            ref_first_exons[gene_id].append((cds_start, cds_end))
                    except: ref_first_exons[gene_id] = [(cds_start, cds_end)]
                transcript_cds[transcript_id]=[]

    print(f"{len(ref_first_exons)} genes imported from the reference GFF")
    return ref_first_exons

def gff_translate(query_gff_file, genome_fasta, ref_gff_file=None, transcript_associations_file=None):
    query_transcripts = parse_gff(query_gff_file)
    
    query_transcript_to_gene, intron_retention_dict = parse_transcript_associations(transcript_associations_file) if os.path.exists(transcript_associations_file) else {}
    ref_first_exons = get_reference_first_exons(ref_gff_file) if os.path.exists(ref_gff_file) else {}
    cds_records, transcript_records, protein_records = extract_cds_and_protein(query_transcripts, genome_fasta, query_transcript_to_gene=query_transcript_to_gene, ref_first_exons=ref_first_exons)
    
    # Output protein information to CSV
    export_protein_summary(protein_records, intron_retention_dict, "protein_summary.txt")

    return cds_records, transcript_records, protein_records

def export_protein_summary(protein_records, intron_retention_dict, output_csv_file):
    summary_data = []
    gene_to_longest_isoform = {}

    # First pass to determine the longest isoform for each gene
    for record in protein_records:
        transcript_id = record.id
        description = record.description
        seq_len = len(record.seq)
        gene_id_match = re.search(r'gene_id:([^ ]+)', description)
        gene_id = gene_id_match.group(1) if gene_id_match else ""

        if gene_id:
            if gene_id not in gene_to_longest_isoform or seq_len > gene_to_longest_isoform[gene_id]['length']:
                gene_to_longest_isoform[gene_id] = {'transcript_id': transcript_id, 'length': seq_len}

    # Second pass to create the summary data
    for record in protein_records:
        transcript_id = record.id
        description = record.description
        prot_seq = record.seq
        aa_start = prot_seq.find('M', 0)
        seq_len = len(prot_seq[aa_start:])
        gene_id_match = re.search(r'gene_id:([^ ]+)', description)
        gene_id = gene_id_match.group(1) if gene_id_match else ""
        nmd_status = "Potential-NMD" if "NMD" in description else "Not-NMD"
        if seq_len<10:
            nmd_status = "NMD"
        intron_retention = intron_retention_dict.get(transcript_id, False)

        longest_isoform_length = gene_to_longest_isoform[gene_id]['length'] if gene_id in gene_to_longest_isoform else "NA"

        summary_data.append([
            gene_id, 
            transcript_id,
            seq_len, 
            nmd_status, 
            intron_retention, 
            longest_isoform_length])

    summary_df = pd.DataFrame(summary_data, columns=[
       "Gene ID",  "Transcript ID", "Protein Length", "NMD Status", "Intron Retention", "Longest Isoform Length"])
    summary_df.to_csv(output_csv_file, index=False, sep='\t')

def find_redundant_proteins(fasta_file, output_file, restrict=None):
    # Dictionary to hold sequences as keys and list of corresponding protein IDs as values
    sequence_dict = defaultdict(list)
    
    # Read the FASTA file and populate the dictionary
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        protein_id = record.id
        if restrict is not None:
            if restrict in protein_id:
                sequence_dict[sequence].append(protein_id)
        else:
            sequence_dict[sequence].append(protein_id)
    
    # Write the output to a two-column text file
    with open(output_file, 'w') as out_f:
        for sequence, protein_ids in sequence_dict.items():
            isoform_ids=[]
            for isoform in protein_ids:
                isoform = isoform.split(';')[0]
                if 'PB.' in isoform:
                    isoform = 'PB'+isoform.split('_PB')[1]
                isoform_ids.append(isoform)
            gene_id = protein_ids[0].split(';')[2]  # Assuming gene_id is the prefix before '|' in protein_id
            gene_id = gene_id.split(':')[1]
            redundant_ids = ','.join(isoform_ids)
            out_f.write(f"{gene_id}\t{redundant_ids}\t{sequence}\n")

    print(f"Output written to {output_file}")

if __name__ == '__main__':

    fasta_file = "protein_sequences.fasta"
    #fasta_file = 'cds_sequences.fasta'
    #fasta_file = 'transcript_sequences.fasta'
    output_file = "redundant_proteins.txt"
    #output_file = "redundant_cDNA.txt"
    #output_file = "redundant_transcripts.txt"
    restrict = 'ENSG00000109689'
    find_redundant_proteins(fasta_file, output_file, restrict=restrict);sys.exit()

    # Example usage
    ref_gff_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/MANE.GRCh38.v0.5.select_ensembl_genomic.gff"
    query_gff_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gff-output/combined.gff"
    genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
    transcript_associations_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gff-output/transcript_associations.txt"

    cds_records, protein_records = gff_translate(
        query_gff_file, 
        genome_fasta, 
        ref_gff_file=ref_gff_file, 
        transcript_associations_file=transcript_associations_file
    )

    # Output results
    with open("cds_sequences.fasta", "w") as cds_file:
        SeqIO.write(cds_records, cds_file, "fasta")

    with open("protein_sequences.fasta", "w") as protein_file:
        SeqIO.write(protein_records, protein_file, "fasta")