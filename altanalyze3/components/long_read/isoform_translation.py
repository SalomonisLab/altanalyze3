import re, sys, os
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

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

def is_nmd(transcript_seq, stop_codon_pos, last_exon_start):
    return stop_codon_pos < last_exon_start - 50

def extract_cds_and_protein(transcripts, genome_fasta, ref_last_exons=None, query_transcript_to_gene={}):
    genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    cds_records = []
    protein_records = []
    for transcript_id, data in tqdm(transcripts.items(), desc="Processing Transcripts"):
        #if transcript_id == 'ENST00000491706.5':
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
            orf_seq, orf_start = get_longest_orf(Seq(transcript_seq))
            if orf_seq:
                cds_seq = transcript_seq[orf_start:orf_start + len(orf_seq) * 3]
                description = f"CDS gene_id:{gene_id}" if gene_id else "CDS"
                cds_records.append(SeqRecord(Seq(cds_seq), id=transcript_id, description=description))
                stop_codon_pos = orf_start + len(orf_seq) * 3
                nmd_annotation = ""
                last_exon_length = exons[-1][2] - exons[-1][1] + 1
                last_exon_start_relative = len(transcript_seq) - last_exon_length
                #print (len(transcript_seq),exons[-1][2],exons[-1][1],last_exon_length,last_exon_start_relative,stop_codon_pos)
                nmd_annotation = " (Potential NMD)" if is_nmd(transcript_seq, stop_codon_pos, last_exon_start_relative) else ""
                protein_description = f"Protein gene_id:{gene_id}{nmd_annotation}" if gene_id else f"Protein{nmd_annotation}"
                protein_records.append(SeqRecord(Seq(orf_seq), id=transcript_id, description=protein_description))
        except KeyError as e:
            print(f"Error: {e}. Check if chromosome identifiers match between GFF and genome FASTA.")
    return cds_records, protein_records

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


def get_reference_last_exons(ref_gff_file):
    ref_last_exons = {}
    with open(ref_gff_file) as file:
        for line in tqdm(file, desc="Parsing Reference GFF"):
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attributes = parts[8]
            match = re.search(r'gene_id[ =]"?([^";]+)"?', attributes)
            if match:
                gene_id = match.group(1).split('.')[0]
                if gene_id not in ref_last_exons:
                    ref_last_exons[gene_id] = int(parts[3]) if parts[6] == '+' else int(parts[4])
                else:
                    if parts[6] == '+':
                        ref_last_exons[gene_id] = max(ref_last_exons[gene_id], int(parts[3]))
                    else:
                        ref_last_exons[gene_id] = min(ref_last_exons[gene_id], int(parts[4]))
    print(f"{len(ref_last_exons)} genes imported from the reference GFF")
    return ref_last_exons

def gff_translate(query_gff_file, genome_fasta, ref_gff_file=None, transcript_associations_file=None):
    query_transcripts = parse_gff(query_gff_file)
    query_transcript_to_gene, intron_retention_dict = parse_transcript_associations(transcript_associations_file) if os.path.exists(transcript_associations_file) else {}
    #ref_last_exons = get_reference_last_exons(ref_gff_file) if os.path.exists(ref_gff_file) else {}
    #query_transcript_to_gene = None
    cds_records, protein_records = extract_cds_and_protein(query_transcripts, genome_fasta, query_transcript_to_gene=query_transcript_to_gene)
    
    # Output protein information to CSV
    export_protein_summary(protein_records, intron_retention_dict, "protein_summary.csv")

    return cds_records, protein_records

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
        seq_len = len(record.seq)
        gene_id_match = re.search(r'gene_id:([^ ]+)', description)
        gene_id = gene_id_match.group(1) if gene_id_match else ""
        nmd_status = "Potential NMD" if "Potential NMD" in description else "Not NMD"
        intron_retention = intron_retention_dict.get(transcript_id, False)

        longest_isoform_length = gene_to_longest_isoform[gene_id]['length'] if gene_id in gene_to_longest_isoform else "NA"

        summary_data.append([
            transcript_id, 
            gene_id, 
            seq_len, 
            nmd_status, 
            intron_retention, 
            longest_isoform_length])

    summary_df = pd.DataFrame(summary_data, columns=[
        "Transcript ID", "Gene ID", "Protein Length", "NMD Status", "Intron Retention", "Longest Isoform Length"])
    summary_df.to_csv(output_csv_file, index=False)


if __name__ == '__main__':
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