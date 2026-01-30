""" This module combines information across GFF files to compare isoform structures, 
annotate exons/transcripts and to nominate the longest exemplar isoforms. Note that
this script is biased against isoforms that are bleeding and internal hybrid."""

import os,sys
import argparse
from collections import defaultdict, OrderedDict
from tqdm import tqdm
import pandas as pd
import collections
import gzip

novelGene = 0
novelGeneLoci = [1,1]
additional_junctions={}

def findNovelSpliceSite(geneID, position, strand):
    exons = geneData[geneID]
    for i, (start, stop, exon) in enumerate(exons):
        """
        if strand == '-' and position>258251:
            print (geneID, position, strand)
            print (exons,start,position,stop,exon)
            pass"""
        if (strand == '+' and position < start) or (strand == '-' and position > start):
            #print ('a');sys.exit()
            if i == 0:
                return f"U1.1_{position}" # Position is upstream of the first exon
            else:
                pass # not in the exon/intron
        elif (strand == '+' and start <= position <= stop) or (strand == '-' and stop <= position <= start):
            #print ('b');sys.exit()
            if 'I' in exon and exon.endswith('.1'): # Position is within the current intron
                #print ('c');sys.exit()
                if (strand == '+' and (position - start < 50)) or (strand == '-' and (start - position < 50)):
                    #if strand == '-':  print(1, exon,position,start, stop, strand,f"{exons[i - 1][2]}_{position}");sys.exit()
                    # Within 50nt of the beginning of the intron - alt SS of prior exon region
                    return f"{exons[i - 1][2]}_{position}"
                elif (strand == '+' and (stop - position < 50)) or (strand == '-' and (position - stop < 50)):
                    #print(1, exon,position,start, stop, strand,f"{exons[i - 1][2]}_{position}")
                    #print ('d');sys.exit()
                    # Within 50nt of the end of the intron - alt SS of next exon region
                    try:
                        return f"{exons[i + 1][2]}_{position}"
                    except:
                        print (position,exons);sys.exit()
                else: # In the middle of the intron
                    #print ('k');sys.exit()
                    return f"{exon}_{position}"
            else:
                #print ('e');sys.exit()
                # Position is within the current exon
                return f"{exon}_{position}"
    #print ('f');sys.exit()
    last_exon_id = exons[-1][2] # Position is after the last exon
    prefix, number = last_exon_id[0], last_exon_id[1:].split('.')
    next_number = str(int(number[0]) + 1)
    return f"U{next_number}.1_{position}"

def rangeCheck(exon_range,novelGeneLoci):
    if exon_range[0] >= novelGeneLoci[0] and exon_range[0] <= novelGeneLoci[1]:
        return True
    elif exon_range[1] >= novelGeneLoci[0] and exon_range[1] <= novelGeneLoci[1]:
        return True
    else:
        return False

def exonAnnotate(chr, exons, strand, transcript_id):
    """Assign exon/intron-level annotations"""
    global novelGene
    global novelGeneLoci
    global additional_junctions

    exon_range = [exons[0][0],exons[-1][-1]]
    ### First determine the genes and initial known assignments
    associations = []
    genes = []
    if strand == '-':
        if exon_range[1]>exon_range[0]: ### Gencode is ordered transcriptomically
            exons.reverse()
        exons = [(b, a) for a, b in exons]
    for (pos1, pos2) in exons:
        if (chr, pos1, strand, 1) in exonCoordinates:
            gene1, exon1, ind1 = exonCoordinates[chr, pos1, strand, 1]
            if gene1 not in genes:
                genes.append(gene1)
        else:
            gene1, exon1, strand1 = None, None, None
        if (chr, pos2, strand, 2) in exonCoordinates:
            gene2, exon2, ind2 = exonCoordinates[chr, pos2, strand, 2]
            if gene2 not in genes:
                genes.append(gene2)
        else:
            gene2, exon2, strand2 = None, None, None
        if gene1 is None and gene2 is None:
            association = 'novel', (pos1,pos2)
        elif gene1 is not None and gene2 is not None:
            if gene1 == gene2:
                if exon1 == exon2:
                    association = 'common', (gene1, exon1)
                else:
                    association = 'mixed', (gene1, exon1, exon2)
            else:
                # Trans-splicing
                association = 'trans', (gene1, gene2, exon1, exon2)
        else:
            if gene1 is None:
                association = 'novel1', (gene2, exon2, pos1)
            else:
                association = 'novel2', (gene1, exon1, pos2)
        associations.append(association)

    ### Determine the spectrum of possible gene IDs (trans-splicing)
    if len(genes) == 0:
        """
        if 'PB.9.48' == transcript_id:
            print(novelGeneLoci,exon_range);sys.exit()
        """
        if rangeCheck(exon_range,novelGeneLoci) == False:
            novelGene += 1
        gene = 'UNK' + str(novelGene)
        novelGeneLoci = exon_range
    elif len(genes) > 1:
        if len(genes) == 3:
            #print('Multi-Trans-splicing', genes)
            pass
        gene = genes[0]
    else:
        gene = genes[0]
    current_gene = gene # used to document current gene for trans-splicing

    ### Annotate novel exon/intron IDs and include in-between regions in the final isoform annotation
    if len(genes) == 0:
        output_list = [f"{start}-{end}" for start, end in exons]
        associations_updated = output_list  # Novel genes don't require exon annotations (likely update - currently position tuple)
        associations_simple = output_list
    else:
        """
        if 'PB.2986.5' == transcript_id:
            print(exons,associations);sys.exit()
        """
        associations_updated = []
        associations_simple = []
        for ind, association in enumerate(associations):
            exonType = association[0]

            if exonType == 'common':
                gene1 = association[1][0]
                exonAnnotation = association[1][1]
                if gene1 != gene or current_gene != gene:
                    exonAnnotation = gene1 + ':' + exonAnnotation
                associations_updated.append(exonAnnotation) # denotes a splice site or end
                associations_simple.append(exonAnnotation)
            elif exonType == 'mixed':
                exonAnnotation1 = association[1][1]
                exonAnnotation2 = association[1][2]
                gene1 = association[1][0]
                start,stop = exons[ind]
                start_index = exonCoordinates[chr, start, strand, 1][-1]
                end_index = exonCoordinates[chr, stop, strand, 2][-1]+1
                if gene1 != gene or current_gene != gene:
                    exonAnnotation1 = gene1 + ':' + exonAnnotation1
                    current_gene = gene1
                associations_updated.append(exonAnnotation1)
                if ind != 0: # don't add the begining of the exon
                    associations_simple.append(exonAnnotation1)
                prior_splice_acceptor = None
                ref_exons = geneData[current_gene]
                # Add exon and intron regions in the middle of a broader exon
                for in_between_region in ref_exons[start_index:end_index]:
                    if current_gene != gene: # trans-splicing
                        associations_updated.append(current_gene+':'+in_between_region[-1])
                    else:
                        associations_updated.append(in_between_region[-1])
                    current_spice_donor = in_between_region[0]
                    if prior_splice_acceptor != None:
                        if transcript_id not in additional_junctions:
                            additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                        else:
                            additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                    prior_splice_acceptor = in_between_region[-2]
                if gene1 != gene or current_gene != gene:
                    exonAnnotation2 = gene1 + ':' + exonAnnotation2
                    current_gene = gene1
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
                if ind != (len(associations)-1): # Don't add the last splice site
                    associations_simple.append(exonAnnotation2)
                """
                if len(ref_exons[start_index:end_index])>1 and 'I2.1' in associations_updated and strand == '-':
                    print(exons)
                    print(associations)
                    print('****',start_index,end_index,ref_exons[start_index:end_index])
                    print(associations_updated);sys.exit()"""
            elif exonType == 'trans':
                gene1 = association[1][0]
                gene2 = association[1][1]
                exonAnnotation1 = association[1][2]
                exonAnnotation2 = association[1][3]
                if gene == gene1:
                    associations_updated.append(exonAnnotation1)
                    associations_simple.append(exonAnnotation1)
                    associations_updated.append(gene2 + ':' + exonAnnotation2)
                    associations_simple.append(gene2 + ':' + exonAnnotation2)
                elif gene == gene2:
                    associations_updated.append(gene1 + ':' + exonAnnotation1)
                    associations_simple.append(gene1 + ':' + exonAnnotation1)
                    associations_updated.append(exonAnnotation2) # denotes a splice site or end
                    associations_simple.append(exonAnnotation2)
                    current_gene = gene2

            elif exonType == 'novel1':
                novel_position = association[1][2]
                exonAnnotation2 = association[1][1]
                gene2 = association[1][0]
                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    #if len(genes)>1: 
                    #print(association,strand,gene2,novel_position,genes);sys.exit()
                    exonAnnotation1 = findNovelSpliceSite(geneID, novel_position, strand)
                    hits.append([exonAnnotation1,geneID])
                    #if len(genes)>1: print(exonAnnotation1,geneID)
                    #if 'U' not in exonAnnotation1:
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation1,geneID = hits[0]
                """if strand == '-' and novel_position>164713:
                    print (exonAnnotation1);sys.exit()"""

                if gene != gene2 or current_gene != gene: # trans-splicing somewhere else in this transcript
                    exonAnnotation1 = gene2+':'+exonAnnotation1
                    current_gene = gene2
                associations_updated.append(exonAnnotation1)
                associations_simple.append(exonAnnotation1)
                # If exons and/or introns span from the begining of a novel exon to a known end, include all regions in between (e.g., U1.1_29934 to E3.1)
                # start_index will not exist in the database (novel end)
                start,stop = exons[ind]
                end_index = exonCoordinates[chr, stop, strand, 2][-1]
                evaluate_between_regions = True
                if 'U' in exonAnnotation1:
                    start_index = 0
                else:
                    try:
                        if 'I' in exonAnnotation1:
                            # If the splice site is in the intron the last portion of the intron is retained
                            # Thus, include the intron splice site coordinates in the search to denote IR
                            start,stop = exonData[geneID,exonAnnotation1.split("_")[0].replace('I','E')]
                            start_index = exonCoordinates[chr, start, strand, 1][-1] + 1
                        else:
                            start,stop = exonData[geneID,exonAnnotation1.split("_")[0]]
                            start_index = exonCoordinates[chr, start, strand, 1][-1]
                    except:
                        evaluate_between_regions = False

                if evaluate_between_regions:
                    ref_exons = geneData[current_gene]
                    prior_splice_acceptor = None
                    # Add exon and intron regions in the middle of a broader exon
                    for in_between_region in ref_exons[start_index:end_index]:
                        if current_gene != gene: # trans-splicing
                            associations_updated.append(current_gene+':'+in_between_region[-1])
                        else:
                            associations_updated.append(in_between_region[-1])
                        current_spice_donor = in_between_region[0]
                        if prior_splice_acceptor != None:
                            if transcript_id not in additional_junctions:
                                additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                            else:
                                additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                        prior_splice_acceptor = in_between_region[-2]

                if gene != gene2 or current_gene != gene: # trans-splicing somewhere else in this transcript
                    exonAnnotation2 = gene2+':'+exonAnnotation2
                    current_gene = gene2
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
                associations_simple.append(exonAnnotation2)
            elif exonType == 'novel2':
                novel_position = association[1][2]
                exonAnnotation1 = association[1][1]
                gene1 = association[1][0]
                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    exonAnnotation2 = findNovelSpliceSite(geneID, novel_position, strand)
                    hits.append([exonAnnotation2,geneID])
                    #if len(genes)>1: print(exonAnnotation1,geneID)
                    #if 'U' not in exonAnnotation1:
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation2,geneID = hits[0]
                """if strand == '-' and novel_position>258251:
                    print (exonAnnotation2);sys.exit()"""

                if gene != gene1 or current_gene != gene: # trans-splicing somewhere else in this transcript
                    exonAnnotation1 = gene1+':'+exonAnnotation1
                    current_gene = gene1
                associations_updated.append(exonAnnotation1)
                associations_simple.append(exonAnnotation1)
                # If exons and/or introns span from the begining of a known exon to the novel end, include all regions in between (e.g., E10.1 to U15.1_1234)
                # end_index will not exist in the database (novel end)
                start,stop = exons[ind]
                start_index = exonCoordinates[chr, start, strand, 1][-1]
                evaluate_between_regions = True
                if 'U' in exonAnnotation2:
                    end_index = None
                else:
                    try:
                        eid = exonAnnotation2.split("_")[0].replace('I','E')
                        start,stop = exonData[geneID,eid]
                        end_index = exonCoordinates[chr, stop, strand, 2][-1]
                        if 'I' in exonAnnotation2 and start_index == end_index: # (E3.1, I3.1_12345) - include the exon-intron boundary (bleeding exon)
                            prefix, number = eid[0], eid[1:].split('.')
                            eid2 =  f"E{str(int(number[0]) + 1)}.1"
                            start,stop = exonData[geneID,eid2]
                            end_index = exonCoordinates[chr, stop, strand, 2][-1]
                    except:
                        # Intronic position is associated with a strange region, such as ENSG00000127483-I15.2
                        evaluate_between_regions = False
                if evaluate_between_regions:
                    ref_exons = geneData[current_gene]
                    prior_splice_acceptor = None
                    # Add exon and intron regions in the middle of a broader exon
                    for in_between_region in ref_exons[start_index:end_index]:
                        if current_gene != gene: # trans-splicing
                            associations_updated.append(current_gene+':'+in_between_region[-1])
                        else:
                            associations_updated.append(in_between_region[-1])
                        current_spice_donor = in_between_region[0]
                        if prior_splice_acceptor != None:
                            if transcript_id not in additional_junctions:
                                additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                            else:
                                additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                        prior_splice_acceptor = in_between_region[-2]

                if gene != geneID or current_gene != gene: 
                    exonAnnotation2 = geneID+':'+exonAnnotation2
                    current_gene = geneID
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
                associations_simple.append(exonAnnotation2)

            elif exonType == 'novel':
                novel_position1 = association[1][0]
                novel_position2 = association[1][1]
                """
                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    exonAnnotation1 = findNovelSpliceSite(geneID, novel_position1, strand)
                    hits.append([exonAnnotation1,geneID])
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation1,geneID = hits[0]
                """
                exonAnnotation1 = findNovelSpliceSite(current_gene, novel_position1, strand)
                if current_gene != gene: 
                    exonAnnotation1 = current_gene+':'+exonAnnotation1
                """
                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    exonAnnotation2 = findNovelSpliceSite(geneID, novel_position2, strand)
                    hits.append([exonAnnotation2,geneID])
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation2,geneID = hits[0]
                """
                exonAnnotation2 = findNovelSpliceSite(current_gene, novel_position2, strand)
                if current_gene != gene: 
                    exonAnnotation2 = current_gene+':'+exonAnnotation2
                ### Need to fill in-between exons
                associations_updated.append(exonAnnotation1)
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
                associations_simple.append(exonAnnotation1)
                associations_simple.append(exonAnnotation2)
    """
    if 'PB.2986.5' == transcript_id:
        print('***',associations)
        print(associations_updated);sys.exit()
    """
    associations_updated = uniqueList(associations_updated)
    associations_simple = uniqueList(associations_simple)

    return gene,associations_updated,associations_simple,genes

def uniqueList(ls):
    return [seen.add(x) or x for seen in (set(),) for x in ls if x not in seen]

def knownIsoform(transcripts):
    search_strings = ["ENST", "NM_", "XM_", "XR_"]
    found = any(any(s in transcript for s in search_strings) for transcript in transcripts)
    if found:
        return 'known'
    else:
        return 'novel'

def selectKnownIsoform(transcripts):
    search_strings = ["ENST", "NM_", "XM_", "XR_"]
    found=[]
    for transcript in transcripts:
        for s in search_strings:
            if transcript.startswith(s):
                found.append(transcript)
    found.sort()
    for transcript in transcripts:
        if transcript not in found:
            found.append(transcript)
    return found


def _split_token(token, default_gene):
    gene_id = default_gene
    core = token
    if ':' in token:
        gene_id, core = token.split(':', 1)
    coord = None
    base = core
    if '_' in core:
        base, coord_str = core.rsplit('_', 1)
        try:
            coord = int(coord_str)
        except ValueError:
            base = core
            coord = None
    return gene_id, base, coord


def _strip_terminal_coords(tokens, default_gene):
    if not tokens:
        return []
    trimmed = list(tokens)
    drop_indices = set()
    for idx in (0, len(trimmed) - 1):
        _, base, coord = _split_token(trimmed[idx], default_gene)
        if coord is not None and base.startswith('E'):
            drop_indices.add(idx)
    if drop_indices:
        return [tok for i, tok in enumerate(trimmed) if i not in drop_indices]
    return trimmed


def _filter_exon_intron_tokens(tokens, default_gene):
    kept_tokens = []
    for tok in tokens:
        gene_id, base, _ = _split_token(tok, default_gene)
        if base.startswith(('E', 'I')):
            if base.startswith('E') and '_' in tok:
                continue
            if gene_id != default_gene and ':' not in tok:
                tok = f"{gene_id}:{base}"
            kept_tokens.append(tok)
    return kept_tokens


def _longest_common_substring_length(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0
    if len(seq_a) > len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    prev = [0] * (len(seq_a) + 1)
    best = 0
    for token in seq_b:
        current = [0]
        for idx, a_token in enumerate(seq_a, start=1):
            if token == a_token:
                val = prev[idx - 1] + 1
            else:
                val = 0
            current.append(val)
            if val > best:
                best = val
        prev = current
    return best


def _substring_similarity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    short = seq_a if len(seq_a) <= len(seq_b) else seq_b
    long = seq_b if short is seq_a else seq_a
    lcs = _longest_common_substring_length(short, long)
    return lcs / float(len(short))


def _cluster_by_substring(items, threshold=0.85):
    clusters = []
    for item in sorted(items, key=lambda x: x['weight'], reverse=True):
        best_idx = None
        best_score = -1.0
        for idx, cluster in enumerate(clusters):
            min_score = 1.0
            for other in cluster['items']:
                score = _substring_similarity(item['cluster_seq'], other['cluster_seq'])
                if score < min_score:
                    min_score = score
                if min_score < threshold:
                    break
            if min_score >= threshold and min_score > best_score:
                best_score = min_score
                best_idx = idx
        if best_idx is not None:
            clusters[best_idx]['items'].append(item)
        else:
            clusters.append({'items': [item]})
    return [cluster['items'] for cluster in clusters]


def _select_cluster_label(group):
    best_label = None
    best_score = None
    for structure in group:
        for iso in structure['items']:
            label = iso.get('isoform_id')
            if not label:
                continue
            tokens = iso.get('tokens_trimmed') or iso.get('tokens') or []
            score = (len(tokens), iso.get('count', 0), len(str(label)))
            if best_score is None or score > best_score:
                best_score = score
                best_label = label
    return best_label


def collapseIsoformsCluster(gene_db, junction_db, gff_organization, mode):
    num_isoforms = 0
    collaped_db = defaultdict(lambda: defaultdict(list))

    for gene, isoforms in tqdm(gene_db.items(), desc="Collapsing Isoforms (cluster)"):
        structure_groups = OrderedDict()
        for isoform in isoforms:
            tokens = [t for t in isoform.split('|') if t]
            if not tokens:
                continue
            support = len(junction_db.get((gene, isoform), []))
            if support == 0:
                support = 1
            key = tuple(tokens)
            tokens_trimmed = _strip_terminal_coords(tokens, gene)
            cluster_seq = _filter_exon_intron_tokens(tokens_trimmed, gene)
            if key not in structure_groups:
                structure_groups[key] = {
                    'structure_key': key,
                    'items': [],
                    'weight': 0,
                    'cluster_seq': cluster_seq
                }
            structure_groups[key]['items'].append({
                'isoform_id': isoform,
                'tokens': tokens,
                'tokens_trimmed': tokens_trimmed,
                'count': support
            })
            structure_groups[key]['weight'] += support

        structure_items = list(structure_groups.values())
        grouped_structures = _cluster_by_substring(structure_items, threshold=0.85)
        for group in grouped_structures:
            label = _select_cluster_label(group)
            if label is None:
                continue
            sub_isoforms = []
            for structure in group:
                for iso in structure['items']:
                    if iso['isoform_id'] != label:
                        sub_isoforms.append(iso['isoform_id'])
            collaped_db[gene][label] = sub_isoforms
        num_isoforms += len(grouped_structures)
    print(num_isoforms, 'collapsed cluster unique isoforms')
    return collaped_db
    
def collapseIsoforms(gene_db, junction_db, gff_organization, mode):
    num_isoforms = 0
    collaped_db = defaultdict(lambda: defaultdict(list))

    def get_transcript_ids(gene, isoforms):
        transcript_ids = defaultdict(lambda: defaultdict(list))
        for iso in isoforms:
            for (file, info) in junction_db[gene, iso]:
                transcript_id = info.split(';')[gff_organization[file][0]].split(gff_organization[file][1])[1]
                transcript_ids[transcript_id] = iso
        return transcript_ids
            
    # Adding a progress bar for the outer loop iterating over genes
    for gene, isoforms in tqdm(gene_db.items(), desc="Collapsing Isoforms"):
        sorted_isoforms = sorted(isoforms, key=len, reverse=True)
        # Initially set every isoform as a key with an empty list
        for isoform in sorted_isoforms:
            collaped_db[gene][isoform] = []
        # Populate the list with sub-isoforms
        for i, isoform in enumerate(sorted_isoforms):
            for shorter_isoform in sorted_isoforms[i+1:]:
                # Check if one isoform is a subset of another
                if shorter_isoform in isoform:
                    if isoform != shorter_isoform:
                        collaped_db[gene][isoform].append(shorter_isoform)
        # Remove super-isoforms that are present as sub-isoforms
        for super_isoform in list(collaped_db[gene].keys()):
            for sub_isoforms in collaped_db[gene].values():
                if super_isoform in sub_isoforms:
                    del collaped_db[gene][super_isoform]
                    break

        num_isoforms += len(collaped_db[gene])
    print(num_isoforms, 'collapsed unique isoforms')

    if mode == 'Ensembl':
        num_isoforms = 0

        super_isoform_collapsed_db = collaped_db
        collaped_db = defaultdict(lambda: defaultdict(list))
        keys_added = set()
        gene_to_isoforms = defaultdict(list)
        for (gene, isoform) in junction_db:
            gene_to_isoforms[gene].append(isoform)

        # Adding progress bar for Ensembl processing
        for gene, isoform_list in tqdm(gene_to_isoforms.items(), desc="Processing Ensembl Isoforms"):
            super_isoforms = super_isoform_collapsed_db[gene]
            if not super_isoforms:
                continue

            iso_to_transcripts = {}
            for isoform in isoform_list:
                transcripts = []
                for (file, info) in junction_db[gene, isoform]:
                    ti, td = gff_organization[file]
                    transcript_id = info.split(';')[ti].split(td)[1]
                    transcripts.append(transcript_id)
                iso_to_transcripts[isoform] = transcripts

            super_to_enst = {}
            sub_to_super = {}
            for super_iso, sub_isoforms in super_isoforms.items():
                combined_iso = [super_iso] + sub_isoforms
                transcript_ids = {}
                for iso in combined_iso:
                    for transcript_id in iso_to_transcripts.get(iso, []):
                        transcript_ids[transcript_id] = iso
                found = next((x for x in transcript_ids if x.startswith('ENST')), False)
                if found:
                    super_to_enst[super_iso] = transcript_ids[found]
                for sub_iso in sub_isoforms:
                    if sub_iso not in sub_to_super:
                        sub_to_super[sub_iso] = super_iso

            for isoform in isoform_list:
                ens_iso = None
                if isoform in super_to_enst:
                    ens_iso = super_to_enst[isoform]
                else:
                    super_iso = sub_to_super.get(isoform)
                    if super_iso and super_iso in super_to_enst:
                        ens_iso = super_to_enst[super_iso]
                if not ens_iso:
                    continue
                if ens_iso == isoform:
                    if isoform not in collaped_db[gene][ens_iso]:
                        collaped_db[gene][ens_iso].append(isoform)
                else:
                    collaped_db[gene][ens_iso].append(isoform)
                keys_added.add((gene, isoform))
                keys_added.add((gene, ens_iso))
        #print ('dd',collaped_db['ENSG00000183072']['E1.1|E1.2|E1.3|E1.4|E3.3|E3.4|E3.5|'])
        # Add isoforms that were not processed earlier
        for (gene, isoform) in junction_db:
            if (gene, isoform) not in keys_added:
                collaped_db[gene][isoform] = []

        # Count the total number of collapsed isoforms
        for gene in collaped_db:
            num_isoforms += len(collaped_db[gene])
        
        print(num_isoforms, 'collapsed Ensembl prioritized unique isoforms')
    #print (collaped_db['ENSG00000183072']['E1.1|E1.2|E1.3|E1.4|E3.3|E3.4|E3.5|'])#;sys.exit()
    return collaped_db


def consolidateLongReadGFFs(directory, exon_reference_dir, mode="collapse"):
    """Only return isoforms with unique splice-site combinations"""
    junction_db = collections.OrderedDict()
    junction_str_db = collections.OrderedDict()
    gene_db = collections.OrderedDict()
    strand_db = collections.OrderedDict()
    trans_spliced_isoforms = collections.OrderedDict()
    global novelGene
    global novelGeneLoci
    global additional_junctions

    # Import preprocessed Ensembl exon information
    importEnsemblGenes(exon_reference_dir)
    
    collapse_isoforms = False
    if isinstance(directory, list):
        # list of gff's from different locations
        files = directory
        directory = os.getcwd() # present working directory
        if len(files)>1:
            collapse_isoforms = True
    elif '.g' in directory.lower():
        # Export exon structure information only
        gff = directory
        directory = os.path.dirname(gff)
        gff = os.path.basename(gff)
        files = [gff]
    else:
        # Combine and report collapsed unique isoform structures 
        files = [file for file in os.listdir(directory) if file.endswith('.gff')]
        files += [file for file in os.listdir(directory) if file.endswith('.gtf')]
        #files = [file for file in os.listdir(directory) if file.endswith('.txt')]
        if len(files)>1:
            collapse_isoforms = True
    
    combined_dir = os.path.join(directory, 'gff-output')
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)
    transcript_associations = os.path.join(combined_dir, 'transcript_associations.txt')
    eo = open(transcript_associations, 'w')

    def getJunctions(exons):
        splice_junctions = []
        if len(exons) == 1:
            splice_junctions = exons  # intron retention or single-exon transcript
        elif strand == '+':
            for i in range(len(exons) - 1):
                splice_junctions.append((exons[i][1], exons[i + 1][0]))
        else:
            for i in range(len(exons) - 1, 0, -1):
                splice_junctions.append((exons[i][0], exons[i - 1][1]))
        return splice_junctions

    def storeJunctions(junctions):
        for coords in junctions:
            uid = chr+':'+str(coords[0])+'-'+str(coords[1])
            if uid in unique_junction_db:
                unique_junction_db[uid].append(file)
            else:
                unique_junction_db[uid] = [file]

    def exon_str(exonIDs):
        if '_' in exonIDs[0]:
            del exonIDs[0]
        if '_' in exonIDs[-1]:
            del exonIDs[-1]
        filtered_exonIDs_str = "|".join(exonIDs)
        return filtered_exonIDs_str,exonIDs

    def process_isoform(chr, strand, info, exons, file):
        transcript_id = info.split(';')[ti].split(td)[1]

        # exonIDs_simple requries further checking as justification for its structure is inconsistent with exonIDs
        gene, exonIDs, exonIDs_simple, genes = exonAnnotate(chr, exons, strand, transcript_id)
        # Restrict the exon string to high confidence exon boundaries
        filtered_exonIDs_str,filtered_exonIDs = exon_str(list(exonIDs))

        """
        if transcript_id == 'NKX2-5|1/3|01H01':
            print (filtered_exonIDs_str,transcript_id,mode);sys.exit()
        """
        if len(genes)>2:
            trans_spliced_isoforms[filtered_exonIDs_str] = []

        if mode == 'collapse':
            eo.write('\t'.join([gene, strand, filtered_exonIDs_str, transcript_id, file]) + '\n')
        else:
            eo.write('\t'.join([gene, strand, "|".join(exonIDs), transcript_id, file]) + '\n')
        #splice_junctions = getJunctions(exons)
        if 'UNK' not in gene:
            #junction_str_db[gene,filtered_exonIDs_str] = tuple(splice_junctions)
            junction_key = filtered_exonIDs_str
            if str(mode).lower().startswith('cluster'):
                junction_key = filtered_exonIDs_str
            try: 
                junction_db[(gene, junction_key)].append((file, info))
            except:
                junction_db[(gene, junction_key)] = [(file, info)]
            try:
                if junction_key not in gene_db[gene]:
                    gene_db[gene].append(junction_key)
            except:
                    gene_db[gene] = [junction_key]
            strand_db[gene] = strand

    gff_organization={}
    for file in files:
        if os.path.exists(file):
            fn = file
            file = os.path.basename(file)
        else:
            fn = os.path.join(directory, file)
        file = file.split('.g')[0]
        firstRow = True
        exons = []
        isoforms = 2
        open_func = gzip.open if fn.endswith('.gz') else open
        with open_func(fn, 'rt') as filepath:
            for line in filepath:
                if line[0] == '#':
                    continue
                data = line.strip()
                t = data.split('\t')
                try: chr, null, type, pos1, pos2, null, strand, null, info = t
                except: 
                    #print (info)
                    continue
                pos1 = int(pos1)
                pos2 = int(pos2)
                if 'chr' not in chr:
                    chr = 'chr'+chr

                if firstRow:
                    # Determine if the transcript is listed first or second
                    if type != 'transcript' and type != 'exon':
                        continue
                    ti = next((i for i, x in enumerate(info.split(';')) if 'transcript_id' in x), -1)
                    if '=' in info:
                        td = '='
                    else:
                        td = '"'
                    gff_organization[file]=ti,td
                    firstRow = False
                else:
                    #print(file,type,chr,pos1,pos2,info);sys.exit()
                    if type == 'transcript':
                        isoforms += 1
                        chr, strand, info = gene_info
                        process_isoform(chr, strand, info, exons, file)
                        exons = []
                    elif type != 'transcript' and type != 'exon':
                        pass
                    else:
                        exons.append((pos1, pos2))
                        gene_info = chr, strand, info
        # for the last isoform in the file
        chr, strand, info = gene_info
        try: 
            process_isoform(chr, strand, info, exons, file)
        except:
            continue
        print(file, '...', isoforms, 'isoforms')

    eo.close()
    print(len(junction_db), 'unique isoforms')

    def sort_isoforms_with_ENST_first(junction_db):
        sorted_junction_db = defaultdict(list)
        for key, value in junction_db.items():
            sorted_value = sorted(value, key=lambda x: 'ENST' not in x[1])
            sorted_junction_db[key] = sorted_value
        return sorted_junction_db

    junction_db = sort_isoforms_with_ENST_first(junction_db)

    if collapse_isoforms == False:
        return transcript_associations
    else:
        """ Most isoforms should be redundant between samples - collapse isoforms based on redundant junctions
        and export a combined minimal GFF and gff-specific associated isoforms to the longest exemplar """ 

        # test case for collapseIsoforms: should give
        # {((1, 2), (3, 4), (5, 6), (7, 8)): [((1, 2), (3, 4), (5, 6))], ((1, 2), (3, 4), (7, 8)): [((1, 2), (3, 4), (7, 8)), ((3, 4), (7, 8))]}
        #a = {'gene1':[((1,2),(3,4),(5,6)),((1,2),(3,4),(5,6),(7,8)),((1,2),(3,4),(7,8)),((3,4),(7,8)),((1,2),(3,4),(7,8))]}
        
        if str(mode).lower().startswith('cluster'):
            super_isoform_db = collapseIsoformsCluster(
                gene_db,
                junction_db,
                gff_organization,
                mode,
            )
        elif mode == 'collapse' or mode == 'Ensembl':
            super_isoform_db = collapseIsoforms(gene_db, junction_db, gff_organization, mode)
        else:
            super_isoform_db = defaultdict(lambda: defaultdict(list))
            for (gene, isoform) in junction_db: 
                super_isoform_db[gene][isoform] = []
        
        isoform_transcript_map = {}
        isoform_pair_map = {}
        for (gene, isoform), entries in junction_db.items():
            transcript_ids = {}
            seen_pairs = set()
            pairs = []
            for (file, info) in entries:
                ti, td = gff_organization[file]
                transcript_id = info.split(';')[ti].split(td)[1]
                transcript_ids[transcript_id] = file
                pair = (transcript_id, file)
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    pairs.append(pair)
            isoform_transcript_map[(gene, isoform)] = transcript_ids
            isoform_pair_map[(gene, isoform)] = pairs

        #print (super_isoform_db['ENSG00000183072']['E1.1|E1.2|E1.3|E1.4|E3.3|E3.4|E3.5|'])

        # Export super- to sub-isoform associations:
        eo = open(os.path.join(combined_dir, 'isoform_links.txt'), 'w')
        #jo = open(os.path.join(combined_dir, 'isoform_junctions.txt'), 'w')
        ao = open(os.path.join(combined_dir, 'isoform_annotations.txt'), 'w')
        isoforms_to_retain = {}
        a=0
        total_genes = len(super_isoform_db)
        # Use tqdm to create a progress bar
        for gene in tqdm(super_isoform_db, total=total_genes, desc="Processing genes"):
            added = set()
            for isoform in super_isoform_db[gene]: 
                # isoform is a tuple of junction coordinates tuples
                strand = strand_db[gene]
                related_transcripts = []
                if 'UNK' not in gene:
                    transcript_ids = isoform_transcript_map.get((gene, isoform), {})
                    if not transcript_ids:
                        continue
                    transcript_list = list(transcript_ids.keys())
                    try:
                        known_candidate = None
                        for t_id in transcript_list:
                            if t_id.startswith(('ENST', 'NM_', 'XM_', 'XR_')):
                                if known_candidate is None or t_id < known_candidate:
                                    known_candidate = t_id
                        ref_super_transcript_id = known_candidate or transcript_list[0]
                    except:
                        print (transcript_ids)
                        print (gene,[isoform])
                        print ('error...');sys.exit()
                    ref_super_file = transcript_ids[ref_super_transcript_id]
                    isoforms_to_retain[(ref_super_file, ref_super_transcript_id)] = [] 
                    related_transcripts.append(ref_super_transcript_id)
                    known_flag = ref_super_transcript_id.startswith(('ENST', 'NM_', 'XM_', 'XR_'))
                    if len(transcript_list)==1:
                        eo.write('\t'.join([gene, ref_super_transcript_id, ref_super_file, '', '']) + '\n')
                        a += 1
                    else: 
                        for t in transcript_list:
                            if t!=ref_super_transcript_id:
                                t_file =  transcript_ids[t]
                                eo.write('\t'.join([gene, ref_super_transcript_id, ref_super_file, t, t_file]) + '\n')
                                if not known_flag and t.startswith(('ENST', 'NM_', 'XM_', 'XR_')):
                                    known_flag = True
                    for sub_isoform in super_isoform_db[gene][isoform]:
                        for (sub_transcript_id, f3) in isoform_pair_map.get((gene, sub_isoform), []):
                            sub_key = (sub_transcript_id, f3)
                            if sub_key not in added:
                                if ref_super_transcript_id !=sub_transcript_id:
                                    eo.write('\t'.join([gene, ref_super_transcript_id, ref_super_file, sub_transcript_id, f3]) + '\n')
                                    added.add(sub_key)
                                    related_transcripts.append(sub_transcript_id)
                                    a += 1
                                    if not known_flag and sub_transcript_id.startswith(('ENST', 'NM_', 'XM_', 'XR_')):
                                        known_flag = True
                    known_isoform = 'known' if known_flag else 'novel'
                    ao.write('\t'.join([gene, ref_super_transcript_id, ref_super_file, ','.join(related_transcripts), known_isoform]) + '\n')

        eo.close()
        ao.close()
        #jo.close()
        print(a, 'isoform pairs')
        print(len(trans_spliced_isoforms), '>2 trans-spliced isoforms')

        # Iterate back through the original GFF files
        combined_gff = os.path.join(combined_dir, 'combined.gff')
        eo = open(combined_gff, 'w')
        for file in files:
            if os.path.exists(file):
                fn = file
                file = os.path.basename(file)
            else:
                fn = os.path.join(directory, file)
            file = file.split('.g')[0]
            firstRow = True
            exons = []
            isoforms = 2
            open_func = gzip.open if fn.endswith('.gz') else open
            with open_func(fn, 'rt') as filepath:
                for line in filepath:
                    data = line.strip()
                    t = data.split('\t')
                    if line[0] =='#':
                        continue
                    try: 
                        chr, a, type, pos1, pos2, b, strand, c, info = t
                    except:
                        continue
                    if 'CDS' != type and 'gene' != type:
                        ti = next((i for i, x in enumerate(info.split(';')) if 'transcript_id' in x), -1)
                        if '=' in info:
                            td = '='
                        else:
                            td = '"'
                        try:
                            transcript_id = info.split(';')[ti].split(td)[1]
                            if (file, transcript_id) in isoforms_to_retain:
                                line = line.replace("PB.", file + '_PB.')
                                eo.write(line)
                        except:
                            pass
        eo.close()
        sorted_gff = combined_gff[:-4]+'-sorted.gff'
        try:
            sort_gff(combined_gff, sorted_gff)
        except:
            print ('sort error')
        return combined_dir

def importEnsemblGenes(exon_file,include_introns=False):
    global exonCoordinates
    global geneData
    global exonData
    exonCoordinates = {}
    geneData = {}
    gene_chr = {}
    exonData = {}
    strandData = {}
    firstRow = True
    prior_gene = None
    with open(exon_file, 'r') as file:
        for line in file:
            data = line.strip()
            t = data.split('\t')
            if firstRow:
                firstRow = False
            else:
                gene, exon, chr, strand, start, stop = t[:6]
                gene_chr[gene] = chr
                start = int(start)
                stop = int(stop)
                if strand == '-':
                    start, stop = stop, start
                exon_info = (start, stop, exon)
                if gene in geneData:
                    geneData[gene].append(exon_info)
                else:
                    geneData[gene] = [exon_info]
                exonData[gene,exon] = start,stop
                strandData[gene] = strand
                prior_gene=gene
    for gene in geneData:
        geneData[gene].sort()
        if strandData[gene] == '-':
            geneData[gene].reverse()
        index=0
        chr = gene_chr[gene]
        for (start, stop, exon) in geneData[gene]:
            if 'E' in exon:
                exonCoordinates[(chr, start, strandData[gene], 1)] = (gene, exon, index)
                exonCoordinates[(chr, stop, strandData[gene], 2)] = (gene, exon, index)
            elif include_introns:
                exonCoordinates[(chr, start, strandData[gene], 1)] = (gene, exon, index)
                exonCoordinates[(chr, stop, strandData[gene], 2)] = (gene, exon, index)
            index+=1
    return exonCoordinates, geneData, strandData
    
def sort_gff(concatenated_gff_file, sorted_gff_file):
    def chromosome_sort_key(chromosome):
        chromosome = chromosome.replace('chr', '').replace('Chr', '')
        if chromosome.isdigit():
            return (0, int(chromosome))
        elif chromosome.upper() == 'Y':
            return (1, 0)
        elif chromosome.upper() == 'X':
            return (1, 1)
        else:
            return (2, chromosome)

    # Read the GFF file and store features with their associated transcript start position
    header_lines = []
    features = []
    transcript_starts = {}
    with open(concatenated_gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                columns = line.strip().split('\t')
                attributes = {attr.split()[0]: attr.split()[1].strip('"') for attr in columns[8].split(';') if attr}
                transcript_id = attributes.get('transcript_id', None)
                if transcript_id:
                    if transcript_id not in transcript_starts:
                        transcript_starts[transcript_id] = int(columns[3])
                    start_position = transcript_starts[transcript_id]
                    features.append((columns[0], start_position, line))  # (chromosome, transcript start, line)
                else:
                    features.append((columns[0], int(columns[3]), line))  # (chromosome, feature start, line)

    # Sort features by chromosome and transcript start position (or feature start position for non-transcript lines)
    features.sort(key=lambda x: (chromosome_sort_key(x[0]), x[1]))

    # Write the sorted GFF file with headers followed by sorted features
    with open(sorted_gff_file, 'w') as sorted_gff:
        sorted_gff.writelines(header_lines)
        for _, _, line in features:
            sorted_gff.write(line)

def convert_alt_chromosome_ids():
    # Load the mapping file from UCSC Genome Browser
    mapping_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.ncbiRefSeq.gtf.gz"
    mapping_df = pd.read_csv(mapping_url, sep='\t', header=None, comment='#')

    # Extract RefSeq ID to chromosome mapping
    mapping = dict(zip(mapping_df[1], mapping_df[0]))

    def convert_gff(input_gff, output_gff):
        with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                    continue
                parts = line.split('\t')
                if parts[0] in mapping:
                    parts[0] = mapping[parts[0]]
                outfile.write('\t'.join(parts))

def reformat_uncoventional_gtf(gtf_file_path,output_file_path):

    # Read the GTF file into a pandas DataFrame
    gtf_df = pd.read_csv(gtf_file_path, sep='\t', header=None, comment='#',
                        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # Replace 'CDS' with 'exon' in the 'feature' column
    gtf_df['feature'] = gtf_df['feature'].replace('CDS', 'exon')

    # Extract transcript IDs
    gtf_df['transcript_id'] = gtf_df['attribute'].str.extract('transcript_id "([^"]+)"')

    # Identify unique transcripts and their first and last exon positions
    transcript_info = gtf_df.groupby('transcript_id').agg({
        'seqname': 'first',
        'source': 'first',
        'strand': 'first',
        'start': 'min',
        'end': 'max',
        'attribute': 'first'
    }).reset_index()

    # Create transcript entries
    transcript_entries = pd.DataFrame({
        'seqname': transcript_info['seqname'],
        'source': transcript_info['source'],
        'feature': 'transcript',
        'start': transcript_info['start'],
        'end': transcript_info['end'],
        'score': '.',
        'strand': transcript_info['strand'],
        'frame': '.',
        'attribute': transcript_info['attribute'].str.replace(r'(\|[^"]+)', '', regex=True)
    })

    # Add transcript_id to transcript_entries to merge with exons
    transcript_entries['transcript_id'] = transcript_info['transcript_id']

    # Initialize a list to hold final entries
    final_entries = []

    # Iterate over each unique transcript ID to gather transcript and exon entries
    for transcript_id in transcript_info['transcript_id']:
        transcript_entry = transcript_entries[transcript_entries['transcript_id'] == transcript_id]
        exon_entries = gtf_df[gtf_df['transcript_id'] == transcript_id]
        final_entries.append(transcript_entry)
        final_entries.append(exon_entries)

    # Combine all entries into a single DataFrame
    final_df = pd.concat(final_entries).reset_index(drop=True)

    # Drop the temporary 'transcript_id' column
    final_df = final_df.drop(columns=['transcript_id'])

    # Save the reformatted DataFrame to a new GTF file
    final_df.to_csv(output_file_path, sep='\t', header=False, index=False, quoting=3)  # quoting=3 is for quote None

    print(f"Reformatted GTF file with transcript entries saved to {output_file_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Combine long-read GFF/GTF files and annotate isoforms.'
    )
    parser.add_argument(
        '--gff-input',
        dest='gff_input',
        required=True,
        help='GFF/GTF file or directory containing GFF/GTF files.'
    )
    parser.add_argument(
        '--gene-model',
        dest='gene_model',
        required=True,
        help='Ensembl exon reference file.'
    )
    parser.add_argument(
        '--mode',
        default='collapse',
        help='Isoform consolidation mode (collapse, Ensembl, cluster, or any other value for raw).'
    )
    args = parser.parse_args()

    consolidateLongReadGFFs(args.gff_input, args.gene_model, mode=args.mode)
