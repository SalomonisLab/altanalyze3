""" This module combines information across GFF files to compare isoform structures, 
annotate exons/transcripts and to nominate the longest exemplar isoforms"""

import os,sys
from collections import defaultdict

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
        if strand == '-':
            genes.reverse()
        if len(genes) == 3:
            #print('Multi-Trans-splicing', genes)
            pass
        gene = genes[0]
    else:
        gene = genes[0]

    ### Annotate novel exon/intron IDs and include in-between regions in the final isoform annotation
    if len(genes) == 0:
        output_list = [f"{start}-{end}" for start, end in exons]
        associations_updated = output_list  # Novel genes don't require exon annotations (likely update - currently position tuple)
    else:
        """
        if 'PB.106.2' == transcript_id:
            print(exons,associations);sys.exit()
        """
        associations_updated = []
        for ind, association in enumerate(associations):
            exonType = association[0]
            if exonType == 'common':
                gene1 = association[1][0]
                exonAnnotation = association[1][1]
                if gene1 != gene:
                    exonAnnotation = gene1 + ':' + exonAnnotation
                associations_updated.append(exonAnnotation) # denotes a splice site or end
            elif exonType == 'mixed':
                exonAnnotation1 = association[1][1]
                exonAnnotation2 = association[1][2]
                gene1 = association[1][0]
                start,stop = exons[ind]
                start_index = exonCoordinates[chr, start, strand, 1][-1]
                end_index = exonCoordinates[chr, stop, strand, 2][-1]+1
                ref_exons = geneData[gene]
                if gene1 != gene:
                    exonAnnotation1 = gene1 + ':' + exonAnnotation1
                associations_updated.append(exonAnnotation1)
                prior_splice_acceptor = None
                # Add exon and intron regions in the middle of a broader exon
                for in_between_region in ref_exons[start_index:end_index]:
                    associations_updated.append(in_between_region[-1])
                    current_spice_donor = in_between_region[0]
                    if prior_splice_acceptor != None:
                        if transcript_id not in additional_junctions:
                            additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                        else:
                            additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                    prior_splice_acceptor = in_between_region[-2]
                if gene1 != gene:
                    exonAnnotation2 = gene1 + ':' + exonAnnotation2
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
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
                    associations_updated.append(gene2 + ':' + exonAnnotation2)
                elif gene == gene2:
                    associations_updated.append(gene1 + ':' + exonAnnotation1)
                    associations_updated.append(exonAnnotation2) # denotes a splice site or end

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

                if gene != gene2 and 'U' not in exonAnnotation1: # trans-splicing somewhere else in this transcript
                    exonAnnotation1 = gene2+':'+exonAnnotation1

                associations_updated.append(exonAnnotation1)
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
                    ref_exons = geneData[gene]
                    prior_splice_acceptor = None
                    # Add exon and intron regions in the middle of a broader exon
                    for in_between_region in ref_exons[start_index:end_index]:
                        associations_updated.append(in_between_region[-1])
                        current_spice_donor = in_between_region[0]
                        if prior_splice_acceptor != None:
                            if transcript_id not in additional_junctions:
                                additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                            else:
                                additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                        prior_splice_acceptor = in_between_region[-2]

                if gene != gene2 and 'U' not in exonAnnotation2: # trans-splicing somewhere else in this transcript
                    exonAnnotation2 = gene2+':'+exonAnnotation2
                associations_updated.append(exonAnnotation2) # denotes a splice site or end

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

                if gene != gene1 and 'U' not in exonAnnotation1: # trans-splicing somewhere else in this transcript
                    exonAnnotation1 = gene1+':'+exonAnnotation1

                associations_updated.append(exonAnnotation1)

                # If exons and/or introns span from the begining of a known exon to the novel end, include all regions in between (e.g., E10.1 to U15.1_1234)
                # end_index will not exist in the database (novel end)
                start,stop = exons[ind]
                start_index = exonCoordinates[chr, start, strand, 1][-1]
                evaluate_between_regions = True
                if 'U' in exonAnnotation2:
                    end_index = None
                else:
                    try:
                        start,stop = exonData[geneID,exonAnnotation2.split("_")[0].replace('I','E')]
                        end_index = exonCoordinates[chr, stop, strand, 2][-1]
                    except:
                        # Intronic position is associated with a strange region, such as ENSG00000127483-I15.2
                        evaluate_between_regions = False

                if evaluate_between_regions:
                    ref_exons = geneData[gene]
                    prior_splice_acceptor = None
                    # Add exon and intron regions in the middle of a broader exon
                    for in_between_region in ref_exons[start_index:end_index]:
                        associations_updated.append(in_between_region[-1])
                        current_spice_donor = in_between_region[0]
                        if prior_splice_acceptor != None:
                            if transcript_id not in additional_junctions:
                                additional_junctions[transcript_id]=[(prior_splice_acceptor,current_spice_donor)]
                            else:
                                additional_junctions[transcript_id].append((prior_splice_acceptor,current_spice_donor))
                        prior_splice_acceptor = in_between_region[-2]

                if gene != geneID and 'U' not in exonAnnotation2: 
                    exonAnnotation2 = geneID+':'+exonAnnotation2
                associations_updated.append(exonAnnotation2) # denotes a splice site or end

            elif exonType == 'novel':
                novel_position1 = association[1][0]
                novel_position2 = association[1][1]
                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    exonAnnotation1 = findNovelSpliceSite(geneID, novel_position1, strand)
                    hits.append([exonAnnotation1,geneID])
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation1,geneID = hits[0]
                if gene != geneID and 'U' not in exonAnnotation1: 
                    exonAnnotation1 = geneID+':'+exonAnnotation1

                hits=[]
                for geneID in genes: ### more than one possible gene the novel exon is in
                    exonAnnotation2 = findNovelSpliceSite(geneID, novel_position2, strand)
                    hits.append([exonAnnotation2,geneID])
                hits.sort() # E and I annotations will be prioritized
                exonAnnotation2,geneID = hits[0]
                if gene != geneID and 'U' not in exonAnnotation2: 
                    exonAnnotation2 = geneID+':'+exonAnnotation2
                ### Need to fill in-between exons
                associations_updated.append(exonAnnotation1)
                associations_updated.append(exonAnnotation2) # denotes a splice site or end
    """
    if 'PB.145.97' == transcript_id:
        print(associations)
        print(associations_updated);sys.exit()
    """
    associations_updated = uniqueList(associations_updated)
    separator = '|'
    associations_updated = separator.join(associations_updated)

    return gene,associations_updated,genes

def uniqueList(ls):
    return [seen.add(x) or x for seen in (set(),) for x in ls if x not in seen]

def collapseIsoforms(gene_db):
    num_isoforms=0
    collaped_db = defaultdict(lambda: defaultdict(list))
    for gene, isoforms in gene_db.items():
        # Sort isoforms by length in descending order
        """
        if num_isoforms==17:
            print (gene,gene_db[gene])"""
        sorted_isoforms = sorted(isoforms, key=len, reverse=True)
        # Initially set every isoform as a key with an empty list
        for isoform in sorted_isoforms:
            collaped_db[gene][isoform] = []
        # Populate the list with sub-isoforms
        for i, isoform in enumerate(sorted_isoforms):
            for shorter_isoform in sorted_isoforms[i+1:]:
                if set(shorter_isoform).issubset(set(isoform)):
                    collaped_db[gene][isoform].append(shorter_isoform)
        # Remove super-isoforms that are present as sub-isoforms
        for super_isoform in list(collaped_db[gene].keys()):
            for sub_isoforms in collaped_db[gene].values():
                if super_isoform in sub_isoforms:
                    del collaped_db[gene][super_isoform]
                    break
        """
        if num_isoforms==17:
            print (collaped_db[gene])
            print (len(collaped_db[gene]),len(gene_db[gene]),num_isoforms)"""
        num_isoforms+= len(collaped_db[gene])
    print (num_isoforms,'collapsed unique isoforms')
    return collaped_db

def consolidateLongReadGFFs(directory, exon_reference_dir):
    """Only return isoforms with unique splice-site combinations"""
    import collections
    junction_db = collections.OrderedDict()
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
        directory = '' # present working directory
        if len(files)>1:
            collapse_isoforms = True
    elif '.g' in directory.lower():
        # Export exon structure information only
        gff = directory
        directory = os.path.dirname(gff)
        gff = os.path.basename(gff)
        combined_dir = os.path.join(directory, 'gff-output')
        os.makedirs(combined_dir, exist_ok=True)
        files = [gff]
    else:
        # Combine and report collapsed unique isoform structures 
        combined_dir = os.path.join(directory, 'gff-output')
        if not os.path.exists(combined_dir):
            os.makedirs(combined_dir)
        files = [file for file in os.listdir(directory) if file.endswith('.gff')]
        files += [file for file in os.listdir(directory) if file.endswith('.gtf')]
        #files = [file for file in os.listdir(directory) if file.endswith('.txt')]
        if len(files)>1:
            collapse_isoforms = True

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

    def process_isoform(chr, strand, info, exons, file):
        transcript_id = info.split(';')[ti].split('"')[1]
        gene, exonIDs, genes = exonAnnotate(chr, exons, strand, transcript_id)
        if len(genes)>2:
            trans_spliced_isoforms[exonIDs] = []
        eo.write('\t'.join([gene, strand, exonIDs, transcript_id, file[:-4]]) + '\n')
        splice_junctions = getJunctions(exons)
        if 'UNK' not in gene:
            try: 
                junction_db[(gene, tuple(splice_junctions))].append((file, info))
            except:
                junction_db[(gene, tuple(splice_junctions))] = [(file, info)]
            try:
                if splice_junctions not in gene_db[gene]:
                    gene_db[gene].append(tuple(splice_junctions))
            except:
                gene_db[gene] = [tuple(splice_junctions)]
            strand_db[gene] = strand

    gff_organization={}
    for file in files:
        #print (file)
        fn = os.path.join(directory, file)
        firstRow = True
        exons = []
        isoforms = 2
        with open(fn, 'r') as filepath:
            for line in filepath:
                data = line.strip()
                t = data.split('\t')
                chr, null, type, pos1, pos2, null, strand, null, info = t
                pos1 = int(pos1)
                pos2 = int(pos2)
                if 'chr' not in chr:
                    chr = 'chr'+chr

                if firstRow:
                    # Determine if the transcript is listed first or second
                    if 'transcript' in info.split(';')[0]: 
                        ti = 0
                    else:
                        ti = 1
                    gff_organization[file]=ti
                    firstRow = False
                else:
                    #print(file,type,chr,pos1,pos2,info);sys.exit()
                    if type == 'transcript':
                        isoforms += 1
                        chr, strand, info = gene_info
                        process_isoform(chr, strand, info, exons, file)
                        exons = []
                    elif type == 'CDS':
                        pass
                    else:
                        exons.append((pos1, pos2))
                        gene_info = chr, strand, info
        # for the last isoform in the file
        chr, strand, info = gene_info
        process_isoform(chr, strand, info, exons, file)
        print(file, '...', isoforms, 'isoforms')

    eo.close()
    print(len(junction_db), 'unique isoforms')

    if collapse_isoforms == False:
        return transcript_associations
    else:
        """ Most isoforms should be redundant between samples - collapse isoforms based on redundant junctions
        and export a combined minimal GFF and gff-specific associated isoforms to the longest exemplar """ 

        # test case for collapseIsoforms: should give
        # {((1, 2), (3, 4), (5, 6), (7, 8)): [((1, 2), (3, 4), (5, 6))], ((1, 2), (3, 4), (7, 8)): [((1, 2), (3, 4), (7, 8)), ((3, 4), (7, 8))]}
        #a = {'gene1':[((1,2),(3,4),(5,6)),((1,2),(3,4),(5,6),(7,8)),((1,2),(3,4),(7,8)),((3,4),(7,8)),((1,2),(3,4),(7,8))]}
        
        super_isoform_db = collapseIsoforms(gene_db)
        
        # Export super- to sub-isoform associations:
        eo = open(os.path.join(combined_dir, 'isoform_links.txt'), 'w')
        jo = open(os.path.join(combined_dir, 'isoform_junctions.txt'), 'w')
        isoforms_to_retain = {}
        a=0
        for gene in super_isoform_db:
            added=[]
            for isoform in super_isoform_db[gene]: 
                # isoform is a tuple of junction coordinates tuples
                strand = strand_db[gene]
                if 'UNK' not in gene:
                    (file,info) = junction_db[gene,isoform][0] # First example of that isoform
                    ref_super_transcript_id = info.split(';')[gff_organization[file]].split('"')[1]
                    isoforms_to_retain[(file,info)] = [] 
                    eo.write('\t'.join([gene, ref_super_transcript_id,file[:-4],'','']) + '\n')
                    a+=1
                    # Get the isoform annotations for the super-isoform
                    for (file, info) in junction_db[gene,isoform][1:]:
                        super_transcript_id = info.split(';')[gff_organization[file]].split('"')[1]
                        if (super_transcript_id,file) not in added:
                            eo.write('\t'.join([gene, ref_super_transcript_id,file[:-4],super_transcript_id,file[:-4]]) + '\n')
                            added.append(super_transcript_id,file)
                            a+=1
                    for sub_isoform in super_isoform_db[gene][isoform]:
                        for (f2, info) in junction_db[gene,sub_isoform]:
                            sub_transcript_id = info.split(';')[gff_organization[file]].split('"')[1]
                            if (sub_transcript_id,f2) not in added:
                                eo.write('\t'.join([gene, ref_super_transcript_id,file[:-4],sub_transcript_id,f2[:-4]]) + '\n')
                                added.append((sub_transcript_id,f2))
                                a+=1
                    # Create an inclusive exon-exon and exon-intron junction object - incorporate in-between-exons/introns
                    isoform_modified = list(isoform)
                    if strand == '-':
                        isoform_modified = [(b, a) for a, b in isoform_modified]
                        if ref_super_transcript_id in additional_junctions:
                            isoform_modified+=list(additional_junctions[ref_super_transcript_id])
                        isoform_modified.sort()
                        isoform_modified.reverse()
                    else:
                        if ref_super_transcript_id in additional_junctions:
                            isoform_modified+=list(additional_junctions[ref_super_transcript_id])
                        isoform_modified.sort()
                    isoform_string = '|'.join([f"{start}-{end}" for start, end in isoform_modified])
                    jo.write('\t'.join([gene, ref_super_transcript_id, file[:-4], isoform_string]) + '\n')

        eo.close()
        jo.close()
        print(a, 'isoform pairs')
        print(len(trans_spliced_isoforms), '>2 trans-spliced isoforms')

        # Iterate back through the original GFF files
        combined_gff = os.path.join(combined_dir, 'combined.gff')
        eo = open(combined_gff, 'w')
        for file in files:
            fn = os.path.join(directory, file)
            firstRow = True
            exons = []
            isoforms = 2
            with open(fn, 'r') as filepath:
                for line in filepath:
                    data = line.strip()
                    t = data.split('\t')
                    chr, a, type, pos1, pos2, b, strand, c, info = t
                    if (file, info) in isoforms_to_retain:
                        if 'CDS' not in line:
                            line = line.replace("PB.", file[:-4] + '_PB.')
                            eo.write(line)
        eo.close()
        sorted_gff = combined_gff[:-4]+'-sorted.gff'
        sort_gff(combined_gff, sorted_gff)
        return combined_dir

def importEnsemblGenes(exon_file):
    global exonCoordinates
    global geneData
    global exonData
    exonCoordinates = {}
    geneData = {}
    exonData = {}
    strandData = {}
    firstRow = True
    prior_gene = None
    index=0
    with open(exon_file, 'r') as file:
        for line in file:
            data = line.strip()
            t = data.split('\t')
            if firstRow:
                firstRow = False
            else:
                gene, exon, chr, strand, start, stop = t[:6]
                start = int(start)
                stop = int(stop)
                if strand == '-':
                    start, stop = stop, start
                if 'E' in exon:
                    exonCoordinates[(chr, start, strand, 1)] = (gene, exon, index)
                    exonCoordinates[(chr, stop, strand, 2)] = (gene, exon, index)
                exon_info = (start, stop, exon)
                if gene!=prior_gene:
                    index = 0
                if gene in geneData:
                    geneData[gene].append(exon_info)
                else:
                    geneData[gene] = [exon_info]
                exonData[gene,exon] = start,stop
                strandData[gene] = strand
                prior_gene=gene
                index+=1
    for gene in geneData:
        geneData[gene].sort()
        if strandData[gene] == '-':
            geneData[gene].reverse()

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

if __name__ == '__main__':
    #concatenated = '/Volumes/salomonis2/LabFiles/Frank-Li/neoantigen/revision/blood/Leucegene_AML/MDS-isoforms/head/combined-GFF-alt/combined.gff'
    #sorted_gff = '/Volumes/salomonis2/LabFiles/Frank-Li/neoantigen/revision/blood/Leucegene_AML/MDS-isoforms/head/combined-GFF-alt/sorted_gff_file.gff'
    #sort_gff(concatenated, sorted_gff);sys.exit()

    exon_reference_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'
    #exon_reference_dir = '/Volumes/salomonis2/software/AltAnalyze-91/AltAnalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'
    exon_reference_dir = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart100/ensembl/Mm/Mm_Ensembl_exon.txt'
    
    gff_input_dir = '/Volumes/salomonis2/LabFiles/Frank-Li/neoantigen/revision/blood/Leucegene_AML/MDS-isoforms/head/'
    gff_input_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory'
    gff_input_dir = '/Users/saljh8/Dropbox/Nanopore/Mouse'
    #gff_input_dir = '/Volumes/salomonis2/LabFiles/Frank-Li/neoantigen/revision/blood/Leucegene_AML/MDS-isoforms/'
    consolidateLongReadGFFs(gff_input_dir,exon_reference_dir)
