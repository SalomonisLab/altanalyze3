import re,sys

def classify_splicing_event(event_id, strand, transcript1, transcript2, coord1, coord2):
    """
    Determines the type of alternative mRNA processing event based on AltAnalyze format.
    Secondarily determines whether the splice/APA site is more proximal (e.g., acceptor vs. donor) 
    or distal in the first junction vs. second. For cassette, ME and RI, inclusion vs. exclusion
    is reporter instead for proximity.

    Parameters:
    - event_id: The event ID containing two splice junctions
    - strand: String indicating forward or reverse strand (+,-)
    - transcript1: A string of exon/intron segments separated by '|'
    - transcript2: A string of exon/intron segments separated by '|'

    Returns:
    - Event type as a string
    - Proximity as a string
    """

    proximity = 'null'
    # Convert transcripts to lists
    if isinstance(transcript1, list):
        pass
    else:
        transcript1 = transcript1.split("|")
        transcript2 = transcript2.split("|")

    # Extract junctions from event ID
    junctions = event_id.split("|")
    junction1, junction2 = junctions[0].split(":")[1], junctions[1].split(":")[1]
    coords1 = list(map(int, coord1.split(":")[1].split("-")))
    coords2 = list(map(int, coord2.split(":")[1].split("-")))
    j1s,j1e = junction1.split("-")
    j2s,j2e = junction2.split("-")
    if strand == "-":
        coords1 = [-x for x in coords1]
        coords2 = [-x for x in coords2]

    # ---- Trans-splicing (two Ensembl genes in one junction) ----
    if junctions[0].count("ENS") > 1 or junctions[1].count("ENS") > 1:
        if junctions[0].count("ENS") > 1 and junctions[1].count("ENS") == 1:
            proximity = 'distal'
        elif junctions[0].count("ENS") == 1 and junctions[1].count("ENS") > 1:
            proximity = 'proximal'
        return "Trans-Splicing", proximity

    # ---- Retained Intron (IR) ----
    sites = [j1s,j1e,j2s,j2e]
    if any("I" in site and "_" not in site for site in sites):
        if "I" in j1s or "I" in j1e:
            proximity = 'inclusion'
        else:
            proximity = 'exclusion'
        return "Intron Retention", proximity
        
    # ---- Alternative Polyadenylation (AltPolyA) ----
    if j1e == transcript1[-1] and j2e == transcript2[-1]: # match to the last transcript exons
        diff = coords1[-1]-coords2[-1]
        if j1s.split(".")[0] == j2s.split(".")[0] and j1e.split(".")[0] == j2e.split(".")[0]: # same exon block
            if diff>0:
                proximity = 'distal'
            else:
                proximity = 'proximal'
            return "AltPolyA", proximity
        elif j1e.split(".")[0] != j2e.split(".")[0]: # the last exon blocks are different
            if diff>0:
                proximity = 'distal'
            else:
                proximity = 'proximal'
            return "Alt C-Terminal Exon", proximity

    # ---- Alternative Promoter (AltPromoter) ----
    if j1s != j2s and (j1s == transcript1[0] or j2s == transcript2[0]): # match to the first transcript exon(s) - not required both match
        diff = coords1[0]-coords2[0]
        if diff>0:
            proximity = 'distal' # This interpretation is only correct when the exon occurs down-stream of the canonical exon 1. If an alternative promoter occurs upstream in the genome it would be considered proximal not distal - unclear which is proper 
        else:
            proximity = 'proximal'
        return "AltPromoter", proximity

    # ---- Alternative 3' or 5' Splice Sites (Alt3' SS / Alt5' SS) ----
    if j1s.split(".")[0] == j2s.split(".")[0] and j1e.split(".")[0] == j2e.split(".")[0]:
        if j1s == j2s:
            diff = coords1[-1]-coords2[-1]
        else:
            diff = coords1[0]-coords2[0]
        if diff>0:
            proximity = 'distal'
        else:
            proximity = 'proximal'
        if j1e.split(".")[1] != j2e.split(".")[1]:
            return "Alt3' SS", proximity
        else:
            return "Alt5' SS", proximity

    # ---- Cassette Exon ----

    if (j1s[1:].split('.')[0] != j2s[1:].split('.')[0]) and (j1e[1:].split('.')[0] == j2e[1:].split('.')[0]):
        # Different exon block detected in the first exon only
        diff = int(j1s[1:].split('.')[0]) - int(j2s[1:].split('.')[0])
    else:
        diff = int(j1e[1:].split('.')[0]) - int(j2e[1:].split('.')[0])

    j1s_ind = transcript1.index(j1s)
    j2s_ind = transcript2.index(j2s)
    j1e_ind = transcript1.index(j1e)
    j2e_ind = transcript2.index(j2e)

    try:
        if diff>0:
            proximity = 'distal'
        else:
            proximity = 'proximal'
        if j1s == j2s:
            if transcript1[j1e_ind+1] == transcript2[j2e_ind+1]: # same distal exon
                return "Mutually-Exclusive Exon", proximity

        if j1e == j2e:
            if transcript1[j1s_ind-1] == transcript2[j2s_ind-1]: # same upstream exon
                return "Mutually-Exclusive Exon", proximity
    except:
        pass

    if diff>0:
        proximity = 'inclusion'
    else:
        proximity = 'exclusion'
    return "Cassette Exon", proximity


if __name__ == '__main__':

    # ---- Test Use cases ----
    examples = [
        # Trans-Splicing
        ("ENSG00000015479:E11.2-E11.3_139322005|ENSG00000015479:E1.14-ENSG00000120727:E3.1", "+",
        "E1.6|E1.7|E1.8|E1.9|E1.10|E1.11|E1.12|E1.13|E1.14|E3.1|E3.2|E3.3|E3.4|E3.5|E3.6|E3.7|E3.8|E3.9|E3.10|E4.2|E5.2|E6.2|E7.1|E8.1|E8.2|E8.3|E8.4|E9.1|E9.2|E9.3|E10.1|E10.2|E11.2|E11.3_139322005",
        "E1.6|E1.7|E1.8|E1.9|E1.10|E1.11|E1.12|E1.13|E1.14|ENSG00000120727:E3.1|ENSG00000120727:E4.1_139364568|ENSG00000120727:E4.2|ENSG00000120727:E5.1|ENSG00000120727:E5.2|ENSG00000120727:E5.3_139368959",
        "chr5:139321994-139322005","chr5:139293805-139363759"),

        # Alternative Polyadenylation (AltPolyA)
        ("ENSG00000055130:E8.1-E8.1_148759573|ENSG00000055130:E8.1-E8.1_148759575", "+",
        "E2.1_148698757|E2.1_148698812|E4.2_148730000|E4.2|E5.1|E5.2|E5.3|E6.1|E7.1|E8.1|E8.1_148759573",
        "E2.1_148698783|E2.1|I2.1|E3.1|E4.1|E4.2|E5.1|E5.2|E5.3|E6.1|E8.1|E8.1_148759575",
        "chr7:148759638-148759573","chr7:148759638-148759575"),

        # Alternative Promoter (AltPromoter)
        ("ENSG00000000457:E1.1-E3.1|ENSG00000000457:E2.2-E3.1", "-",
        "E1.1|E3.1|E3.2|E4.1|E5.1|E6.1|E7.1|E8.2|E9.1|E10.1|E12.1|E13.1|E13.2|E16.1",
        "E2.1|E2.2|E3.1|E3.2|E4.1|E5.1|E6.1|E7.1|E8.2|E9.1|E10.1|E12.1|E13.1|E13.2|E15.1|E16.1|E16.2|E17.1",
        "chr1:169894007-169888890","chr1:169893788-169888890"),

        # Alternative 3' Splice Site (Alt3' SS)
        ("ENSG00000197555:E63.1-E64.1|ENSG00000197555:E63.1-E64.2", "+",
        "E31.2|E45.1|E45.2|E45.3|E45.4|E45.5|E47.2|E48.1|E49.2|E50.1|E51.1|E51.2|E52.1|E53.1|E54.1|E55.1|E56.2|E57.1|E57.3|E58.1|E58.2|E59.1|E60.1|E61.1|E61.2|E62.1|E63.1|E64.1|E64.2|E65.1|E65.2|E65.3|E65.4|E65.5",
        "E44.1|E45.1|E45.2|E45.3|E45.4|E45.5|E47.2|E48.1|E49.2|E50.1|E51.1|E51.2|E52.1|E53.1|E54.1|E55.1|E56.2|E57.3|E58.1|E58.2|E59.1|E60.1|E61.1|E61.2|E62.1|E63.1|E64.2|E65.1|E65.2|E65.3",
        "chr14:71735391-71738241","chr14:71735391-71738244"),

        # Retained Intron (RI)
        ("ENSG00000244952:E1.1-E2.1|ENSG00000244952:E1.1-I1.1", "+",
        "E1.1|E2.1|E3.1|E4.1|E5.1",
        "E1.1|I1.1|E3.1|E4.1|E5.1",
        "chr15:32613998-32614109","chr15:32613998-32613999"),

        # Cassette Exon
        ("ENSG00000075539:E91.1-E93.1|ENSG00000075539:E90.2-E93.1", "-",
        "E63.1_48547613|E63.1|E63.2|I63.1|E64.1|E64.2|I64.1|E65.1|E65.2|E65.3|E67.1|E69.1|E71.1|E71.2|E72.1|E72.2|E73.1|E76.1|E77.1|E79.1|E79.2|E83.1|E84.1|I84.1|E85.1|E88.3|E88.4|E90.1|E90.2|I90.1|E91.1|E93.1|E93.2_48515071",
        "E2.2_48780279|E2.2|E2.3|E2.4|E2.5|E2.6|E5.1|E7.1|E14.1|E14.2|E14.3|E16.1|E18.2|E18.3|E19.1|E20.1|E21.1|E22.1|E23.2|E24.1|E25.1|E26.1|E27.1|E28.1|E29.1|E30.1|E31.1|E32.1|E33.1|E34.1|E35.1|E36.1|E37.1|E38.1|E40.1|E41.1|E44.1|E45.1|E45.2|E45.3|E46.1|E48.1|E49.1|E52.1|E52.2|E53.1|E56.1|E56.2|E58.3|E58.4|E60.1|E61.2|E61.3|E62.1|E62.2|E63.1|E65.2|E65.3|E67.1|E69.1|E71.1|E71.2|E72.1|E72.2|E73.1|E76.1|E77.1|E79.1|E79.2|E83.1|E84.1|E85.1|E88.3|E88.4|E90.1|E90.2|E93.1|E93.2|E96.1|E98.1|E99.1|E101.2|E104.1|E106.1|E106.2|E107.1|E108.1|E108.2|E108.3|E108.4|E108.5|E108.6|U109.1_48497357",
        "chr4:48520524-48515275","chr4:48521048-48515275"),

        # Mutually Exclusive Exons
        ("ENSG00000114861:E55.1-E56.1|ENSG00000114861:E54.1-E56.1", "-",
        "E17.1|E17.2|E24.1|E32.1|E32.2|E33.1|E35.1|E36.1|E45.1|E49.1|E50.1|E51.1|E52.1|E52.2|E53.1|E55.1|E56.1|E58.2|E59.1|E59.2|E59.3|E59.4",
        "E1.1_71583978|E1.1|E1.2|E3.2|E5.1|E7.1|E12.1|E12.2|E17.1|E17.2|E24.1|E32.1|E32.2|E33.1|E35.1|E36.1|E45.1|E49.1|E50.1|E51.1|E52.1|E52.2|E53.1|E54.1|E56.1|E58.2|E59.1|E59.2|E59.3|E59.4|E59.5|E59.6|E59.7",
        "chr3:70972011-70970805","chr3:70972555-70970805"),

        # Alt-C terminus
        ("ENSG00000129450:E6.1-E8.1|ENSG00000129450:E6.1-E9.1", "+",
        "U1.1_51124906|E1.1|E2.1|E3.1|E4.1|E5.1|E6.1|E8.1",
        "E1.1|E2.1|E3.1|E4.1|E5.1|E6.1|E9.1",
        "chr19:51128510-51129891","chr19:51128510-51135962"),

        # Alternative 5' Splice Site (Alt5' SS)
        ("ENSG00000000457:E10.1-E12.1|ENSG00000000457:E10.1_169864458-E12.1", "-",
        "E1.1|E3.1|E3.2|E4.1|E5.1|E6.1|E7.1|E8.2|E9.1|E10.1|E12.1|E13.1|E13.2|E16.1",
        "E7.1_169870345|E7.1|E8.2|E9.1|E10.1|E10.1_169864458|E12.1|E13.1|E13.2|E15.1|E16.1|E16.2|I16.1|I16.1_169854048",
        "chr1:169864369-169862797","chr1:169864458-169862797"),

    ]

    # ---- Run classification and print results ----
    for event_id, strand, t1, t2, c1, c2 in examples:
        print(f"{event_id}: {classify_splicing_event(event_id, strand, t1, t2, c1, c2)}")
