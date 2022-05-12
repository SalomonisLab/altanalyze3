
import requests
import xml.etree.ElementTree as ET


def getProteinCoordinates():
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "chromosome_name" value = "6"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_peptide_id" />
		<Attribute name = "ensembl_exon_id" />
		<Attribute name = "genomic_coding_start" />
		<Attribute name = "genomic_coding_end" />
		<Attribute name = "cds_start" />
		<Attribute name = "cds_end" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "exon_chrom_start" />
		<Attribute name = "exon_chrom_end" />
	</Dataset>
</Query>
    """
    headers = {'Content-Type': 'application/xml'} # set what your server accepts

    response = requests.get('http://www.ensembl.org/biomart/martview/231b5fcfe421fcfb07745b3b39ef18d0', data=xml, headers=headers).text
    with open("response.txt", "w") as f:
        f.write(response)

getProteinCoordinates()