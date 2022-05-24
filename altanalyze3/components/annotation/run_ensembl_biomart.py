
import requests
# import pandas as pd

import xmltodict
from xml.etree import ElementTree
from io import StringIO
"""
This is a generalized python module for getting data from Ensemble using Biomart server.
"""

DEFAULT_HOST = 'http://www.biomart.org'
DEFAULT_PATH = '/biomart/martservice'
DEFAULT_PORT = 80
DEFAULT_SCHEMA = 'default'
writePath = '/Users/sin9gp/altanalyze3/src/annotation'

def query(
                attributes=None,
                only_unique=True,
                use_attr_names=False,
                dtypes = None,
                name = "hsapiens_gene_ensembl"
                ):
            """Queries the dataset to retrieve the contained data.
            Args:
                attributes (list[str]): Names of attributes to fetch in query.
                    Attribute names must correspond to valid attributes. See
                    the attributes property for a list of valid attributes.

                only_unique (bool): Whether to return only rows containing
                    unique values (True) or to include duplicate rows (False).
                use_attr_names (bool): Whether to use the attribute names
                    as column names in the result (True) or the attribute
                    display names (False).
                dtypes (dict[str,any]): Dictionary of attributes --> data types
                    to describe to pandas how the columns should be handled
            Returns:
                pandas.DataFrame: DataFrame containing the query results.
            """

            root = ElementTree.Element('Query')
            root.set('virtualSchemaName', 'default')
            root.set('formatter', 'TSV')
            root.set('header', '1')
            root.set('uniqueRows', str(int(only_unique)))
            root.set('datasetConfigVersion', '0.6')

            # Add dataset element.
            dataset = ElementTree.SubElement(root, 'Dataset')
            dataset.set('name', name)
            dataset.set('interface', 'default')

            # Default to default attributes if none requested.
            # if attributes is None:
            #     attributes = list(default_attributes.keys())
            attributes = ["ensembl_exon_id", "ensembl_peptide_id", "transcript_start","transcript_end","interpro","interpro_short_description","interpro_start","interpro_end"]
            # Add attribute elements.
            for name in range(0,len(attributes)-1):
                try:
                    attr = attributes[name]
                    add_attr_node(dataset, attr)
                except KeyError:
                    raise Exception(
                        'Unknown attribute {}, check dataset attributes '
                        'for a list of valid attributes.'.format(name))
            
            # Fetch response.
            response = get(query=ElementTree.tostring(root))

            # Raise exception if an error occurred.
            if 'Query ERROR' in response.text:
                raise Exception(response.text)

            # Parse results into a DataFrame.
            try:
                with open(writePath, 'a') as f:
                    result = pd.read_csv(StringIO(response.text), sep='\t', dtype=dtypes)
                    f.write(result)

                result = pd.read_csv(StringIO(response.text), sep='\t', dtype=dtypes)
            # Type error is raised of a data type is not understood by pandas
            except TypeError as err:
                raise ValueError("Non valid data type is used in dtypes")

            if use_attr_names:
                # Rename columns with attribute names instead of display names.
                column_map = {
                    attributes[attr].display_name: attr
                    for attr in attributes
                }
                result.rename(columns=column_map, inplace=True)

            return result


def add_attr_node(root, attr):
        attr_el = ElementTree.SubElement(root, 'Attribute')
        attr_el.set('name', attr)

@property
def url(self):
    """Url used to connect to the biomart service."""
    return '{}:{}{}'.format(DEFAULT_HOST, DEFAULT_PORT, DEFAULT_PATH)

def get(**params):
            """Performs get request to the biomart service.
            Args:
                **params (dict of str: any): Arbitrary keyword arguments, which
                    are added as parameters to the get request to biomart.
            Returns:
                requests.models.Response: Response from biomart for the request.
            """
          
            r = requests.get(url, params=params)
            print("biomart request is made", r) 
            r.raise_for_status()
       


query(attributes=["ensembl_exon_id", "ensembl_peptide_id", "transcript_start","transcript_end","interpro","interpro_short_description","interpro_start","interpro_end"])
              