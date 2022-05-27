"""
This is a generalized python module for getting data from Ensemble using Biomart server.
"""

from __future__ import absolute_import, division, print_function
import requests

from future.utils import native_str
from builtins import *
from xml.etree import ElementTree
import pandas as pd
from io import StringIO
from xml.etree.ElementTree import fromstring as xml_from_string


DEFAULT_HOST = 'http://www.biomart.org'
DEFAULT_PATH = '/biomart/martservice'
DEFAULT_PORT = 80
DEFAULT_SCHEMA = 'default'


class ServerBase(object):
    """Base class that handles requests to the biomart server.
    Attributes:
        host (str): Host to connect to for the biomart service.
        path (str): Path to the biomart service on the host.
        port (str): Port to connect to on the host.
        url (str): Url used to connect to the biomart service.
        use_cache (bool): Whether to cache requests to biomart.
    """

    def __init__(self, host=None, path=None, port=None):
        """ServerBase constructor.
        Args:
            host (str): Url of host to connect to.
            path (str): Path on the host to access to the biomart service.
            port (int): Port to use for the connection.
            use_cache (bool): Whether to cache requests.
        """
        # Use defaults if arg is None.
        host = host or DEFAULT_HOST
        path = path or DEFAULT_PATH
        port = port or DEFAULT_PORT

        # Add http prefix and remove trailing slash.
        host = self._add_http_prefix(host)
        host = self._remove_trailing_slash(host)

        # Ensure path starts with slash.
        if not path.startswith('/'):
            path = '/' + path

        self._host = host
        self._path = path
        self._port = port

    @property
    def host(self):
        """Host to connect to for the biomart service."""
        return self._host

    @property
    def path(self):
        """Path to the biomart service on the host."""
        return self._path

    @property
    def port(self):
        """Port to connect to on the host."""
        return self._port

    @property
    def url(self):
        """Url used to connect to the biomart service."""
        return '{}:{}{}'.format(self._host, self._port, self._path)

    @staticmethod
    def _add_http_prefix(url, prefix='http://'):
        if not url.startswith('http://') or url.startswith('https://'):
            url = prefix + url
        return url

    @staticmethod
    def _remove_trailing_slash(url):
        if url.endswith('/'):
            url = url[:-1]
        return url

    def get(self, **params):
        """Performs get request to the biomart service.
        Args:
            **params (dict of str: any): Arbitrary keyword arguments, which
                are added as parameters to the get request to biomart.
        Returns:
            requests.models.Response: Response from biomart for the request.
        """

        r = requests.get(self.url, params=params)
        r.raise_for_status()
        return r


class BiomartException(Exception):
    """Basic exception class for biomart exceptions."""
    pass


class Dataset(ServerBase):
    """Class representing a biomart dataset.
    This class is responsible for handling queries to biomart
    datasets. Queries can select a subset of attributes and can be filtered
    using any available filters. A list of valid attributes is available in
    the attributes property. If no attributes are given, a set of default
    attributes is used. A list of valid filters is available in the filters
    property. The type of value that can be specified for a given filter
    depends on the filter as some filters accept single values, whilst others
    can take lists of values.
    Args:
        name (str): Id of the dataset.
        display_name (str): Display name of the dataset.
        host (str): Url of host to connect to.
        path (str): Path on the host to access to the biomart service.
        port (int): Port to use for the connection.
        use_cache (bool): Whether to cache requests.
        virtual_schema (str): The virtual schema of the dataset.
    Examples:
        Directly connecting to a dataset:
            >>> dataset = Dataset(name='hsapiens_gene_ensembl',
            >>>                   host='http://www.ensembl.org')
        Querying the dataset:
            >>> dataset.query(attributes=['ensembl_gene_id',
            >>>                           'external_gene_name'],
            >>>               filters={'chromosome_name': ['1','2']})
        Listing available attributes:
            >>> dataset.attributes
            >>> dataset.list_attributes()
        Listing available filters:
            >>> dataset.filters
            >>> dataset.list_filters()
    """

    def __init__(self,
                 name,
                 display_name='',
                 host=None,
                 path=None,
                 port=None,
                 virtual_schema=DEFAULT_SCHEMA):
        super().__init__(host=host, path=path, port=port)

        self._name = name
        self._display_name = display_name
        self._virtual_schema = virtual_schema
        self._filters = None
        self._attributes = None
        self._default_attributes = None
        self._datatype = None

    @property
    def name(self):
        """Name of the dataset (used as dataset id)."""
        return self._name

    @property
    def display_name(self):
        """Display name of the dataset."""
        return self._display_name

    @property
    def filters(self):
        """List of filters available for the dataset."""
        if self._filters is None:
            self._filters, self._attributes = self._fetch_configuration()
        return self._filters

    @property
    def attributes(self):
        """List of attributes available for the dataset (cached)."""
        if self._attributes is None:
            self._filters, self._attributes = self._fetch_configuration()
        return self._attributes

    @property
    def default_attributes(self):
        """List of default attributes for the dataset."""
        if self._default_attributes is None:
            self._default_attributes = {
                name: attr
                for name, attr in self.attributes.items()
                if attr.default is True
            }
        return self._default_attributes

    def list_attributes(self):
        """Lists available attributes in a readable DataFrame format.
        Returns:
            pd.DataFrame: Frame listing available attributes.
        """

        def _row_gen(attributes):
            for attr in attributes.values():
                yield (attr.name, attr.display_name, attr.description)

        return pd.DataFrame.from_records(
            _row_gen(self.attributes),
            columns=['name', 'display_name', 'description'])

    def list_filters(self):
        """Lists available filters in a readable DataFrame format.
        Returns:
            pd.DataFrame: Frame listing available filters.
        """

        def _row_gen(attributes):
            for attr in attributes.values():
                yield (attr.name, attr.type, attr.description)

        return pd.DataFrame.from_records(
            _row_gen(self.filters), columns=['name', 'type', 'description'])

    def _fetch_configuration(self):
        # Get datasets using biomart.
        response = self.get(type='configuration', dataset=self._name)

        # Check response for problems.
        if 'Problem retrieving configuration' in response.text:
            raise BiomartException('Failed to retrieve dataset configuration, '
                                   'check the dataset name and schema.')

        # Get filters and attributes from xml.
        xml = ElementTree.fromstring(response.content)

        filters = {f.name: f for f in self._filters_from_xml(xml)}
        attributes = {a.name: a for a in self._attributes_from_xml(xml)}

        return filters, attributes

    @staticmethod
    def _filters_from_xml(xml):
        for node in xml.iter('FilterDescription'):
            attrib = node.attrib
            yield Filter(
                name=attrib['internalName'], type=attrib.get('type', ''))

    @staticmethod
    def _attributes_from_xml(xml):
        for page_index, page in enumerate(xml.iter('AttributePage')):
            for desc in page.iter('AttributeDescription'):
                attrib = desc.attrib

                # Default attributes can only be from the first page.
                default = (page_index == 0 and
                           attrib.get('default', '') == 'true')

                yield Attribute(
                    name=attrib['internalName'],
                    display_name=attrib.get('displayName', ''),
                    description=attrib.get('description', ''),
                    default=default)

    def query(self,
              attributes=None,
              filters=None,
              only_unique=True,
              use_attr_names=False,
              dtypes=None,
              datatype=None
              ):
        """Queries the dataset to retrieve the contained data.
        Args:
            attributes (list[str]): Names of attributes to fetch in query.
                Attribute names must correspond to valid attributes. See
                the attributes property for a list of valid attributes.
            filters (dict[str,any]): Dictionary of filters --> values
                to filter the dataset by. Filter names and values must
                correspond to valid filters and filter values. See the
                filters property for a list of valid filters.
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

        # Example query from Ensembl biomart:
        #
        # <?xml version="1.0" encoding="UTF-8"?>
        # <!DOCTYPE Query>
        # <Query  virtualSchemaName = "default" formatter = "TSV" header = "0"
        #  uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
        #   <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        #       <Filter name = "chromosome_name" value = "1,2"/>
        #       <Filter name = "end" value = "10000000"/>
        #       <Filter name = "start" value = "1"/>
        #       <Attribute name = "ensembl_gene_id" />
        #       <Attribute name = "ensembl_transcript_id" />
        #   </Dataset>
        # </Query>

        # Setup query element.
        root = ElementTree.Element('Query')
        root.set('virtualSchemaName', self._virtual_schema)
        root.set('formatter', 'TSV')
        root.set('header', '1')
        root.set('uniqueRows', native_str(int(only_unique)))
        root.set('datasetConfigVersion', '0.6')

        # Add dataset element.
        dataset = ElementTree.SubElement(root, 'Dataset')
        dataset.set('name', self.name)
        dataset.set('interface', 'default')

        # Default to default attributes if none requested.
        if attributes is None:
            attributes = list(self.default_attributes.keys())

        # Add attribute elements.
        for name in attributes:
            try:
                attr = self.attributes[name]
                self._add_attr_node(dataset, attr)
            except KeyError:
                raise BiomartException(
                    'Unknown attribute {}, check dataset attributes '
                    'for a list of valid attributes.'.format(name))

        if filters is not None:
            # Add filter elements.
            for name, value in filters.items():
                try:
                    filter_ = self.filters[name]
                    self._add_filter_node(dataset, filter_, value)
                except KeyError:
                    raise BiomartException(
                        'Unknown filter {}, check dataset filters '
                        'for a list of valid filters.'.format(name))

        # Fetch response.
        response = self.get(query=ElementTree.tostring(root))

        # Raise exception if an error occurred.
        if 'Query ERROR' in response.text:
            raise BiomartException(response.text)

        # Parse results into a DataFrame.
        try:
            result = pd.read_csv(StringIO(response.text),
                                 sep='\t', dtype=dtypes)
            if (datatype == "protein_coordinates"):
                result.to_csv(
                    'Hs_ProteinCoordinates_build_100_38.csv', sep='\t')
            elif(datatype == "protein_feature"):
                result.to_csv(
                    'Hs_ProteinFeatures_build_100_38.csv', sep='\t')
        # Type error is raised of a data type is not understood by pandas
        except TypeError as err:
            raise ValueError("Non valid data type is used in dtypes")

        if use_attr_names:
            # Rename columns with attribute names instead of display names.
            column_map = {
                self.attributes[attr].display_name: attr
                for attr in attributes
            }
            result.rename(columns=column_map, inplace=True)

        return result

    @staticmethod
    def _add_attr_node(root, attr):
        attr_el = ElementTree.SubElement(root, 'Attribute')
        attr_el.set('name', attr.name)

    @staticmethod
    def _add_filter_node(root, filter_, value):
        """Adds filter xml node to root."""
        filter_el = ElementTree.SubElement(root, 'Filter')
        filter_el.set('name', filter_.name)

        # Set filter value depending on type.
        if filter_.type == 'boolean':
            # Boolean case.
            if value is True or value.lower() in {'included', 'only'}:
                filter_el.set('excluded', '0')
            elif value is False or value.lower() == 'excluded':
                filter_el.set('excluded', '1')
            else:
                raise ValueError('Invalid value for boolean filter ({})'
                                 .format(value))
        elif isinstance(value, list) or isinstance(value, tuple):
            # List case.
            filter_el.set('value', ','.join(map(str, value)))
        else:
            # Default case.
            filter_el.set('value', str(value))

    def __repr__(self):
        return ('<biomart.Dataset name={!r}, display_name={!r}>'
                .format(self._name, self._display_name))


class Attribute(object):
    """Biomart dataset attribute.
    Attributes:
        name (str): Attribute name.
        display_name (str): Attribute display name.
        description (str): Attribute description.
    """

    def __init__(self, name, display_name='', description='', default=False):
        """Attribute constructor.
        Args:
            name (str): Attribute name.
            display_name (str): Attribute display name.
            description (str): Attribute description.
            default (bool): Whether the attribute is a default
                attribute of the corresponding datasets.
        """
        self._name = name
        self._display_name = display_name
        self._description = description
        self._default = default

    @property
    def name(self):
        """Name of the attribute."""
        return self._name

    @property
    def display_name(self):
        """Display name of the attribute."""
        return self._display_name

    @property
    def description(self):
        """Description of the attribute."""
        return self._description

    @property
    def default(self):
        """Whether this is a default attribute."""
        return self._default

    def __repr__(self):
        return (('<biomart.Attribute name={!r},'
                 ' display_name={!r}, description={!r}>')
                .format(self._name, self._display_name, self._description))


class Filter(object):
    """Biomart dataset filter.
    Attributes:
        name (str): Filter name.
        type (str): Type of the filter (boolean, int, etc.).
        description (str): Filter description.
    """

    def __init__(self, name, type, description=''):
        """ Filter constructor.
        Args:
            name (str): Filter name.
            type (str): Type of the filter (boolean, int, etc.).
            description (str): Filter description.
        """
        self._name = name
        self._type = type
        self._description = description

    @property
    def name(self):
        """Filter name."""
        return self._name

    @property
    def type(self):
        """Filter type."""
        return self._type

    @property
    def description(self):
        """Filter description."""
        return self._description

    def __repr__(self):
        return ('<biomart.Filter name={!r}, type={!r}>'
                .format(self.name, self.type))


dataset = Dataset(name='apolyacanthus_gene_ensembl',
                  host='http://www.ensembl.org')

# Protein Coordinates
dataset.query(attributes=["ensembl_transcript_id", "ensembl_exon_id", "ensembl_peptide_id", "start_position",
              "end_position", "transcript_start", "transcript_end", "cdd", "cdd_start", "cdd_end"], datatype='protein_coordinates')


# Protein Features
dataset.query(attributes=["ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id_version",
              "interpro", "interpro_description", "interpro_start", "interpro_end", "cdd", "cdd_start", "cdd_end"], datatype='protein_feature')
