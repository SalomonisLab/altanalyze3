"""
This is a generalized python module for getting data from Ensemble using Biomart server.
"""

import requests
from xml.etree import ElementTree
# import pandas as pd
from io import StringIO
from xml.etree.ElementTree import fromstring as xml_from_string

DEFAULT_HOST = 'http://www.biomart.org'
DEFAULT_PATH = '/biomart/martservice'
DEFAULT_PORT = 80
DEFAULT_SCHEMA = 'default'

ensemble_server = 'http://www.ensembl.org'
species = 'hsapiens_gene_ensembl'

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


class Server(ServerBase):
    """Class representing a biomart server.
    Typically used as main entry point to the biomart server. Provides
    functionality for listing and loading the marts that are available
    on the server.
    Args:
        host (str): Url of host to connect to.
        path (str): Path on the host to access to the biomart service.
        port (int): Port to use for the connection.
        use_cache (bool): Whether to cache requests.
    Examples:
        Connecting to a server and listing available marts:
            >>> server = Server(host='http://www.ensembl.org')
            >>> server.list_marts()
        Retrieving a mart:
            >>> mart = server['ENSEMBL_MART_ENSEMBL']
    """

    _MART_XML_MAP = {
        'name': 'name',
        'database_name': 'database',
        'display_name': 'displayName',
        'host': 'host',
        'path': 'path',
        'virtual_schema': 'serverVirtualSchema'
    }

    def __init__(self, host=None, path=None, port=None, use_cache=True):
        super().__init__(host=host, path=path, port=port, use_cache=use_cache)
        self._marts = None

    def __getitem__(self, name):
        return self.marts[name]

    @property
    def marts(self):
        """List of available marts."""
        if self._marts is None:
            self._marts = self._fetch_marts()
        return self._marts

    def list_marts(self):
        """Lists available marts in a readable DataFrame format.
        Returns:
            pd.DataFrame: Frame listing available marts.
        """

        def _row_gen(attributes):
            for attr in attributes.values():
                yield (attr.name, attr.display_name)

        return pd.DataFrame.from_records(
            _row_gen(self.marts), columns=['name', 'display_name'])

    def _fetch_marts(self):
        response = self.get(type='registry')

        xml = xml_from_string(response.content)
        marts = [
            self._mart_from_xml(child)
            for child in xml.findall('MartURLLocation')
        ]

        return {m.name: m for m in marts}

    def _mart_from_xml(self, node):
        params = {k: node.attrib[v] for k, v in self._MART_XML_MAP.items()}
        params['extra_params'] = {
            k: v
            for k, v in node.attrib.items()
            if k not in set(self._MART_XML_MAP.values())
        }
        return Mart(use_cache=self.use_cache, **params)

    def __repr__(self):
        return ('<biomart.Server host={!r}, path={!r}, port={!r}>'
                .format(self.host, self.path, self.port))


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
                    use_cache=True,
                    virtual_schema=DEFAULT_SCHEMA):
            super().__init__(host=host, path=path, port=port, use_cache=use_cache)

            self._name = name
            self._display_name = display_name
            self._virtual_schema = virtual_schema

            self._filters = None
            self._attributes = None
            self._default_attributes = None


    def _fetch_configuration(self):
            # Get datasets using biomart.
            response = ServerBase.get(type='configuration', dataset=self._name)

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
                dtypes = None
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

            root = ElementTree.Element('Query')
            root.set('virtualSchemaName', self._virtual_schema)
            root.set('formatter', 'TSV')
            root.set('header', '1')
            root.set('uniqueRows', str(int(only_unique)))
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
            print(response)

            # Raise exception if an error occurred.
            if 'Query ERROR' in response.text:
                raise BiomartException(response.text)

            # Parse results into a DataFrame.
            try:
                result = pd.read_csv(StringIO(response.text), sep='\t', dtype=dtypes)
            #Type error is raised of a data type is not understood by pandas
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




dataset = Dataset(name = species,host = ensemble_server)
dataset.query(attributes = ["ensembl_exon_id", "ensembl_peptide_id", "transcript_start","transcript_end","interpro","interpro_short_description","interpro_start","interpro_end"])

