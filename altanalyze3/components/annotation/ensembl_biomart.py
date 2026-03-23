"""
Generalized module for querying Ensembl BioMart, including support for
archived Ensembl versions and bulk export of transcript/protein annotations.
"""

import os
import re
import gzip
import logging
import requests
from xml.etree import ElementTree
import pandas as pd
from io import StringIO
from urllib.parse import urljoin

logger = logging.getLogger(__name__)

DEFAULT_HOST = 'http://www.biomart.org'
DEFAULT_PATH = '/biomart/martservice'
DEFAULT_PORT = 80
DEFAULT_SCHEMA = 'default'

# Ensembl release version -> archive host mapping
ENSEMBL_ARCHIVE_HOSTS = {
    75:  'http://feb2014.archive.ensembl.org',
    80:  'http://may2015.archive.ensembl.org',
    85:  'http://jul2016.archive.ensembl.org',
    90:  'http://aug2017.archive.ensembl.org',
    95:  'http://jan2019.archive.ensembl.org',
    100: 'http://apr2020.archive.ensembl.org',
    101: 'http://aug2020.archive.ensembl.org',
    102: 'http://nov2020.archive.ensembl.org',
    103: 'http://feb2021.archive.ensembl.org',
    104: 'http://may2021.archive.ensembl.org',
    105: 'http://dec2021.archive.ensembl.org',
    106: 'http://apr2022.archive.ensembl.org',
    107: 'http://jul2022.archive.ensembl.org',
    108: 'http://oct2022.archive.ensembl.org',
    109: 'http://feb2023.archive.ensembl.org',
    110: 'http://jul2023.archive.ensembl.org',
    111: 'http://jan2024.archive.ensembl.org',
    112: 'http://may2024.archive.ensembl.org',
}

# Human-readable chromosome list used to chunk large sequence queries
HUMAN_CHROMOSOMES = [str(c) for c in range(1, 23)] + ['X', 'Y', 'MT']

ENSEMBL_FTP_BASE = 'https://ftp.ensembl.org/pub'
ENSEMBL_SPECIES_ALIASES = {
    'hsapiens': 'homo_sapiens',
    'mmusculus': 'mus_musculus',
    'drerio': 'danio_rerio',
    'rnorvegicus': 'rattus_norvegicus',
}
ENSEMBL_REQUEST_TIMEOUT = 120


def _normalize_species_ftp_name(species):
    value = str(species or '').strip().lower()
    if not value:
        raise ValueError("species must be provided")
    if value in ENSEMBL_SPECIES_ALIASES:
        return ENSEMBL_SPECIES_ALIASES[value]
    if '_' in value:
        return value
    return value


def _species_filename_prefix(species_ftp_name):
    parts = species_ftp_name.split('_', 1)
    if len(parts) == 1:
        return parts[0].capitalize()
    return f"{parts[0].capitalize()}_{parts[1]}"


def _strip_ensembl_version(identifier):
    value = str(identifier or '').strip()
    if not value:
        return ''
    return value.split('.', 1)[0]


def _http_get(url, session=None, stream=False):
    requester = session or requests
    response = requester.get(url, timeout=ENSEMBL_REQUEST_TIMEOUT, stream=stream)
    response.raise_for_status()
    return response


def _list_ftp_filenames(url, session=None):
    response = _http_get(url, session=session)
    names = set()
    for href in re.findall(r'href="([^"]+)"', response.text):
        candidate = os.path.basename(href.strip().rstrip('/'))
        if not candidate or candidate in {'.', '..'} or candidate.startswith('?'):
            continue
        names.add(candidate)
    if not names:
        raise RuntimeError(f"Could not list FTP directory contents: {url}")
    return sorted(names)


def _resolve_release_asset_url(ensembl_version, relative_dir, filename_pattern,
                               session=None):
    dir_url = f"{ENSEMBL_FTP_BASE}/release-{ensembl_version}/{relative_dir.rstrip('/')}/"
    filenames = _list_ftp_filenames(dir_url, session=session)
    matches = [name for name in filenames if re.match(filename_pattern, name)]
    if not matches:
        raise RuntimeError(
            f"No Ensembl release file matched {filename_pattern!r} under {dir_url}"
        )
    if len(matches) > 1:
        logger.info("Multiple matching files found under %s, selecting %s", dir_url, matches[0])
    return urljoin(dir_url, matches[0])


def _download_to_path(url, output_path, session=None):
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        logger.info("Using existing Ensembl asset: %s", output_path)
        return output_path

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    logger.info("Downloading %s -> %s", url, output_path)
    response = _http_get(url, session=session, stream=True)
    with open(output_path, 'wb') as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    return output_path


def _iter_fasta_records_gz(path):
    header = None
    seq_chunks = []
    with gzip.open(path, 'rt') as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, ''.join(seq_chunks)


def _parse_ensembl_fasta_header(header):
    tokens = str(header or '').strip().split()
    if not tokens:
        return {}
    parsed = {'sequence_id': _strip_ensembl_version(tokens[0])}
    for token in tokens[1:]:
        if ':' not in token:
            continue
        key, value = token.split(':', 1)
        parsed[key] = value
    if 'gene' in parsed:
        parsed['gene'] = _strip_ensembl_version(parsed['gene'])
    if 'transcript' in parsed:
        parsed['transcript'] = _strip_ensembl_version(parsed['transcript'])
    return parsed


def _parse_cdna_fasta(cdna_fasta_path):
    rows = []
    metadata = {}
    for header, sequence in _iter_fasta_records_gz(cdna_fasta_path):
        parsed = _parse_ensembl_fasta_header(header)
        transcript_id = parsed.get('transcript') or parsed.get('sequence_id')
        gene_id = parsed.get('gene', '')
        transcript_biotype = parsed.get('transcript_biotype', '')
        if not transcript_id or not sequence:
            continue
        rows.append({
            'ensembl_transcript_id': transcript_id,
            'ensembl_gene_id': gene_id,
            'cdna': sequence,
        })
        metadata[transcript_id] = {
            'ensembl_gene_id': gene_id,
            'transcript_biotype': transcript_biotype,
            'gene_symbol': parsed.get('gene_symbol', ''),
            'gene_biotype': parsed.get('gene_biotype', ''),
            'description': parsed.get('description', ''),
        }
    df = pd.DataFrame(rows, columns=[
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'cdna',
    ])
    if not df.empty:
        df = df.drop_duplicates(subset=['ensembl_transcript_id']).reset_index(drop=True)
    return df, metadata


def _parse_pep_fasta(pep_fasta_path):
    rows = []
    for header, sequence in _iter_fasta_records_gz(pep_fasta_path):
        parsed = _parse_ensembl_fasta_header(header)
        peptide_id = parsed.get('sequence_id')
        transcript_id = parsed.get('transcript', '')
        gene_id = parsed.get('gene', '')
        if not peptide_id or not sequence:
            continue
        rows.append({
            'ensembl_peptide_id': peptide_id,
            'ensembl_transcript_id': transcript_id,
            'ensembl_gene_id': gene_id,
            'transcript_biotype': parsed.get('transcript_biotype', ''),
            'gene_biotype': parsed.get('gene_biotype', ''),
            'peptide': sequence,
        })
    df = pd.DataFrame(rows, columns=[
        'ensembl_peptide_id',
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'transcript_biotype',
        'gene_biotype',
        'peptide',
    ])
    if not df.empty:
        df = df.drop_duplicates(subset=['ensembl_peptide_id']).reset_index(drop=True)
    return df


def _parse_gtf_attributes(field):
    return {
        key: value
        for key, value in re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)";', str(field or ''))
    }


def _parse_transcript_annotations_from_gtf(gtf_path):
    rows = []
    with gzip.open(gtf_path, 'rt') as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9 or parts[2] != 'transcript':
                continue
            attrs = _parse_gtf_attributes(parts[8])
            transcript_id = _strip_ensembl_version(attrs.get('transcript_id'))
            gene_id = _strip_ensembl_version(attrs.get('gene_id'))
            if not transcript_id:
                continue
            transcript_biotype = (
                attrs.get('transcript_biotype')
                or attrs.get('transcript_type')
                or ''
            )
            gene_biotype = attrs.get('gene_biotype') or attrs.get('gene_type') or ''
            rows.append({
                'ensembl_gene_id': gene_id,
                'ensembl_transcript_id': transcript_id,
                'gene_symbol': attrs.get('gene_name', ''),
                'gene_biotype': gene_biotype,
                'transcript_biotype': transcript_biotype,
                'chromosome_name': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
            })
    df = pd.DataFrame(rows, columns=[
        'ensembl_gene_id',
        'ensembl_transcript_id',
        'gene_symbol',
        'gene_biotype',
        'transcript_biotype',
        'chromosome_name',
        'start',
        'end',
        'strand',
    ])
    if not df.empty:
        df = df.drop_duplicates(subset=['ensembl_transcript_id']).reset_index(drop=True)
    return df


def get_ensembl_archive_host(version):
    """Return the archive host URL for a given Ensembl release version.

    Args:
        version (int): Ensembl release number (e.g. 100).

    Returns:
        str: Archive host URL.

    Raises:
        ValueError: If the version is not in the known archive mapping.
    """
    if version not in ENSEMBL_ARCHIVE_HOSTS:
        raise ValueError(
            f"Ensembl version {version} not in archive mapping. "
            f"Known versions: {sorted(ENSEMBL_ARCHIVE_HOSTS.keys())}. "
            f"Pass host= directly for unlisted versions."
        )
    return ENSEMBL_ARCHIVE_HOSTS[version]


def get_ensembl_dataset(species='hsapiens', ensembl_version=None, host=None):
    """Create a Dataset connected to a specific Ensembl version.

    Args:
        species (str): Species prefix (default 'hsapiens').
        ensembl_version (int or None): Ensembl release number. If None, uses
            the current release at ensembl.org.
        host (str or None): Explicit host URL. Overrides ensembl_version.

    Returns:
        Dataset: A connected Dataset instance.
    """
    if host is None:
        if ensembl_version is not None:
            host = get_ensembl_archive_host(ensembl_version)
        else:
            host = 'http://www.ensembl.org'
    return Dataset(
        name=f'{species}_gene_ensembl',
        host=host,
        path='/biomart/martservice',
        virtual_schema='default',
    )


class ServerBase(object):
    """Base class that handles requests to the biomart server.

    Attributes:
        host (str): Host to connect to for the biomart service.
        path (str): Path to the biomart service on the host.
        port (str): Port to connect to on the host.
        url (str): Url used to connect to the biomart service.
    """

    def __init__(self, host=None, path=None, port=None):
        host = host or DEFAULT_HOST
        path = path or DEFAULT_PATH
        port = port or DEFAULT_PORT

        host = self._add_http_prefix(host)
        host = self._remove_trailing_slash(host)

        if not path.startswith('/'):
            path = '/' + path

        self._host = host
        self._path = path
        self._port = port

    @property
    def host(self):
        return self._host

    @property
    def path(self):
        return self._path

    @property
    def port(self):
        return self._port

    @property
    def url(self):
        return '{}:{}{}'.format(self._host, self._port, self._path)

    @staticmethod
    def _add_http_prefix(url, prefix='http://'):
        if not url.startswith('http://') and not url.startswith('https://'):
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
            **params: Arbitrary keyword arguments added as request parameters.

        Returns:
            requests.models.Response: Response from biomart.
        """
        r = requests.get(self.url, params=params)
        r.raise_for_status()
        return r


class BiomartException(Exception):
    """Basic exception class for biomart exceptions."""
    pass


class Dataset(ServerBase):
    """Class representing a biomart dataset.

    Handles queries to biomart datasets. Queries can select a subset of
    attributes and can be filtered using any available filters.

    Args:
        name (str): Id of the dataset.
        display_name (str): Display name of the dataset.
        host (str): Url of host to connect to.
        path (str): Path on the host to access to the biomart service.
        port (int): Port to use for the connection.
        virtual_schema (str): The virtual schema of the dataset.

    Examples:
        Directly connecting to a dataset:
            >>> dataset = Dataset(name='hsapiens_gene_ensembl',
            ...                   host='http://www.ensembl.org')
        Querying the dataset:
            >>> dataset.query(attributes=['ensembl_gene_id',
            ...                           'external_gene_name'],
            ...               filters={'chromosome_name': ['1','2']})
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

    @property
    def name(self):
        return self._name

    @property
    def display_name(self):
        return self._display_name

    @property
    def filters(self):
        if self._filters is None:
            self._filters, self._attributes = self._fetch_configuration()
        return self._filters

    @property
    def attributes(self):
        if self._attributes is None:
            self._filters, self._attributes = self._fetch_configuration()
        return self._attributes

    @property
    def default_attributes(self):
        if self._default_attributes is None:
            self._default_attributes = {
                name: attr
                for name, attr in self.attributes.items()
                if attr.default is True
            }
        return self._default_attributes

    def list_attributes(self):
        """Lists available attributes in a readable DataFrame format."""
        def _row_gen(attributes):
            for attr in attributes.values():
                yield (attr.name, attr.display_name, attr.description)
        return pd.DataFrame.from_records(
            _row_gen(self.attributes),
            columns=['name', 'display_name', 'description'])

    def list_filters(self):
        """Lists available filters in a readable DataFrame format."""
        def _row_gen(attributes):
            for attr in attributes.values():
                yield (attr.name, attr.type, attr.description)
        return pd.DataFrame.from_records(
            _row_gen(self.filters), columns=['name', 'type', 'description'])

    def _fetch_configuration(self):
        response = self.get(type='configuration', dataset=self._name)
        if 'Problem retrieving configuration' in response.text:
            raise BiomartException('Failed to retrieve dataset configuration, '
                                   'check the dataset name and schema.')
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
              dtypes=None):
        """Queries the dataset to retrieve the contained data.

        Args:
            attributes (list[str]): Names of attributes to fetch.
            filters (dict[str, any]): Filters to apply to the query.
            only_unique (bool): Whether to return only unique rows.
            use_attr_names (bool): Use attribute names as column names
                instead of display names.
            dtypes (dict[str, any]): Column data types for pandas.

        Returns:
            pandas.DataFrame: Query results.
        """
        root = ElementTree.Element('Query')
        root.set('virtualSchemaName', self._virtual_schema)
        root.set('formatter', 'TSV')
        root.set('header', '1')
        root.set('uniqueRows', str(int(only_unique)))
        root.set('datasetConfigVersion', '0.6')

        dataset = ElementTree.SubElement(root, 'Dataset')
        dataset.set('name', self.name)
        dataset.set('interface', 'default')

        if attributes is None:
            attributes = list(self.default_attributes.keys())

        for name in attributes:
            try:
                attr = self.attributes[name]
                self._add_attr_node(dataset, attr)
            except KeyError:
                raise BiomartException(
                    'Unknown attribute {}, check dataset attributes '
                    'for a list of valid attributes.'.format(name))

        if filters is not None:
            for name, value in filters.items():
                try:
                    filter_ = self.filters[name]
                    self._add_filter_node(dataset, filter_, value)
                except KeyError:
                    raise BiomartException(
                        'Unknown filter {}, check dataset filters '
                        'for a list of valid filters.'.format(name))

        response = self.get(query=ElementTree.tostring(root))

        if 'Query ERROR' in response.text:
            raise BiomartException(response.text)

        try:
            result = pd.read_csv(StringIO(response.text), sep='\t', dtype=dtypes)
        except TypeError:
            raise ValueError("Invalid data type in dtypes argument")

        if use_attr_names:
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
        filter_el = ElementTree.SubElement(root, 'Filter')
        filter_el.set('name', filter_.name)
        if filter_.type == 'boolean':
            if value is True or (isinstance(value, str) and
                                 value.lower() in {'included', 'only'}):
                filter_el.set('excluded', '0')
            elif value is False or (isinstance(value, str) and
                                    value.lower() == 'excluded'):
                filter_el.set('excluded', '1')
            else:
                raise ValueError(
                    'Invalid value for boolean filter ({})'.format(value))
        elif isinstance(value, (list, tuple)):
            filter_el.set('value', ','.join(map(str, value)))
        else:
            filter_el.set('value', str(value))

    def __repr__(self):
        return ('<biomart.Dataset name={!r}, display_name={!r}>'
                .format(self._name, self._display_name))


class Attribute(object):
    """Biomart dataset attribute."""

    def __init__(self, name, display_name='', description='', default=False):
        self._name = name
        self._display_name = display_name
        self._description = description
        self._default = default

    @property
    def name(self):
        return self._name

    @property
    def display_name(self):
        return self._display_name

    @property
    def description(self):
        return self._description

    @property
    def default(self):
        return self._default

    def __repr__(self):
        return ('<biomart.Attribute name={!r}, display_name={!r}, '
                'description={!r}>').format(
                    self._name, self._display_name, self._description)


class Filter(object):
    """Biomart dataset filter."""

    def __init__(self, name, type, description=''):
        self._name = name
        self._type = type
        self._description = description

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def description(self):
        return self._description

    def __repr__(self):
        return '<biomart.Filter name={!r}, type={!r}>'.format(
            self.name, self.type)


# ---------------------------------------------------------------------------
# Bulk export helpers
# ---------------------------------------------------------------------------

def export_ensembl_annotations_from_release_files(
    ensembl_version,
    output_dir,
    species='hsapiens',
    species_ftp_name=None,
    session=None,
):
    """Export Ensembl transcript annotations and sequences from release files.

    This path avoids BioMart entirely. It downloads the release-specific
    cDNA FASTA, peptide FASTA, and GTF from the Ensembl FTP archive and
    parses them locally.
    """
    species_ftp_name = _normalize_species_ftp_name(species_ftp_name or species)
    prefix = _species_filename_prefix(species_ftp_name)
    release_cache_dir = os.path.join(output_dir, f'ensembl_release_{ensembl_version}_cache')
    os.makedirs(release_cache_dir, exist_ok=True)

    gtf_url = _resolve_release_asset_url(
        ensembl_version,
        f'gtf/{species_ftp_name}',
        rf'^{re.escape(prefix)}\..+\.{ensembl_version}\.gtf\.gz$',
        session=session,
    )
    cdna_url = _resolve_release_asset_url(
        ensembl_version,
        f'fasta/{species_ftp_name}/cdna',
        rf'^{re.escape(prefix)}\..+\.cdna\.all\.fa\.gz$',
        session=session,
    )
    pep_url = _resolve_release_asset_url(
        ensembl_version,
        f'fasta/{species_ftp_name}/pep',
        rf'^{re.escape(prefix)}\..+\.pep\.all\.fa\.gz$',
        session=session,
    )

    gtf_path = _download_to_path(gtf_url, os.path.join(release_cache_dir, os.path.basename(gtf_url)), session=session)
    cdna_path = _download_to_path(cdna_url, os.path.join(release_cache_dir, os.path.basename(cdna_url)), session=session)
    pep_path = _download_to_path(pep_url, os.path.join(release_cache_dir, os.path.basename(pep_url)), session=session)

    logger.info("Parsing Ensembl cDNA FASTA: %s", cdna_path)
    cdna_df, cdna_meta = _parse_cdna_fasta(cdna_path)

    logger.info("Parsing Ensembl peptide FASTA: %s", pep_path)
    pep_df = _parse_pep_fasta(pep_path)

    logger.info("Parsing Ensembl transcript annotations from GTF: %s", gtf_path)
    transcript_annotations_df = _parse_transcript_annotations_from_gtf(gtf_path)

    enst_to_ensp_df = pep_df[['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'transcript_biotype']].copy()
    if not enst_to_ensp_df.empty:
        missing_gene = enst_to_ensp_df['ensembl_gene_id'].eq('')
        missing_biotype = enst_to_ensp_df['transcript_biotype'].eq('')
        if missing_gene.any() or missing_biotype.any():
            meta_df = pd.DataFrame.from_dict(cdna_meta, orient='index').reset_index()
            meta_df.rename(columns={'index': 'ensembl_transcript_id'}, inplace=True)
            enst_to_ensp_df = enst_to_ensp_df.merge(
                meta_df[['ensembl_transcript_id', 'ensembl_gene_id', 'transcript_biotype']],
                on='ensembl_transcript_id',
                how='left',
                suffixes=('', '_cdna'),
            )
            enst_to_ensp_df['ensembl_gene_id'] = enst_to_ensp_df['ensembl_gene_id'].mask(
                enst_to_ensp_df['ensembl_gene_id'].eq(''),
                enst_to_ensp_df['ensembl_gene_id_cdna'],
            )
            enst_to_ensp_df['transcript_biotype'] = enst_to_ensp_df['transcript_biotype'].mask(
                enst_to_ensp_df['transcript_biotype'].eq(''),
                enst_to_ensp_df['transcript_biotype_cdna'],
            )
            enst_to_ensp_df.drop(
                columns=['ensembl_gene_id_cdna', 'transcript_biotype_cdna'],
                inplace=True,
            )
        enst_to_ensp_df = enst_to_ensp_df[
            enst_to_ensp_df['ensembl_peptide_id'].astype(str).str.strip().astype(bool)
        ]
        enst_to_ensp_df = enst_to_ensp_df.drop_duplicates().reset_index(drop=True)

    mrna_df = cdna_df.copy()
    if not mrna_df.empty:
        mrna_df = mrna_df.dropna(subset=['cdna'])
        mrna_df = mrna_df[mrna_df['cdna'].str.strip().astype(bool)]
        mrna_df = mrna_df.drop_duplicates(subset=['ensembl_transcript_id']).reset_index(drop=True)

    protein_df = pep_df[['ensembl_transcript_id', 'ensembl_peptide_id', 'peptide']].copy()
    if not protein_df.empty:
        protein_df = protein_df.dropna(subset=['peptide'])
        protein_df = protein_df[protein_df['peptide'].str.strip().astype(bool)]
        protein_df = protein_df.drop_duplicates(subset=['ensembl_peptide_id']).reset_index(drop=True)

    os.makedirs(output_dir, exist_ok=True)
    v = ensembl_version
    results = {}

    transcript_annotations_path = os.path.join(output_dir, f'ENST_transcript_annotations_build_{v}.tsv')
    transcript_annotations_df.to_csv(transcript_annotations_path, sep='\t', index=False)
    results['transcript_annotations'] = (transcript_annotations_path, transcript_annotations_df)

    enst_ensp_path = os.path.join(output_dir, f'ENST_to_ENSP_build_{v}.tsv')
    enst_to_ensp_df.to_csv(enst_ensp_path, sep='\t', index=False)
    results['enst_to_ensp'] = (enst_ensp_path, enst_to_ensp_df)

    mrna_path = os.path.join(output_dir, f'ENST_mRNA_sequences_build_{v}.tsv')
    mrna_df.to_csv(mrna_path, sep='\t', index=False)
    results['mrna_sequences'] = (mrna_path, mrna_df)

    prot_path = os.path.join(output_dir, f'ENSP_protein_sequences_build_{v}.tsv')
    protein_df.to_csv(prot_path, sep='\t', index=False)
    results['protein_sequences'] = (prot_path, protein_df)

    logger.info("Static Ensembl release export complete -> %s", output_dir)
    return results

def _query_by_chromosome(dataset, attributes, chromosomes=None,
                         use_attr_names=True):
    """Query a dataset chromosome-by-chromosome to avoid BioMart timeouts.

    Large queries (especially those fetching sequences) frequently time out
    when requesting all chromosomes at once.  This helper iterates over
    chromosomes, concatenates the partial results, and logs progress.

    Args:
        dataset (Dataset): Connected Dataset instance.
        attributes (list[str]): BioMart attribute names.
        chromosomes (list[str] or None): Chromosomes to iterate over.
            Defaults to HUMAN_CHROMOSOMES.
        use_attr_names (bool): Use internal attribute names as columns.

    Returns:
        pandas.DataFrame: Concatenated results across all chromosomes.
    """
    if chromosomes is None:
        chromosomes = HUMAN_CHROMOSOMES

    frames = []
    for chrom in chromosomes:
        logger.info("  querying chromosome %s ...", chrom)
        try:
            df = dataset.query(
                attributes=attributes,
                filters={'chromosome_name': chrom},
                use_attr_names=use_attr_names,
            )
            if not df.empty:
                frames.append(df)
        except (BiomartException, requests.RequestException) as exc:
            logger.warning("  chromosome %s failed: %s", chrom, exc)
    if not frames:
        return pd.DataFrame(columns=attributes)
    return pd.concat(frames, ignore_index=True)


def export_enst_to_ensp(dataset, output_path, chromosomes=None):
    """Export ENST -> ENSP (transcript to protein) associations as TSV.

    Args:
        dataset (Dataset): Connected Ensembl BioMart Dataset.
        output_path (str): Destination file path for the TSV.
        chromosomes (list[str] or None): Chromosomes to query.

    Returns:
        pandas.DataFrame: The exported data.
    """
    logger.info("Exporting ENST -> ENSP associations to %s", output_path)
    attributes = [
        'ensembl_gene_id',
        'ensembl_transcript_id',
        'ensembl_peptide_id',
        'transcript_biotype',
    ]
    df = _query_by_chromosome(dataset, attributes, chromosomes)
    # Drop rows with no peptide ID (non-coding transcripts)
    df = df.dropna(subset=['ensembl_peptide_id'])
    df = df[df['ensembl_peptide_id'].str.strip().astype(bool)]
    df = df.drop_duplicates().reset_index(drop=True)
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info("  wrote %d ENST->ENSP associations", len(df))
    return df


def export_mrna_sequences(dataset, output_path, chromosomes=None):
    """Export ENST -> mRNA (cDNA) sequences as TSV.

    Columns: ensembl_transcript_id, ensembl_gene_id, cdna

    Args:
        dataset (Dataset): Connected Ensembl BioMart Dataset.
        output_path (str): Destination file path for the TSV.
        chromosomes (list[str] or None): Chromosomes to query.

    Returns:
        pandas.DataFrame: The exported data.
    """
    logger.info("Exporting mRNA sequences to %s", output_path)
    attributes = [
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'cdna',
    ]
    df = _query_by_chromosome(dataset, attributes, chromosomes)
    df = df.dropna(subset=['cdna'])
    df = df[df['cdna'].str.strip().astype(bool)]
    df = df.drop_duplicates(subset=['ensembl_transcript_id']).reset_index(drop=True)
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info("  wrote %d mRNA sequences", len(df))
    return df


def export_protein_sequences(dataset, output_path, chromosomes=None):
    """Export ENST/ENSP -> protein (peptide) sequences as TSV.

    Columns: ensembl_transcript_id, ensembl_peptide_id, peptide

    Args:
        dataset (Dataset): Connected Ensembl BioMart Dataset.
        output_path (str): Destination file path for the TSV.
        chromosomes (list[str] or None): Chromosomes to query.

    Returns:
        pandas.DataFrame: The exported data.
    """
    logger.info("Exporting protein sequences to %s", output_path)
    attributes = [
        'ensembl_transcript_id',
        'ensembl_peptide_id',
        'peptide',
    ]
    df = _query_by_chromosome(dataset, attributes, chromosomes)
    df = df.dropna(subset=['peptide'])
    df = df[df['peptide'].str.strip().astype(bool)]
    df = df.drop_duplicates(subset=['ensembl_peptide_id']).reset_index(drop=True)
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info("  wrote %d protein sequences", len(df))
    return df


def export_ensembl_annotations(ensembl_version, output_dir, species='hsapiens',
                               host=None, chromosomes=None, source='ftp',
                               species_ftp_name=None, session=None):
    """Export all three annotation files for a given Ensembl version.

    Produces:
        - {output_dir}/ENST_transcript_annotations_build_{version}.tsv
        - {output_dir}/ENST_to_ENSP_build_{version}.tsv
        - {output_dir}/ENST_mRNA_sequences_build_{version}.tsv
        - {output_dir}/ENSP_protein_sequences_build_{version}.tsv

    Args:
        ensembl_version (int): Ensembl release number (e.g. 100).
        output_dir (str): Directory to write output files.
        species (str): Species prefix (default 'hsapiens').
        host (str or None): Explicit host URL, overrides ensembl_version.
        chromosomes (list[str] or None): Chromosomes to query. Defaults
            to human autosomes + X, Y, MT.
        source (str): One of 'ftp', 'auto', or 'biomart'. Default 'ftp'
            uses release-specific static Ensembl files directly. 'auto'
            prefers FTP and falls back to BioMart if needed.
        species_ftp_name (str or None): Explicit Ensembl FTP species directory
            name, e.g. 'homo_sapiens'. Defaults to a normalized version of
            species.
        session (requests.Session or None): Optional requests session for
            static-file downloads.

    Returns:
        dict: Mapping of output type -> (file_path, DataFrame).

    Example:
        >>> from components.annotation.ensembl_biomart import export_ensembl_annotations
        >>> results = export_ensembl_annotations(100, output_dir='./ensembl_export')
    """
    if source not in {'auto', 'ftp', 'biomart'}:
        raise ValueError("source must be one of: 'auto', 'ftp', 'biomart'")

    if source in {'auto', 'ftp'} and host is None:
        try:
            return export_ensembl_annotations_from_release_files(
                ensembl_version=ensembl_version,
                output_dir=output_dir,
                species=species,
                species_ftp_name=species_ftp_name,
                session=session,
            )
        except Exception as exc:
            if source == 'ftp':
                raise
            logger.warning(
                "Static Ensembl release export failed for release %s: %s. "
                "Falling back to BioMart.",
                ensembl_version,
                exc,
            )

    dataset = get_ensembl_dataset(
        species=species, ensembl_version=ensembl_version, host=host)

    os.makedirs(output_dir, exist_ok=True)
    v = ensembl_version

    results = {}

    enst_ensp_path = os.path.join(output_dir, f'ENST_to_ENSP_build_{v}.tsv')
    results['enst_to_ensp'] = (
        enst_ensp_path,
        export_enst_to_ensp(dataset, enst_ensp_path, chromosomes),
    )

    mrna_path = os.path.join(output_dir, f'ENST_mRNA_sequences_build_{v}.tsv')
    results['mrna_sequences'] = (
        mrna_path,
        export_mrna_sequences(dataset, mrna_path, chromosomes),
    )

    prot_path = os.path.join(output_dir, f'ENSP_protein_sequences_build_{v}.tsv')
    results['protein_sequences'] = (
        prot_path,
        export_protein_sequences(dataset, prot_path, chromosomes),
    )

    logger.info("All exports complete -> %s", output_dir)
    return results
