import tarfile
import pooch
from pathlib import Path
import anndata as ad


def fetch_geo_raw():
    """
    Download and extract the raw GEO tarball (~698 MB).
    Returns path to extracted folder of per-sample .h5 files.
    """
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207633&format=file"
    path = pooch.retrieve(
        url=url,
        fname="GSE207633_RAW.tar",
        path=pooch.os_cache("pyudon"),
        known_hash=None,
    )
    extract_dir = Path(path).with_suffix("")  # remove .tar suffix
    if not extract_dir.exists():
        with tarfile.open(path, "r") as tar:
            tar.extractall(extract_dir)
    return extract_dir


def fetch_processed_h5ad():
    """
    Download and load the pre-built .h5ad (~2.5 GB, hosted on Zenodo).
    Returns an AnnData object.
    """
    url = "https://zenodo.org/records/17194840/files/sjia_adata_raw_nometadata.h5ad"
    path = pooch.retrieve(
        url=url,
        fname="sjia_adata_raw_nometadata.h5ad",
        path=pooch.os_cache("pyudon"),
        known_hash=None,
    )
    return ad.read_h5ad(path)

