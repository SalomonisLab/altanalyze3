"""Build bundled fastCNV gene-coordinate resources from Ensembl GTF files."""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Optional, Sequence
from urllib.parse import urljoin
from urllib.request import urlopen, urlretrieve

import pandas as pd


LOGGER = logging.getLogger(__name__)

ENSEMBL_CURRENT_GTF = "https://ftp.ensembl.org/pub/current_gtf/"
SPECIES = {
    "human": {
        "directory": "homo_sapiens/",
        "pattern": r"Homo_sapiens\.GRCh38\.(\d+)\.chr\.gtf\.gz",
        "template": "Homo_sapiens.GRCh38.{release}.chr.gtf.gz",
        "output": "Hs_Ensembl_GRCh38_genes.tsv",
        "alias_output": "Hs_Ensembl_hg38_genes.tsv",
        "assembly": "GRCh38",
    },
    "mouse": {
        "directory": "mus_musculus/",
        "pattern": r"Mus_musculus\.GRCm39\.(\d+)\.chr\.gtf\.gz",
        "template": "Mus_musculus.GRCm39.{release}.chr.gtf.gz",
        "output": "Mm_Ensembl_GRCm39_genes.tsv",
        "alias_output": None,
        "assembly": "GRCm39",
    },
}


def _extract_attr(attrs: str, key: str) -> str:
    match = re.search(rf'{re.escape(key)} "([^"]+)"', attrs)
    return match.group(1) if match else ""


def _latest_gtf_url(species_config: Dict[str, str]) -> tuple[str, str]:
    directory_url = urljoin(ENSEMBL_CURRENT_GTF, species_config["directory"])
    with urlopen(directory_url, timeout=60) as handle:
        html = handle.read().decode("utf-8", errors="replace")
    matches = re.findall(species_config["pattern"], html)
    if not matches:
        raise RuntimeError(f"No matching GTF found at {directory_url}")
    release = str(max(int(match) for match in matches))
    filename = species_config["template"].format(release=release)
    return urljoin(directory_url, filename), release


def parse_gtf_genes(gtf_gz: Path, *, species: str, assembly: str, release: str) -> pd.DataFrame:
    records = []
    with gzip.open(gtf_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            gene_id = _extract_attr(attrs, "gene_id")
            gene_name = _extract_attr(attrs, "gene_name") or gene_id
            gene_biotype = _extract_attr(attrs, "gene_biotype") or _extract_attr(attrs, "gene_type")
            if not gene_id:
                continue
            records.append(
                {
                    "gene": gene_name,
                    "gene_id": gene_id,
                    "chr": chrom if str(chrom).startswith("chr") else f"chr{chrom}",
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_biotype": gene_biotype,
                    "species": species,
                    "assembly": assembly,
                    "ensembl_release": release,
                    "source": source,
                }
            )
    df = pd.DataFrame.from_records(records)
    if df.empty:
        raise RuntimeError(f"No gene records parsed from {gtf_gz}")
    df = df.sort_values(["chr", "start", "end", "gene"], kind="mergesort").drop_duplicates("gene", keep="first")
    return df


def build_resource(species: str, output_dir: Path, keep_gtf: bool = False) -> Path:
    config = SPECIES[species]
    output_dir.mkdir(parents=True, exist_ok=True)
    url, release = _latest_gtf_url(config)
    LOGGER.info("Downloading %s Ensembl release %s GTF: %s", species, release, url)
    with tempfile.TemporaryDirectory() as temp_dir:
        gtf_path = Path(temp_dir) / Path(url).name
        urlretrieve(url, gtf_path)
        df = parse_gtf_genes(gtf_path, species=species, assembly=config["assembly"], release=release)
        output_path = output_dir / config["output"]
        df.to_csv(output_path, sep="\t", index=False)
        alias_output = config.get("alias_output")
        if alias_output:
            shutil.copyfile(output_path, output_dir / alias_output)
        if keep_gtf:
            shutil.copyfile(gtf_path, output_dir / gtf_path.name)
    LOGGER.info("Wrote %s", output_path)
    return output_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build fastCNV Ensembl gene-coordinate resources.")
    parser.add_argument("--species", choices=sorted(SPECIES) + ["all"], default="all")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "resources",
        help="Directory for generated coordinate TSVs.",
    )
    parser.add_argument("--keep-gtf", action="store_true", help="Keep downloaded GTF files in the output directory.")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s | %(message)s")
    selected = sorted(SPECIES) if args.species == "all" else [args.species]
    for species in selected:
        build_resource(species, args.output_dir, keep_gtf=args.keep_gtf)


if __name__ == "__main__":
    main()
