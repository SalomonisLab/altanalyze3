"""Curated RNA query aliases; native dataset identifiers always take precedence.

The CD206/MRC1 mapping was supplied with the viewer issue on 2026-10-25.
This module performs no network access or dataset transformation.
"""
import re


RNA_GENE_ALIASES = {"CD206": "MRC1"}


def rna_alias_key(value):
    return re.sub(r"[\s_.-]+", "", str(value or "").strip()).upper()
