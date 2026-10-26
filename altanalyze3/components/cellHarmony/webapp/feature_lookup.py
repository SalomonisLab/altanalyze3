"""Resolve supplied RNA aliases without renaming or removing matrix features."""
from __future__ import annotations

import re

import pandas as pd
from altanalyze3.components.annotation.gene_aliases import RNA_GENE_ALIASES


# RNA query names only; predicted modalities have independent feature identities.
# Exact dataset IDs and supplied annotations always take precedence.
class AmbiguousFeatureError(ValueError):
    """The query needs an exact native feature ID, not a server error."""


def token(value) -> str:
    text = str(value or "").strip().upper()
    text = re.sub(r"^(ENS[A-Z]*G\d+)\.\d+$", r"\1", text)
    return re.sub(r"[\s_.-]+", "", text)


class FeatureLookup:
    def __init__(self, names, var=None, *, rna_aliases=True):
        self.rna_aliases = rna_aliases
        self.names = [str(name) for name in names]
        self.rows = {name: i for i, name in enumerate(self.names)}
        self.normalized = {}
        self.aliases = {}
        for name in self.names:
            self._add(self.normalized, token(name), name)
        if var is not None:
            for column in ("gene_symbols", "gene_symbol", "feature_name", "Gene", "Symbol"):
                if column not in var:
                    continue
                for name, value in zip(self.names, var[column].astype(object)):
                    if pd.isna(value) or not str(value).strip():
                        continue
                    alias = str(value).strip()
                    self._add(self.aliases, token(alias), name)
        # Suggestions retain every primary ID and add unambiguous supplied
        # aliases. No feature is collapsed when symbols are shared.
        aliases = []
        if var is not None:
            for column in ("gene_symbols", "gene_symbol", "feature_name", "Gene", "Symbol"):
                if column in var:
                    aliases.extend(str(v).strip() for v in var[column].astype(object)
                                   if not pd.isna(v) and str(v).strip() and self.aliases.get(token(v)))
        self.suggestions = list(dict.fromkeys(self.names + aliases))
        if rna_aliases:
            for alias in RNA_GENE_ALIASES:
                try:
                    if self.resolve(alias) is not None and alias not in self.suggestions:
                        self.suggestions.append(alias)
                except AmbiguousFeatureError:
                    pass  # All native IDs remain available; do not offer an ambiguous alias.
        # One readable suggestion per feature keeps the browser's datalist from
        # doubling in size. Exact primary IDs and other supplied aliases remain
        # accepted query keys. Shared symbols display the distinct primary IDs.
        self.display_names = list(self.names)
        if var is not None:
            for column in ("gene_symbols", "gene_symbol", "feature_name", "Gene", "Symbol"):
                if column not in var:
                    continue
                for i, value in enumerate(var[column].astype(object)):
                    if pd.isna(value) or self.display_names[i] != self.names[i]:
                        continue
                    alias = str(value).strip()
                    if alias and not re.fullmatch(r"ENS[A-Z]*G\d+(?:\.\d+)?", alias.upper()) and self.aliases.get(token(alias)) == self.names[i]:
                        self.display_names[i] = alias

    @staticmethod
    def _add(index, key, name):
        if key in index and index[key] != name:
            index[key] = None
        else:
            index[key] = name

    def resolve(self, requested):
        requested = str(requested or "").strip()
        if requested in self.rows:
            return requested
        key = token(requested)
        for index in (self.normalized, self.aliases):
            if key in index:
                if index[key] is None:
                    raise AmbiguousFeatureError(f"Feature alias '{requested}' maps to multiple primary IDs. Select an exact feature ID.")
                return index[key]
        canonical = RNA_GENE_ALIASES.get(key) if self.rna_aliases else None
        if canonical:
            for index in (self.normalized, self.aliases):
                target = token(canonical)
                if target in index:
                    if index[target] is None:
                        raise AmbiguousFeatureError(f"Feature alias '{requested}' maps to multiple primary IDs. Select an exact feature ID.")
                    return index[target]
        return None


def expression_lookup(cache):
    cached = cache.get("feature_lookup")
    if cached is None:
        var = getattr(cache.get("adata"), "var", None) if cache.get("modality", "rna") == "rna" else None
        cached = FeatureLookup(cache["var_names"], var, rna_aliases=cache.get("modality", "rna") == "rna")
        cached = cache.setdefault("feature_lookup", cached)
    return cached
