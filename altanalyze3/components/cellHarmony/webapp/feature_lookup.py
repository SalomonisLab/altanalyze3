"""Resolve supplied RNA aliases without renaming or removing matrix features."""
from __future__ import annotations

import re

import pandas as pd


def token(value) -> str:
    text = str(value or "").strip().upper()
    text = re.sub(r"^(ENS[A-Z]*G\d+)\.\d+$", r"\1", text)
    return re.sub(r"[\s_.-]+", "", text)


class FeatureLookup:
    def __init__(self, names, var=None):
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
                    raise ValueError(f"Feature alias '{requested}' maps to multiple primary IDs. Select an exact feature ID.")
                return index[key]
        return None


def expression_lookup(cache):
    cached = cache.get("feature_lookup")
    if cached is None:
        var = getattr(cache.get("adata"), "var", None) if cache.get("modality", "rna") == "rna" else None
        cached = FeatureLookup(cache["var_names"], var)
        cached = cache.setdefault("feature_lookup", cached)
    return cached
