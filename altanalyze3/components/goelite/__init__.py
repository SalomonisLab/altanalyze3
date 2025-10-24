"""
AltAnalyze3 GO-Elite component.

Provides a modernised GO-Elite workflow including ontology parsing,
enrichment scoring, and GO DAG prioritisation logic.
"""

from .runner import GOEliteRunner, EnrichmentSettings
from .parser import build_tree_with_annotations, ParsedGO
from .structures import GOTree, GOTermNode, EnrichmentResult
from .prio import prioritize_terms, PrioritizationSettings
from .resources import (
    DEFAULT_CACHE_ENV,
    AVAILABLE_SPECIES,
    load_cached_resources,
    prepare_species_resources,
    resolve_cache_dir,
)

__all__ = [
    "GOEliteRunner",
    "EnrichmentSettings",
    "PrioritizationSettings",
    "build_tree_with_annotations",
    "ParsedGO",
    "GOTree",
    "GOTermNode",
    "EnrichmentResult",
    "prioritize_terms",
    "prepare_species_resources",
    "load_cached_resources",
    "resolve_cache_dir",
    "DEFAULT_CACHE_ENV",
    "AVAILABLE_SPECIES",
]
