"""iso2function — isoform-resolved molecular-interaction & functional annotation for AltAnalyze3.

Attaches definitively measured, isoform-resolved interaction/function data (PPI, PDI, activation,
condensate behavior) to long-read isoforms, anchored on the canonical Ensembl91 structure string, then
builds Cytoscape-ready interaction networks and runs isoform-focused gene-set enrichment to interpret
the functional consequences of isoform switches / expression differences.

See GOALS.md for the full design, identifier model, phased plan, and constraints.

Layers:
  ingest/     parse raw interaction tables -> tidy long-form edges
  crosswalk/  reconcile clone IDs <-> structure string <-> observed isoform / ENST / ENSP
  network/    assemble interaction graph + Cytoscape export + isoform-switch differential
  enrichment/ isoform-focused gene-set enrichment over interaction partners / targets
"""

__version__ = "0.0.1"
