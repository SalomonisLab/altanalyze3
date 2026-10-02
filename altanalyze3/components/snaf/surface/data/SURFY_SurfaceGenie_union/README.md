# SURFY + SurfaceGenie union: the default SNAF-B surfaceome

SNAF-B uses this database when `--surface_db` is not given (`surface_db.DEFAULT_SURFACE_DB`).
Nathan made it the default on 2026-10-01. Pass `--surface_db alt91` for the legacy
Alt91_db list of 2,810 genes.

## Contents

| File | What it holds |
|---|---|
| `surface_genes.txt` | the 4,009 genes of the union, with `has_reference` per gene |
| `surface_reference.fasta` | 7,598 reference protein sequences for the 3,584 genes that have one |
| `surface_db_params.json` | build provenance and counts |
| `sources/` | the two source tables, the union script, the union table and its report |

## How it was built

1. `sources/build_surface_union.py` joined SURFY (`sources/SURFY-Ensembl.txt`, 2,760 Ensembl genes)
   and SurfaceGenie (`sources/SurfaceGenie.txt`, 3,141 Ensembl genes) on 2026-08-16. The union holds
   4,009 genes: 1,892 in both, 868 SURFY only and 1,249 SurfaceGenie only. No threshold was
   applied. `sources/surface_union_report.txt` gives every count.
2. `altanalyze3 snaf-build-surface-db` (mode `replace`) built this folder from
   `sources/surface_union_all.txt` on 2026-08-16 18:21. 2,467 references came from the Alt91_db
   UniProt FASTA and 1,117 from EnsMart100 UniProt. 425 of the 4,009 genes have no reference protein,
   so SNAF-B cannot score them.

The files are byte-identical copies of the database the pediatric-AML SNAF-B v3 run used
(`/Users/saljh8/Dropbox/SNAF/Pedatric-AML/snaf_b_2026-08-15/surface_db/UNION`) and of its sources
(`/Users/saljh8/Dropbox/SNAF/SurfaceProtein-Database`).

Sources: SURFY, Bausch-Fluck et al., PNAS 2018 (in silico human surfaceome); SurfaceGenie,
Waas et al., Bioinformatics 2020.
