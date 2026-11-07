# Combined seven-source rsID dataset

This folder combines the normalized AlphaFold, JASPAR, ClinVar, GWAS Catalog,
gnomAD, GTR, and PheGenI generator outputs without making cross-source biological
claims. There is exactly one top-level record per supplied rsID.

- Input rsIDs: 20
- Present sources: alpha_fold, jaspar, clinvar, gwas_catalog, gnomad, gtr, phegeni
- Missing sources: none

`combined_data.json` is the authoritative nested representation.
`combined_summary.tsv` has one row per rsID and keeps every source object in a
JSON-valued column, alongside useful counts. `combined_data.sqlite` mirrors
those rows with one JSON column per source. `run_manifest.json` records all
consumed files, SHA-256 hashes, embedded source manifests, and source-reported
errors.

One-to-many records are nested by their actual relation keys: genomic variant,
ClinVar record, GWAS/PheGenI source row, gnomAD variant/callset, or GTR test.
They are never flattened into a Cartesian product. Source-specific cautions
and evidence scopes remain authoritative.
