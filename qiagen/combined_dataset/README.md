# Combined eight-source rsID dataset

This folder combines the normalized AlphaFold, JASPAR, ClinVar, GWAS Catalog,
gnomAD, GTR, PheGenI, and PharmCAT generator outputs without making cross-source
biological claims. There is exactly one top-level record per supplied rsID.

- Input rsIDs: 20
- Present sources: alpha_fold, jaspar, clinvar, gwas_catalog, gnomad, gtr, phegeni, pharmcat
- Missing sources: none

`combined_data.json` is the authoritative nested representation.
`combined_summary.tsv` has one row per rsID and keeps every source object in a
JSON-valued column, alongside useful counts. `combined_data.sqlite` mirrors
those rows with one JSON column per source and one shared dbSNP gene/consequence
column. `run_manifest.json` records all consumed files, SHA-256 hashes,
embedded source manifests, and source-reported errors.

Every present source contains a validated identical copy of
`dbsnp_gene_consequences.json`. Its `genes` and `functional_consequences` are
independent record-level aggregates from dbSNP ESummary; they must not be read
as gene-to-consequence, allele-specific, or transcript-specific pairs.

One-to-many records are nested by their actual relation keys: genomic variant,
ClinVar record, GWAS/PheGenI source row, gnomAD variant/callset, GTR test, or
PharmCAT reference position. They are never flattened into a Cartesian
product. Source-specific cautions and evidence scopes remain authoritative.

The PharmCAT source reports only exact membership in the pinned official
GRCh38 reference-position catalog. An rsID-only panel contains no sample
genotype, phase, copy number, or quality, so this dataset does not call star
alleles, haplotypes, diplotypes, phenotypes, drug response, or prescribing
recommendations.
