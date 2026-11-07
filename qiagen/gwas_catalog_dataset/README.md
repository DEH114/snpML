# GWAS Catalog dataset for the QIAGEN 20-rsID panel

This folder is a reproducible snapshot of literature-curated top associations
from the NHGRI-EBI GWAS Catalog for the exact 20 public test rsIDs. It contains
no patient genotypes.

- Metadata-reported current release: 2026-10-29
- Archive SHA-256: d24750b05a62aafee399616b068286054d2fae4713d53003725e565910bd9a23
- Input rsIDs: 20
- rsIDs with one or more associations: 15
- Normalized rsID-association relationships: 3318
- Losslessly retained source rows: 3318
- Multi-variant/haplotype-or-interaction relationships: 1

## Files

`rsid_associations.tsv` is the analysis-friendly table with one relationship
between a panel rsID and a source association row. `matched_source_records.json`
retains every matching official row and all original strings without numeric
coercion; its TSV copy adds only provenance and match columns. A source row can
produce more than one relationship when it contains multiple panel rsIDs.

`panel_status.tsv` always contains all 20 inputs in their original order.
`dbsnp_resolution.tsv` records current identifiers used during matching.
`gwas_catalog_metadata.json`, `archive_provenance.json`, and `run_manifest.json`
capture metadata and file provenance. The metadata release is recorded
separately and is not asserted to identify a user-supplied or previously
cached archive; the archive hash is authoritative. `gwas_catalog.sqlite` mirrors the main
tables. `errors.json` records non-fatal lookup warnings.

`association_scope=multi_variant_haplotype_or_interaction` means the official
`SNPS` field contains multiple rsIDs. Such a set-level association must not be
assigned to one component SNP alone. Reported risk/effect alleles are retained
as published; they are not assumed to be forward-strand GRCh38 alternate
alleles and are not harmonized here.

GWAS association is not proof of causality, a diagnosis, or an individual risk
estimate. This archive contains curated top associations, not every result from
full genome-wide summary statistics.

Official sources:
- https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/
- https://www.ebi.ac.uk/gwas/docs/file-downloads/
- https://www.ebi.ac.uk/gwas/docs/fileheaders/
