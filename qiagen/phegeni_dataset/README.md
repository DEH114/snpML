# PheGenI dataset for a 20-rsID panel

This folder is a reproducible, public association-data snapshot from NCBI's
Phenotype-Genotype Integrator (PheGenI). It contains no patient-level data.

- Source table SHA-256: f782842b5444c07a323dff287fb647b5bdbc59e2a60d35684347e13e8cd21722
- Input rsIDs: 20
- rsIDs with one or more associations: 12
- Normalized rsID-association relationships: 87
- Losslessly retained source rows: 87
- NHGRI-source relationships: 82
- dbGaP-source relationships: 5

## Files

`phegeni_associations.tsv` contains one normalized relationship between a panel
rsID and a matching PheGenI row. `matched_source_records.json` retains every
matching source row, all original columns, and all original strings without
numeric coercion. `panel_status.tsv` includes every input in original order,
including rsIDs with no match. `dbsnp_resolution.tsv` records current IDs used
for matching. `source_provenance.json` and `run_manifest.json` record hashes,
retrieval details, and counts. `phegeni.sqlite` mirrors the normalized tables.

Matching is exact: the complete numeric `SNP rs` value must equal the numeric
part of the supplied rsID or its current dbSNP ID. Substring, list, range, and
free-text matches are never inferred.

PheGenI is an integrator, not an independent seventh experiment. Its table
combines NHGRI GWAS Catalog records with selected top results from dbGaP, so
some NHGRI records can overlap the separate GWAS Catalog dataset. Do not add
counts across those datasets as if they were independent evidence. A reported
association is not proof of causality, a diagnosis, or an individual risk
estimate. PheGenI's bulk file does not expose an explicit release identifier;
the recorded SHA-256 is the authoritative snapshot identifier.

Official source:
https://www.ncbi.nlm.nih.gov/gap/phegeni
