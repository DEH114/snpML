# ClinVar dataset for the QIAGEN 20-rsID panel

This folder contains an NCBI dbSNP-to-ClinVar snapshot for the exact 20 public
test rsIDs. It contains no patient genotypes.

- Input rsIDs: 20
- Merged input rsIDs resolved by dbSNP: 1
- Explicit rsID-to-ClinVar links: 42
- Direct ClinVar records: 23
- Related haplotype/genotype records: 19
- Linked sets without returned component xrefs: 2
- Input rsIDs with at least one direct record: 16

## Files

`dbsnp_resolution.tsv` maps each supplied ID to its current dbSNP ID. The
supplied ID remains the primary panel key. `rsid_clinvar_links.tsv` is the
auditable link table from each panel ID to every linked ClinVar Variation ID.
`clinvar_records.tsv` contains direct single-variant records.
`related_set_records.tsv` contains multi-variant records in which the rsID is
only one component; do not assign their set-level classifications to that SNP
alone. `panel_status.tsv` keeps all 20 inputs in original order. JSON copies
preserve nested values, `raw/` preserves API responses, and `clinvar.sqlite`
contains `dbsnp_resolution`, `links`, `records`, `panel_status`, `traits`, and
`metadata` tables.

`traits.tsv` expands the trait labels returned in each VCV-level ESummary
record. Its classification fields are VCV aggregates repeated per trait; they
are not condition-specific RCV assertion records.

`allele_role=identity_reference` marks a linked ClinVar record whose canonical
SPDI has identical deleted and inserted sequence. Such a record is retained
for provenance but is not a changed alternate allele. `variant_keys` are only
created for changed single-base GRCh38 alleles.

Some explicit ELink relationships lead to a haplotype/genotype record for
which ClinVar ESummary does not return component objects. Those records remain
in `related_set_records` with `record_scope=linked_set_without_component_xref`
and `exact_current_dbsnp_xref=false`; this preserves the link without treating
its set-level classification as an allele-level interpretation.

ClinVar reports submitted interpretations and review status. Germline,
clinical-impact, and oncogenicity classifications remain separate. These data
are not a diagnosis, a personal risk estimate, or an independently verified
truth set. Interpret every classification with its exact allele, condition,
review status, and evaluation date.

Sources: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/ and
https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/
