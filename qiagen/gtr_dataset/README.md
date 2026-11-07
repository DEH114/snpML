# GTR registered-test matches for a 20-rsID panel

This snapshot searches NCBI's Genetic Testing Registry (GTR) for the supplied
public rsIDs. It contains registered laboratory test summaries, not patient
records or variant classifications.

- Input rsIDs: 20
- rsIDs with at least one indexed test match: 3
- rsID-to-test relationships: 11
- Unique registered tests: 11

`panel_status.tsv` contains all inputs in order. `gtr_test_matches.tsv`
preserves which exact quoted rsID search returned each test. `gtr_tests.tsv`
contains deduplicated test summaries; `conditions.tsv` and
`variant_targets.tsv` expand selected nested fields. JSON copies preserve
nested values, `raw/` preserves Entrez responses, and `gtr.sqlite` contains the
normalized tables.

GTR does not expose a dedicated rsID field or dbSNP-to-GTR link. Therefore the
relationship basis is an exact quoted `[WORD]` Entrez index match. An rsID may
not appear literally in ESummary even when the search index links the test.
The test may assay a gene, panel, or equivalent HGVS target rather than only
that SNP.

GTR content is submitted by laboratories and is not independently verified by
NIH. A match is not a pathogenicity assertion, FDA approval, recommendation,
or proof of clinical validity. Absence of a match is not evidence that testing
is unavailable.

Official access documentation:
https://www.ncbi.nlm.nih.gov/gtr/docs/maintenance_use/
