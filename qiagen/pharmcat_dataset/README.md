# PharmCAT reference-position dataset

This folder records which loci in an ordered 20-rsID
panel occur in the official PharmCAT 3.4.0 GRCh38
reference-position VCF.

- Source VCF SHA-256: 63f49b9004b6a5ca921deebe8ec59f8e45f6990543e07b7890c0d6201f5e69fc
- Input rsIDs: 20
- rsIDs with one or more reference-position matches: 2
- Normalized position relationships: 2
- Losslessly retained source rows: 2

`pharmcat_position_matches.tsv` contains exact supplied/current-rsID matches and
the source `PX` pharmacogene and `POI` flag. `matched_source_records.json`
preserves the complete original UTF-8 VCF source line, its exact line ending,
all fields, and a SHA-256 of the original bytes. `panel_status.tsv` contains
every input in original order. `pharmcat.sqlite` mirrors these normalized
tables, and the provenance and run manifests identify the pinned source.

## Scope limitation

This is a **locus-only reference-position index**, not a PharmCAT clinical
report. The `PharmCAT` sample column in the official reference VCF is a template
and is retained only for source fidelity; this generator never interprets its
`GT` value. An rsID list provides no patient genotype, phase, copy number, or
quality evidence. Consequently this dataset does not infer star alleles,
diplotypes, metabolizer phenotypes, drug response, or prescribing
recommendations. A missing position is not evidence of no pharmacogenomic
importance.

Official sources:
- https://github.com/PharmGKB/PharmCAT/releases/tag/v3.4.0
- https://github.com/PharmGKB/PharmCAT
