# gnomAD dataset for a 20-rsID panel

This folder is an allele-frequency snapshot from the official gnomAD GraphQL
API for the supplied public rsIDs. dbSNP was used first to resolve merged
identifiers. It contains no patient genotypes.

- Input rsIDs: 20
- Merged input rsIDs resolved by dbSNP: 1
- gnomAD variant links: 28
- Unique gnomAD variants: 28
- Input rsIDs with at least one gnomAD variant: 20

## Files

`dbsnp_resolution.tsv` preserves each supplied ID and its current dbSNP ID.
`rsid_gnomad_links.tsv` records the verified rsID-to-gnomAD-allele mapping.
`gnomad_variants.tsv` contains allele identity, while
`callset_frequencies.tsv` keeps exome, genome, and joint totals separate.
`population_frequencies.tsv` expands the population rows and calculates AF as
AC/AN. JSON copies preserve the same tables, `panel_status.tsv` contains all
inputs in original order, `raw/` preserves GraphQL and dbSNP responses, and
`source_provenance.json` records retrieval metadata. `gnomad.sqlite` is only
created for a run with no lookup or identity-verification errors.

Population rows overlap: ancestry totals, XX/XY subdivisions, overall XX/XY,
and some HGDP/1000 Genomes rows are not independent groups. Never sum these
rows. For X/Y variants, retain hemizygote counts and do not infer AN as twice
the sample count. Mitochondrial variants use a different gnomAD API schema and
are not supported by this normalizer.

An absent gnomAD record is not evidence that an allele is benign or absent in
all humans. Frequencies, filters, flags, coverage, ancestry composition, and
the exact normalized allele must be considered together.

Official sources:
- https://github.com/broadinstitute/gnomad-browser/blob/main/graphql-api/README.md
- https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1
