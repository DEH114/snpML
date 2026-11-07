# Seven rsID Dataset Generators

This guide introduces the seven dataset-generation scripts in this toolkit:

1. `generate_alphafold_dataset.py`
2. `generate_jaspar_dataset.py`
3. `generate_clinvar_dataset.py`
4. `generate_gwas_catalog_dataset.py`
5. `generate_phegeni_dataset.py`
6. `generate_gnomad_dataset.py`
7. `generate_gtr_dataset.py`

Together they annotate an ordered list of public rsIDs with protein, regulatory,
clinical, association, population-frequency, and registered-test evidence. The
input is a list of loci, not an individual's genotype, phased haplotype, or
phenotype. The scripts do not produce a diagnosis or personal risk estimate.

## 1. Generator overview

| Script | Main resources | Dataset purpose | Default output folder |
| --- | --- | --- | --- |
| [`generate_alphafold_dataset.py`](generate_alphafold_dataset.py) | Ensembl variation/VEP, AlphaMissense, UniProt, AlphaFold DB | Resolve GRCh38 SNV alleles and annotate missense transcripts, protein context, and existing reference-protein structures | `results/qiagen_20_snp_dataset/alpha_fold_dataset/` |
| [`generate_jaspar_dataset.py`](generate_jaspar_dataset.py) | Ensembl variation/sequence and a local JASPAR PFM collection | Compare reference and alternate allele motif scores in the same sequence window | `results/qiagen_20_snp_dataset/jaspar_dataset/` |
| [`generate_clinvar_dataset.py`](generate_clinvar_dataset.py) | NCBI dbSNP and ClinVar | Resolve merged identifiers, follow explicit ClinVar links, and keep direct variant records separate from related haplotype/genotype records | `results/qiagen_20_snp_dataset/clinvar_dataset/` |
| [`generate_gwas_catalog_dataset.py`](generate_gwas_catalog_dataset.py) | NCBI dbSNP and the official NHGRI-EBI GWAS Catalog association archive | Retain exact rsID association relationships and lossless official source rows | `results/qiagen_20_snp_dataset/gwas_catalog_dataset/` |
| [`generate_phegeni_dataset.py`](generate_phegeni_dataset.py) | NCBI dbSNP and the official NCBI PheGenI association table | Retain exact rsID relationships to PheGenI's integrated phenotype-genotype associations, including their reported source | `results/qiagen_20_snp_dataset/phegeni_dataset/` |
| [`generate_gnomad_dataset.py`](generate_gnomad_dataset.py) | NCBI dbSNP and the gnomAD v4 GraphQL API | Verify rsID-to-allele mappings and export callset and population frequencies | `results/qiagen_20_snp_dataset/gnomad_dataset/` |
| [`generate_gtr_dataset.py`](generate_gtr_dataset.py) | NCBI dbSNP and Genetic Testing Registry Entrez records | Find registered laboratory tests indexed by an exact rsID search | `results/qiagen_20_snp_dataset/gtr_dataset/` |

## 2. Requirements

Use Python 3.11 or newer and install the project dependencies from the toolkit
root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

The online workflow accesses public services. No personal API key is required.
The JASPAR generator also requires a local JASPAR-format PFM collection; its
default path is `data/vertebrates.jaspar`.

## 3. Input contract

Every generator accepts either a tab-separated file with an `rsid` column:

```text
rank	rsid
1	rs334
2	rs429358
```

or a plain file containing one identifier per line:

```text
rs334
rs429358
```

Identifiers must be unique lowercase values matching `rs<number>`. Input order
is preserved throughout the status tables and combined dataset.

The supplied `rsid` remains the durable panel key. When dbSNP reports a merged
identifier, the scripts record `current_rsid` and `was_merged` separately; they
do not silently replace the supplied ID.

Direct calls to a generator require the original QIAGEN 20-rsID panel in its
original order by default. Add `--allow-custom-panel` for another validated
panel. The all-datasets runner enables custom-panel mode automatically.

## 4. Recommended: generate all seven datasets

Run the exact QIAGEN panel:

```powershell
.\.venv\Scripts\python.exe generate_seven_datasets.py `
  --input examples\qiagen_20_rsids.tsv `
  --out results\qiagen_20_snp_dataset
```

Run a new panel in a fresh output folder:

```powershell
.\.venv\Scripts\python.exe generate_seven_datasets.py `
  --input path\to\new_rsids.tsv `
  --out results\new_rsid_panel
```

[`generate_seven_datasets.py`](generate_seven_datasets.py) validates and
normalizes the input, invokes all seven scripts, captures their logs, and runs
[`combine_seven_datasets.py`](combine_seven_datasets.py) after all seven
sources succeed. `generate_all_datasets.py` remains a generic compatibility
entry point for the same seven-source workflow.

The seven-source runner also performs one shared dbSNP enrichment step. It
normalizes the dbSNP ESummary `genes` and `fxn_class` values once, then writes
the same byte-identical `dbsnp_gene_consequences.json` into every one of the
seven dataset folders before combination. This avoids seven independent dbSNP
lookups while giving every source folder the same ordered, panel-complete
annotation.

The execution order is:

1. ClinVar
2. GWAS Catalog
3. gnomAD
4. GTR
5. AlphaFold/protein
6. JASPAR
7. PheGenI
8. Combined dataset

Sequential execution is intentional. The scripts share a response cache and
use public services with rate limits. GWAS Catalog can also reuse ClinVar's
validated `dbsnp_resolution.json` instead of resolving the same panel again.

The suite root contains:

```text
<out>/
|-- panel_input.tsv
|-- all_datasets_run_manifest.json
|-- logs/
|-- alpha_fold_dataset/
|-- jaspar_dataset/
|-- clinvar_dataset/
|-- gwas_catalog_dataset/
|-- phegeni_dataset/
|-- gnomad_dataset/
|-- gtr_dataset/
`-- combined_dataset/
```

`all_datasets_run_manifest.json` records the input hash, exact child commands,
timestamps, durations, exit codes, source status, and combination status.
`logs/<step>.stdout.log` and `logs/<step>.stderr.log` capture each subprocess.

By default, one failed source does not stop the remaining generators, but the
suite exits nonzero and skips combination. Useful controls are:

- `--fail-fast`: stop launching later generators after the first failure.
- `--allow-partial-combine`: combine complete available source folders, while
  retaining an overall failed suite status.
- `--skip-combine`: run only the seven source generators.

Use a fresh output root for a different panel, especially for partial runs, so
older source files cannot be confused with the new run.

## 5. Run the generators independently

The default 20-rsID panel can be regenerated with:

```powershell
.\.venv\Scripts\python.exe generate_alphafold_dataset.py
.\.venv\Scripts\python.exe generate_jaspar_dataset.py
.\.venv\Scripts\python.exe generate_clinvar_dataset.py
.\.venv\Scripts\python.exe generate_gwas_catalog_dataset.py
.\.venv\Scripts\python.exe generate_phegeni_dataset.py
.\.venv\Scripts\python.exe generate_gnomad_dataset.py
.\.venv\Scripts\python.exe generate_gtr_dataset.py
```

For a custom list, specify the input, output, and opt-in flag explicitly. For
example:

```powershell
.\.venv\Scripts\python.exe generate_clinvar_dataset.py `
  --input path\to\new_rsids.tsv `
  --out results\new_rsid_panel\clinvar_dataset `
  --allow-custom-panel
```

All seven scripts share these options:

| Option | Meaning |
| --- | --- |
| `--input PATH` | rsID TSV or one-rsID-per-line input |
| `--out PATH` | Source-specific output directory |
| `--cache PATH` | Persistent response/download cache; default `.cache` |
| `--allow-custom-panel` | Accept a unique validated panel other than the original 20 IDs |
| `--offline` | Use cached content only; a cache miss is an explicit failure |
| `--refresh` | Refresh cached public data; incompatible with `--offline` |

Run `python <script> --help` for every source-specific option.

## 6. AlphaFold/protein dataset

The AlphaFold generator resolves each rsID into all supported non-reference
GRCh38 SNV alleles. It then runs VEP for the exact allele and retains all
missense transcript annotations. Optional enrichment adds AlphaMissense,
UniProt features, and an existing AlphaFold reference-protein model.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dataset_summary.json/.tsv`: one ordered row per input rsID.
- `resolution_details.json`: Ensembl resolution evidence for every input.
- `resolved_snvs.tsv`: one row per supported GRCh38 SNV allele.
- `unsupported_or_unresolved.tsv`: indels, ambiguous mappings, and unresolved
  inputs that could not enter the SNV branch.
- `variant_status.json`: one status object per resolved allele.
- `protein_results.json/.tsv`: one row per missense transcript annotation.
- `vep/`: raw allele-specific VEP responses.
- `structures/`: downloaded reference-protein PDB files when available.
- `dataset_manifest.json`, `run_manifest.json`, and `errors.json`.

The join key between allele and transcript records is `variant_key`, formatted
as `GRCh38:<chromosome>:<1-based-position>:<REF>><ALT>`.

Source-specific controls include `--all-transcripts`, `--vep-only`,
`--no-alphamissense`, and `--annotate-other/--no-annotate-other`.

AlphaFold files are existing reference-protein models, not mutant structure
predictions. pLDDT measures local model confidence, not pathogenicity or
variant effect. AlphaMissense is a precomputed model annotation, not a clinical
classification.

## 7. JASPAR dataset

The JASPAR generator uses the same conservative GRCh38 SNV preparation logic.
For each supported allele, it retrieves one reference sequence window, creates
the corresponding alternate sequence, and compares motif scores on the same
window and strand against the selected PFM collection.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dataset_summary.json/.tsv` and `resolution_details.json`.
- `resolved_snvs.tsv` and `unsupported_or_unresolved.tsv`.
- `variant_status.json`: one row per scanned allele.
- `motif_results.json/.tsv`: qualifying allele-matrix-window-strand hits.
- `motif_profiles.json`: the matrix definitions referenced by `matrix_id`.
- `sequences/`: auditable reference and alternate sequence contexts.
- A portable copy of the selected PFM file.
- `dataset_manifest.json`, `run_manifest.json`, and `errors.json`.

Key options and defaults are `--pfm-file data/vertebrates.jaspar`,
`--flank 40`, `--threshold 0.8`, `--min-delta 0.1`, and
`--pseudocount 0.5`.

Motif-score changes are computational hypotheses. They are not p-values, do
not prove transcription-factor binding or regulatory activity, and do not
identify a target gene. No qualifying hit does not exclude regulatory effects.

## 8. ClinVar dataset

The ClinVar generator first resolves current and merged identifiers through
dbSNP, follows explicit dbSNP-to-ClinVar links, and retrieves linked ClinVar
Variation summaries. It verifies current dbSNP cross-references and separates
direct variant records from haplotype, genotype, and other set-level records.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dbsnp_resolution.json/.tsv`: supplied-to-current identifier mapping.
- `rsid_clinvar_links.json/.tsv`: auditable rsID-to-ClinVar relationships.
- `clinvar_records.json/.tsv`: direct variant records.
- `related_set_records.json/.tsv`: related multi-variant records.
- `traits.json/.tsv`: trait and classification rows joined by rsID, ClinVar UID,
  and record scope.
- `panel_status.json/.tsv`: one row per input, including no-match states.
- `raw/`, `clinvar.sqlite`, `run_manifest.json`, and `errors.json`.

`--max-records-per-rsid` defaults to 1,000 and fails rather than silently
truncating an unexpectedly large result set.

ClinVar classifications are submitted, allele-specific, condition-specific,
and review-status-dependent. A set-level classification must not be assigned
to one component SNP. Absence from ClinVar is not evidence of benignity.

## 9. GWAS Catalog dataset

The GWAS generator streams the official ontology-annotated association ZIP
without extracting it. It matches exact rsID tokens against both supplied and
resolved current identifiers. Every matched official row is retained without
numeric coercion, and a separate normalized table records its relationship to
each panel rsID.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dbsnp_resolution.json/.tsv` and `panel_status.json/.tsv`.
- `rsid_associations.json/.tsv`: one panel-rsID-to-source-row relationship.
- `matched_source_records.json/.tsv`: lossless official source rows.
- `archive_provenance.json`: archive member, size, and SHA-256 evidence.
- `gwas_catalog_metadata.json`: optional release metadata.
- `gwas_catalog.sqlite`, `run_manifest.json`, and `errors.json`.

Relationships join to source rows through `source_row_number`; the
`source_record_sha256` value provides an additional content-integrity check.

Use `--archive PATH --no-metadata` for a fully local official ZIP. The default
download and uncompressed safety limits are 1 GB and 4 GB. An existing
`--dbsnp-resolution` file can be supplied explicitly; otherwise the script
first checks the sibling ClinVar dataset and then queries NCBI.

GWAS associations are statistical literature associations, not proof of
causality or personal risk. Reported risk alleles are preserved as published
and are not assumed to be the forward-strand GRCh38 ALT allele. Multi-rsID
associations remain set-level.

## 10. PheGenI dataset

The PheGenI generator filters the official NCBI Phenotype-Genotype Integrator
bulk TSV by exact supplied and resolved rsID. PheGenI integrates association
evidence from the NHGRI GWAS Catalog with selected, most-significant dbGaP
summary associations and NCBI context. This makes it a distinct integrated
view, but not an independent replication of GWAS Catalog evidence.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dbsnp_resolution.json/.tsv` and `panel_status.json/.tsv`.
- `phegeni_associations.json/.tsv`: normalized panel-rsID-to-source-row
  relationships.
- `matched_source_records.json/.tsv`: lossless matched bulk-table rows.
- `source_provenance.json`: bulk-table URL, retrieval time, size, and
  SHA-256 evidence.
- `phegeni.sqlite`, `run_manifest.json`, and `errors.json`.

Use `--table PATH` to filter a previously downloaded official TSV. The
generator preserves the PheGenI `Source` field so NHGRI-derived and dbGaP-
derived rows remain distinguishable. Chromosome and location values are
source-reported: the bulk export does not declare a coordinate assembly, so
the toolkit does not label those positions GRCh37 or GRCh38 and does not join
them to allele-specific coordinates without separate validation.

PheGenI reports association evidence, not causal variants, diagnoses, or
personal risk. Its selected dbGaP rows are summary results rather than all
study p-values or participant-level data. See the official
[PheGenI page](https://www.ncbi.nlm.nih.gov/gap/phegeni) and
[FAQ](https://www.ncbi.nlm.nih.gov/books/NBK570253/).

## 11. gnomAD dataset

The gnomAD generator resolves current dbSNP identifiers, searches the official
gnomAD v4 GraphQL service for every represented allele, and fetches each
variant detail record. A variant is accepted only when its returned dbSNP
cross-reference verifies the resolved current rsID.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dbsnp_resolution.json/.tsv` and `panel_status.json/.tsv`.
- `rsid_gnomad_links.json/.tsv`: verified rsID-to-variant relationships.
- `gnomad_variants.json/.tsv`: normalized allele identity.
- `callset_frequencies.json/.tsv`: separate exome, genome, and joint totals.
- `population_frequencies.json/.tsv`: population and sex-stratified rows.
- `source_provenance.json`, `raw/`, `run_manifest.json`, and `errors.json`.
- `gnomad.sqlite`, created only for a run with no lookup or identity errors.

The relationship chain is:

```text
input rsID
  -> rsid_gnomad_links.variant_id
  -> gnomad_variants.variant_id
  -> callset_frequencies.(variant_id, callset)
  -> population_frequencies.(variant_id, callset, population_id)
```

Population rows overlap. Exome, genome, joint, ancestry, and sex strata must
not be summed. Frequency is population-level evidence, not a patient genotype
or disease-risk estimate. Absence from gnomAD is not evidence of benignity.

## 12. GTR dataset

The GTR generator resolves supplied/current rsIDs and performs an exact quoted
`[WORD]` Entrez search. It then downloads and deduplicates registered-test
summaries while retaining the query-to-test relationship.

Important outputs include:

- `dbsnp_gene_consequences.json`: the ordered, shared dbSNP record-level
  gene/consequence annotation added by `generate_seven_datasets.py`.
- `dbsnp_resolution.json/.tsv` and `panel_status.json/.tsv`.
- `gtr_test_matches.json/.tsv`: input-rsID-to-test relationships.
- `gtr_tests.json/.tsv`: one normalized row per unique registered test.
- `conditions.json/.tsv` and `variant_targets.json/.tsv`.
- `raw/`, `gtr.sqlite`, `run_manifest.json`, and `errors.json`.

`--max-records-per-rsid` defaults to 1,000.

GTR is a laboratory-test registry, not a variant-classification database. An
index match is not a pathogenicity assertion, FDA approval, recommendation,
clinical-validity finding, or patient result. It also does not prove that the
registered test targets only the matched SNP.

## 13. Common output and audit rules

- `rank` and `rsid` preserve the supplied panel identity and order.
- `current_rsid` and `was_merged` record dbSNP identity changes separately.
- `panel_status.json/.tsv` or `dataset_summary.json/.tsv` preserves one ordered
  status row per input, including zero-hit and unresolved inputs.
- JSON is authoritative for nested values. TSV cells containing lists or
  objects are JSON-encoded strings.
- `errors.json` records row-, request-, or validation-level problems.
- Run and dataset manifests record parameters, counts, timestamps, source
  provenance, and cautions.
- Raw-response folders make retrieval and normalization auditable.
- A nonzero exit can still leave useful partial output. Always inspect the
  manifest, status table, and `errors.json` before using a dataset.

### Shared dbSNP gene/consequence annotation

Each of the seven source folders contains
`dbsnp_gene_consequences.json`. The seven copies are deliberately
byte-identical and contain one ordered record per supplied rsID. Each record
keeps `genes` and `functional_consequences` as separate arrays, along with the
supplied/current rsID identity and dbSNP provenance.

These arrays are record-level aggregates from dbSNP ESummary. They do **not**
state which consequence belongs to which gene, allele, or transcript. In
particular, consumers must not construct a gene-by-consequence Cartesian
product or describe the values as exact gene:consequence pairs. Use an
allele- and transcript-specific source such as VEP when that precision is
required.

The eight primary enriched JSON files are:

1. `alpha_fold_dataset/dbsnp_gene_consequences.json`
2. `jaspar_dataset/dbsnp_gene_consequences.json`
3. `clinvar_dataset/dbsnp_gene_consequences.json`
4. `gwas_catalog_dataset/dbsnp_gene_consequences.json`
5. `phegeni_dataset/dbsnp_gene_consequences.json`
6. `gnomad_dataset/dbsnp_gene_consequences.json`
7. `gtr_dataset/dbsnp_gene_consequences.json`
8. `combined_dataset/combined_data.json`

Run `generate_seven_datasets.py` to create all eight as part of the normal
workflow. Directly running only one source generator does not perform this
suite-level shared enrichment.

## 14. Cache, refresh, and offline reproduction

Complete the first run online with a persistent cache. The same panel can then
be regenerated without network access:

```powershell
.\.venv\Scripts\python.exe generate_seven_datasets.py `
  --input path\to\new_rsids.tsv `
  --out results\new_rsid_panel_offline `
  --cache .cache `
  --offline
```

Offline cache misses fail explicitly and never silently access the network.
Use `--refresh` when a new public-data snapshot is desired. Public databases
change over time, so preserve the run manifests and hashes with any analysis.

GWAS Catalog can instead use a local official ZIP, and PheGenI can use a local
official TSV through `--table`. Other sources require their cached API
responses for a completely offline run.

## 15. Combined dataset

Combine seven existing normalized source folders without any network request:

```powershell
.\.venv\Scripts\python.exe combine_seven_datasets.py `
  --input path\to\new_rsids.tsv `
  --dataset-root results\new_rsid_panel
```

`combined_dataset/combined_data.json` contains one top-level object per input
rsID with these source sections:

```text
rank
rsid
identity
dbsnp_gene_consequence
alpha_fold
jaspar
clinvar
gwas_catalog
gnomad
gtr
phegeni
```

The top-level `dbsnp_gene_consequence` object is the canonical combined copy
for that rsID. Its `genes` and `functional_consequences` remain separate
record-level aggregates; it does not assert exact gene/consequence pairing.

One-to-many data remain nested by allele, transcript, ClinVar record, GWAS or
PheGenI source row, gnomAD callset/population, motif hit, or GTR test. This
avoids inventing false Cartesian combinations between unrelated source
records.

The combined folder also contains `combined_summary.tsv`,
`combined_data.sqlite`, `run_manifest.json`, and a generated `README.md`.

## 16. Current validated 20-rsID snapshot

The cached 20-panel run completed successfully on 2026-11-07. At that snapshot:

| Dataset | Generated content |
| --- | --- |
| AlphaFold/protein | 34 resolved SNV alleles, 377 missense-transcript rows, and 10 reference PDB files |
| JASPAR | 1,019 motif profiles and 1,644 qualifying motif windows |
| ClinVar | 42 explicit links: 23 direct and 19 related-set records |
| GWAS Catalog | 3,318 normalized association relationships |
| PheGenI | 87 exact matched associations across 12 inputs: 82 NHGRI-derived and 5 dbGaP-derived rows |
| gnomAD | 28 verified allele links and 6,347 population/sex-stratified rows |
| GTR | 11 rsID-to-registered-test relationships |

The final combined output preserves the 20 ordered rsID records. Check the
current combined run manifest, SQLite integrity result, and local test run
before publishing a refreshed snapshot. Counts can change after refreshing
any public source.

## 17. Interpretation limits

These datasets integrate different evidence types; they are not interchangeable
labels:

- AlphaMissense is a prediction, while AlphaFold describes reference-protein
  structure confidence.
- JASPAR reports motif-score hypotheses, not experimentally measured binding.
- ClinVar reports submitted clinical interpretations tied to exact alleles and
  conditions.
- GWAS Catalog reports statistical associations that may tag rather than cause
  a trait.
- PheGenI integrates NHGRI and selected dbGaP association evidence; overlap
  with GWAS Catalog must not be counted as independent support.
- gnomAD reports aggregate population frequencies.
- GTR reports registered tests, not classifications or recommendations.

A missing record from any source does not establish benignity or absence of
biological or clinical relevance. Variant- and condition-specific clinical
interpretation requires qualified review of the underlying evidence.
