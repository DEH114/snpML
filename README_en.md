## NEW — QIAGEN BioHack 2026 20-rsID batch runner

This release adds `challenge_batch.py` and `examples/qiagen_20_rsids.tsv`. Start with **[`QIAGEN_CHALLENGE_QUICKSTART_EN.md`](QIAGEN_CHALLENGE_QUICKSTART_EN.md)** before using the original branch commands. It supports rsID-only batch input, multiallelic splitting, indel/merged-ID reporting, and original-order coverage summaries.

# BioHack SNP Functional Interpretation Toolkit

Two independent analysis modules for **QIAGEN BioTech Question 1**. This toolkit is **not** a complete SNP-to-disease diagnostic system.

- **Protein branch:** explicitly specified genomic allele → VEP → missense transcripts → precomputed AlphaMissense annotations → UniProt functional features → variant location and pLDDT in an existing AlphaFold reference-protein model.
- **JASPAR branch:** explicitly specified genomic allele → reference/alternate DNA sequences → differences in motif-match scores for the same genomic window and strand.
- **Dataset utilities:** download ClinVar or GWAS Catalog data and construct small demonstration inputs.

**Scope: human GRCh38, 1-based forward-strand genomic coordinates, and single-nucleotide substitutions involving A/C/G/T.** Indels, CNVs, haplotypes, and automatic liftover from GRCh37 are not supported. For rsID-only input, the scripts expand each rsID into all valid non-reference SNV alleles; this does **not** infer which allele a person carries.

## 1. Installation and Quick Start

Use a dedicated Python 3.11 or newer environment. Run the following commands from the extracted `biohack_snp_toolkit` directory:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python protein_branch.py --input examples/protein_input.tsv --out results/protein
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar
```

The analysis scripts use public APIs. They do not require Nebius, Brev, Base44, paid API keys, or a GPU.

The first online run requires access to Ensembl, UniProt, and AlphaFold DB. The JASPAR command above uses the single CEBPA motif included in the package, so it does not need to contact the JASPAR API. It still needs Ensembl to retrieve the reference sequence.

The protein example, `rs334 / ALT=A`, uses the genomic forward strand. The corresponding ClinVar record describes GRCh38 `11:5227002 T>A` and HBB `p.Glu7Val`. Do not interpret the transcript-level change `c.20A>T` as genomic forward-strand `A>T`.

The regulatory example, `rs12740374 / ALT=T`, is intended to check C/EBP motif changes. Actual scores depend on the selected matrix, sequence, and thresholds. The code does not hard-code a required positive result.

**Testing status recorded in the original package:** 27 local unit and mocked-workflow tests passed. The execution environment did not have network access, so live-API end-to-end testing was not completed and no real-SNP analysis results were generated. `SYNTHETIC` test data are not biological evidence. The bundled examples are input panels, not complete database downloads; use the download commands for larger datasets.

## 2. Input Formats

### Option A: rsID and genomic forward-strand ALT — recommended

```text
rsid	alt
rs334	A
```

### Option B: explicit coordinates and alleles

```text
chrom	pos	ref	alt	assembly	rsid
11	5227002	T	A	GRCh38	rs334
```

Use **tab-separated** files. Ready-to-use TSV examples are included in `examples/`.

For coordinate-based input, `chrom`, `pos`, `ref`, `alt`, and `assembly` are the authoritative variant identifiers. The optional `rsid` is only a label; the script does not separately verify that it matches those coordinates. Use Option A when automatic rsID resolution is needed.

### Option C: rsID only

Provide a file with an `rsid` column, or one rsID per line. The scripts explicitly expand different ALT alleles. Do not collapse all resulting alleles into one overall “pathogenicity” label for the rsID.

## 3. Protein Branch

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein
```

### Processing steps

1. Resolve the variant to a GRCh38 primary-chromosome location and verify REF. Ambiguous mappings, invalid alleles, and reference mismatches produce explicit errors.
2. Query VEP using the exact coordinates and ALT allele, retaining all missense transcripts. AlphaMissense parsing supports both nested REST JSON and flat fields. Missing values are not replaced with zero.
3. By default, select one missense transcript per gene for the slower UniProt and structure queries, prioritizing MANE Select and then canonical transcripts. Other missense transcripts remain in the results table. Use `--all-transcripts` to perform detailed annotation for every missense transcript.
4. Map functional features only when the **Ensembl translated sequence exactly matches the full-length UniProt sequence**, and verify the reference amino acid. An identical numeric position in different isoforms is not automatically treated as the same residue.
5. Retrieve the existing PDB URL returned by the AlphaFold API rather than hard-coding a model version such as `model_v4` or `model_v6`. Accept only full-length matching models, verify the complete PDB chain sequence, and use the actual chain and residue sequence number (`resSeq`) to locate the residue.
6. Return local pLDDT, overlapping functional annotations, the reference-model download URL, and a local PDB file. Missing models or functional information retain explicit status fields; the script does not guess missing results.

### Output files

| File | Contents |
|---|---|
| `protein_results.tsv` | One row per genomic allele × missense transcript; suitable for a frontend. |
| `protein_results.json` | Full nested functional annotations, protein-level disease background, and structure records. |
| `variant_status.json` | Variant-level status, including no missense consequence, failure, and success. |
| `vep/*.json` | Raw VEP responses. |
| `structures/*.pdb` | Existing AlphaFold models that passed full-sequence verification. |
| `errors.json` | Recorded failures and their causes. |
| `run_manifest.json` | Request sources, retrieval timestamps, and content hashes. |

### How to interpret the results

- `vep_am_pathogenicity` is a **precomputed AlphaMissense prediction returned by VEP**. It is neither a newly trained model from this toolkit nor a patient's probability of developing a disease.
- `vep_am_class` is a model label, **not** a ClinVar classification. VEP's default AlphaMissense matching does not guarantee that the prediction's source transcript is the same as the displayed transcript. The output explicitly records this limitation; do not describe it as an exact isoform-specific prediction.
- `plddt` measures **local confidence in the reference structure**, not the effect of the variant.
- `protein_level_disease_annotations` provide protein-level background only. They must not be treated as disease evidence for the specific SNP.
- The script **does not predict mutant structures, RMSD, ΔΔG, clinical phenotypes, drug responses, or disease risk**.

If Ensembl returns an AlphaMissense plugin error, explicitly disable AlphaMissense and rerun. The other annotations can still be useful:

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein_basic --no-alphamissense
```

## 4. JASPAR Branch

Start with the real, versioned CEBPA profile included in the package:

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar
```

The 40 counts in this position frequency matrix (PFM) come from the official JASPAR MA0102.5 page. Source information and CC-BY-4.0 attribution are recorded in the accompanying `.source.json` file.

Alternatively, select current CORE vertebrate profiles online by transcription-factor name:

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --motif-names CEBPA CEBPB --out results/jaspar_online
```

Or request a specific matrix version:

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --matrix-ids MA0102.5 --out results/jaspar_versioned
```

These commands screen only the specified transcription factors. **They do not scan the entire JASPAR collection.**

For a broader screen, use the JASPAR downloads page and select **CORE → Vertebrates → non-redundant → Single batch file → JASPAR**. Pass the downloaded file to the script:

```bash
python jaspar_branch.py --input data/gwas_ldl.tsv --pfm-file data/vertebrates.jaspar --out results/jaspar_full
```

JASPAR is a **motif database**, not a SNP dataset. Vertebrate matrices are not necessarily all derived from humans; check the motif's species, TF family, and relevant tissue context.

### Scoring definitions

- Add a pseudocount of **0.5 per base per position**. Use background probabilities of 0.25 for each of A, C, G, and T.
- Compute the position-specific scoring matrix (PSSM) as `log2(position probability / background probability)`.
- For each motif, compare REF and ALT only within the **same genomic window overlapping the SNP**, testing both strands.
- `delta_log2_odds = ALT_score - REF_score`.
- `relative_score = (score - matrix_min) / (matrix_max - matrix_min)`.
- By default, retain windows with `max(REF relative score, ALT relative score) >= 0.8` and `abs(delta relative score) >= 0.1`. A gain or loss means crossing the manually specified match threshold.

**The 0.8 and 0.1 thresholds are demonstration filters, not calibrated statistical significance thresholds.** The script does not compute p-values or FDR. Scores should not be directly treated as comparable effect sizes across different motifs.

Use tissue context, chromatin accessibility, TF expression, and eQTL or experimental results to evaluate candidate mechanisms. No motif hit does not establish the absence of a regulatory effect.

### Output files and interpretation

The script produces `motif_results.tsv`, `motif_results.json`, `variant_status.json`, `motif_profiles.json`, `sequences/*.json`, error records, and a provenance manifest.

`target_gene` is deliberately left empty. The script does not treat the motif's TF name as the regulated target gene and does not automatically choose the nearest gene.

The protein and JASPAR scripts are independent modules. The JASPAR script does not itself establish that an input variant lies in an active regulatory element.

## 5. Downloading and Building Small Datasets

### ClinVar: the first choice for the protein branch

For reproducibility, the downloader is pinned to the GRCh38 VCF snapshot dated **2026-11-01**, recorded as verified in the original documentation, rather than silently switching to a newer release on each run.

```bash
python datasets.py download --source clinvar --out data/clinvar_GRCh38.vcf.gz
python datasets.py clinvar --vcf data/clinvar_GRCh38.vcf.gz --missense-only --n-per-class 25 --out data/clinvar_demo.tsv
python protein_branch.py --input data/clinvar_demo.tsv --out results/clinvar_demo
```

The selection criteria are:

- GRCh38 single-nucleotide substitutions from biallelic records.
- At least one missense consequence.
- Non-conflicting pathogenic/likely pathogenic or benign/likely benign classifications.
- Review status indicating agreement among multiple submitters, expert-panel review, or practice guidelines.

Reservoir sampling is applied across the complete file rather than taking its first few dozen rows, avoiding that source of chromosome-selection bias. If too few records meet the criteria, the output reports the actual number selected.

The ClinVar VCF `ID` field is a **Variation ID**, not an rsID. The code reads rsIDs from `INFO/RS` and uses explicit coordinates when an rsID is unavailable. Multiallelic VCF records are excluded to avoid incorrectly assigning an aggregated classification to one ALT allele.

This panel is intended for development, evidence retrieval, and demonstration. Do not supply ClinVar labels to a summarization model and then use those same labels to claim independent predictive accuracy. The filtered panel is also not a random population sample.

### GWAS Catalog: an entry point for traits and publications

```bash
python datasets.py download --source gwas --out data/gwas_associations.zip
python datasets.py gwas --input data/gwas_associations.zip --trait "LDL|low density lipoprotein" --n 20 --out data/gwas_ldl.tsv
```

The script selects inputs only from association records containing a single rsID. It saves the complete original rows—including trait, risk allele, p-value, and PMID—to `.gwas_records.json`. It does not split haplotype or SNP × SNP interaction associations.

The input TSV contains only rsIDs. Do not make allele-specific interpretations until the GWAS risk allele has been aligned with the current GRCh38 REF/ALT representation.

A GWAS lead SNP is not necessarily causal and is not guaranteed to produce a detectable motif change. Start with the two individual demonstration SNPs before expanding to dozens of GWAS SNPs. Do not treat the motif-hit rate as disease-prediction accuracy.

### eQTL Catalogue and 1000 Genomes

See [DATASETS_AND_SOURCES_EN.md](DATASETS_AND_SOURCES_EN.md). eQTL Catalogue can provide independent, tissue-relevant candidate target-gene evidence; 1000 Genomes can provide population-frequency context.

The current package does not implement complete eQTL or 1000 Genomes integration. This avoids incorporating gene links or risk-allele directions into the output before they have been checked.

## 6. Caching, Reproducibility, and Troubleshooting

Successful requests are automatically cached in `.cache`. Subsequent runs with the same parameters can use `--offline`. A cache miss produces an explicit error rather than silently making a network request. `--refresh` requests updated data and cannot be combined with `--offline`.

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein_offline --offline
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar_offline --offline
```

Before the live demonstration, complete a real online run, inspect `variant_status.json` and `errors.json`, and then verify an offline rerun.

An empty TSV does not automatically mean “no biological effect.” Possible explanations include no missense consequence, no motif passing the thresholds, a service error, or failure to obtain an exact sequence mapping.

The package does not contain full ClinVar/GWAS databases, real patient data, or API credentials. Do not commit `.cache`, large datasets, private data, or credentials to a public repository. Retain attribution to the source databases, motif matrices, and original publications when publishing the project.

## 7. Interface for the Frontend or Summarization Module

Use `variant_key = GRCh38:chrom:pos:REF>ALT` as the primary join key, and retain transcript identifiers and provenance.

At minimum, separate the information supplied to a summarization model into these categories:

1. **Existing database or literature evidence:** ClinVar, GWAS, or eQTL records, with the condition/trait, allele, tissue, and publication specified.
2. **Computational predictions:** AlphaMissense or JASPAR results, with the model/matrix and scope of applicability stated.
3. **Structural and protein background:** reference-structure location, functional features, and protein-level disease annotations.
4. **Missing, not applicable, or failed requests:** preserve these as distinct statuses.

The summarization prompt must not promote computational predictions or structural/protein background into validated variant-to-phenotype conclusions.

The code can be incorporated into a Streamlit application or backend service. The current package does **not** include a web UI, automated PubMed retrieval, or LLM-generated summaries.
