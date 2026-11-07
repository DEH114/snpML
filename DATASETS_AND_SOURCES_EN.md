# Public SNP Datasets and Interface Sources

**Source-check date recorded in the original documentation: November 6, 2026.** This English version preserves those source records; it does not represent a new live verification.

The recommendations below are based on the user-provided Question 1 screenshot: integrate public data and original publications for a list of SNPs to explain function, frequency, and disease associations. The complete competition data and intellectual-property rules were not provided. Accordingly, this document assesses relevance to the challenge, not compliance with every competition rule.

## Recommended Dataset Combination

| Resource | Best use in this project | Suggested download or selection | Limitations |
|---|---|---|---|
| ClinVar | Protein-coding variants, submitted disease classifications, and review status. | GRCh38 VCF; select 25 pathogenic/likely pathogenic and 25 benign/likely benign missense SNVs with stronger review status. | Records are not universally undisputed ground truth. Retain condition and review evidence. Do not treat retrieved labels as an independent predictive evaluation. |
| GWAS Catalog | Statistical SNP–trait associations and PMIDs for original studies. | Ontology-annotated all-associations TSV ZIP; start with 20 LDL/low-density-lipoprotein-related rsIDs. | A lead SNP is not necessarily causal. The risk allele is not necessarily ALT in the current reference assembly. |
| eQTL Catalogue | Candidate target-gene, tissue, and effect-direction evidence for the regulatory branch. | Download summary statistics or fine-mapped credible sets by study and tissue. | An eQTL or credible set alone does not prove a functional mechanism. Check alleles, tissue, and assembly/version. |
| 1000 Genomes / IGSR | Population-frequency context and population-variation inputs. | A GRCh38 call set, or population allele frequencies for target variants retrieved through Ensembl. | Do not infer that each sample is disease-free. This is not a disease-label dataset. Some historical call sets use GRCh37. |

## Official Download Sources and Licensing

### ClinVar

- [Download guide](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/)
- [GRCh38 VCF directory](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)
- [Pinned snapshot used by this toolkit](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261101.vcf.gz)
- [Rolling VCF file](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz)
- [TSV/FTP documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
- [Use and attribution policy](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
- [Public-data overview](https://www.ncbi.nlm.nih.gov/clinvar/intro/)

ClinVar is publicly accessible and downloadable. Cite ClinVar and the relevant submissions or original publications. It is not a stand-alone diagnostic tool. This package does not download or process personal genomes.

### GWAS Catalog

- [Official download documentation](https://www.ebi.ac.uk/gwas/docs/file-downloads/)
- [Association ZIP used by this toolkit](https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip)
- [Full summary statistics](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
- [About and terms of use](https://www.ebi.ac.uk/gwas/docs/about/)
- [API documentation](https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/)

General Catalog data are governed by the EMBL-EBI Services Terms of Use. Summary statistics are generally CC0, but individual studies may have exceptions that must be checked. **Do not describe all GWAS Catalog data as CC0.**

At the source-check date recorded above, the supported REST API was v2; v1 had been retired in August 2026. This toolkit uses the official ZIP download rather than depending on the old API.

### eQTL Catalogue

- [Data access and download documentation](https://www.ebi.ac.uk/eqtl/Data_access/)
- [Study list](https://www.ebi.ac.uk/eqtl/Studies/)
- [FTP root directory](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/)
- [License information](https://www.ebi.ac.uk/eqtl/License/)

Data are licensed under CC-BY-4.0, and code under Apache-2.0. **The original source check recorded that the legacy REST API was deprecated and unavailable.** Do not rely on outdated API tutorials. Prefer downloading study/tissue-specific files using the links in the official metadata.

ALT is the statistical effect allele; the assembly and individual records still need to be checked. The official guidance recommends avoiding frequent remote tabix requests, which may trigger protective restrictions.

### 1000 Genomes / IGSR

- [Data portal](https://www.internationalgenome.org/data/)
- [Choosing variant call files](https://www.internationalgenome.org/faq/how-do-i-find-specific-variant-files/)
- [Single-variant and frequency queries](https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/)
- [Ensembl variation API](https://rest.ensembl.org/documentation/info/variation_id)

Select a GRCh38 call set. Different collections have their own data-reuse policies; do not assume they are all CC0.

For a project annotating a few dozen SNPs, downloading whole-chromosome VCFs for all samples is unnecessary. The Ensembl variation endpoint can return population allele frequencies with `pops=1`. Present frequencies by population, specific allele, and source/version rather than merging them into one universal frequency.

## Branch Interfaces and Model Documentation

- [Ensembl VEP exact-genomic-variant endpoint](https://rest.ensembl.org/documentation/info/vep_region_get)
- [Ensembl VEP parameters](https://rest.ensembl.org/documentation/info/vep_id_get)
- [Ensembl sequence endpoint](https://rest.ensembl.org/documentation/info/sequence_region)
- [Official AlphaMissense VEP plugin source](https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/main/AlphaMissense.pm)
- [UniProt REST documentation](https://www.uniprot.org/api-documentation)
- [AlphaFold API documentation](https://alphafold.ebi.ac.uk/api-docs)
- [Example response used for the original API-schema check](https://alphafold.ebi.ac.uk/api/prediction/P68871)
- [AlphaFold FAQ](https://alphafold.ebi.ac.uk/faq)
- [JASPAR API documentation](https://jaspar.elixir.no/api/v1/docs/)
- [JASPAR downloads](https://jaspar.elixir.no/downloads/)
- [Bundled CEBPA PFM: MA0102.5](https://jaspar.elixir.no/matrix/MA0102.5/)

JASPAR uses CC-BY-4.0. Preserve attribution to the database and the specific matrix version.

## Sources for the Two Demonstration Cases

### rs334

Source: [ClinVar Variation ID 15333](https://www.ncbi.nlm.nih.gov/clinvar/variation/15333/).

GRCh38 `NC_000011.10:g.5227002T>A`; HBB `p.Glu7Val`. Traditional mature-protein numbering also describes this substitution as E6V. The code uses the full reference-protein sequence and does not mix numbering conventions.

### rs12740374

Source: Musunuru et al., *From noncoding variant to phenotype via SORT1 at the 1p13 cholesterol locus*, Nature (2010). [Original publication](https://www.nature.com/articles/nature09266).

The study provides evidence for the C/EBP motif and liver SORT1 mechanism. The publication is an independent source for interpretation, not a conclusion automatically established by this toolkit's motif score.

These established cases are suitable for checking mappings and demonstrating the workflow. They do not constitute an independent model benchmark. For the hackathon, prioritize reporting provenance, missing information, workflow reliability, and the evidence level supporting each interpretation.
