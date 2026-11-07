# QIAGEN BioHack 2026 — 20-rsID Challenge Input

This update uses the **20 rsIDs supplied in the QIAGEN Bioinformatics Challenge screenshot** and runs the existing protein and JASPAR pipelines without assuming that an rsID uniquely identifies a single alternate allele.

> Scope: human genome **GRCh38**. The protein and motif-scoring branches support **single-nucleotide substitutions only**. Indels, merged IDs, unresolved IDs, and ambiguous mappings are reported instead of silently being scored as SNPs. The variant list is not a sample genotype file. It does not support personal risk prediction or APOE haplotype calls.

## 1. Files added/changed

| File | Purpose |
| --- | --- |
| `examples/qiagen_20_rsids.tsv` | Exact 20 IDs in the organizer's original order. No invented ALT alleles. |
| `examples/qiagen_20_rsids.txt` | Same inputs, one ID per line. |
| `examples/qiagen_20_rsids.csv` | Excel-friendly copy. |
| `challenge_batch.py` | **New:** resolve all IDs, expand valid SNV ALT alleles, run optional branches, summarize per rsID, flag unsupported variants. |
| `protein_branch.py` | **Updated:** optional `--vep-only` for faster annotation without UniProt/AlphaFold downloads. |
| `tests/test_challenge_batch.py` | **New:** tests for list fidelity, multi-allelic SNPs, indels, ambiguous mappings, mixed alleles and branch integration. |

The original entrypoints and example workflows remain available.

## 2. Install

```bash
# Run from the biohack_snp_toolkit directory
python -m pip install -r requirements.txt
```

The public Ensembl and UniProt/AlphaFold endpoints do not require a personal API key, but live runs require Internet connectivity.

## 3. Recommended workflow

### A. First resolve the 20 rsIDs

```bash
python challenge_batch.py --out results/qiagen --prepare-only --annotate-other
```

Inspect these outputs before spending time on structural analysis:

- `results/qiagen/challenge_summary.tsv` — **one row per provided rsID**, in the original order. The most convenient file for the team.
- `results/qiagen/resolved_snvs.tsv` — **one row per supported SNV ALT allele**; this may have more rows than the original 20-rsID input.
- `results/qiagen/unsupported_or_unresolved.tsv` — unsupported indel alleles, ambiguous coordinates, failed lookups or mixed-allele entries needing manual review.
- `results/qiagen/resolution_details.json` — full mapping evidence and errors.
- `results/qiagen/other_variants_vep/` — optional raw VEP consequence information for unresolved/non-SNV entries; **rsID level only**, not an allele-specific effect prediction.

### B. Fast hackathon run: VEP + AlphaMissense + a limited JASPAR pilot

```bash
python challenge_batch.py --out results/qiagen --vep-only --annotate-other
```

This invokes both branches but avoids expensive protein and 3D-structure fetches. **By default the JASPAR branch loads only the bundled CEBPA motif.** That is intended as a reproducible proof of concept, not a comprehensive TF scan. A zero-hit result here does **not** imply absence of regulatory activity.

### C. Full protein context and full JASPAR motif set

Download a JASPAR **CORE / Vertebrates / non-redundant / Single batch file / JASPAR format** from:

https://jaspar.elixir.no/downloads/

Save that file under `data/vertebrates.jaspar` (or any other local filename), then run:

```bash
python challenge_batch.py --out results/qiagen_full \
  --pfm-file data/vertebrates.jaspar --annotate-other
```

This attempts UniProt functional annotations, preexisting AlphaFold **reference** structures and pLDDT for supported *missense* transcript annotations, plus a broader regulatory motif scan. Live data requests and structure downloads can be slow; start with the fast run or a small subset first.

### D. Smoke-test a subset

```bash
python challenge_batch.py --max-rsids 3 --out results/first_three --vep-only
python challenge_batch.py --only-rsids rs334,rs429358,rs7412 \
  --out results/apoe_hbb --vep-only
```

The second command demonstrates annotation for the selected IDs; **it cannot call APOE ε2/ε3/ε4 genotypes** because this input contains neither individual genotypes nor haplotype phase.

### E. Run each branch separately after resolution

```bash
python protein_branch.py --input results/qiagen/resolved_snvs.tsv \
  --out results/protein_only --vep-only

python jaspar_branch.py --input results/qiagen/resolved_snvs.tsv \
  --pfm-file data/vertebrates.jaspar --out results/jaspar_only
```

Windows PowerShell: put each command on one line (or replace the backslash line continuations as appropriate).

## 4. Output fields for integration

### `challenge_summary.tsv` (one record per original rsID)

- `rank`, `rsid`: exact original input order and ID.
- `resolution_status`: `resolved_snv`, `resolved_snv_with_unhandled_alleles`, `unsupported_or_ambiguous`, or `lookup_failed`.
- `n_snv_alt_alleles`, `variant_keys`: how many single-base substitutions were successfully normalized, and their `GRCh38:chrom:position:REF>ALT` keys.
- `has_non_snv_alleles`, `var_class`, `most_severe_consequence`: initial triage information from Ensembl (the consequence is rsID-level metadata, not an exact ALT clinical assessment).
- `protein_status`, `n_missense_transcripts`: VEP annotation coverage. Zero may mean a nonmissense mutation, not that the variant has no effect.
- `jaspar_status`, `n_motif_hits`: motif pilot/broad-scan coverage. A zero is conditional on the motif collection and heuristic thresholds.
- `special_note`, `reason`: scientifically important cautions and error explanations.

### Per-ALT outputs

`protein/protein_results.tsv` and `jaspar/motif_results.tsv` use `variant_key` to match the **same specific ALT allele** across branches. `protein/variant_status.json` and `jaspar/variant_status.json` also preserve cases with zero hits or annotation errors. These are distinct from a clinical phenotype or disease diagnosis.

## 5. Specific issues in the 20-rsID list

1. **rs121913529 (KRAS)** is multiallelic. The runner expands resolvable SNV ALT alleles, rather than arbitrarily assigning one amino-acid substitution to the rsID.
2. **rs113993960 (CFTR p.Phe508del)** is a deletion allele. It belongs in a VEP/ClinVar/literature explanation, but is not eligible for SNV-only AlphaMissense or single-base JASPAR score changes.
3. **rs80357713** has been reported as a historical rsID merged into `rs80357783` and associated with an insertion. The code keeps the supplied ID untouched, flags lookup/mapping issues, and does **not** automatically replace an identifier without verifying the actual genomic allele.
4. **rs429358 and rs7412 (APOE)** define common APOE protein haplotypes together. A list of untyped rsIDs provides no genotype or phase; do not predict a person's APOE haplotype or disease risk from the list alone.

References: [Ensembl variation REST](https://rest.ensembl.org/documentation/info/variation_id), [VEP REST](https://rest.ensembl.org/documentation/info/vep_id_get), [JASPAR downloads](https://jaspar.elixir.no/downloads/), [CFTR p.Phe508del](https://pmc.ncbi.nlm.nih.gov/articles/PMC10454170/), [APOE genotype combinations](https://www.ncbi.nlm.nih.gov/books/NBK605938/table/lecanemab.T.common_apoe_allele_minimum_v/).

## 6. Test / validation status

- **34 local unit and synthetic integration tests passed** at packaging time.
- Tests include all 20 IDs exactly as transcribed, multi-ALT splitting, indel/non-SNV rejection with recorded statuses, ambiguous mapping handling, and summary generation.
- **Real Ensembl, JASPAR, UniProt and AlphaFold API end-to-end calls were not tested in the packaging environment** (no outbound Internet). API drift, response differences or per-ID special cases may need adjustment on a connected machine.
- No real annotations or predicted phenotypes are bundled or implied by the example IDs.

## 7. Important limits and privacy

All 20 IDs are public test **locus identifiers**, not genotypes from a patient. The pipeline enumerates possible alternate alleles at those loci, not observed alleles in a person. AlphaMissense is a model annotation rather than a variant-specific clinical diagnosis; AlphaFold pLDDT measures reference-structure model confidence, not pathogenicity; JASPAR motif score differences are hypotheses, not measured TF binding. Clinical interpretation requires variant- and condition-specific ClinVar, GWAS, functional and literature evidence with source references.

Use pre-reviewed public data in demonstrations and verify terms before publishing data or outputs. Keep token/API credentials out of Git repositories.
