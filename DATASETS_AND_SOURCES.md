# 公开 SNP 数据与接口来源

核对日期：2026-11-06。以下依据是用户提供的 Question 1 截图：从 SNP 列表整合公开数据与原始文献，解释功能、频率及疾病关联。未提供完整竞赛数据/IP规则，所以这里判断的是题目匹配度，不是对全部比赛条款的保证。

## 推荐数据组合

| 资源 | 最适合本项目的用途 | 建议下载/选择 | 限制 |
|---|---|---|---|
| ClinVar | 蛋白编码变异、已提交疾病分类、审核状态 | GRCh38 VCF；筛选 25 个 P/LP 与 25 个 B/LB 的高审核 missense SNV | 记录不都是无争议真值；保留疾病和审核证据；不要把检索标签当独立预测评估 |
| GWAS Catalog | SNP—性状统计关联与原始研究 PMID | all associations ontology-annotated TSV ZIP；先选 LDL/low density lipoprotein 相关 20 个 rsID | lead SNP 不等于因果位点；risk allele 不一定是当前参考基因组的 ALT |
| eQTL Catalogue | 调控支路的候选靶基因、组织、方向证据 | 按研究/组织下载 summary statistics 或 fine-mapped credible sets | eQTL/credible set 不是单独的功能机制证明；必须核对等位基因、组织与版本 |
| 1000 Genomes / IGSR | 群体背景频率及正常人群变异输入 | GRCh38 call set，或从 Ensembl 获取目标变异的 population allele frequencies | 不能推断每个样本“无病”；不是疾病标签集；部分历史 call set 使用 GRCh37 |

## 官方下载入口与许可

### ClinVar

- 下载指南：https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/
- GRCh38 VCF 目录：https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
- 本工具固定的已核实快照：https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261101.vcf.gz
- 当前滚动文件：https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
- TSV 说明：https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/
- 使用/署名政策：https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- 公开数据说明：https://www.ncbi.nlm.nih.gov/clinvar/intro/

ClinVar 可公开访问并下载；引用 ClinVar 和相应提交/论文来源。不是直接诊断工具。本包不下载或处理个人基因组。

### GWAS Catalog

- 官方下载说明：https://www.ebi.ac.uk/gwas/docs/file-downloads/
- 本工具使用的 ZIP：https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip
- 完整 summary statistics：https://www.ebi.ac.uk/gwas/downloads/summary-statistics
- 使用政策：https://www.ebi.ac.uk/gwas/docs/about/
- API 文档：https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/

一般 Catalog 数据适用 EMBL-EBI Services Terms of Use；summary statistics 通常为 CC0，个别研究例外必须检查。**不能把所有 Catalog 数据统一写成 CC0。** 当前 REST API 为 v2，v1 于 2026 年 8 月退役。本工具通过官方 ZIP 避免依赖旧 API。

### eQTL Catalogue

- 数据/下载说明：https://www.ebi.ac.uk/eqtl/Data_access/
- 研究列表：https://www.ebi.ac.uk/eqtl/Studies/
- 数据 FTP 根目录：https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/
- 许可：https://www.ebi.ac.uk/eqtl/License/

数据 CC-BY-4.0，代码 Apache-2.0。**官方目前明确旧 REST API 已弃用且不可用**；不要套用过时 API 教程。优先按官方 metadata 中的链接下载具体研究/组织的文件。ALT 为统计效应等位基因；仍需核对 assembly 和具体记录。官方建议避免频繁远程 tabix 请求，以免触发防护。

### 1000 Genomes / IGSR

- 数据入口：https://www.internationalgenome.org/data/
- 如何选择 call set：https://www.internationalgenome.org/faq/how-do-i-find-specific-variant-files/
- 单变异/频率查询说明：https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/
- Ensembl variation API：https://rest.ensembl.org/documentation/info/variation_id

选择 GRCh38 的 call set；不同 collection 有各自 data reuse policy，不能统一当作 CC0。本项目是注释几十个 SNP，不必下载所有样本的整个染色体 VCF。Ensembl variation endpoint 的 `pops=1` 可取得 population allele frequencies；需要按 population、具体 allele 和来源版本展示，而非混成一个统一频率。

## 支路接口与模型说明

- Ensembl VEP exact genomic variant endpoint：https://rest.ensembl.org/documentation/info/vep_region_get
- Ensembl VEP 参数：https://rest.ensembl.org/documentation/info/vep_id_get
- Ensembl 序列接口：https://rest.ensembl.org/documentation/info/sequence_region
- AlphaMissense 官方 VEP 插件源码：https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/main/AlphaMissense.pm
- UniProt REST 文档：https://www.uniprot.org/api-documentation
- AlphaFold API：https://alphafold.ebi.ac.uk/api-docs
- 实际 API schema 核对样例：https://alphafold.ebi.ac.uk/api/prediction/P68871
- AlphaFold 说明：https://alphafold.ebi.ac.uk/faq
- JASPAR API：https://jaspar.elixir.no/api/v1/docs/
- JASPAR 下载：https://jaspar.elixir.no/downloads/
- 随包 CEBPA PFM：https://jaspar.elixir.no/matrix/MA0102.5/
- JASPAR 使用 CC-BY-4.0，保留数据库和矩阵版本署名。

## 两个演示案例的来源

- **rs334**：ClinVar https://www.ncbi.nlm.nih.gov/clinvar/variation/15333/ 。GRCh38 `NC_000011.10:g.5227002T>A`，HBB `p.Glu7Val`。成熟蛋白的传统编号也写作 E6V；本代码依据完整参考蛋白序列，不混用编号。
- **rs12740374**：Musunuru et al., *From noncoding variant to phenotype via SORT1 at the 1p13 cholesterol locus*, Nature (2010). https://www.nature.com/articles/nature09266 。研究提供 C/EBP motif 和肝脏 SORT1 机制证据。该论文是独立解释来源，不是本工具 motif 评分自动证明的结论。

经典案例适合检查映射和展示，不构成独立模型 benchmark。对本次黑客松，优先报告信息来源、缺失情况、流程可靠性和解释的证据层级。
