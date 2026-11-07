## 更新：QIAGEN BioHack 2026 出题方 20 个 rsID

新版增加 `challenge_batch.py` 和 `examples/qiagen_20_rsids.tsv`。建议先阅读 **[`QIAGEN_CHALLENGE_使用说明.md`](QIAGEN_CHALLENGE_使用说明.md)**，再使用下方原有的单支路命令。新增多等位基因展开、indel/旧 rsID 问题报告和按输入顺序汇总结果。

# BioHack SNP 功能解释工具包

面向 QIAGEN BioTech Question 1 的两个独立分析模块；不是完整的 SNP→疾病诊断系统。

- **蛋白支路**：明确的基因组等位基因 → VEP → 错义转录本 → AlphaMissense 预计算注释 → UniProt 功能位点 → 现有 AlphaFold 参考蛋白模型中的位置与 pLDDT。
- **JASPAR 支路**：明确的基因组等位基因 → 参考/替代 DNA 序列 → 同一窗口、同一方向的 motif 匹配差异。
- **数据工具**：下载 ClinVar / GWAS Catalog，构建小型演示输入。

**范围：人类 GRCh38、1-based 正链坐标、单碱基替换 A/C/G/T。** 不支持 indel、CNV、单倍型、从 GRCh37 自动 liftover。rsID-only 输入会拆出该 rsID 的所有合法非参考 SNV 等位基因；这不是推断用户携带哪一个等位基因。

## 1. 安装与最快运行路径

建议使用独立 Python 3.11 或更新环境。从解压后的 `biohack_snp_toolkit` 目录执行：

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python protein_branch.py --input examples/protein_input.tsv --out results/protein
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar
```

这些命令使用公开 API，不需要 Nebius、Brev、Base44 或付费 API key，也不需要 GPU。
第一次运行需要网络访问 Ensembl、UniProt 和 AlphaFold DB。第二条支路的上述命令使用随包附带的单个 JASPAR CEBPA motif，不需要访问 JASPAR API；仍需 Ensembl 获取参考序列。

蛋白示例 `rs334 / ALT=A` 使用基因组正链；ClinVar 对应 GRCh38 `11:5227002 T>A`、HBB p.Glu7Val。不能把转录本上的 `c.20A>T` 当成基因组正链 `A>T`。
调控示例 `rs12740374 / ALT=T` 用于检查 C/EBP motif 变化；实际评分由所选矩阵、序列和阈值决定，代码不硬编码“必须检出”的结果。

**测试范围：27 项本地单元/模拟工作流测试已通过。执行环境无法联网，因此没有完成真实 API 端到端测试，未生成真实 SNP 分析结果。测试中的 SYNTHETIC 数据不是生物学证据。** 随包示例是输入面板，不是全量数据库；大型数据需要运行下载命令。

## 2. 输入格式

方式 A：明确 rsID 和基因组正链 ALT（推荐）。

```text
rsid	alt
rs334	A
```

方式 B：明确坐标和等位基因。

```text
chrom	pos	ref	alt	assembly	rsid
11	5227002	T	A	GRCh38	rs334
```

真实文件应使用制表符。`examples/` 中已提供可直接使用的 TSV。坐标输入以 `chrom,pos,ref,alt,assembly` 为权威标识；可选 rsid 只是标签，不会再次证明该 rsID 与坐标一致。需要自动解析 rsID 时使用方式 A。

方式 C：只有 rsID 列或每行一个 rsID。脚本会明确拆出不同 ALT，不能把所有结果汇总为一个 rsID 的“统一致病性”。

## 3. 蛋白支路

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein
```

实现细节：

1. 解析 GRCh38 主染色体位置，核对 REF；歧义多映射、非法等位基因和 REF 不匹配均报错。
2. VEP 按具体坐标与 ALT 请求，保留所有错义转录本。AlphaMissense 支持嵌套 REST JSON 和平铺字段，不把缺失值改成 0。
3. 默认每个基因在错义转录本中优先 MANE Select，再 canonical，选择一个进行较慢的 UniProt/结构查询。其他转录本仍在结果表中。`--all-transcripts` 可对所有错义转录本做深度注释。
4. 只在 **Ensembl 翻译序列与 UniProt 全长序列完全一致** 时映射功能位点；同时核对参考氨基酸。不自动把不同 isoform 的同一个数字位置视为同一残基。
5. 读取 AlphaFold API 返回的现有 PDB 地址，不硬编码 model_v4/v6；只接收全长匹配模型。再按完整 PDB 链序列核对，用实际 chain/resSeq 定位。
6. 输出局部 pLDDT、重叠功能注释、参考模型下载地址和本地 PDB。模型/功能信息缺失会保留状态，不猜测结果。

输出：

- `protein_results.tsv`：每行一个基因组等位基因×错义转录本，适合交给前端。
- `protein_results.json`：完整嵌套功能、蛋白层面疾病背景及结构记录。
- `variant_status.json`：包括无错义后果、失败与成功状态。
- `vep/*.json`：原始 VEP 输出。
- `structures/*.pdb`：通过全序列核对的现有 AlphaFold 模型。
- `errors.json` / `run_manifest.json`：失败原因、请求来源、检索时间及内容哈希。

结果解释：

- `vep_am_pathogenicity` 是 VEP 返回的 AlphaMissense 预计算预测，不是本代码的新模型，也不是患者患病概率。
- `vep_am_class` 是模型标签，不是 ClinVar 分类。VEP 的默认 AlphaMissense 匹配不保证预测来源转录本与当前展示转录本相同；输出中明确标记这一范围限制。不能将其称为精确 isoform 特异预测。
- `plddt` 是参考结构的局部预测置信度，不是变异效应。
- `protein_level_disease_annotations` 仅为蛋白层面背景，不能变成这个具体 SNP 的疾病证据。
- 脚本**没有预测突变结构、RMSD、ΔΔG、临床表型、药物反应或疾病风险**。

若 Ensembl 返回 AlphaMissense 插件错误，可显式关闭重跑，其他注释仍有用途：

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein_basic --no-alphamissense
```

## 4. JASPAR 支路

优先先跑随包的一个真实、带版本的 CEBPA profile：

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar
```

该 PFM 的 40 个计数来自 JASPAR 官方 MA0102.5 页面，来源与 CC-BY-4.0 署名记录在同目录的 `.source.json` 文件。

也可联网按 TF 名称选择当前 CORE vertebrate profiles：

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --motif-names CEBPA CEBPB --out results/jaspar_online
```

或者指定版本：

```bash
python jaspar_branch.py --input examples/regulatory_input.tsv --matrix-ids MA0102.5 --out results/jaspar_versioned
```

以上都只是指定 TF 的筛选，**不是整个 JASPAR 的全量扫描**。要扩展，在 JASPAR 下载页选择 CORE → Vertebrates → non-redundant → Single batch file → JASPAR，并传入下载后的文件：

```bash
python jaspar_branch.py --input data/gwas_ldl.tsv --pfm-file data/vertebrates.jaspar --out results/jaspar_full
```

JASPAR 是 motif 数据库，不是 SNP 数据集。脊椎动物矩阵也不等于全部来自人类；请检查 motif 的物种、TF 家族和组织背景。

计算定义：

- 每碱基/每位置加 0.5 pseudocount，背景 A/C/G/T 各 0.25。
- PSSM = log2(位置概率 / 背景概率)。
- 对每个 motif，仅比较**覆盖 SNP 的同一基因组窗口**在 REF/ALT 下的分数；两条链均测试。
- `delta_log2_odds = ALT_score - REF_score`。
- `relative_score = (score - matrix_min) / (matrix_max - matrix_min)`。
- 默认保留 max(REF relative, ALT relative) ≥0.8 且 |Δrelative| ≥0.1 的窗口；gain/loss 指跨过这个人为设置的阈值。

**0.8 和 0.1 是演示筛选参数，不是校准的显著性阈值。没有计算 p 值/FDR，不宜跨 motif 把分数直接当作统一效应量排序。** 使用组织背景、开放染色质、TF 表达、eQTL/实验结果验证候选机制。无命中不表示无调控效应。

输出：`motif_results.tsv/json`、`variant_status.json`、`motif_profiles.json`、`sequences/*.json`、错误和来源清单。`target_gene` 故意保留空值；本代码不会把 motif 的 TF 名字当成被调控的靶基因，也不会默认选择最近基因。

这两个脚本是独立模块；JASPAR 脚本不自动证明输入位点处于有效调控元件内。

## 5. 下载和构建小数据集

### ClinVar：蛋白支路优先

脚本固定到已核实存在的 2026-11-01 GRCh38 VCF，以便复现，而不是每次无记录地更换到新版本。

```bash
python datasets.py download --source clinvar --out data/clinvar_GRCh38.vcf.gz
python datasets.py clinvar --vcf data/clinvar_GRCh38.vcf.gz --missense-only --n-per-class 25 --out data/clinvar_demo.tsv
python protein_branch.py --input data/clinvar_demo.tsv --out results/clinvar_demo
```

筛选：GRCh38、双等位记录中的单碱基替换、至少一个 missense 后果、无冲突的病理/可能病理或良性/可能良性标签、多个提交者一致/专家组/指南审核状态。对完整文件做 reservoir sampling，不取前几十行以免染色体偏倚。若符合条件不足，输出真实数量。

ClinVar VCF 的 `ID` 是 Variation ID，不是 rsID；代码从 `INFO/RS` 取 rsID，缺少时使用明确坐标。过滤掉多等位 VCF 记录，避免把聚合分类错误赋给某一个 ALT。

这套面板用于开发、证据检索与展示。不能先把 ClinVar 标签输入摘要模型，再拿这些标签宣称独立预测准确率；也不能把通过这套筛选得到的样本当成人群随机样本。

### GWAS Catalog：性状/文献入口

```bash
python datasets.py download --source gwas --out data/gwas_associations.zip
python datasets.py gwas --input data/gwas_associations.zip --trait "LDL|low density lipoprotein" --n 20 --out data/gwas_ldl.tsv
```

脚本只从单 rsID 关联记录抽取输入，原始 trait、risk allele、p-value、PMID 等整行保存到 `.gwas_records.json`。不拆解单倍型或 SNP×SNP 交互关联。输出仅是 rsID 列；将 GWAS risk allele 与当前 GRCh38 REF/ALT 对齐前，不能做等位基因特异的解释。

GWAS lead SNP 不保证就是因果变异，也不保证能检测到 motif 变化。建议先做两个经典单 SNP，再扩大到几十个 GWAS SNP；不要用 motif 命中率作为疾病预测准确率。

### eQTL Catalogue 与 1000 Genomes

见 `DATASETS_AND_SOURCES.md`。前者适合提供独立的组织相关候选靶基因证据，后者适合补充群体频率。当前包不实现 eQTL/1000 Genomes 的完整数据连接，避免把基因/风险等位基因方向未经核对就混入输出。

## 6. 缓存、可复现性和故障排查

首次成功请求自动缓存到 `.cache`。随后相同参数可使用 `--offline`；缓存未命中会明确报错，不偷偷联网。`--refresh` 会重新请求更新数据，不可与 `--offline` 同时使用。

```bash
python protein_branch.py --input examples/protein_input.tsv --out results/protein_offline --offline
python jaspar_branch.py --input examples/regulatory_input.tsv --pfm-file examples/CEBPA_MA0102.5.jaspar --out results/jaspar_offline --offline
```

现场前先真实运行一次、检查 `variant_status.json` 和 `errors.json`，再验证离线重跑。空 TSV 不自动等于没有生物效应，可能是无错义后果、没有 motif 达到阈值、服务出错或无法严格映射。

代码不包含全量 ClinVar/GWAS、任何真实患者数据或 API credentials。不要提交 `.cache`、大型数据、隐私数据或密钥到公开仓库。公开项目时保留源数据库、矩阵和原始论文署名。

## 7. 给前端/摘要模块的接口

以 `variant_key = GRCh38:chrom:pos:REF>ALT` 为主连接键，保留 transcript 和来源。给摘要模型传递的内容至少分为：

- 已有数据库/文献证据（ClinVar/GWAS/eQTL，标明条件、等位基因、组织及文献）；
- 计算预测（AlphaMissense/JASPAR，标明模型/矩阵及适用范围）；
- 结构和蛋白背景（参考结构位置、功能位点、蛋白层面疾病背景）；
- 未找到/不适用/请求失败。

提示词禁止把后两类升级成已验证的变异→表型结论。代码能够嵌入 Streamlit/服务端，但当前交付不含 Web UI、自动 PubMed 检索或 LLM 生成摘要。
