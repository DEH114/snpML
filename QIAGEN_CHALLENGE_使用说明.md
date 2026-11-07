# QIAGEN BioHack 2026 — 20 个 rsID 的使用说明

这版保留出题方截图中的 **全部 20 个 rsID 和原始顺序**，不提前猜测 REF/ALT 或疾病结论。支持人类 **GRCh38**。

## 最快启动

在 `biohack_snp_toolkit` 文件夹内运行：

```bash
python -m pip install -r requirements.txt

# 第一步：只解析 rsID → GRCh38 REF/ALT，发现多等位基因及 indel
python challenge_batch.py --out results/qiagen --prepare-only --annotate-other

# 第二步：批量跑 VEP/AlphaMissense + JASPAR 简单示例
python challenge_batch.py --out results/qiagen --vep-only --annotate-other

# 第三步：运行蛋白注释/参考结构 + 全量 JASPAR（先下载 JASPAR 文件）
python challenge_batch.py --out results/qiagen_full --pfm-file data/vertebrates.jaspar --annotate-other
```

注意：第二步默认只使用压缩包附带的 **CEBPA 一个 JASPAR motif**，为了确保随时能做最小演示；不能因此判断其他转录因子的结合不受影响。需要到 <https://jaspar.elixir.no/downloads/> 下载 `CORE / Vertebrates / non-redundant / Single batch file / JASPAR` 才能扩大筛选范围。

如果时间紧，可以先针对三个 rsID 测试：

```bash
python challenge_batch.py --only-rsids rs334,rs429358,rs7412 --vep-only --out results/quick_demo
```

## 对前端最重要的文件

| 路径 | 内容 |
| --- | --- |
| `results/qiagen/challenge_summary.tsv` | **每个原始 rsID 一行，永远按原顺序**：是否解析、支路状态、命中数、错误原因。 |
| `results/qiagen/resolved_snvs.tsv` | **每个可支持的具体 ALT 等位基因一行**：可以直接分别传给 Protein / JASPAR。 |
| `results/qiagen/unsupported_or_unresolved.tsv` | insertion、deletion、无法确认的旧 rsID、多重定位等，不会静默丢弃。 |
| `results/qiagen/protein/protein_results.tsv` | missense 注释、AlphaMissense、UniProt 和可用 AlphaFold 参考结构信息。 |
| `results/qiagen/jaspar/motif_results.tsv` | 参考和替代序列的同位点 motif 匹配差异。 |
| `results/qiagen/other_variants_vep/` | 启用 `--annotate-other` 后非 SNV/无法解析记录的 rsID 级 VEP 后备注释，不是特定 ALT 的病理预测。 |

## 这 20 个 rsID 的特殊注意事项

- **rs121913529（KRAS）**：一个 rsID 可以有多个单核苷酸 ALT。代码会分别分析，不随机只选择一种替换。
- **rs113993960（CFTR）**：经典的 p.Phe508del 属于删除变异，不应送入仅适用于单碱基替换的 AlphaMissense 或单碱基 JASPAR 分支。
- **rs80357713**：公开资料提示旧 ID 合并到 rs80357783，且涉及插入变异。代码保留原 ID，记录解析/不支持状态，不自动认定两者是同一个具体等位基因。
- **rs429358 和 rs7412（APOE）**：必须结合两个位点和个体基因型/相位信息才能判断常见 APOE 单倍型。本测试表只有 rsID，不能推断个人 APOE 基因型。

## 校验与限制

- 本地使用**合成 API 响应**和临时文件运行了 **34 项单元与流程测试，全部通过**；并核对原始 20 个 rsID 的顺序。
- 打包环境不能进行真实外部 API 全流程验证；因此**没有声称已成功运行 20 个真实 SNP 的注释**。
- 结果不是临床结论。AlphaMissense 是计算预测；AlphaFold pLDDT 是参考模型置信度；JASPAR 命中是候选调控机制。具体疾病关系应整合 ClinVar、GWAS 和原始文献，并保留来源。

**英文版（适合直接发给队友）**：`QIAGEN_CHALLENGE_QUICKSTART_EN.md`。
