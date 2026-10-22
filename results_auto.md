# ClinVar-VEP 干净基准 v0.1 · 结果

> 生成方式：`python evaluate.py --all`

| Dataset | Scorer | Test pos. rate | Random AUPRC | **AUROC** | **AUPRC** | AUPRC/rand | s |
|---|---|---|---|---|---|---|---|
| vep_dataset.tsv | random | 0.1110 | 0.1110 | 0.4891 | 0.1116 | 1.01x | 0.0 |
| vep_dataset.tsv | gc | 0.1110 | 0.1110 | 0.5349 | 0.1221 | 1.10x | 0.1 |
| vep_dataset.tsv | kmer_ref | 0.1110 | 0.1110 | 0.5372 | 0.1265 | 1.14x | 1.2 |
| vep_dataset.tsv | kmer_diff | 0.1110 | 0.1110 | 0.6633 | 0.2388 | 2.15x | 4.1 |
| vep_dataset.tsv | kmer_diff+shuffle | 0.1110 | 0.1110 | 0.5107 | 0.1261 | 1.14x | 3.9 |
| vep_dataset_annot.tsv | consequence_rule | 0.1110 | 0.1110 | 0.6802 | 0.4925 | 4.44x | 0.0 |
| vep_notrunc.tsv | kmer_diff | 0.0740 | 0.0740 | 0.6781 | 0.2313 | 3.13x | 3.8 |
| vep_noprot.tsv | kmer_diff | 0.0918 | 0.0918 | 0.8070 | 0.3521 | 3.83x | 1.3 |

## 读法

- **AUROC** 与正类比例无关，跨数据集可比；**AUPRC 绝对值不可**（它的随机基线 = 正类比例）
- **AUPRC/rand** = AUPRC ÷ random baseline; safer than the absolute value
- `kmer_diff+shuffle` 是阴性对照，应回到 AUROC ≈ 0.5
- `consequence_rule` (a heuristic using only the HGVS consequence label, no sequence) is a **contamination indicator**: the higher it scores on a dataset, the more protein-level consequence that dataset contains
