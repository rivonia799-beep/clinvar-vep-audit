# Results — ClinVar VEP Benchmark Audit

- Raw machine-generated matrix: [`results_auto.md`](results_auto.md) (`python evaluate.py --all`)
- This file: curated tables and reading notes.

## 1. Main matrix (test split)

| Dataset | Scorer | Pos. rate | Random AUPRC | **AUROC** | **AUPRC** | AUPRC / random |
|---|---|---|---|---|---|---|
| `vep_dataset.tsv` (unstratified) | random | 0.1110 | 0.1110 | 0.4891 | 0.1116 | 1.01× |
| `vep_dataset.tsv` | GC content | 0.1110 | 0.1110 | 0.5349 | 0.1221 | 1.10× |
| `vep_dataset.tsv` | 4-mer context (ALT unused) | 0.1110 | 0.1110 | 0.5372 | 0.1265 | 1.14× |
| `vep_dataset.tsv` | **4-mer difference** | 0.1110 | 0.1110 | **0.6633** | **0.2388** | **2.15×** |
| `vep_dataset.tsv` | 4-mer difference + shuffled labels | 0.1110 | 0.1110 | 0.5107 | 0.1261 | 1.14× |
| `vep_dataset_annot.tsv` | **consequence-rule baseline** | 0.1110 | 0.1110 | **0.6802** | **0.4925** | **4.44×** |
| `vep_notrunc.tsv` | 4-mer difference | 0.0740 | 0.0740 | 0.6781 | 0.2313 | 3.13× |
| `vep_noprot.tsv` | 4-mer difference | 0.0918 | 0.0918 | 0.8070 | 0.3521 | 3.83× |
| `vep_noprot.tsv` | **Carbon-500M probe** | 0.0918 | 0.0918 | **0.8502** | **0.5852** | **6.37×** |
| `vep_dataset.tsv` (unstratified) | Carbon-500M probe | 0.1110 | 0.1110 | 0.7110 | 0.4701 | 4.23× |

## 2. How to read this table

- **AUROC is comparable across datasets** — its random baseline is always 0.500.
- **AUPRC absolute values are not** — its random baseline *is* the positive rate.
  Compare **AUPRC / random** instead, or report both.
- Two built-in sanity checks:
  - `random` → AUPRC 0.1116 ≈ positive rate 0.1110 → **metric implementation is correct**
  - `4-mer difference + shuffled labels` → AUROC 0.5107 → **no leakage in the pipeline**

## 3. Molecular-consequence composition of the test split

| Consequence class | Positives (n = 222) | Negatives (n = 1,778) |
|---|---|---|
| Non-protein (intronic / splice / UTR) | 64 (28.8%) | 633 (35.6%) |
| Synonymous | 1 (0.5%) | 953 (53.6%) |
| Missense | 77 (34.7%) | 192 (10.8%) |
| **Nonsense (stop-gained)** | **80 (36.0%)** | **0 (0.0%)** |

This is the basis of **finding 1**. The 0.0% on the negative side is not a coincidence —
it follows from the label definition: in ClinVar a protein-truncating variant is essentially
never classified as benign.

**Note:** the dataset is SNV-only, so *no frameshift variants exist in it at all* (verified:
0 out of 276 protein-truncating variants carry a `fs` annotation). The correct description of
this class is **nonsense / stop-gained**.

Two further observations:

- **Synonymous variants dominate the negative class** (53.6% of negatives vs 0.5% of positives).
  They are excluded from the non-protein subset used in the reported results, because a synonymous
  change can still carry exonic splicing signals.
- **Missense is enriched in positives** (34.7% vs 10.8%). An earlier version of this analysis folded
  synonymous variants into the missense class, which inverted this ratio.

## 4. Effect decomposition

| Factor | How it was isolated | AUROC change | Magnitude |
|---|---|---|---|
| Task stratification | same model (4-mer LR), different subsets | 0.663 → 0.807 | **+0.144** |
| Model scale | same subset, 4-mer LR → Carbon-500M probe | 0.807 → 0.850 | +0.043 |

Scaling gain measured on both subsets:

| Subset | 4-mer LR | Carbon-500M | Δ |
|---|---|---|---|
| Unstratified | 0.663 | 0.711 | **+0.048** |
| Non-protein | 0.807 | 0.850 | **+0.043** |

The two gains are nearly identical → model scale behaves like an **additive** improvement that
is independent of task contamination. **Contamination cannot be compensated by scaling.**

## 5. Reference point from the literature (not directly comparable)

Avsec et al., *Nature* **649**, 1206–1218 (2026) reports, for the `clinvar_splice_site_region`
task with metric `auprc_max_abs_track_aggregation`:

| | AUPRC |
|---|---|
| Random baseline | 0.115 |
| Pangolin | 0.554 |
| AlphaGenome | 0.570 |

⚠️ **Not comparable** with the values above: the task subset, negative-sampling procedure and
evaluated models all differ. Task subset alone changes difficulty substantially — our
non-protein subset is an easier task by construction.

Note in particular that the Carbon-500M probe reaches AUPRC 0.5852, numerically above the
reported 0.570. **This is not evidence that the probe is better**; it reflects a different and
easier task. Do not present it as a comparison.

## 6. Terminology

| Development term | Reporting term used here |
|---|---|
| "smoking-gun baseline" | consequence-rule baseline |
| "dirty task" | unstratified set |
| "cleaned task" | molecule-consequence stratified subset |
| "dummy rule" | rule-based baseline |
