# Development Log — ClinVar VEP Benchmark Audit

> Internal record of how this audit was carried out, including the wrong turns.
> Kept because the *reasoning path* is part of the result.

## 1. Motivation

While running a DNA language model on variant effect prediction (VEP), I noticed that
reported metrics on ClinVar-style evaluation sets did not match my intuition about what a
sequence-only model can and cannot do. Two specific doubts:

1. **Is the task measuring sequence modelling at all?** ClinVar labels are assigned by human
   curators using evidence that includes protein-level consequences.
2. **Is the metric comparable across subsets?** AUPRC is prevalence-dependent.

I decided to quantify both before trusting any number.

## 2. What was built

```
ClinVar variant_summary (GRCh38)
    └─ filter: autosomes, SNV, review status >= 1 star, P/LP vs B/LB
       └─ 1,423,685 variants
          └─ chromosome-wise split (test: 3,6,9,12,16,18,19,21)
             └─ sample 4,000 / 2,000, fetch ±100 bp from UCSC hg38
                └─ featurise + fit logistic regression (pure NumPy)
                   └─ evaluate: AUROC / AUPRC / random baseline
```

## 3. Timeline and wrong turns

| Step | Action | What it changed |
|---|---|---|
| 1 | Built the dataset, ran a 4-mer baseline | AUROC 0.663 / AUPRC 0.239 |
| 2 | Added two control scorers (GC content; reference-window k-mers with ALT unused) | Both ≈ 0.53 — i.e. context alone carries almost no usable signal |
| 3 | Annotated every variant with its molecular consequence | Found the composition asymmetry (below) |
| 4 | Re-ran on two stratified subsets | AUROC 0.663 → 0.678 → 0.807 |
| 5 | Swapped the model for Carbon-500M | +0.043 AUROC, **same size on both subsets** |
| 6 | Wrote the attribution | 3-factor decomposition completed |

### Wrong turns (kept for honesty)

1. **First hypothesis was wrong — and so was my first correction.** I initially assumed
   *missense* variants were the main contaminant, then "corrected" myself using a classifier that
   silently folded **synonymous** variants (HGVS `p.Arg41=`) into missense. With that fixed, the
   four classes are: non-protein 28.8%/35.6%, synonymous 0.5%/53.6%, missense 34.7%/10.8%,
   nonsense 36.0%/0.0% (positives/negatives). **Synonymous variants dominate the negatives**, and
   the protein-level contaminant is **nonsense (stop-gained)**. Lesson: a classifier bug can
   survive two rounds of reasoning and appear to "confirm" a story.
2. **First metric framing was wrong.** I claimed the benchmark results were "not comparable"
   because the metrics differed. In fact the AUPRC values *are* aligned; the correct statement
   is that they can be computed side by side but **cannot be causally attributed** (three
   variables differ simultaneously).
3. **One consequence-classification bug is known and unfixed** (see README, known limitations):
   synonymous variants written as `p.Arg41=` are counted as missense. Its direction makes the
   non-protein subset *smaller* than it should be, i.e. the bias is conservative.

## 4. Result

**Finding 1 — contamination.** 36.0% of positives are stop-gain/frameshift vs 0.0% of
negatives. A consequence-rule baseline using **no sequence features** reaches AUPRC 0.4925,
above the trained 4-mer model (0.2388).

**Finding 2 — attribution.** Task stratification: **+0.144 AUROC**. Model scale: **+0.043**,
and the same on both subsets → scaling cannot compensate for contamination.

**Finding 3 — metric validity.** AUPRC's random baseline is the positive rate; absolute AUPRC
is not cross-dataset comparable. Report AUROC alongside.

## 5. What I would do next

1. **Matched negative sampling** (distance-to-splice-site binning, same-gene matching) —
   removes the gene-identity shortcut.
2. **Real consequence annotation** (VEP / SnpEff) instead of the HGVS string heuristic.
3. **Held-out validation split** — the current numbers use fixed hyperparameters.
4. **A non-neural but domain-informed baseline** (e.g. existing splice predictors) — currently
   the strongest sequence model present is a 500M LM linear probe.

## 6. Notes on terminology

During development I used informal shorthand ("smoking-gun baseline", "dirty task").
All reporting artefacts use the standard terms instead — see `TERMINOLOGY.md`.
