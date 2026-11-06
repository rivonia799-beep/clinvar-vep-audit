# AI4Science / DNA Sequence Modeling Portfolio Notes

This note is a short guide for AI4Science, DNA sequence modeling, and algorithm engineering
roles. It summarizes what the repository demonstrates and, just as importantly, what it does
not claim.

## 30-second summary

I built a reproducible audit of a fixed-window DNA variant-effect benchmark. The pipeline
starts from 9.3M ClinVar rows, filters 1.42M high-confidence autosomal SNVs, constructs 6,000
reference windows of 201 bp, performs a chromosome-wise train/test split, and compares:

- GC content and 4-mer frequency baselines
- 4-mer difference features with a from-scratch NumPy logistic regression
- a consequence-only rule baseline
- a frozen Carbon-500M DNA language model with a linear probe

The main finding is that the unstratified test set contains a strong protein-level shortcut:
36.0% of positives are nonsense/stop-gained variants, while negatives contain 0.0%. Since the
model only sees a fixed window of naked genomic DNA, the premature-stop consequence is not
recoverable from the input. A zero-parameter consequence rule reaches AUPRC 0.4925, higher than
the 4-mer model (0.2388) and the Carbon-500M probe on the same unstratified task (0.4701).

## What this repository demonstrates

### 1. DNA sequence representation

- fixed-window DNA inputs centered on the variant (`±100 bp`, 201 bp total)
- k-mer frequency features and ALT−REF difference features
- frozen pretrained DNA language model embeddings using mean pooling
- direct comparison between classical features and a 500M-parameter DNA model

Relevant files:

- `build_step2_seq.py`
- `evaluate.py`
- `carbon_probe.py`

### 2. Model evaluation and benchmark design

- chromosome-wise hard split to reduce sequence homology leakage
- AUROC and AUPRC reported together
- random baseline and shuffled-label negative control
- consequence-only rule baseline
- explicit analysis of dataset shift and class imbalance

Relevant files:

- `evaluate.py`
- `results.md`
- `REPORT_ONE_PAGER.md`

### 3. Algorithm engineering

- pure NumPy implementation of logistic regression, sigmoid, AUROC, and average precision
- PyTorch/Transformers inference path for Carbon-500M embeddings
- frozen data hashes and rebuild instructions
- Dockerfile with a reproducible default command
- CLI paths for individual scorers and full baseline matrix

Relevant files:

- `Dockerfile`
- `requirements.txt`
- `build.py`
- `evaluate.py`
- `carbon_probe.py`

### 4. Scientific diagnosis rather than leaderboard chasing

The project does not ask only whether a larger model wins. It asks whether the benchmark is
measuring the intended capability. The metric attribution separates two effects:

- task stratification: AUROC 0.663 → 0.807 (+0.144)
- model scale: 4-mer LR → Carbon-500M probe on the stratified task: 0.807 → 0.850 (+0.043)

This is the part most directly transferable to AI-driven sequence design: define the target
function, control the data distribution, and verify that evaluation actually measures the
intended biological capability.

## Relevance to an mRNA / sequence-design team

| Job requirement | Evidence in this repository | Level |
|---|---|---|
| DNA sequence representation | 201 bp windows, k-mer and frozen DNA-LM embeddings | Strong |
| Model evaluation | chromosome split, AUROC/AUPRC, shuffled-label control | Strong |
| Data QC and pipeline engineering | 9.3M-row filtering, 1.42M final variants, frozen hashes | Strong |
| PyTorch / pretrained model use | Carbon-500M hidden-state extraction and linear probe | Working knowledge |
| Algorithm implementation | from-scratch NumPy logistic regression and metrics | Strong |
| Reproducibility | CLI scripts, Dockerfile, rebuild instructions | Strong |
| mRNA generation / UTR design | not implemented here | Gap |
| VAE / Diffusion / autoregressive sequence generation | not implemented here | Gap |
| Multi-objective sequence optimization | not implemented here | Gap |

## Honest scope statement

This repository is a benchmark audit and DNA sequence-model probing project. It is not a
generative mRNA design system and does not claim experience with 5′UTR/CDS/3′UTR optimization,
VAE, diffusion, or autoregressive sequence generation.

The relevant transferable skills are sequence representation, controlled evaluation, data
pipeline construction, and reproducible engineering. Those are the parts I would bring into a
sequence-design team first.

