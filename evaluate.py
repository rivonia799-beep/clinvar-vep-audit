# -*- coding: utf-8 -*-
"""ClinVar VEP Benchmark Audit -- evaluation script (pure numpy, no sklearn).

用法:
    python evaluate.py --scorer kmer_diff
    python evaluate.py --scorer dummy_trunc --data data/vep_dataset_annot.tsv
    python evaluate.py --all              # 跑完整矩阵并写出 results.md

设计原则（见 README §2）:
    * 强制同时报 AUROC / AUPRC / 倍数   —— AUPRC 绝对值会随正类比例变
    * 内置"照妖镜"基线 dummy_trunc      —— 检验基准本身是否干净
"""
import argparse, csv, os, sys, time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
B = "ACGT"
KMERS = [a + b + c + d for a in B for b in B for c in B for d in B]
KIDX = {k: i for i, k in enumerate(KMERS)}
FLANK = 100
TRUNCATING = {"终止获得(Nonsense)", "移码(Frameshift)"}   # 蛋白截短型

SCORERS = ["random", "gc", "kmer_ref", "kmer_diff", "consequence_rule"]


# ---------------- 特征 ----------------
def kmer_freq(seq, k=4):
    v = np.zeros(4 ** k, dtype=np.float64)
    n = len(seq) - k + 1
    for i in range(n):
        j = KIDX.get(seq[i:i + k])
        if j is not None:
            v[j] += 1.0
    return v / max(n, 1)


def featurize(rows, mode, seed=0):
    """返回 (X, fit_needed)"""
    if mode == "random":
        rng = np.random.default_rng(seed)
        return rng.random((len(rows), 1)), False
    if mode == "consequence_rule":
        X = np.array([[1.0 if r.get("cons") in TRUNCATING else 0.0] for r in rows])
        return X, False
    if mode == "gc":
        return np.array([[(s.count("G") + s.count("C")) / len(s)] for s in (r["seq"] for r in rows)]), True
    F = []
    for r in rows:
        ref = r["seq"]
        alt = ref[:FLANK] + r["alt"] + ref[FLANK + 1:]
        fr, fa = kmer_freq(ref), kmer_freq(alt)
        if mode == "kmer_ref":
            F.append(fr)
        elif mode == "kmer_diff":
            F.append(fa - fr)
    return np.vstack(F), True


# ---------------- 指标 ----------------
def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    order = np.argsort(s, kind="mergesort")
    ss = s[order]
    ranks = np.empty(len(s))
    i, r = 0, 1.0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        ranks[order[i:j + 1]] = (r + (r + (j - i))) / 2.0
        r += (j - i + 1); i = j + 1
    npos = int((y == 1).sum()); nneg = len(y) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    return (ranks[y == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)


def average_precision(y, s):
    order = np.argsort(-np.asarray(s), kind="mergesort")
    y = np.asarray(y)[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    prec = tp / np.maximum(tp + fp, 1e-12)
    rec = tp / tp[-1]
    return float(np.sum(np.diff(np.concatenate([[0.0], rec])) * prec))


def fit_lr(X, y, lr=0.5, iters=2000, l2=1e-4):
    n, d = X.shape
    w = np.zeros(d); b = 0.0
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ w + b, -30, 30)))
        g = p - y
        w -= lr * (X.T @ g / n + l2 * w)
        b -= lr * g.mean()
    return w, b


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------- 单次评测 ----------------
def run_one(data_path, scorer, shuffle=False, seed=0, quiet=False):
    rows = read_tsv(data_path)
    tr = [r for r in rows if r["split"] == "train"]
    te = [r for r in rows if r["split"] == "test"]
    ytr = np.array([int(r["label"]) for r in tr], dtype=np.float64)
    yte = np.array([int(r["label"]) for r in te], dtype=np.float64)
    if shuffle:
        ytr = np.random.default_rng(0).permutation(ytr)

    Xtr, need_fit = featurize(tr, scorer, seed)
    Xte, _ = featurize(te, scorer, seed)
    if need_fit:
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
        Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
        t0 = time.time()
        w, b = fit_lr(Xtr, ytr)
        cost = time.time() - t0
        s = Xte @ w + b
    else:
        cost = 0.0
        s = Xte[:, 0]

    prev = float(yte.mean())
    a, p = auc(yte, s), average_precision(yte, s)
    return dict(data=os.path.basename(data_path), scorer=scorer + ("+shuffle" if shuffle else ""),
                n_tr=len(tr), n_te=len(te), pos_tr=int(ytr.sum()), pos_te=int(yte.sum()),
                prev=prev, auroc=a, auprc=p, ratio=p / prev if prev else float("nan"), sec=cost)


def show(r):
    print(f"  {r['scorer']:16s} | 训练 {r['n_tr']:5d} 测试 {r['n_te']:5d} (pos {r['pos_te']:4d}) | "
          f"正类率 {r['prev']:.4f} | AUROC {r['auroc']:.4f} | AUPRC {r['auprc']:.4f} | "
          f"AUPRC/rand {r['ratio']:5.2f}x | {r['sec']:.1f}s")


MATRIX = [
    ("data/vep_dataset.tsv", "random", False),
    ("data/vep_dataset.tsv", "gc", False),
    ("data/vep_dataset.tsv", "kmer_ref", False),
    ("data/vep_dataset.tsv", "kmer_diff", False),
    ("data/vep_dataset.tsv", "kmer_diff", True),
    ("data/vep_dataset_annot.tsv", "consequence_rule", False),
    ("data/vep_notrunc.tsv", "kmer_diff", False),
    ("data/vep_noprot.tsv", "kmer_diff", False),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/vep_dataset.tsv")
    ap.add_argument("--scorer", default="kmer_diff", choices=SCORERS)
    ap.add_argument("--shuffle", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    if not a.all:
        p = a.data if os.path.isabs(a.data) else os.path.join(HERE, a.data)
        r = run_one(p, a.scorer, a.shuffle)
        print("=" * 100)
        print(f"数据集 {r['data']}")
        print(f"随机基线 (random)   AUROC 0.5000   AUPRC {r['prev']:.4f}")
        print("-" * 100)
        show(r)
        print("=" * 100)
        return

    print("=" * 100)
    print("ClinVar VEP Benchmark Audit -- full matrix")
    print("=" * 100)
    res, cur = [], None
    for data, scorer, sh in MATRIX:
        if data != cur:
            cur = data
            print(f"\n[{data}]")
        p = os.path.join(HERE, data)
        try:
            r = run_one(p, scorer, sh)
        except Exception as e:
            print(f"  {scorer:16s} | 跳过：{type(e).__name__} {e}")
            continue
        show(r); res.append(r)

    out = os.path.join(HERE, "results_auto.md")   # 自动生成的矩阵；手写分析见 results.md
    with open(out, "w", encoding="utf-8") as f:
        f.write("# ClinVar-VEP 干净基准 v0.1 · 结果\n\n")
        f.write("> 生成方式：`python evaluate.py --all`\n\n")
        f.write("| Dataset | Scorer | Test pos. rate | Random AUPRC | **AUROC** | **AUPRC** | AUPRC/rand | s |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in res:
            f.write(f"| {r['data']} | {r['scorer']} | {r['prev']:.4f} | {r['prev']:.4f} | "
                    f"{r['auroc']:.4f} | {r['auprc']:.4f} | {r['ratio']:.2f}x | {r['sec']:.1f} |\n")
        f.write("\n## 读法\n\n")
        f.write("- **AUROC** 与正类比例无关，跨数据集可比；**AUPRC 绝对值不可**（它的随机基线 = 正类比例）\n")
        f.write("- **AUPRC/rand** = AUPRC ÷ random baseline; safer than the absolute value\n")
        f.write("- `kmer_diff+shuffle` 是阴性对照，应回到 AUROC ≈ 0.5\n")
        f.write("- `consequence_rule` (a heuristic using only the HGVS consequence label, no sequence) is a **contamination indicator**: the higher it scores on a dataset, the more protein-level consequence that dataset contains\n")
    print(f"\n已写出 {out}")


if __name__ == "__main__":
    main()
