# -*- coding: utf-8 -*-
"""Carbon-500M linear probe -- swap the model, keep task/data/split/metric fixed.

目的（对应开发记录里的"三因素分解"）:
    只动【模型能力】，任务 / 数据 / 切分 / 指标全不变。
    对照: 4-mer diff 在 data/vep_noprot.tsv 上 = AUROC 0.807 / AUPRC 0.352

用法:
    python carbon_probe.py --max 40            # 先用 40 条冒烟（约 1 分钟）
    python carbon_probe.py                     # 全量（默认 data/vep_noprot.tsv，约 4198 次前向）
    python carbon_probe.py --data data/vep_dataset.tsv     # 原始池子（较慢）

注意: 第一次跑要加载 500M 权重，约 1 分钟；之后会缓存嵌入。
"""
import os, sys, time, argparse, csv
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

WORK = r"C:\Users\12495\Documents\Codex\2026-09-20\fen\work"
os.environ.setdefault("HF_MODULES_CACHE", os.path.join(WORK, "hf_modules"))
CARBON = os.path.join(WORK, "carbon-500m")
QWEN_TOK = os.path.join(WORK, "qwen3-tokenizer")

HERE = os.path.dirname(os.path.abspath(__file__))
FLANK = 100


def load_carbon():
    import torch
    from transformers import AutoTokenizer, AutoModelForCausalLM
    _orig = AutoTokenizer.from_pretrained

    def _patched(name, *a, **kw):
        if name == "Qwen/Qwen3-4B-Base":
            name = QWEN_TOK
        return _orig(name, *a, **kw)

    AutoTokenizer.from_pretrained = _patched
    tok = _orig(CARBON, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(CARBON, dtype=torch.float32)
    model.eval()
    return tok, model


def embed(tok, model, seqs, bs=8, tag=""):
    """最后一层、去掉首尾特殊 token 后的 mean pooling。返回 (n, H)。"""
    import torch
    out = []
    t0 = time.time()
    with torch.no_grad():
        for i in range(0, len(seqs), bs):
            chunk = seqs[i:i + bs]
            enc = tok(["<dna>" + s for s in chunk], return_tensors="pt",
                      padding=True, add_special_tokens=False)
            h = model(input_ids=enc["input_ids"], output_hidden_states=True).hidden_states[-1]
            h = h.numpy()
            for j in range(len(chunk)):
                n = int(enc["attention_mask"][j].sum())
                out.append(h[j, 1:max(n - 1, 2)].mean(0))
            if (i // bs) % 20 == 0:
                print(f"    {tag} {i+len(chunk)}/{len(seqs)}  {time.time()-t0:.0f}s", flush=True)
    return np.array(out, dtype=np.float32)


def fit_lr(X, y, lr=0.5, iters=2000, l2=1e-4):
    mu, sd = X.mean(0), X.std(0) + 1e-9
    A = (X - mu) / sd
    n, d = A.shape
    w = np.zeros(d); b = 0.0
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(A @ w + b, -30, 30)))
        g = p - y
        w -= lr * (A.T @ g / n + l2 * w)
        b -= lr * g.mean()
    return w, b, mu, sd


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    o = np.argsort(s, kind="mergesort"); sp = s[o]
    r = np.empty(len(s)); i = 0; base = 1.0
    while i < len(sp):
        j = i
        while j + 1 < len(sp) and sp[j + 1] == sp[i]:
            j += 1
        r[o[i:j + 1]] = (base + (base + (j - i))) / 2.0
        base += (j - i + 1); i = j + 1
    n1 = int((y == 1).sum()); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def ap(y, s):
    o = np.argsort(-np.asarray(s), kind="mergesort")
    yy = np.asarray(y)[o]
    tp = np.cumsum(yy); fp = np.cumsum(1 - yy)
    prec = tp / np.maximum(tp + fp, 1e-12); rec = tp / tp[-1]
    return float(np.sum(np.diff(np.concatenate([[0.0], rec])) * prec))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--data", default="data/vep_noprot.tsv")
    ap_.add_argument("--mode", default="diff", choices=["diff", "ref", "both"])
    ap_.add_argument("--bs", type=int, default=8)
    ap_.add_argument("--max", type=int, default=0, help="只取前 N 条（冒烟用）")
    ap_.add_argument("--cache", default=None)
    a = ap_.parse_args()

    path = a.data if os.path.isabs(a.data) else os.path.join(HERE, a.data)
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    if a.max:
        # 文件里 train 块在前、test 块在后 -> 前后各取一半，保证两个 split 都有
        half = a.max // 2
        rows = rows[:half] + rows[-(a.max - half):]
    tag = f"carbon_{os.path.basename(a.data)}"
    cache = a.cache or os.path.join(HERE, f"_cache_{tag}.npz")

    refs = [r["seq"] for r in rows]
    alts = [r["seq"][:FLANK] + r["alt"] + r["seq"][FLANK + 1:] for r in rows]

    Mr = Ma = None
    if os.path.exists(cache):
        # ★ 用 with 读，读完立刻关文件句柄（Windows 上不关就删不掉/盖不了）
        with np.load(cache) as z:
            Mr, Ma = z["Mr"], z["Ma"]
        if Mr.shape[0] != len(rows):     # 缓存条数和当前数据对不上 -> 忽略，末尾会覆盖写
            print(f"[!] 缓存 {Mr.shape[0]} 条 != 数据 {len(rows)} 条，忽略缓存，重新提取")
            Mr = Ma = None
        else:
            print(f"用缓存 {os.path.basename(cache)}  形状 {Mr.shape}")
    if Mr is None:
        print(f"数据 {os.path.basename(a.data)}  {len(rows)} 条  -> 需要 {2*len(rows)} 次前向")
        t0 = time.time()
        tok, model = load_carbon()
        print(f"  模型加载完成 {time.time()-t0:.0f}s")
        Mr = embed(tok, model, refs, bs=a.bs, tag="REF")
        Ma = embed(tok, model, alts, bs=a.bs, tag="ALT")
        np.savez_compressed(cache, Mr=Mr, Ma=Ma)
        print(f"  嵌入已缓存 -> {os.path.basename(cache)}")

    split = np.array([r["split"] for r in rows])
    y = np.array([int(r["label"]) for r in rows], dtype=np.float64)
    if a.mode == "diff":
        X = Ma - Mr
    elif a.mode == "ref":
        X = Mr
    else:
        X = np.concatenate([Mr, Ma - Mr], axis=1)

    tr, te = split == "train", split == "test"
    w, b, mu, sd = fit_lr(X[tr], y[tr])
    s = ((X[te] - mu) / sd) @ w + b
    yte = y[te]
    prev = float(yte.mean())
    A, P = auc(yte, s), ap(yte, s)

    print("=" * 84)
    print(f"Carbon-500M 探针   mode={a.mode}   训练 {int(tr.sum())}  测试 {int(te.sum())} (pos {int(yte.sum())})")
    print(f"随机基线 (random)   AUROC 0.5000   AUPRC {prev:.4f}")
    print(f"  Carbon-500M probe   AUROC {A:.4f}   AUPRC {P:.4f}   AUPRC/rand {P/prev:.2f}x")
    print(f"  对照: 4-mer logistic regression -- see results.md (same dataset)")
    print("=" * 84)


if __name__ == "__main__":
    main()
