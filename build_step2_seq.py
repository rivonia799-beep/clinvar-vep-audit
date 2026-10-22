# -*- coding: utf-8 -*-
"""抽样 + 从 UCSC 取参考序列，构造 VEP 数据集。

用法:
    python fetch_seq.py --n_train 4000 --n_test 2000 --flank 100 --workers 12

输出 vep_dataset.tsv:
    chrom pos ref alt label gene split seq   (seq 长度 = 2*flank+1, 变异位在正中间)
"""
import argparse, csv, json, os, random, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {"User-Agent": "Mozilla/5.0"}
TEST_CHR = {"3", "6", "9", "12", "16", "18", "19", "21"}


def sample(path, n, split, seed):
    """按自然正负比分层抽样。"""
    pos, neg = [], []
    with open(path, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            is_test = r["chrom"] in TEST_CHR
            if (split == "test") != is_test:
                continue
            (pos if r["label"] == "1" else neg).append(r)
    rng = random.Random(seed)
    rng.shuffle(pos); rng.shuffle(neg)
    p_rate = len(pos) / (len(pos) + len(neg))
    n_pos = int(round(n * p_rate)); n_neg = n - n_pos
    out = pos[:n_pos] + neg[:n_neg]
    rng.shuffle(out)
    print(f"  [{split}] 池子 pos={len(pos):,} neg={len(neg):,} 正类率={p_rate:.4f}"
          f" -> 抽 {n_pos} pos + {n_neg} neg = {len(out)}")
    return out, p_rate


def fetch(chrom, pos, flank, tries=3):
    start = pos - 1 - flank
    end = pos + flank
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38"
           f";chrom=chr{chrom};start={start};end={end}")
    for k in range(tries):
        try:
            r = OP.open(urllib.request.Request(url, headers=UA), timeout=30)
            j = json.loads(r.read().decode())
            return j.get("dna", "").upper()
        except Exception:
            time.sleep(0.5 * (k + 1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(DATA, "clinvar_all.tsv"))
    ap.add_argument("--n_train", type=int, default=4000)
    ap.add_argument("--n_test", type=int, default=2000)
    ap.add_argument("--flank", type=int, default=100)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    tr, pr_tr = sample(a.src, a.n_train, "train", a.seed)
    te, pr_te = sample(a.src, a.n_test, "test", a.seed + 1)
    jobs = [(r, "train") for r in tr] + [(r, "test") for r in te]
    print(f"共 {len(jobs)} 条，开始取序列（{a.workers} 线程）...")

    t0 = time.time()
    done = 0
    rows = []
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch, r["chrom"], int(r["pos"]), a.flank): (r, sp) for r, sp in jobs}
        for fu in futs:
            r, sp = futs[fu]
            seq = fu.result()
            done += 1
            if done % 500 == 0:
                el = time.time() - t0
                print(f"   {done}/{len(jobs)}  {el:.0f}s  ({el/done*1000:.0f} ms/条)", flush=True)
            if not seq or len(seq) != 2 * a.flank + 1:
                continue
            if seq[a.flank] != r["ref"]:
                continue          # 参考碱基对不上 -> 丢弃
            rows.append([r["chrom"], r["pos"], r["ref"], r["alt"], r["label"], r["gene"], sp, seq])

    el = time.time() - t0
    print(f"\n取到 {len(rows)}/{len(jobs)} 条可用序列，用时 {el:.0f}s")
    out = os.path.join(DATA, "vep_dataset.tsv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "label", "gene", "split", "seq"])
        w.writerows(rows)
    import collections
    c = collections.Counter((r[6], r[4]) for r in rows)
    print("写出:", out)
    for k in sorted(c):
        print("   ", k, c[k])


if __name__ == "__main__":
    main()
