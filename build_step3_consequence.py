# -*- coding: utf-8 -*-
"""Annotate the 6,000 sampled variants with their molecular consequence, and
derive two stratified subsets.

Outputs (written into data/):
    vep_dataset_annot.tsv   6,000 rows + `cons` column
    vep_notrunc.tsv         drops protein-truncating variants
    vep_noprot.tsv          keeps only non-protein, non-synonymous
                            (intronic / splice / UTR / regulatory).

Synonymous variants are labelled separately (`同义(Synonymous)`) rather than
being folded into missense, but they are NOT part of the non-protein subset used
in the reported results -- a synonymous change can still carry exonic splicing
signals. Including them is left to a future version.

Note on terminology
-------------------
This dataset contains **single-nucleotide variants (SNVs) only**. An SNV can
create a premature stop codon (nonsense / stop-gained) but **cannot by itself
cause a frameshift** -- frameshifts require an insertion or deletion. We
therefore report two separate classes and assert that the frameshift class is
empty for this data.

Usage:
    python build_step3_consequence.py
"""
import urllib.request, gzip, io, sys, os, csv, collections

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
COLS = ["chrom", "pos", "ref", "alt", "label", "gene", "split", "seq"]

NONPROT = "非蛋白质层"
NONSENSE = "终止获得(Nonsense)"
FRAMESHIFT = "移码(Frameshift)"
INFRAME = "框内插入缺失(Inframe-indel)"
MISSENSE = "错义(Missense)"
SYNONYMOUS = "同义(Synonymous)"
TRUNCATING = {NONSENSE, FRAMESHIFT}


def classify(name):
    """Classify by HGVS protein consequence. String heuristic, not a real annotator."""
    if "(p." not in name:
        return NONPROT
    p = name.split("(p.")[-1].rstrip(")")
    if "fs" in p:
        return FRAMESHIFT
    if "Ter" in p or "*" in p:
        return NONSENSE
    if "del" in p or "ins" in p or "dup" in p:
        return INFRAME
    if "=" in p:                       # synonymous, e.g. p.Arg41=
        return SYNONYMOUS
    return MISSENSE


rows = list(csv.DictReader(open(os.path.join(HERE, "vep_dataset.tsv"), encoding="utf-8"), delimiter="\t"))
K = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in rows}
print("样本数:", len(K))

op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
raw = op.open(urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"}), timeout=600)
gz = gzip.GzipFile(fileobj=io.BufferedReader(raw, 1 << 20))
rd = csv.reader(io.TextIOWrapper(gz, encoding="utf-8", errors="ignore"), delimiter="\t")
h = next(rd); h[0] = h[0].lstrip("#")
i = {c: n for n, c in enumerate(h)}

hit, n = {}, 0
for row in rd:
    n += 1
    if n % 2000000 == 0:
        print(f"  扫 {n:,} 行，命中 {len(hit)}", flush=True)
    if row[i["Assembly"]] != "GRCh38" or row[i["Type"]] != "single nucleotide variant":
        continue
    k = (row[i["Chromosome"]], row[i["PositionVCF"]],
         row[i["ReferenceAlleleVCF"]].upper(), row[i["AlternateAlleleVCF"]].upper())
    if k in K and k not in hit:
        hit[k] = classify(row[i["Name"]])
raw.close()
print(f"扫完 {n:,} 行，匹配 {len(hit)}/{len(K)}")

for r in rows:
    r["cons"] = hit.get((r["chrom"], r["pos"], r["ref"], r["alt"]), "未匹配")


def write(path, sel):
    with open(os.path.join(HERE, path), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            if sel(r):
                w.writerow({c: r[c] for c in COLS})
    print(f"  写出 {path}")


with open(os.path.join(HERE, "vep_dataset_annot.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS + ["cons"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("  写出 vep_dataset_annot.tsv")

write("vep_notrunc.tsv", lambda r: r["cons"] not in TRUNCATING)
write("vep_noprot.tsv", lambda r: r["cons"] == NONPROT)

cnt = collections.Counter(r["cons"] for r in rows)
print("\n后果分类计数:")
for k2, v2 in cnt.most_common():
    print(f"   {k2:32s} {v2}")

# --- self-check: an SNV-only set must contain no frameshift variants ---
n_fs = cnt.get(FRAMESHIFT, 0)
print(f"\n[自检] 移码变异数量 = {n_fs}（SNV 数据集应为 0）")
assert n_fs == 0, "SNV 数据集中出现移码，分类或筛选逻辑有误"

for name, sel in [("原始 6000", lambda r: True),
                  ("去掉蛋白截断", lambda r: r["cons"] not in TRUNCATING),
                  ("只留非蛋白质层", lambda r: r["cons"] == NONPROT)]:
    sub = [r for r in rows if sel(r)]
    c = collections.Counter((r["split"], r["label"]) for r in sub)
    print(f"[{name}] 训练 {c[('train','0')]+c[('train','1')]} (pos {c[('train','1')]}) | "
          f"测试 {c[('test','0')]+c[('test','1')]} (pos {c[('test','1')]}) | "
          f"测试正类率 {c[('test','1')]/max(c[('test','0')]+c[('test','1')],1):.4f}")
