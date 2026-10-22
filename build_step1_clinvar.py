# -*- coding: utf-8 -*-
"""从 ClinVar variant_summary 流式筛出 VEP 二分类数据集（不落地 455MB 原始文件）。

筛法（尽量对齐 AlphaGenome 论文 Table 4 的做法）：
  Assembly == GRCh38 ; 常染色体 1..22 ; Type == single nucleotide variant
  ReviewStatus != "no assertion criteria provided"      （即 >=1 星）
  P/LP -> label 1 ; B/LB -> label 0
  PIL 去重: (chrom, pos, ref, alt)

切分沿用论文 Table 4 的染色体划分：
  test  = 3, 6, 9, 12, 16, 18, 19, 21   ;   train = 其余
"""
import urllib.request, gzip, io, sys, os, csv, collections

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

POS_SIG = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
NEG_SIG = {"Benign", "Likely benign", "Benign/Likely benign"}
TEST_CHR = {"3", "6", "9", "12", "16", "18", "19", "21"}

op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
raw = op.open(req, timeout=600)
gz = gzip.GzipFile(fileobj=io.BufferedReader(raw, 1 << 20))

reader = csv.reader(io.TextIOWrapper(gz, encoding="utf-8", errors="ignore"), delimiter="\t")
header = next(reader)
header[0] = header[0].lstrip("#")
idx = {c: i for i, c in enumerate(header)}
need = ["Assembly", "Type", "Chromosome", "PositionVCF", "ReferenceAlleleVCF",
        "AlternateAlleleVCF", "ClinicalSignificance", "ReviewStatus", "GeneSymbol"]
for c in need:
    assert c in idx, c

seen = set()
rows = []
n = 0
AUTOS = {str(i) for i in range(1, 23)}
for row in reader:
    n += 1
    if n % 1000000 == 0:
        print(f"  scanned {n:,} rows, kept {len(rows):,}", flush=True)
    if len(row) <= idx["ReviewStatus"]:
        continue
    if row[idx["Assembly"]] != "GRCh38":
        continue
    if row[idx["Type"]] != "single nucleotide variant":
        continue
    c = row[idx["Chromosome"]]
    if c not in AUTOS:
        continue
    if row[idx["ReviewStatus"]].strip().lower().startswith("no assertion"):
        continue
    sig = row[idx["ClinicalSignificance"]].strip()
    if sig in POS_SIG:
        y = 1
    elif sig in NEG_SIG:
        y = 0
    else:
        continue
    ref = row[idx["ReferenceAlleleVCF"]].strip().upper()
    alt = row[idx["AlternateAlleleVCF"]].strip().upper()
    pos = row[idx["PositionVCF"]].strip()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        continue
    if not pos.isdigit():
        continue
    key = (c, pos, ref, alt)
    if key in seen:
        continue
    seen.add(key)
    rows.append((c, int(pos), ref, alt, y, row[idx["GeneSymbol"]], sig))

raw.close()
print(f"\nDONE scanned={n:,}  kept={len(rows):,}")

cnt = collections.Counter((("test" if r[0] in TEST_CHR else "train"), r[4]) for r in rows)
print("按 split / label 计数：")
for k in sorted(cnt):
    print("   ", k, cnt[k])
print("染色体分布(前10):", collections.Counter(r[0] for r in rows).most_common(10))

out = os.path.join(OUT, "clinvar_all.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["chrom", "pos", "ref", "alt", "label", "gene", "clinsig"])
    w.writerows(rows)
print("已写出:", out)
