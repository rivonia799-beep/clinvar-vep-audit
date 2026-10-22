# -*- coding: utf-8 -*-
"""一键重建 data/ 下的全部数据（需要联网，约 30 分钟）。

用法:
    python build.py --dry-run      # 只打印将要执行什么
    python build.py --all
    python build.py --step 2       # 只重跑第 2 步

三个 build_stepN.py 直接把结果写进 data/。
"""
import argparse, os, shutil, subprocess, sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PY = sys.executable

STEPS = {
    1: ("build_step1_clinvar.py", [], ["clinvar_all.tsv"], "流式筛 ClinVar（下 455 MB，约 8 分钟）"),
    2: ("build_step2_seq.py", ["--n_train", "4000", "--n_test", "2000", "--workers", "14"],
        ["vep_dataset.tsv"], "抽样 + 从 UCSC 取序列（约 17 分钟）"),
    3: ("build_step3_consequence.py", [],
        ["vep_dataset_annot.tsv", "vep_notrunc.tsv", "vep_noprot.tsv"], "标注分子后果 + 生成干净子集（下 455 MB，约 8 分钟）"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--step", type=int, choices=[1, 2, 3])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    steps = [1, 2, 3] if (a.all or a.dry_run) else ([a.step] if a.step else [])
    if not steps:
        ap.print_help(); return

    os.makedirs(DATA, exist_ok=True)
    for n in steps:
        script, args, outs, desc = STEPS[n]
        print(f"\n=== 第 {n} 步 · {desc} ===")
        print(f"    {script} {' '.join(args)}")
        if a.dry_run:
            print(f"    -> 产出 {outs}，然后搬到 data/"); continue
        r = subprocess.run([PY, os.path.join(HERE, script)] + args, cwd=HERE)
        if r.returncode != 0:
            print(f"!!! 第 {n} 步失败（退出码 {r.returncode}）"); raise SystemExit(1)
        for f in outs:
            dst = os.path.join(DATA, f)
            print(f"    data/{f} " + ("已更新" if os.path.exists(dst) else "!!! 未生成"))
    print("\n完成。")


if __name__ == "__main__":
    main()
