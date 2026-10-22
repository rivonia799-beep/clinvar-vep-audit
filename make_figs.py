# -*- coding: utf-8 -*-
"""Generate report figures as SVG -- no third-party dependencies (stdlib only).

Usage:
    python make_figs.py

Outputs:
    fig_A_auroc.svg
    fig_B_consequence_rule.svg

SVG can be embedded directly in Markdown (GitHub renders it) and inserted into
Word (Insert > Pictures). To get PNG, open the SVG in a browser and screenshot,
or use any online SVG->PNG converter.
"""
import os

W, H = 900, 560
L, R, T, B = 120, 850, 80, 470          # plot area
FS = 'font-family="Helvetica, Arial, Segoe UI, sans-serif"'


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def head():
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
            f'viewBox="0 0 {W} {H}">\n<rect width="{W}" height="{H}" fill="#ffffff"/>\n')


def ymap(v, ymin, ymax):
    return B - (v - ymin) / (ymax - ymin) * (B - T)


def yaxis(ymin, ymax, step, ylabel, ticks=6):
    """Returns (svg_string, list_of_value_boundaries)"""
    out = []
    # axis lines
    out.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{B}" stroke="#333" stroke-width="1.5"/>')
    out.append(f'<line x1="{L}" y1="{B}" x2="{R}" y2="{B}" stroke="#333" stroke-width="1.5"/>')
    n = int(round((ymax - ymin) / step))
    for i in range(n + 1):
        v = ymin + i * step
        y = ymap(v, ymin, ymax)
        out.append(f'<line x1="{L-6}" y1="{y:.1f}" x2="{L}" y2="{y:.1f}" stroke="#333"/>')
        out.append(f'<text x="{L-12}" y="{y+4:.1f}" font-size="12" {FS} fill="#333" '
                   f'text-anchor="end">{v:.2f}</text>')
    out.append(f'<text x="30" y="{(T+B)/2:.0f}" font-size="13" {FS} fill="#333" '
               f'text-anchor="middle" transform="rotate(-90 30 {(T+B)/2:.0f})">{esc(ylabel)}</text>')
    return "\n".join(out)


def dashed(v, ymin, ymax, label):
    y = ymap(v, ymin, ymax)
    return (f'<line x1="{L}" y1="{y:.1f}" x2="{R}" y2="{y:.1f}" stroke="#888" '
            f'stroke-width="1.2" stroke-dasharray="6,4"/>\n'
            f'<text x="{R}" y="{T-16}" font-size="11" {FS} fill="#666" '
            f'text-anchor="end">{esc(label)}</text>')


def fig_a(path):
    ymin, ymax, step = 0.4, 0.95, 0.1
    groups = ["Unstratified\n(n = 2,000)", "Non-protein subset\n(n = 697)"]
    series = [("4-mer logistic regression", [0.663, 0.807], "#8ab4f8"),
              ("Carbon-500M linear probe", [0.711, 0.850], "#1a73e8")]

    s = [head()]
    s.append(f'<text x="{(L+R)/2:.0f}" y="34" font-size="16" {FS} fill="#111" '
             f'text-anchor="middle">Effect of task stratification vs. model scale</text>')
    s.append(yaxis(ymin, ymax, step, "AUROC (test split)"))
    s.append(dashed(0.5, ymin, ymax, "random baseline (0.500)"))

    gw = (R - L) / len(groups)
    bw = gw * 0.30
    for gi, gname in enumerate(groups):
        cx = L + gw * (gi + 0.5)
        for si, (sname, vals, color) in enumerate(series):
            v = vals[gi]
            off = (si - 0.5) * (bw + 6)
            x = cx + off - bw / 2
            y = ymap(v, ymin, ymax)
            s.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" '
                     f'height="{(B-y):.1f}" fill="{color}" stroke="#ffffff"/>')
            s.append(f'<text x="{x+bw/2:.1f}" y="{y-7:.1f}" font-size="12" {FS} '
                     f'fill="#111" text-anchor="middle" font-weight="bold">{v:.3f}</text>')
        for li, line in enumerate(gname.split("\n")):
            s.append(f'<text x="{cx:.1f}" y="{B+22+li*16:.0f}" font-size="12" {FS} '
                     f'fill="#333" text-anchor="middle">{esc(line)}</text>')

    # legend
    lx, ly = L + 8, T + 10
    for i, (sname, _, color) in enumerate(series):
        y = ly + i * 22
        s.append(f'<rect x="{lx}" y="{y-10}" width="14" height="12" fill="{color}"/>')
        s.append(f'<text x="{lx+20}" y="{y}" font-size="12" {FS} fill="#333">{esc(sname)}</text>')
    s.append("</svg>\n")
    open(path, "w", encoding="utf-8").write("\n".join(s))
    print("wrote", path)


def fig_b(path):
    ymin, ymax, step = 0.0, 0.6, 0.1
    bars = [("Consequence-rule\nbaseline\n(0 params)", 0.4925, "#d93025"),
            ("Carbon-500M\nlinear probe\n(500M params)", 0.4701, "#1a73e8"),
            ("4-mer logistic\nregression\n(256-d)", 0.2388, "#8ab4f8")]

    s = [head()]
    s.append(f'<text x="{(L+R)/2:.0f}" y="34" font-size="16" {FS} fill="#111" '
             f'text-anchor="middle">A consequence rule outperforms a DNA language model</text>')
    s.append(yaxis(ymin, ymax, step, "AUPRC (test split, unstratified)"))
    s.append(dashed(0.111, ymin, ymax, "random baseline (0.111)"))

    gw = (R - L) / len(bars)
    bw = gw * 0.45
    for gi, (name, v, color) in enumerate(bars):
        cx = L + gw * (gi + 0.5)
        x = cx - bw / 2
        y = ymap(v, ymin, ymax)
        s.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" '
                 f'height="{(B-y):.1f}" fill="{color}"/>')
        s.append(f'<text x="{cx:.1f}" y="{y-8:.1f}" font-size="13" {FS} '
                 f'fill="#111" text-anchor="middle" font-weight="bold">{v:.4f}</text>')
        for li, line in enumerate(name.split("\n")):
            s.append(f'<text x="{cx:.1f}" y="{B+22+li*16:.0f}" font-size="12" {FS} '
                     f'fill="#333" text-anchor="middle">{esc(line)}</text>')
    s.append("</svg>\n")
    open(path, "w", encoding="utf-8").write("\n".join(s))
    print("wrote", path)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fig_a(os.path.join(here, "fig_A_auroc.svg"))
    fig_b(os.path.join(here, "fig_B_consequence_rule.svg"))
