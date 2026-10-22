# -*- coding: utf-8 -*-
"""Render the two report figures as PNG (Pillow). Supersampled 2x for crisp print.

Run with a Python that has Pillow (e.g. the Codex bundled runtime):
    python make_figs_png.py
"""
import os
from PIL import Image, ImageDraw, ImageFont

S = 2                                   # supersample
W, H = 900 * S, 560 * S
L, R, T, B = 130 * S, 850 * S, 90 * S, 470 * S

FONT_CANDIDATES = [
    r"C:\Windows\Fonts\arial.ttf",
    r"C:\Windows\Fonts\segoeui.ttf",
    r"C:\Windows\Fonts\msyh.ttc",
]


def font(size):
    for p in FONT_CANDIDATES:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size * S)
            except Exception:
                pass
    return ImageFont.load_default()


F_TITLE = font(17)
F_LABEL = font(12)
F_TICK = font(11)
F_VAL = font(12)
F_VALB = font(13)


def ymap(v, ymin, ymax):
    return B - (v - ymin) / (ymax - ymin) * (B - T)


def axes(dr, ymin, ymax, step, ylabel):
    dr.line([(L, T), (L, B)], fill=(51, 51, 51), width=2 * S)
    dr.line([(L, B), (R, B)], fill=(51, 51, 51), width=2 * S)
    n = int(round((ymax - ymin) / step))
    for i in range(n + 1):
        v = ymin + i * step
        y = ymap(v, ymin, ymax)
        dr.line([(L - 7 * S, y), (L, y)], fill=(51, 51, 51), width=S)
        dr.text((L - 12 * S, y), f"{v:.2f}", font=F_TICK, fill=(51, 51, 51), anchor="rm")
    lab = Image.new("RGBA", (300 * S, 24 * S), (255, 255, 255, 0))
    ImageDraw.Draw(lab).text((150 * S, 12 * S), ylabel, font=F_LABEL,
                             fill=(51, 51, 51), anchor="mm")
    lab = lab.rotate(90, expand=True)
    return lab, ((T + B) / 2)


def dashed(dr, v, ymin, ymax, label):
    y = ymap(v, ymin, ymax)
    x = L
    while x < R:
        dr.line([(x, y), (min(x + 8 * S, R), y)], fill=(140, 140, 140), width=S)
        x += 14 * S
    dr.text((R, T - 16 * S), label, font=F_TICK, fill=(110, 110, 110), anchor="rs")


def fig_a(path):
    img = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(img)
    ymin, ymax, step = 0.40, 0.95, 0.10
    dr.text(((L + R) / 2, 40 * S), "Effect of task stratification vs. model scale",
            font=F_TITLE, fill=(17, 17, 17), anchor="mm")
    lab, cy = axes(dr, ymin, ymax, step, "AUROC (test split)")
    img.paste(lab, (14 * S, int(cy - lab.height / 2)), lab)
    dashed(dr, 0.5, ymin, ymax, "random baseline (0.500)")

    groups = [("Unstratified (n = 2,000)", 0.663, 0.711),
              ("Non-protein subset (n = 697)", 0.807, 0.850)]
    series = [("4-mer logistic regression", (138, 180, 248)),
              ("Carbon-500M linear probe", (26, 115, 232))]
    gw = (R - L) / len(groups)
    bw = gw * 0.30
    for gi, (gname, v_kmer, v_carbon) in enumerate(groups):
        cx = L + gw * (gi + 0.5)
        for si, (sname, color, v) in enumerate(zip([s[0] for s in series],
                                                   [s[1] for s in series],
                                                   [v_kmer, v_carbon])):
            off = (si - 0.5) * (bw + 7 * S)
            x0 = cx + off - bw / 2
            y = ymap(v, ymin, ymax)
            dr.rectangle([x0, y, x0 + bw, B], fill=color)
            dr.text((x0 + bw / 2, y - 9 * S), f"{v:.3f}", font=F_VAL,
                    fill=(17, 17, 17), anchor="mb")
        dr.text((cx, B + 26 * S), gname, font=F_LABEL, fill=(51, 51, 51), anchor="ma")

    for i, (sname, color) in enumerate(series):
        y = T + 8 * S + i * 26 * S
        dr.rectangle([L + 10 * S, y, L + 24 * S, y + 13 * S], fill=color)
        dr.text((L + 32 * S, y + 6 * S), sname, font=F_LABEL, fill=(51, 51, 51), anchor="lm")
    img.save(path, dpi=(200, 200))
    print("wrote", path)


def fig_b(path):
    img = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(img)
    ymin, ymax, step = 0.0, 0.6, 0.1
    dr.text(((L + R) / 2, 40 * S), "A consequence rule outperforms a DNA language model",
            font=F_TITLE, fill=(17, 17, 17), anchor="mm")
    lab, cy = axes(dr, ymin, ymax, step, "AUPRC (unstratified test split)")
    img.paste(lab, (14 * S, int(cy - lab.height / 2)), lab)
    dashed(dr, 0.111, ymin, ymax, "random baseline (0.111)")

    bars = [("Consequence-rule baseline\n(0 parameters)", 0.4925, (217, 48, 37)),
            ("Carbon-500M linear probe\n(500M parameters)", 0.4701, (26, 115, 232)),
            ("4-mer logistic regression\n(256 features)", 0.2388, (138, 180, 248))]
    gw = (R - L) / len(bars)
    bw = gw * 0.42
    for gi, (name, v, color) in enumerate(bars):
        cx = L + gw * (gi + 0.5)
        x0 = cx - bw / 2
        y = ymap(v, ymin, ymax)
        dr.rectangle([x0, y, x0 + bw, B], fill=color)
        dr.text((cx, y - 10 * S), f"{v:.4f}", font=F_VALB, fill=(17, 17, 17), anchor="mb")
        for li, line in enumerate(name.split("\n")):
            dr.text((cx, B + 26 * S + li * 18 * S), line, font=F_LABEL,
                    fill=(51, 51, 51), anchor="ma")
    img.save(path, dpi=(200, 200))
    print("wrote", path)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fig_a(os.path.join(here, "fig_A_auroc.png"))
    fig_b(os.path.join(here, "fig_B_consequence_rule.png"))
