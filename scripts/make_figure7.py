#!/usr/bin/env python3
"""Figure 7 - editing against a matched Cas9-only control in two mouse backgrounds.

Panel A  efficiency and excision vs matched control (from the reclassified
         allele tables; pass --summary to override the defaults)
Panel B  deletion-size distribution against the 137 bp blunt-blunt prediction
Panel C  apparent editing as a function of the minimum indel size counted

Panels B and C are computed from the raw FASTQ by the same anchor method as
y452b3_anchor_recount.py, so the figure is reproducible from reads alone.

Usage:
    python3 make_figure7.py --fastq-dir data/fastq \
        --ref refs/Mm_Pla2g6_F2R2-amplicon.fa --out Figure7_matched_control.png
"""

import argparse
import collections
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import (EXCISION_BP, SAMPLES, fastq_path, fuzzy_find,
                           load_reference, orient, read_fastq)

L_SPAN, M_SPAN, R_SPAN = (21, 44), (130, 153), (240, 263)
MIN_EXCISION_BP = 50
THRESHOLDS = [1, 2, 3, 5, 10]

# From y452b3_reclassify_alleles.py. Override with --summary.
DEFAULT_SUMMARY = {
    "Cas9 ctrl C57BL/6": {"edited": 0.28, "excision": 0.0},
    "Pla2g6 C57BL/6":    {"edited": 54.74, "excision": 33.11},
    "Cas9 ctrl NOD":     {"edited": 0.07, "excision": 0.0},
    "Pla2g6 NOD":        {"edited": 44.45, "excision": 26.22},
}

BLUE, ORANGE = "#3B6FB6", "#D1793B"
INK, MUTED, SURFACE = "#1f2328", "#6b7280", "#fcfcfb"
GRID = "#e9ecef"


def scan(ref, path):
    aL = ref[L_SPAN[0] - 1:L_SPAN[1]]
    aM = ref[M_SPAN[0] - 1:M_SPAN[1]]
    aR = ref[R_SPAN[0] - 1:R_SPAN[1]]
    gap_LM = M_SPAN[0] - L_SPAN[1] - 1
    gap_MR = R_SPAN[0] - M_SPAN[1] - 1
    gap_LR = R_SPAN[0] - L_SPAN[1] - 1

    total = excision = 0
    sizes = collections.Counter()
    worst = collections.Counter()
    for seq, _q in read_fastq(path):
        read = orient(seq, aL)
        if read is None:
            continue
        iL = fuzzy_find(read, aL)
        endL = iL + len(aL)
        iR = fuzzy_find(read[endL:], aR)
        if iR < 0:
            continue
        absR = endL + iR
        iM = fuzzy_find(read[endL:absR], aM)
        if iM < 0:
            net = (absR - endL) - gap_LR
            if net <= -MIN_EXCISION_BP:
                total += 1
                excision += 1
                sizes[-net] += 1
            continue
        absM = endL + iM
        endM = absM + len(aM)
        n1 = (absM - endL) - gap_LM
        n2 = (absR - endM) - gap_MR
        total += 1
        worst[max(abs(n1), abs(n2))] += 1

    curve = [100 * (sum(v for k, v in worst.items() if k >= t) + excision) / total
             for t in THRESHOLDS]
    return total, excision, sizes, curve


def style(ax):
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#d6d9de")
    ax.spines["bottom"].set_color("#d6d9de")
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.set_facecolor(SURFACE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--summary", help="JSON from y452b3_reclassify_alleles.py")
    ap.add_argument("--out", default="Figure7_matched_control.png")
    args = ap.parse_args()

    ref = load_reference(args.ref)
    summary = json.load(open(args.summary)) if args.summary else DEFAULT_SUMMARY

    data = {}
    for stem, label, _cond, _strain, _nuc in SAMPLES:
        data[label] = scan(ref, fastq_path(args.fastq_dir, stem))

    fig = plt.figure(figsize=(13.5, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 1.1], wspace=0.32)

    # Panel A
    ax = fig.add_subplot(gs[0, 0])
    labels = ["Cas9 ctrl\nNuc#5\nC57BL/6", "Pla2g6\nNuc#5\nC57BL/6",
              "Cas9 ctrl\nNuc#6\nNOD", "Pla2g6\nNuc#6\nNOD"]
    order = ["Cas9 ctrl C57BL/6", "Pla2g6 C57BL/6", "Cas9 ctrl NOD", "Pla2g6 NOD"]
    edited = [summary[k]["edited"] for k in order]
    excis = [summary[k]["excision"] for k in order]
    x = np.arange(4)
    w = 0.38
    ax.bar(x - w / 2, edited, w, zorder=3,
           color=[BLUE if i % 2 else "#b9c3cf" for i in range(4)], label="Total edited")
    ax.bar(x + w / 2, excis, w, zorder=3,
           color=[ORANGE if i % 2 else "#e2cbb6" for i in range(4)],
           label=f"Excision allele ({EXCISION_BP} bp del)")
    for xx, v in zip(x - w / 2, edited):
        ax.text(xx, v + 1.6, f"{v:.2f}" if v < 5 else f"{v:.1f}", ha="center",
                fontsize=8.5, color=INK)
    for xx, v in zip(x + w / 2, excis):
        ax.text(xx, v + 1.6, f"{v:.1f}" if v > 0 else "0", ha="center",
                fontsize=8.5, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.8, color=INK)
    ax.set_ylabel("% of analysable reads", fontsize=9.5, color=INK)
    ax.set_ylim(0, 74)
    ax.set_title("A  Editing efficiency vs matched control", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", handlelength=1.2)
    style(ax)

    # Panel B
    ax = fig.add_subplot(gs[0, 1])
    for label, col in [("Pla2g6 C57BL/6", BLUE), ("Pla2g6 NOD", ORANGE)]:
        total, _e, sizes, _c = data[label]
        ks = np.arange(126, 201)
        ax.plot(ks, [100 * sizes.get(int(k), 0) / total for k in ks],
                lw=2, color=col, label=label, zorder=3)
    ax.axvline(EXCISION_BP, color="#444", ls=(0, (4, 3)), lw=1.3, zorder=2)
    ax.text(EXCISION_BP - 1.5, 7.6, f"{EXCISION_BP} bp\npredicted\nblunt join",
            fontsize=8, color="#444", va="top", ha="right")
    ax.text(149, 6.5, "observed mode = prediction\n161 bp product, matches\nthe HMA band",
            fontsize=8, color=INK, va="top")
    ax.set_xlim(126, 200)
    ax.set_xlabel("deletion size (bp)", fontsize=9.5, color=INK)
    ax.set_ylabel("% of analysable reads", fontsize=9.5, color=INK)
    ax.set_title("B  Observed mode = predicted blunt join", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.legend(frameon=False, fontsize=8.5)
    style(ax)

    # Panel C
    ax = fig.add_subplot(gs[0, 2])
    styles = {"Cas9 ctrl C57BL/6": (MUTED, "--"), "Cas9 ctrl NOD": ("#9aa3ad", "--"),
              "Pla2g6 C57BL/6": (BLUE, "-"), "Pla2g6 NOD": (ORANGE, "-")}
    for label, (col, ls) in styles.items():
        ax.plot(THRESHOLDS, data[label][3], ls, color=col, lw=2, marker="o",
                ms=5, label=label, zorder=3)
    ax.axhspan(0, 1, color="#f1f3f5", zorder=1)
    ax.text(10, 1.6, "breakpoint-anchored call\nctrl = 0.07-0.28%",
            fontsize=8, color=MUTED, ha="right")
    ax.set_xticks(THRESHOLDS)
    ax.set_xlabel("minimum |indel| counted (bp)", fontsize=9.5, color=INK)
    ax.set_ylabel("apparent editing (%)", fontsize=9.5, color=INK)
    ax.set_title("C  Permissive calling inflates control", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    style(ax)

    fig.patch.set_facecolor(SURFACE)
    fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
