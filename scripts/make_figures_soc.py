#!/usr/bin/env python3
"""Main-text Figures 4, 5 and 6, regenerated from the control-matched dataset (SO-C).

Every panel that scores editing carries its own Cas9-only control, which is the
point of the revision: no main-text number rests on an uncontrolled dataset.

Usage:
    python3 make_figures_soc.py --spectrum so_c_spectrum.json \
        --fastq-dir data/fastq --ref refs/Mm_Pla2g6_F2R2-amplicon.fa --outdir .
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import (ANCHOR_LEN, ANCHOR_OFFSET, CUT1, CUT2, NO_CUT_WINDOWS,
                           SAMPLES, fastq_path, load_reference, read_fastq)
from y452b3_error_floor import score_site

BLUE, ORANGE = "#3B6FB6", "#D1793B"
GREY, LGREY = "#6b7280", "#9aa3ad"
INK, MUTED, SURFACE, GRID = "#1f2328", "#6b7280", "#fcfcfb", "#e9ecef"

DUAL = ["Pla2g6 C57BL/6", "Pla2g6 NOD"]
CTRL = ["Cas9 ctrl C57BL/6", "Cas9 ctrl NOD"]
SHORT = {"Pla2g6 C57BL/6": "SO-C\nC57BL/6", "Pla2g6 NOD": "SO-C\nNOD",
         "Cas9 ctrl C57BL/6": "ctrl\nC57BL/6", "Cas9 ctrl NOD": "ctrl\nNOD"}

# Four analyses of the same II reads (preprint Section 2.8).
II_ESTIMATES = [("CR2\nas run", 1.51), ("CR2\nindels only", 1.36),
                ("CR2 re-quant.\nat cut sites", 1.80), ("anchor\nmethod", 0.91)]
II_NOCUT_ANY = 0.0039   # Illumina no-cut-window indel rate, %


def style(ax, legend=False):
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#d6d9de")
    ax.spines["bottom"].set_color("#d6d9de")
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.set_facecolor(SURFACE)
    if legend:
        ax.legend(frameon=False, fontsize=8)


def finish(fig, path):
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    print("wrote", path)
    plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def figure4(spec, outdir):
    fig = plt.figure(figsize=(13.5, 8.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26)

    # A: allele classes at two thresholds
    ax = fig.add_subplot(gs[0, 0])
    classes = ["unedited", "cut1_only", "cut2_only", "both_independent", "excision"]
    nice = ["unedited", "cut 1 only", "cut 2 only", "both,\nindependent", "excision"]
    order = DUAL + CTRL
    x = np.arange(len(classes))
    w = 0.2
    colors = {DUAL[0]: BLUE, DUAL[1]: ORANGE, CTRL[0]: GREY, CTRL[1]: LGREY}
    for i, lab in enumerate(order):
        e = spec[lab]
        n = e["scorable"]
        vals = [100 * e["classes"][c] / n for c in classes]
        vals3 = [100 * e["classes_min3bp"][c] / n for c in classes]
        pos = x + (i - 1.5) * w
        ax.bar(pos, vals, w * 0.92, color=colors[lab], zorder=3,
               label=lab if i < 2 else lab + " (control)")
        ax.plot(pos, vals3, ls="none", marker="_", ms=9, mew=2, color=INK, zorder=4)
    ax.plot([], [], ls="none", marker="_", ms=9, mew=2, color=INK,
            label="same, counting only\nindels $\\geq$3 bp")
    ax.set_xticks(x)
    ax.set_xticklabels(nice, fontsize=8.5, color=INK)
    ax.set_ylabel("% of scorable reads", fontsize=9.5, color=INK)
    ax.set_title("A  Allele classes, with matched controls", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.set_ylim(0, 128)
    style(ax)
    ax.legend(frameon=False, fontsize=7.4, ncol=2, loc="upper center")

    # B, C: per-site indel spectra with control trace
    for panel, (key, site, title) in enumerate(
            [("site1_spectrum", CUT1, "B  Indel spectrum at cut 1 (73|74)"),
             ("site2_spectrum", CUT2, "C  Indel spectrum at cut 2 (210|211)")]):
        ax = fig.add_subplot(gs[0, 1] if panel == 0 else gs[1, 0])
        ks = np.array([k for k in range(-20, 11) if k != 0])  # unmodified excluded
        for lab, col in [(DUAL[0], BLUE), (DUAL[1], ORANGE)]:
            e = spec[lab][key]
            v = [100 * e["hist"].get(str(int(k)), 0) / e["n"] for k in ks]
            ax.bar(ks + (0.2 if lab == DUAL[1] else -0.2), v, 0.4, color=col,
                   zorder=3, label=lab)
        ctrl = spec[CTRL[0]][key]
        cv = [100 * ctrl["hist"].get(str(int(k)), 0) / ctrl["n"] for k in ks]
        ax.step(ks, cv, where="mid", color=INK, lw=1.4, zorder=5,
                label="Cas9-only control")
        ax.set_xlim(-20.5, 10.5)
        shown = [100 * spec[l][key]["hist"].get(str(int(k)), 0) / spec[l][key]["n"]
                 for l in DUAL for k in ks] + cv
        ax.set_ylim(0, max(shown) * 1.28)
        ax.set_xlabel("indel size (bp, deletions negative)", fontsize=9.5, color=INK)
        ax.set_ylabel("% of reads scored at this site", fontsize=9.5, color=INK)
        ax.annotate("unmodified reads excluded", xy=(0.02, 0.92), xycoords="axes fraction",
                    fontsize=7.6, color=MUTED)
        ax.set_title(title, fontsize=10.5, color=INK, loc="left", pad=10)
        style(ax)
        ax.legend(frameon=False, fontsize=8)

    # D: reading frame
    ax = fig.add_subplot(gs[1, 1])
    labs = [SHORT[l] for l in order]
    fs = [100 * spec[l]["frame"]["frameshift"] / spec[l]["scorable"] for l in order]
    inf = [100 * spec[l]["frame"]["in_frame"] / spec[l]["scorable"] for l in order]
    xx = np.arange(len(order))
    ax.bar(xx, fs, 0.55, color=[colors[l] for l in order], zorder=3, label="frameshift")
    ax.bar(xx, inf, 0.55, bottom=fs, color="#cfd6de", zorder=3, label="in-frame")
    for i, (a, b) in enumerate(zip(fs, inf)):
        ax.text(i, a / 2, f"{a:.0f}", ha="center", va="center", fontsize=8.5, color="white")
        ax.text(i, a + b / 2, f"{b:.0f}", ha="center", va="center", fontsize=8.5, color=INK)
    ax.set_xticks(xx)
    ax.set_xticklabels(labs, fontsize=8.5, color=INK)
    ax.set_ylabel("% of scorable reads", fontsize=9.5, color=INK)
    ax.set_title("D  Edited alleles by reading frame", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    style(ax, legend=True)
    finish(fig, os.path.join(outdir, "Figure4_spectrum_soc.pdf"))


# ---------------------------------------------------------------- Figure 5
def figure5(spec, ref, fastq_dir, outdir):
    fig = plt.figure(figsize=(13.5, 4.3))
    gs = fig.add_gridspec(1, 3, wspace=0.34)

    # A: background by position - sham windows vs real cut sites
    ax = fig.add_subplot(gs[0, 0])
    sites = [(f"{w}", w) for w in NO_CUT_WINDOWS] + [("cut 1", CUT1), ("cut 2", CUT2)]
    for lab, col in [(CTRL[0], GREY), (CTRL[1], LGREY)]:
        stem = [s for s, l, *_ in SAMPLES if l == lab][0]
        reads = [s for s, _q in read_fastq(fastq_path(fastq_dir, stem))]
        vals = [score_site(ref, reads, site)[1] for _n, site in sites]
        ax.plot(range(len(sites)), vals, marker="o", ms=6, lw=2, color=col,
                label=lab, zorder=3)
    ax.axhline(II_NOCUT_ANY, color=BLUE, ls=(0, (4, 3)), lw=1.6, zorder=2)
    ax.text(len(sites) - 1, II_NOCUT_ANY * 1.4, "II no-cut windows\n0.0039%",
            fontsize=7.8, color=BLUE, ha="right")
    ax.set_yscale("log")
    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels([n for n, _ in sites], fontsize=8.5, color=INK)
    ax.set_xlabel("amplicon position scored", fontsize=9.5, color=INK)
    ax.set_ylabel("apparent indel rate (%)", fontsize=9.5, color=INK)
    ax.set_title("A  Background depends on where you look", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    style(ax, legend=True)

    # B: clustering, with the control floor
    ax = fig.add_subplot(gs[0, 1])
    rows = [("II", 0.680, 45.9, BLUE), ("SO-C\nC57BL/6", spec[DUAL[0]]["clustering"]["p_site2"],
             spec[DUAL[0]]["clustering"]["p_site2_given_site1"], ORANGE),
            ("SO-C\nNOD", spec[DUAL[1]]["clustering"]["p_site2"],
             spec[DUAL[1]]["clustering"]["p_site2_given_site1"], "#8a5a2b"),
            ("control\nfloor", spec[CTRL[0]]["clustering"]["p_site2"],
             spec[CTRL[0]]["clustering"]["p_site2_given_site1"], GREY)]
    xx = np.arange(len(rows))
    ax.bar(xx - 0.19, [r[1] for r in rows], 0.36, color="#cfd6de", zorder=3,
           label="P(site 2)")
    ax.bar(xx + 0.19, [r[2] for r in rows], 0.36, color=[r[3] for r in rows],
           zorder=3, label="P(site 2 | site 1)")
    ax.set_yscale("log")
    ax.set_xticks(xx)
    ax.set_xticklabels([r[0] for r in rows], fontsize=8, color=INK)
    ax.set_ylabel("probability (%)", fontsize=9.5, color=INK)
    ax.set_title("B  Cutting is clustered, above a real floor", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.set_ylim(0.3, 400)
    style(ax)
    ax.legend(frameon=False, fontsize=7.6, loc="upper left", ncol=2)

    # C: II estimate range
    ax = fig.add_subplot(gs[0, 2])
    names = [n for n, _ in II_ESTIMATES]
    vals = [v for _n, v in II_ESTIMATES]
    ax.bar(range(len(vals)), vals, 0.6, color=BLUE, zorder=3)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.04, f"{v:.2f}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(names, fontsize=7.8, color=INK)
    ax.set_ylabel("core editing estimate (%)", fontsize=9.5, color=INK)
    ax.set_ylim(0, 2.1)
    ax.set_title("C  One II dataset, four defensible answers", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    style(ax)
    finish(fig, os.path.join(outdir, "Figure5_validation_soc.pdf"))


# ---------------------------------------------------------------- Figure 6
def figure6(spec, outdir):
    fig = plt.figure(figsize=(13.5, 4.3))
    gs = fig.add_gridspec(1, 3, wspace=0.32)
    rar = spec[DUAL[0]]["rarefaction"]
    depths = sorted(int(d) for d in rar)

    # A: class frequency rarefaction
    ax = fig.add_subplot(gs[0, 0])
    series = [("unedited", "#7b8794"), ("excision", ORANGE),
              ("cut1_only", BLUE), ("cut2_only", "#8a5a2b"),
              ("both_independent", "#b0b7c3")]
    for key, col in series:
        m = np.array([rar[str(d)][key]["mean"] for d in depths])
        sd = np.array([rar[str(d)][key]["sd"] for d in depths])
        ax.plot(depths, m, lw=2, color=col, label=key.replace("_", " "), zorder=3)
        ax.fill_between(depths, m - sd, m + sd, color=col, alpha=0.18, zorder=2)
    ax.axvline(1000, color="#444", ls=(0, (4, 3)), lw=1.2, zorder=2)
    ax.set_xlabel("reads subsampled", fontsize=9.5, color=INK)
    ax.set_ylabel("% of scorable reads", fontsize=9.5, color=INK)
    ax.set_title("A  Spectrum saturates by ~1,000 reads", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    ax.set_ylim(0, 52)
    style(ax)
    ax.legend(frameon=False, fontsize=7.4, ncol=3, loc="upper center")

    # B: distinct alleles
    ax = fig.add_subplot(gs[0, 1])
    for lab, col in [(DUAL[0], BLUE), (DUAL[1], ORANGE)]:
        r = spec[lab]["rarefaction"]
        ds = sorted(int(d) for d in r)
        m = [r[str(d)]["distinct_alleles"]["mean"] for d in ds]
        sd = [r[str(d)]["distinct_alleles"]["sd"] for d in ds]
        ax.errorbar(ds, m, yerr=sd, lw=2, color=col, marker="o", ms=4,
                    capsize=2, label=lab, zorder=3)
    ax.set_xlabel("reads subsampled", fontsize=9.5, color=INK)
    ax.set_ylabel("distinct indel sizes recovered", fontsize=9.5, color=INK)
    ax.set_title("B  The allele tail is not saturated", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    style(ax, legend=True)

    # C: detection floor
    ax = fig.add_subplot(gs[0, 2])
    ns = {"SO-C\n(~1,900 reads)": 1949, "II\n(401,563 reads)": 401563}
    for k, (name, n) in enumerate(ns.items()):
        for j, reads_needed in enumerate([3, 10]):
            ax.bar(k + (j - 0.5) * 0.36, 100 * reads_needed / n, 0.34,
                   color=[ORANGE, BLUE][j], zorder=3,
                   label=f"{reads_needed} supporting reads" if k == 0 else None)
    ax.set_yscale("log")
    ax.set_xticks(range(len(ns)))
    ax.set_xticklabels(list(ns), fontsize=8.5, color=INK)
    ax.set_ylabel("lowest detectable allele freq. (%)", fontsize=9.5, color=INK)
    ax.set_title("C  Depth sets the detection floor", fontsize=10.5,
                 color=INK, loc="left", pad=10)
    style(ax, legend=True)
    finish(fig, os.path.join(outdir, "Figure6_depth_soc.pdf"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spectrum", required=True)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    spec = json.load(open(args.spectrum))
    ref = load_reference(args.ref)
    figure4(spec, args.outdir)
    figure5(spec, ref, args.fastq_dir, args.outdir)
    figure6(spec, args.outdir)


if __name__ == "__main__":
    main()
