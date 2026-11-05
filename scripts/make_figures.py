#!/usr/bin/env python3
"""Report figures. Palette: validated categorical slots 1-3 (blue/orange/aqua)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

C = {"g1": "#2a78d6", "g2": "#eb6834", "g3": "#1baf7a"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e3e2df"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlecolor": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
})


def style(ax):
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)


def fig_ranking(projects, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.2))
    for ax, (name, d, g) in zip(axes, projects):
        g = g.sort_values("core_indel_pct")
        for i, (_, r) in enumerate(g.iterrows()):
            lo, hi = (float(x) for x in r["core_ci95"].split("-"))
            ax.barh(i, r["core_indel_pct"], height=.5, color=C[r["guide"]],
                    zorder=3)
            ax.plot([lo, hi], [i, i], color="white", lw=3.0, zorder=4)
            ax.plot([lo, hi], [i, i], color=INK2, lw=1.2, zorder=5)
            ax.plot([r["background_pct"]], [i], marker="|", ms=11, mew=1.8,
                    color=INK, zorder=6)
            ax.text(hi + 2.5, i, "%.1f%%   adj %.1f%%"
                    % (r["core_indel_pct"], r["adjusted_mutagenicity_pct"]),
                    va="center", fontsize=8.5, color=INK)
        ax.set_yticks(range(len(g)))
        ax.set_yticklabels(["%s  (%s)" % (r["guide"], r["sample"])
                            for _, r in g.iterrows()])
        ax.set_xlim(0, 118)
        ax.set_xticks([0, 20, 40, 60, 80, 100])
        ax.set_xlabel("core indel % of analyzable reads")
        ax.set_title(name, loc="left", fontweight="bold")
        style(ax)
    axes[0].plot([], [], marker="|", ls="none", ms=11, mew=1.8, color=INK,
                 label="internal background")
    axes[0].plot([], [], color=INK2, lw=1.2, label="95% Wilson CI")
    axes[0].legend(frameon=False, fontsize=7.5, loc="lower right", ncol=1)
    fig.suptitle("Guide mutagenicity, ranked  ·  bar = raw core indel rate, "
                 "adj = after subtracting the internal background",
                 x=0.008, ha="left", fontsize=9.5, color=INK2)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_profiles(d, out, title):
    ss = list(d["samples"].items())
    fig, axes = plt.subplots(len(ss), 1, figsize=(9.6, 1.55 * len(ss)),
                             sharex=True)
    for ax, (sname, v) in zip(axes, ss):
        xs = sorted(int(k) for k in v["bp_profile"])
        ys = [100 * v["bp_profile"][str(x)] / v["N"] for x in xs]
        ax.fill_between(xs, ys, color=C[v["guide"]], alpha=.18, zorder=2)
        ax.plot(xs, ys, color=C[v["guide"]], lw=2, zorder=3)
        for g, c in d["cuts"].items():
            own = g == v["guide"]
            ax.axvline(c, color=C[g], ls="-" if own else ":",
                       lw=1.6 if own else 1.0, alpha=.9 if own else .55, zorder=4)
            ax.text(c, ax.get_ylim()[1] * .02, " " + g, fontsize=7.5,
                    color=C[g], rotation=90, va="bottom")
        ax.set_ylabel("% reads", fontsize=8)
        ax.text(.005, .88, "%s  ·  injected %s  ·  N=%d"
                % (sname, v["guide"], v["N"]), transform=ax.transAxes,
                fontsize=8.5, fontweight="bold", color=INK)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color=GRID, lw=.7)
    axes[-1].set_xlabel("amplicon position (1-based)  ·  breakpoint within ±%d bp"
                        % d["window"])
    fig.suptitle(title, x=0.008, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_spectrum(d, out, title):
    ss = list(d["samples"].items())
    fig, axes = plt.subplots(1, len(ss), figsize=(9.6, 2.7), sharey=False)
    for ax, (sname, v) in zip(axes, ss):
        tot = v["on_target_reads"] or 1
        xs, ys = [], []
        for k, n in sorted(v["del_sizes"].items(), key=lambda kv: int(kv[0])):
            xs.append(-int(k))
            ys.append(100 * n / tot)
        for k, n in sorted(v["ins_sizes"].items(), key=lambda kv: int(kv[0])):
            xs.append(int(k))
            ys.append(100 * n / tot)
        ax.bar(xs, ys, width=.85, color=C[v["guide"]], zorder=3)
        ax.axvline(0, color=INK2, lw=.8, zorder=2)
        ax.set_xlim(-32, 32)
        ax.set_title("%s  (%s)" % (sname, v["guide"]), loc="left", fontsize=9,
                     fontweight="bold")
        ax.set_xlabel("← deletion    size (bp)    insertion →", fontsize=8)
        ax.grid(axis="y", color=GRID, lw=.7)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("% of on-target edited reads", fontsize=8)
    fig.suptitle(title, x=0.008, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_frame(projects, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 2.6))
    for ax, (name, d, _) in zip(axes, projects):
        labels, fs, inf = [], [], []
        for sname, v in d["samples"].items():
            t = sum(v["net_on"].values()) or 1
            f = sum(n for k, n in v["net_on"].items() if int(k) % 3)
            labels.append("%s\n%s" % (v["guide"], sname))
            fs.append(100 * f / t)
            inf.append(100 * (t - f) / t)
        x = range(len(labels))
        ax.bar(x, fs, width=.55, color="#2a78d6", label="frameshift", zorder=3)
        ax.bar(x, inf, width=.55, bottom=[a + 0.6 for a in fs], color="#eb6834",
               label="in-frame", zorder=3)
        for i, a in enumerate(fs):
            ax.text(i, a / 2, "%.0f%%" % a, ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, 105)
        ax.set_ylabel("% of on-target edited reads", fontsize=8)
        ax.set_title(name, loc="left", fontweight="bold")
        ax.grid(axis="y", color=GRID, lw=.7)
        ax.set_axisbelow(True)
    axes[0].legend(frameon=False, fontsize=8, loc="lower left", ncol=2)
    fig.suptitle("Predicted frame consequence of on-target alleles",
                 x=0.008, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_crosstalk(d, out):
    v = d["samples"]["bend33"]
    fig, ax = plt.subplots(figsize=(9.6, 2.6))
    xs = sorted(int(k) for k in v["bp_positions"])
    ys = [100 * v["bp_positions"][str(x)] / v["N"] for x in xs]
    ax.bar(xs, ys, width=1.0, color="#52514e", zorder=3)
    for g, c in d["cuts"].items():
        ax.axvline(c, color=C[g], lw=1.8, zorder=4)
        ax.text(c, max(ys) * .95, " " + g, color=C[g], fontsize=9,
                fontweight="bold")
    ax.annotate("expected: g3 only", xy=(113, max(ys) * .62),
                xytext=(70, max(ys) * .8), fontsize=8.5, color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=.9))
    ax.annotate("unexpected second cluster:\nbend32 (g2) allele repertoire",
                xy=(134, max(ys) * .8), xytext=(150, max(ys) * .72),
                fontsize=8.5, color="#e34948",
                arrowprops=dict(arrowstyle="->", color="#e34948", lw=1.1))
    ax.set_xlabel("amplicon position (1-based)")
    ax.set_ylabel("% reads with a breakpoint here", fontsize=8)
    ax.set_title("bend33 (guide 3) carries a second, foreign breakpoint cluster",
                 loc="left", fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    projects = []
    for p, label in (("bend3", "bend3  (Dr_bend3_Ex3, 252 bp)"),
                     ("dlx3b", "dlx3b  (Dr_dlx3b_Ex1, 377 bp)")):
        d = json.load(open("%s/report_data.json" % p))
        g = pd.read_csv("%s/03_summary/guide_efficiency.tsv" % p, sep="\t")
        projects.append((label, d, g))
        os.makedirs("%s/04_report" % p, exist_ok=True)
        fig_profiles(d, "%s/04_report/fig2_positional_profile.png" % p,
                     "Breakpoint-anchored positional profile — %s" % label)
        fig_spectrum(d, "%s/04_report/fig3_indel_spectrum.png" % p,
                     "Indel size spectrum at the injected guide's cut — %s" % label)
    fig_ranking(projects, "figures/fig1_ranking.png")
    fig_frame(projects, "figures/fig4_frame.png")
    fig_crosstalk(projects[0][1], "figures/fig5_crosstalk.png")
    print("figures written")
