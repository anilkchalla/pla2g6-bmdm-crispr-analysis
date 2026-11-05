#!/usr/bin/env python3
"""
Sham-cut-site positional scan — internal negative control when no matched
untreated sample exists.

WHY THIS EXISTS
---------------
CRISPResso2's `Effect_vector_deletion.txt` reports, for each reference position,
how many reads have that position deleted. A single 53 bp deletion therefore
contributes to 53 positions. Any region swept by large deletions is inflated,
and a guide producing mostly small discrete indels looks weak next to a guide
producing large ones. Using that vector to compare two guides is length-biased
and can invert the ranking.

This script instead counts, for every candidate nick position in the amplicon,
the number of reads with an indel BREAKPOINT within +/-W. Each read is counted
once regardless of indel size. Every non-guide position acts as a sham cut site,
so the amplicon-wide median is an empirical background built from the sample
itself.

A real cut site appears as a localized peak that decays to the median within
~15-20 bp on both flanks. Without a matched control this is suggestive, not
proof: report results as putative (protocol section 1).

Usage:
  python3 positional_scan.py --crispresso-dir <dir> --config config.yaml \
      --outdir 03_summary [--window 5]
"""
import argparse
import os
import statistics
import zipfile

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crispresso-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="03_summary")
    ap.add_argument("--amplicon", default="Reference")
    ap.add_argument("--window", type=int, default=5)
    a = ap.parse_args()

    import analyze_alleles as A
    cfg = A.load_cfg(a.config)
    amp = cfg["_amp"]
    p = cfg.get("params", {})
    lo = p.get("exclude_bp_from_left", 15) + 1
    hi = len(amp) - p.get("exclude_bp_from_right", 15)
    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]

    df = A.load_alleles(a.crispresso_dir)
    if "Reference_Name" in df.columns:
        df = df[df["Reference_Name"] == a.amplicon]

    events, N = [], 0
    for _, r in df.iterrows():
        d, i, _ = A.parse_events(r["Aligned_Sequence"], r["Reference_Sequence"])
        d = [(s, e) for s, e in d if e >= lo and s <= hi]
        i = [(l, rr, q) for l, rr, q in i if lo <= l <= hi]
        n = int(r["#Reads"])
        N += n
        if d or i:
            events.append((d, i, n))

    W = a.window

    def anchored(c):
        x, y = c - W, c + W - 1
        t = 0
        for d, i, n in events:
            if (any((x <= s <= y) or (x <= e + 1 <= y) for s, e in d)
                    or any(x <= l <= y or x <= rr <= y for l, rr, _ in i)):
                t += n
        return t

    positions = list(range(lo + W, hi - W))
    prof = {c: anchored(c) for c in positions}
    med = statistics.median(prof.values())

    rows = []
    ranked = sorted(prof.items(), key=lambda kv: -kv[1])
    rank_of = {c: r for r, (c, _) in enumerate(ranked, start=1)}
    for c in positions:
        v = prof[c]
        g = next((names[k] for k, cc in enumerate(cuts) if cc == c), "")
        rows.append({
            "position_1based": c, "guide": g,
            "anchored_reads": v,
            "pct_of_analyzable": round(100 * v / N, 4),
            "fold_over_median": round(v / med, 2) if med else float("nan"),
            "rank": rank_of[c],
        })
    out = pd.DataFrame(rows)
    os.makedirs(a.outdir, exist_ok=True)
    out.to_csv(os.path.join(a.outdir, "positional_scan.tsv"), sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(list(prof), [100 * v / N for v in prof.values()],
            color="#2c7fb8", lw=1.4)
    ax.axhline(100 * med / N, color="k", ls="--", lw=1,
               label="amplicon median (sham background) %.4f%%" % (100 * med / N))
    for c, nm in zip(cuts, names):
        ax.axvline(c, color="#c1392b", ls=":", lw=1.6)
        ax.text(c, ax.get_ylim()[1] * .9, " " + nm, color="#c1392b", fontsize=9)
    ax.set_xlabel("candidate nick position (1-based)")
    ax.set_ylabel("%% reads with indel breakpoint within +/-%d bp" % W)
    ax.set_title("Breakpoint-anchored positional scan (length-unbiased)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "fig5_positional_scan.png"), dpi=200)
    plt.close(fig)

    print("N = %d | callable %d-%d | window +/-%d" % (N, lo, hi, W))
    print("sham background (amplicon median) = %d reads (%.4f%%)"
          % (med, 100 * med / N))
    print()
    for c, nm in zip(cuts, names):
        v = prof.get(c, 0)
        print("%-14s pos %3d : %6d reads  %.4f%%  %.1fx median  rank %d/%d"
              % (nm, c, v, 100 * v / N, v / med, rank_of.get(c, -1), len(positions)))
    print()
    print("top 10 positions overall:")
    for c, v in ranked[:10]:
        tag = next((" <== %s" % names[k] for k, cc in enumerate(cuts)
                    if abs(cc - c) <= 2), "")
        print("  pos %3d  %6d  %.4f%%  %.1fx%s" % (c, v, 100 * v / N, v / med, tag))
    print()
    print("Interpretation: a real cut site is a localized peak decaying to the")
    print("median within ~15-20 bp. Without a matched untreated control this is")
    print("suggestive, not proof - report as putative.")
    print("wrote positional_scan.tsv + fig5_positional_scan.png to %s" % a.outdir)


if __name__ == "__main__":
    main()
