#!/usr/bin/env python3
"""
Long-read (Oxford Nanopore / Plasmidsaurus) CRISPR amplicon analyzer.
Framework v2.0 section 16.

WHY NOT CRISPResso2 DIRECTLY
----------------------------
CRISPResso2's global Needleman-Wunsch alignment and its default gap penalties
are tuned for short, low-indel-error Illumina reads. Nanopore reads are long and
indel-error-prone, and the biologically interesting alleles here are 100-300 bp
excisions. This module uses k-mer anchor chaining (the approach minimap2 takes)
which handles large structural deletions natively and is not confounded by the
1-2 bp homopolymer indels that dominate the ONT error profile.

METHOD
------
1. Index every k-mer of the reference that occurs exactly once (repeat-safe).
2. For each read, and its reverse complement, collect anchors (read_pos, ref_pos).
   Keep the orientation with more anchors.
3. Sort anchors by read position and track the diagonal (ref_pos - read_pos).
   A jump of +n in the diagonal is an n bp deletion; -n is an n bp insertion.
4. Merge adjacent anchors on the same diagonal; emit indels above --min-indel.
5. Attribute each breakpoint to a guide cut site by +/- window proximity, and
   classify reads against the predicted paired-guide excision products.

Homopolymer-adjacent 1-2 bp indels are flagged (ONT's dominant error mode) and
excluded from the primary editing call by default.

Usage:
  python3 ont_amplicon.py --fastq sample.fastq --config config.yaml \
      --outdir 03_summary [--control-fastq control.fastq]
"""
import argparse
import gzip
import math
import os
from collections import Counter, defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
K = 14
MIN_INDEL = 3          # bp; below this ONT error dominates
HOMOPOLYMER_MIN = 4    # run length that makes a short indel suspect


def rc(s):
    return s.translate(COMP)[::-1]


def read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            if seq:
                yield seq


def index_unique_kmers(ref, k=K):
    seen = Counter()
    for i in range(len(ref) - k + 1):
        seen[ref[i:i + k]] += 1
    return {km: i for i, km in
            ((i, ref[i:i + k]) for i in range(len(ref) - k + 1))
            if seen[km] == 1}


def homopolymer_runs(ref, minlen=HOMOPOLYMER_MIN):
    """Set of reference positions (1-based) inside a homopolymer run."""
    out = set()
    i = 0
    while i < len(ref):
        j = i
        while j < len(ref) and ref[j] == ref[i]:
            j += 1
        if j - i >= minlen:
            out.update(range(i + 1, j + 1))
        i = j
    return out


def call_indels(read, idx, k=K, min_indel=MIN_INDEL):
    """
    Return (indels, n_anchors, strand) where indels is a list of
    (kind, ref_start_1based, ref_end_1based, size); kind in {"D","I"}.
    """
    best = None
    for strand, s in (("+", read), ("-", rc(read))):
        anchors = []
        for i in range(0, len(s) - k + 1):
            p = idx.get(s[i:i + k])
            if p is not None:
                anchors.append((i, p))
        if best is None or len(anchors) > len(best[1]):
            best = (strand, anchors)
    strand, anchors = best
    if len(anchors) < 5:
        return [], len(anchors), strand

    # keep a colinear, monotonically increasing chain
    chain = []
    for rp, fp in anchors:
        if chain and fp <= chain[-1][1]:
            continue
        chain.append((rp, fp))

    indels = []
    for (r1, f1), (r2, f2) in zip(chain, chain[1:]):
        dr, df = r2 - r1, f2 - f1
        delta = df - dr
        if delta >= min_indel:                      # deletion in the read
            start = f1 + k                          # 1-based first deleted base
            indels.append(("D", start, start + delta - 1, delta))
        elif -delta >= min_indel:                   # insertion into the read
            size = -delta
            pos = f1 + k
            indels.append(("I", pos, pos, size))
    return indels, len(anchors), strand


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def wstr(k, n):
    lo, hi = wilson(k, n)
    return "%.3f-%.3f" % (100 * lo, 100 * hi)


def analyse(fastq, ref, cuts, names, window, pair_tol, min_indel, hp_positions,
            callable_lo=1, callable_hi=10**9, min_anchors=5):
    idx = index_unique_kmers(ref)
    sig_counts = Counter()
    struct = Counter()
    guide_hits = Counter()
    lengths = []
    n_reads = n_analyzable = 0
    hp_flagged = 0
    n_core = 0          # reads with >=1 breakpoint within +/-window of a cut
    n_any = 0           # reads with >=1 qualifying indel anywhere callable

    pair_defs = []
    for i in range(len(cuts)):
        for j in range(i + 1, len(cuts)):
            c1, c2 = sorted((cuts[i], cuts[j]))
            pair_defs.append((names[i], names[j], c1, c2, c2 - c1))

    # Guard: pairing tolerance must not exceed half the closest inter-cut
    # distance, or adjacent cut sites become mutually assignable.
    if len(cuts) > 1:
        min_gap = min(abs(b - a) for a in cuts for b in cuts if b > a)
        cap = max(3, min_gap // 2 - 1)
        if pair_tol > cap:
            print("  NOTE: pairing tolerance reduced %d -> %d bp "
                  "(closest cut sites are %d bp apart)" % (pair_tol, cap, min_gap))
            pair_tol = cap

    for read in read_fastq(fastq):
        n_reads += 1
        lengths.append(len(read))
        indels, na, strand = call_indels(read, idx, min_indel=min_indel)
        if na < min_anchors:
            continue
        n_analyzable += 1

        keep = []
        for kind, s, e, size in indels:
            if size <= 2 and (s in hp_positions or e in hp_positions):
                hp_flagged += 1
                continue
            # drop events in primer-derived / unreliable read-end sequence
            if e < callable_lo or s > callable_hi:
                continue
            keep.append((kind, s, e, size))

        if not keep:
            sig_counts["WT"] += 1
            struct["WT"] += 1
            continue

        sig = "|".join("%s%d:%d-%d" % (kd, sz, s, e) for kd, s, e, sz in
                       sorted(keep, key=lambda x: x[1]))
        sig_counts[sig] += 1
        n_any += 1

        # guide attribution by breakpoint proximity (one count per guide per read)
        hit_guides = set()
        for kd, s, e, sz in keep:
            for gi, c in enumerate(cuts):
                if (abs(s - c) <= window) or (abs(e + 1 - c) <= window):
                    hit_guides.add(names[gi])
        for g in hit_guides:
            guide_hits[g] += 1
        if hit_guides:
            n_core += 1

        # Structural classification against predicted excision products.
        # BEST match, not first match: when two cut sites are closer together
        # than the pairing tolerance (Ephx2 g1/g3 are 19 bp apart), a first-match
        # rule assigns excisions to whichever pair happens to be tested first.
        # Score each candidate by total endpoint distance and take the minimum.
        label, best_score = "other_indel", None
        for gA, gB, c1, c2, span in pair_defs:
            for kd, s, e, sz in keep:
                if kd != "D":
                    continue
                d1, d2 = abs(s - c1), abs(e + 1 - c2)
                if d1 <= pair_tol and d2 <= pair_tol:
                    score = d1 + d2
                    if best_score is None or score < best_score:
                        best_score, label = score, \
                            "excision_%s_x_%s(%dbp)" % (gA, gB, span)
        struct[label] += 1

    return {"n_reads": n_reads, "n_analyzable": n_analyzable,
            "signatures": sig_counts, "structural": struct,
            "guide_hits": guide_hits, "lengths": lengths,
            "homopolymer_filtered": hp_flagged,
            "n_core": n_core, "n_any": n_any}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--control-fastq", default=None)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="03_summary")
    ap.add_argument("--window", type=int, default=10,
                    help="+/- bp for breakpoint attribution (ONT default 10)")
    ap.add_argument("--pair-tolerance", type=int, default=25)
    ap.add_argument("--min-indel", type=int, default=MIN_INDEL)
    a = ap.parse_args()

    import analyze_alleles as A
    cfg = A.load_cfg(a.config)
    ref = cfg["_amp"]
    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]
    hp = homopolymer_runs(ref)
    os.makedirs(a.outdir, exist_ok=True)

    p = cfg.get("params", {})
    c_lo = p.get("exclude_bp_from_left", 20) + 1
    c_hi = len(ref) - p.get("exclude_bp_from_right", 20)
    ed = analyse(a.fastq, ref, cuts, names, a.window, a.pair_tolerance,
                 a.min_indel, hp, c_lo, c_hi)
    ct = (analyse(a.control_fastq, ref, cuts, names, a.window, a.pair_tolerance,
                  a.min_indel, hp, c_lo, c_hi) if a.control_fastq else None)
    print("callable interval: %d-%d" % (c_lo, c_hi))

    Ne, Nc = ed["n_analyzable"], (ct["n_analyzable"] if ct else 0)

    print("=" * 76)
    print("LONG-READ AMPLICON ANALYSIS  (reference %d bp, k=%d)" % (len(ref), K))
    print("=" * 76)
    print("edited  : %s" % os.path.basename(a.fastq))
    print("  reads %d | analyzable %d (%.1f%%) | homopolymer indels filtered %d"
          % (ed["n_reads"], Ne, 100 * Ne / ed["n_reads"], ed["homopolymer_filtered"]))
    if ct:
        print("control : %s" % os.path.basename(a.control_fastq))
        print("  reads %d | analyzable %d (%.1f%%)"
              % (ct["n_reads"], Nc, 100 * Nc / ct["n_reads"]))
    print()

    # ---- structural allele table -----------------------------------------
    rows = []
    keys = sorted(set(list(ed["structural"]) + (list(ct["structural"]) if ct else [])),
                  key=lambda x: -ed["structural"].get(x, 0))
    for k in keys:
        e = ed["structural"].get(k, 0)
        c = ct["structural"].get(k, 0) if ct else 0
        pe, pc = 100 * e / Ne, (100 * c / Nc if Nc else float("nan"))
        rows.append({
            "structural_class": k,
            "edited_reads": e, "edited_pct": round(pe, 4),
            "edited_ci95": wstr(e, Ne),
            "control_reads": c if ct else "",
            "control_pct": round(pc, 4) if ct else "",
            "absolute_excess_pct": round(pe - pc, 4) if ct else "",
        })
    st = pd.DataFrame(rows)
    st.to_csv(os.path.join(a.outdir, "structural_alleles.tsv"), sep="\t", index=False)
    print("STRUCTURAL ALLELES")
    print(st.to_string(index=False))
    print()

    # ---- guide-level -------------------------------------------------------
    gr = []
    for i, g in enumerate(names):
        e = ed["guide_hits"].get(g, 0)
        c = ct["guide_hits"].get(g, 0) if ct else 0
        pe, pc = 100 * e / Ne, (100 * c / Nc if Nc else float("nan"))
        gr.append({"guide": g, "cut_site_1based": cuts[i],
                   "edited_reads_with_breakpoint": e,
                   "edited_pct": round(pe, 4), "ci95": wstr(e, Ne),
                   "control_pct": round(pc, 4) if ct else "",
                   "absolute_excess_pct": round(pe - pc, 4) if ct else ""})
    gd = pd.DataFrame(gr)
    gd.to_csv(os.path.join(a.outdir, "guide_efficiency_ont.tsv"), sep="\t", index=False)
    print("GUIDE BREAKPOINT ATTRIBUTION (+/-%d bp)" % a.window)
    print(gd.to_string(index=False))
    print()

    # ---- allele signatures -------------------------------------------------
    sr = []
    for sig, n in ed["signatures"].most_common(200):
        c = ct["signatures"].get(sig, 0) if ct else 0
        sr.append({"allele_signature": sig, "edited_reads": n,
                   "edited_pct": round(100 * n / Ne, 4),
                   "ci95": wstr(n, Ne),
                   "control_reads": c if ct else "",
                   "control_pct": round(100 * c / Nc, 4) if Nc else ""})
    pd.DataFrame(sr).to_csv(os.path.join(a.outdir, "allele_signatures_ont.tsv"),
                            sep="\t", index=False)

    # ---- figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 3.6))
    bins = range(0, max(ed["lengths"]) + 25, 10)
    ax.hist(ed["lengths"], bins=bins, alpha=.65, label="edited", color="#c1392b")
    if ct:
        ax.hist(ct["lengths"], bins=bins, alpha=.55, label="control", color="#2c7fb8")
    ax.axvline(len(ref), color="k", ls="--", lw=1, label="WT %d bp" % len(ref))
    for gA, gB, c1, c2, span in [(names[i], names[j], sorted((cuts[i], cuts[j]))[0],
                                  sorted((cuts[i], cuts[j]))[1],
                                  abs(cuts[j] - cuts[i]))
                                 for i in range(len(cuts))
                                 for j in range(i + 1, len(cuts))]:
        L = len(ref) - span
        ax.axvline(L, color="#e67e22", ls=":", lw=1.4)
        ax.text(L, ax.get_ylim()[1] * .8, " %s×%s\n %dbp" % (gA[:2], gB[:2], L),
                fontsize=7, color="#e67e22")
    ax.set_xlim(0, len(ref) * 1.4)
    ax.set_xlabel("read length (bp)")
    ax.set_ylabel("reads")
    ax.set_title("Read-length spectrum vs predicted excision products")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "fig6_ont_length_spectrum.png"), dpi=200)
    plt.close(fig)

    print("PRIMARY EDITING METRICS")
    summ = []
    for lbl, key in [("core indel efficiency (breakpoint within +/-%d of a cut)"
                      % a.window, "n_core"),
                     ("any qualifying indel in callable interval", "n_any")]:
        e = ed[key]
        pe = 100 * e / Ne
        row = {"metric": lbl, "edited_reads": e, "edited_pct": round(pe, 3),
               "ci95": wstr(e, Ne)}
        if ct:
            c = ct[key]
            pc = 100 * c / Nc
            row["control_pct"] = round(pc, 3)
            row["absolute_excess_pct"] = round(pe - pc, 3)
            row["normalized_adjusted_pct"] = round(
                100 * max(0.0, (pe / 100 - pc / 100) / (1 - pc / 100)), 3)
        summ.append(row)
    sm = pd.DataFrame(summ)
    sm.to_csv(os.path.join(a.outdir, "sample_metrics_ont.tsv"), sep="\t", index=False)
    print(sm.to_string(index=False))
    print()
    print("wrote 3 tables + 1 figure to %s" % a.outdir)


if __name__ == "__main__":
    main()
