#!/usr/bin/env python3
"""
Substitution audit (protocol section 8).

Substitutions are reported but NEVER counted toward nuclease editing efficiency.
Run this on the second CRISPResso2 pass (quantification_window_size 0, without
--ignore_substitutions).

Reporting rules applied (section 8.1). A noticeable SNP candidate must satisfy:
  - frequency >= 1.0% of analyzable reads AND >= 20 supporting reads
  - outside primer-derived sequence (callable interval)
  - not confined to the first/last unreliable bases
  - control frequency compared when a matched control is supplied

Classification labels (section 8.2):
  baseline_polymorphism | noticeable_SNP_candidate | likely_technical_artifact
  | indel_associated_substitution

Usage:
  python3 snp_audit.py --crispresso-dir <audit pass dir> --config config.yaml \
      --outdir 03_summary [--control-dir <control audit dir>]
"""
import argparse
import os
from collections import defaultdict

import pandas as pd

MIN_PCT = 1.0
MIN_READS = 20
END_UNRELIABLE = 5     # bp at each end of the callable interval treated as unreliable


def collect(crispresso_dir, config, amplicon, callable_lo, callable_hi):
    import analyze_alleles as A
    cfg = A.load_cfg(config)
    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]
    df = A.load_alleles(crispresso_dir)
    if "Reference_Name" in df.columns:
        df = df[df["Reference_Name"] == amplicon]

    pos_counts = defaultdict(int)      # (pos, ref, alt) -> reads
    on_indel = defaultdict(int)
    total = 0
    for _, r in df.iterrows():
        d, i, s = A.parse_events(r["Aligned_Sequence"], r["Reference_Sequence"])
        n = int(r["#Reads"])
        total += n
        d = [(x, y) for x, y in d if y >= callable_lo and x <= callable_hi]
        i = [(l, rr, q) for l, rr, q in i if callable_lo <= l <= callable_hi]
        has_indel = bool(d or i)
        for p, ref, alt in s:
            if callable_lo <= p <= callable_hi:
                pos_counts[(p, ref, alt)] += n
                if has_indel:
                    on_indel[(p, ref, alt)] += n
    return pos_counts, on_indel, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crispresso-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="03_summary")
    ap.add_argument("--amplicon", default="Reference")
    ap.add_argument("--control-dir", default=None)
    ap.add_argument("--min-pct", type=float, default=MIN_PCT)
    ap.add_argument("--min-reads", type=int, default=MIN_READS)
    a = ap.parse_args()

    import analyze_alleles as A
    cfg = A.load_cfg(a.config)
    p = cfg.get("params", {})
    amp = cfg["_amp"]
    lo = p.get("exclude_bp_from_left", 15) + 1
    hi = len(amp) - p.get("exclude_bp_from_right", 15)

    counts, on_indel, N = collect(a.crispresso_dir, a.config, a.amplicon, lo, hi)
    ctrl, _, Nc = ({}, {}, 0)
    if a.control_dir:
        ctrl, _, Nc = collect(a.control_dir, a.config, a.amplicon, lo, hi)

    rows = []
    for (pos, ref, alt), k in sorted(counts.items(), key=lambda x: -x[1]):
        pct = 100 * k / N
        cpct = (100 * ctrl.get((pos, ref, alt), 0) / Nc) if Nc else None
        near_end = (pos < lo + END_UNRELIABLE) or (pos > hi - END_UNRELIABLE)
        meets = (pct >= a.min_pct and k >= a.min_reads)

        if cpct is not None and cpct > 0 and pct <= 2 * cpct:
            label = "baseline_polymorphism"
        elif near_end:
            label = "likely_technical_artifact"
        elif not meets:
            label = "below_reporting_threshold"
        else:
            label = "noticeable_SNP_candidate"

        rows.append({
            "amplicon_position_1based": pos,
            "ref>alt": "%s>%s" % (ref, alt),
            "edited_reads": k,
            "edited_pct": round(pct, 4),
            "control_pct": (round(cpct, 4) if cpct is not None else "no_control"),
            "on_indel_haplotype_reads": on_indel.get((pos, ref, alt), 0),
            "within_unreliable_end": near_end,
            "meets_reporting_threshold": meets,
            "classification": label,
        })

    df = pd.DataFrame(rows)
    os.makedirs(a.outdir, exist_ok=True)
    dst = os.path.join(a.outdir, "snps_not_editing.tsv")
    df.to_csv(dst, sep="\t", index=False)

    rep = df[df["classification"] == "noticeable_SNP_candidate"]
    print("analyzable reads N = %d | callable interval %d-%d" % (N, lo, hi))
    print("matched control    : %s" % (a.control_dir or "NONE"))
    print("distinct substitution positions observed : %d" % len(df))
    print("noticeable SNP candidates (>=%.1f%% and >=%d reads) : %d"
          % (a.min_pct, a.min_reads, len(rep)))
    if len(rep):
        print()
        print(rep.head(25).to_string(index=False))
    print()
    print("SUBSTITUTIONS ARE NOT COUNTED TOWARD EDITING EFFICIENCY.")
    print("wrote %s" % dst)


if __name__ == "__main__":
    main()
