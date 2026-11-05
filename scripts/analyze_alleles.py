#!/usr/bin/env python3
"""
Primary allele analysis (framework v2.0).

Supersedes summarize_alleles.py. Implements the protocol requirements:
  - substitutions excluded from the editing numerator
  - alleles collapsed by normalized indel signature; haplotypes kept separately
  - stable event-signature naming (D7:123-129, I3:123_124:ACT, PAIR_D137:g1-g2:74-210)
  - guide attribution: unique / shared / paired-guide excision / distal
  - Wilson 95% confidence intervals on every frequency
  - evidence tiers 1/2/3 with predeclared support thresholds
  - exact-WT, window-unmodified, and SNP-only reported separately

Coordinate convention: CRISPResso2 0-based internally; ALL OUTPUT IS 1-BASED.

Usage:
  python3 analyze_alleles.py --crispresso-dir <dir> --config <config.yaml> \
      --outdir 03_summary [--control-dir <control CRISPResso dir>]
"""
import argparse
import json
import math
import os
import re
import zipfile
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Predeclared thresholds (protocol section 7.2). Declared before inspecting data.
MAIN_MIN_READS = 10
MAIN_MIN_PCT = 0.1
APPENDIX_MIN_READS = 3
CORE_WINDOW = 5        # +/- bp around each cut
DISTAL_WINDOW = 20     # +/- bp sensitivity region
PAIR_TOLERANCE = 20    # +/- bp endpoint tolerance for paired excision


# ---------------------------------------------------------------- statistics
def wilson(k, n, z=1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def wilson_str(k, n):
    lo, hi = wilson(k, n)
    return "%.4f-%.4f" % (100 * lo, 100 * hi)


# ------------------------------------------------------------- indel parsing
def parse_events(aligned, reference):
    """
    Walk a CRISPResso2 pairwise alignment and extract indel events in
    1-BASED reference coordinates.

    Returns (deletions, insertions, substitutions) where
      deletions     = [(start1, end1)]                inclusive reference interval
      insertions    = [(left1, right1, seq)]          between two reference bases
      substitutions = [(pos1, ref_base, alt_base)]
    """
    dels, ins, subs = [], [], []
    refpos = 0            # 1-based position of the last consumed reference base
    i, n = 0, len(reference)
    while i < n:
        r, a = reference[i], aligned[i]
        if r == "-":                                   # insertion into reference
            j = i
            while j < n and reference[j] == "-":
                j += 1
            seq = aligned[i:j].replace("-", "")
            if seq:
                ins.append((refpos, refpos + 1, seq))
            i = j
            continue
        refpos += 1
        if a == "-":                                   # deletion from reference
            j = i
            while j < n and aligned[j] == "-" and reference[j] != "-":
                j += 1
            length = j - i
            dels.append((refpos, refpos + length - 1))
            refpos += length - 1
            i = j
            continue
        if a != r and a != "N" and r != "N":
            subs.append((refpos, r, a))
        i += 1
    return dels, ins, subs


def signature(dels, ins):
    """Stable, sortable event signature. Empty string when no indel."""
    parts = []
    for s, e in sorted(dels):
        parts.append("D%d:%d-%d" % (e - s + 1, s, e))
    for l, r, seq in sorted(ins):
        parts.append("I%d:%d_%d:%s" % (len(seq), l, r, seq))
    return "|".join(parts)


# ------------------------------------------------------------- attribution
def attribute(dels, ins, cuts, names, window, pair_tol):
    """
    Return (guides_hit, class_label, pair_label, guides_endpoint_anchored).

    Cuts are 1-based nick positions: the nick lies between cut-1 and cut.

    Two attribution strengths are computed:
      SPAN     - the protocol's rule (section 6.2): a deletion qualifies when its
                 interval intersects the guide window.
      ENDPOINT - stricter: a deletion *breakpoint* lies within the window.

    A large deletion arising at one cut can sweep across a second, distant cut
    and be credited to it under the SPAN rule alone. Both are reported so that
    inflation of a weak guide's rate by unrelated large deletions is visible.
    """
    hit, anchored = set(), set()
    for gi, c in enumerate(cuts):
        lo, hi = c - window, c + window - 1
        for s, e in dels:
            if s <= hi and e >= lo:                       # SPAN
                hit.add(names[gi])
            if (lo <= s <= hi) or (lo <= e + 1 <= hi):     # ENDPOINT
                anchored.add(names[gi])
        for l, r, _ in ins:
            if lo <= l <= hi or lo <= r <= hi:
                hit.add(names[gi])
                anchored.add(names[gi])

    # paired-guide excision: one deletion spanning between two intended cuts
    pair = None
    if len(cuts) >= 2:
        for a in range(len(cuts)):
            for b in range(a + 1, len(cuts)):
                c1, c2 = sorted((cuts[a], cuts[b]))
                for s, e in dels:
                    if (abs(s - c1) <= pair_tol and abs(e - (c2 - 1)) <= pair_tol
                            and (e - s + 1) >= 0.5 * (c2 - c1)):
                        pair = "PAIR_D%d:%s-%s:%d-%d" % (
                            e - s + 1, names[a], names[b], s, e)
                        hit.update([names[a], names[b]])
                        anchored.update([names[a], names[b]])
    if pair:
        return hit, "paired_guide_excision", pair, anchored
    if not hit:
        return hit, "distal_or_background_indel", None, anchored
    if len(hit) == 1:
        return hit, "guide_unique", None, anchored
    return hit, "guide_shared", None, anchored


# ------------------------------------------------------------------- loaders
def load_alleles(d):
    z = os.path.join(d, "Alleles_frequency_table.zip")
    if os.path.exists(z):
        with zipfile.ZipFile(z) as zf:
            with zf.open(zf.namelist()[0]) as fh:
                return pd.read_csv(fh, sep="\t")
    p = os.path.join(d, "Alleles_frequency_table.txt")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t")
    raise FileNotFoundError("no allele frequency table in %s" % d)


def load_cfg(path):
    import yaml
    c = yaml.safe_load(open(path))
    c["_amp"] = "".join(c["amplicon"]["sequence"].upper().split())
    return c


def crispresso_version(d):
    p = os.path.join(d, "CRISPResso2_info.json")
    if os.path.exists(p):
        try:
            j = json.load(open(p))
            return j.get("running_info", {}).get("version", "unknown")
        except Exception:
            pass
    return "unknown"


# ---------------------------------------------------------------------- main
def build_table(df, cuts, names, window, pair_tol, callable_lo, callable_hi):
    """Per-haplotype event extraction and classification."""
    rows = []
    for _, r in df.iterrows():
        dels, ins, subs = parse_events(r["Aligned_Sequence"], r["Reference_Sequence"])
        # restrict to the callable interval (primer-derived bases excluded)
        dels = [(s, e) for s, e in dels if e >= callable_lo and s <= callable_hi]
        ins = [(l, rr, sq) for l, rr, sq in ins if callable_lo <= l <= callable_hi]
        subs = [(p, a, b) for p, a, b in subs if callable_lo <= p <= callable_hi]

        sig = signature(dels, ins)
        has_indel = bool(sig)
        exact_wt = (r["Aligned_Sequence"] == r["Reference_Sequence"])

        if not has_indel:
            cls = "exact_reference_WT" if exact_wt else (
                "SNP_only_unmodified" if subs else "window_unmodified_other")
            aid = {"exact_reference_WT": "WT",
                   "SNP_only_unmodified": "SNP_ONLY",
                   "window_unmodified_other": "UNMOD_OTHER"}[cls]
            guides, attr, pair, anchored = set(), "none", None, set()
        else:
            guides, attr, pair, anchored = attribute(
                dels, ins, cuts, names, window, pair_tol)
            cls = ("paired_guide_excision" if attr == "paired_guide_excision"
                   else "CRISPR_associated_indel" if guides
                   else "distal_or_background_indel")
            aid = pair or sig

        rows.append({
            "allele_id": aid,
            "event_signature": sig,
            "classification": cls,
            "guide_attribution": ",".join(sorted(guides)) if guides else "",
            "guide_endpoint_anchored": ",".join(sorted(anchored)) if anchored else "",
            "attribution_class": attr,
            "n_del_events": len(dels), "n_ins_events": len(ins),
            "bp_deleted": sum(e - s + 1 for s, e in dels),
            "bp_inserted": sum(len(sq) for _, _, sq in ins),
            "n_substitutions": len(subs),
            "substitutions": ";".join("%d%s>%s" % (p, a, b) for p, a, b in subs[:20]),
            "exact_wt": exact_wt,
            "reads": int(r["#Reads"]),
            "aligned_sequence": r["Aligned_Sequence"],
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crispresso-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="03_summary")
    ap.add_argument("--amplicon", default="Reference")
    ap.add_argument("--control-dir", default=None,
                    help="matched control CRISPResso2 dir; enables background correction")
    ap.add_argument("--window", type=int, default=CORE_WINDOW)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    cfg = load_cfg(a.config)
    amp = cfg["_amp"]
    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]
    p = cfg.get("params", {})
    excl_l = p.get("exclude_bp_from_left", 15)
    excl_r = p.get("exclude_bp_from_right", 15)
    callable_lo, callable_hi = excl_l + 1, len(amp) - excl_r

    df = load_alleles(a.crispresso_dir)
    if "Reference_Name" in df.columns:
        df = df[df["Reference_Name"] == a.amplicon]
    tbl = build_table(df, cuts, names, a.window, PAIR_TOLERANCE,
                      callable_lo, callable_hi)
    N = int(tbl["reads"].sum())      # analyzable reads for this amplicon

    # ---------------- collapse haplotypes into edit alleles ----------------
    edit = (tbl.groupby(["allele_id", "event_signature", "classification",
                         "guide_attribution", "guide_endpoint_anchored",
                         "attribution_class", "bp_deleted", "bp_inserted"],
                        dropna=False)
               .agg(reads=("reads", "sum"),
                    n_haplotypes=("reads", "size"))
               .reset_index())
    edit["pct_analyzable"] = (100 * edit["reads"] / N).round(4)
    edit["ci95_pct"] = [wilson_str(k, N) for k in edit["reads"]]
    edit["indel_size"] = edit["bp_inserted"] - edit["bp_deleted"]

    def tier(r):
        if r["classification"] in ("exact_reference_WT", "SNP_only_unmodified",
                                   "window_unmodified_other"):
            return "reference"
        if r["classification"] == "distal_or_background_indel":
            return "distal"
        meets_main = (r["reads"] >= MAIN_MIN_READS
                      and r["pct_analyzable"] >= MAIN_MIN_PCT)
        if not meets_main:
            return "tier3_putative" if r["reads"] < APPENDIX_MIN_READS else "tier2_supported"
        # Tier 1 requires control or replicate corroboration
        return "tier2_supported" if a.control_dir is None else "tier1_candidate"

    edit["evidence_tier"] = edit.apply(tier, axis=1)
    edit = edit.sort_values("reads", ascending=False)

    # ---------------- sample-level metrics ---------------------------------
    def frac(mask):
        k = int(tbl.loc[mask, "reads"].sum())
        return k, 100 * k / N, wilson_str(k, N)

    is_indel = tbl["event_signature"] != ""
    core = tbl["attribution_class"].isin(
        ["guide_unique", "guide_shared", "paired_guide_excision"])

    metrics = []
    for label, mask in [
        ("exact_reference_WT", tbl["classification"] == "exact_reference_WT"),
        ("SNP_only", tbl["classification"] == "SNP_only_unmodified"),
        ("window_unmodified (no qualifying indel)", ~core),
        ("raw_core_indel_efficiency", core),
        ("paired_guide_excision", tbl["attribution_class"] == "paired_guide_excision"),
        ("distal_or_background_indel",
         is_indel & (tbl["attribution_class"] == "distal_or_background_indel")),
        ("any_indel_anywhere", is_indel),
    ]:
        k, pct, ci = frac(mask)
        metrics.append({"metric": label, "reads": k, "pct_of_analyzable": round(pct, 4),
                        "ci95_pct": ci, "denominator_N": N})

    wt_reads = int(tbl.loc[tbl["classification"] == "exact_reference_WT", "reads"].sum())
    mut_reads = int(tbl.loc[core, "reads"].sum())
    metrics.append({
        "metric": "mutant_to_exactWT_ratio",
        "reads": mut_reads,
        "pct_of_analyzable": (round(mut_reads / wt_reads, 6) if wt_reads
                              else "undefined"),
        "ci95_pct": "", "denominator_N": wt_reads})
    met = pd.DataFrame(metrics)

    # ---------------- guide-level metrics ----------------------------------
    def _has(col, g):
        return tbl[col].fillna("").apply(lambda s: g in s.split(",") if s else False)

    grows = []
    for g in names:
        hit = _has("guide_attribution", g)                # SPAN rule (protocol)
        anch = _has("guide_endpoint_anchored", g)         # ENDPOINT rule (stricter)
        uniq = hit & (tbl["attribution_class"] == "guide_unique")
        shar = hit & (tbl["attribution_class"].isin(["guide_shared",
                                                     "paired_guide_excision"]))
        k = int(tbl.loc[hit, "reads"].sum())
        ka = int(tbl.loc[anch, "reads"].sum())
        c = cuts[names.index(g)]
        grows.append({
            "guide": g,
            "cut_site_1based": c,
            "core_window_1based": "%d-%d" % (c - a.window, c + a.window - 1),
            "marginal_span_reads": k,
            "marginal_span_pct": round(100 * k / N, 4),
            "marginal_span_ci95": wilson_str(k, N),
            "marginal_endpoint_reads": ka,
            "marginal_endpoint_pct": round(100 * ka / N, 4),
            "marginal_endpoint_ci95": wilson_str(ka, N),
            "span_only_inflation_pct": round(100 * (k - ka) / N, 4),
            "unique_pct": round(100 * int(tbl.loc[uniq, "reads"].sum()) / N, 4),
            "shared_pct": round(100 * int(tbl.loc[shar, "reads"].sum()) / N, 4),
            "insertion_pct": round(100 * int(
                tbl.loc[hit & (tbl["bp_inserted"] > 0), "reads"].sum()) / N, 4),
            "deletion_pct": round(100 * int(
                tbl.loc[hit & (tbl["bp_deleted"] > 0), "reads"].sum()) / N, 4),
        })
    gd = pd.DataFrame(grows)

    def band(x):
        return ("Low" if x < 5 else "Modest" if x < 20
                else "Moderate" if x < 50 else "High")
    # Descriptive band uses the stricter endpoint-anchored rate.
    gd["descriptive_band"] = gd["marginal_endpoint_pct"].apply(band)
    gd["note"] = gd["marginal_endpoint_pct"].apply(
        lambda x: "below 1% - consider redesign" if x < 1 else "")

    # ---------------- window sensitivity analysis ---------------------------
    # Protocol 5.4: report a +/-1 bp estimate for direct comparability with the
    # CRISPResso2 default, and a wider window to show threshold dependence.
    sens = []
    for w in (1, 5, 20):
        t = build_table(df, cuts, names, w, PAIR_TOLERANCE, callable_lo, callable_hi)
        m = t["attribution_class"].isin(
            ["guide_unique", "guide_shared", "paired_guide_excision"])
        k = int(t.loc[m, "reads"].sum())
        row = {"core_window_bp": "+/-%d" % w, "reads": k,
               "pct_of_analyzable": round(100 * k / N, 4),
               "ci95_pct": wilson_str(k, N)}
        for g in names:
            hm = t["guide_attribution"].fillna("").apply(
                lambda s: g in s.split(",") if s else False)
            row[g + "_pct"] = round(100 * int(t.loc[hm, "reads"].sum()) / N, 4)
        sens.append(row)
    sens_df = pd.DataFrame(sens)
    sens_df.to_csv(os.path.join(a.outdir, "window_sensitivity.tsv"),
                   sep="\t", index=False)

    # ---------------- optional control comparison ---------------------------
    if a.control_dir:
        from control_compare import compare
        cmp_df = compare(edit, a.control_dir, a.config, a.amplicon, N,
                         a.window, callable_lo, callable_hi)
        cmp_df.to_csv(os.path.join(a.outdir, "control_comparison.tsv"),
                      sep="\t", index=False)

    # ---------------- write deliverables ------------------------------------
    o = a.outdir
    main_tbl = edit[(edit["reads"] >= MAIN_MIN_READS)
                    & (edit["pct_analyzable"] >= MAIN_MIN_PCT)]
    lowfreq = edit[((edit["reads"] < MAIN_MIN_READS)
                    | (edit["pct_analyzable"] < MAIN_MIN_PCT))
                   & (edit["reads"] >= APPENDIX_MIN_READS)]

    main_tbl.to_csv(os.path.join(o, "alleles_indel_collapsed.tsv"), sep="\t", index=False)
    edit.to_csv(os.path.join(o, "alleles_indel_collapsed_all.tsv"), sep="\t", index=False)
    lowfreq.to_csv(os.path.join(o, "low_frequency_candidates.tsv"), sep="\t", index=False)
    tbl.drop(columns=["aligned_sequence"]).to_csv(
        os.path.join(o, "alleles_full_haplotype.tsv"), sep="\t", index=False)
    met.to_csv(os.path.join(o, "sample_metrics.tsv"), sep="\t", index=False)
    gd.to_csv(os.path.join(o, "guide_efficiency.tsv"), sep="\t", index=False)

    # ---------------- figure: indel spectrum by attribution -----------------
    ind = edit[edit["event_signature"] != ""]
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for cls, col in [("guide_unique", "#2c7fb8"), ("guide_shared", "#7fcdbb"),
                     ("paired_guide_excision", "#e67e22"),
                     ("distal_or_background_indel", "#bdc3c7")]:
        s = ind[ind["attribution_class"] == cls]
        if len(s):
            ax.scatter(s["indel_size"], s["reads"], s=18, alpha=.75,
                       color=col, label=cls)
    ax.set_yscale("log")
    ax.set_xlabel("indel size (bp; negative = deletion)")
    ax.set_ylabel("reads (log)")
    ax.set_title("Edit alleles by guide attribution")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(o, "fig4_alleles_by_attribution.png"), dpi=200)
    plt.close(fig)

    # ---------------- console summary ---------------------------------------
    print("CRISPResso2 version : %s" % crispresso_version(a.crispresso_dir))
    print("coordinates         : 1-based, inclusive")
    print("callable interval   : %d-%d (excluded %d left / %d right)"
          % (callable_lo, callable_hi, excl_l, excl_r))
    print("analyzable reads N  : %d" % N)
    print("matched control     : %s" % (a.control_dir or "NONE - results are PUTATIVE"))
    print()
    print(met.to_string(index=False))
    print()
    print(gd.to_string(index=False))
    print()
    print("window sensitivity (span rule):")
    print(sens_df.to_string(index=False))
    print()
    n_edit_main = int((main_tbl["event_signature"] != "").sum())
    n_edit_all = int((edit["event_signature"] != "").sum())
    print("distinct edit alleles observed : %d" % n_edit_all)
    print("  meeting main-report threshold (>=%d reads AND >=%.2f%% of N=%d, "
          "i.e. >=%d reads) : %d"
          % (MAIN_MIN_READS, MAIN_MIN_PCT, N,
             math.ceil(MAIN_MIN_PCT * N / 100), n_edit_main))
    print("  low-frequency candidates (>=%d reads) : %d"
          % (APPENDIX_MIN_READS, len(lowfreq)))
    if n_edit_main == 0 and n_edit_all > 0:
        print("  NOTE: no single edit allele reaches the main-report threshold.")
        print("        Editing is distributed across many low-frequency alleles")
        print("        rather than concentrated in a dominant repair outcome.")
    print()
    print("wrote 7 tables + 1 figure to %s" % o)


if __name__ == "__main__":
    main()
