#!/usr/bin/env python3
"""
Matched-control comparison (protocol section 7.4).

For each normalized indel signature, tests enrichment in the edited sample over
a matched untreated/mock control with a one-sided Fisher exact test, then applies
Benjamini-Hochberg FDR control across all tested alleles.

Also computes:
  - absolute excess              p_edited - p_control
  - normalized adjusted efficiency  max(0, (p_e - p_c) / (1 - p_c))
  - the control's upper 95% binomial bound as an empirical background limit

CAVEAT, enforced in the output: without UMIs, PCR-amplified reads are not
independent observations. These p-values are screening evidence. Biological
replicates are the unit of inference.

Importable as a module (used by analyze_alleles.py --control-dir) or runnable:
  python3 control_compare.py --edited-dir <dir> --control-dir <dir> \
      --config config.yaml --outdir 03_summary
"""
import argparse
import os

import pandas as pd

try:
    from scipy.stats import fisher_exact
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def bh_adjust(pvals):
    """Benjamini-Hochberg step-up FDR adjustment."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        val = min(prev, pvals[i] * n / k)
        adj[i] = val
        prev = val
    return adj


def wilson_upper(k, n, z=1.96):
    import math
    if n == 0:
        return float("nan")
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return min(1.0, (c + h) / d)


def _signature_counts(crispresso_dir, config, amplicon, window,
                      callable_lo, callable_hi):
    """Read a CRISPResso2 dir and return {signature: reads}, total_reads."""
    import analyze_alleles as A
    cfg = A.load_cfg(config)
    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]
    df = A.load_alleles(crispresso_dir)
    if "Reference_Name" in df.columns:
        df = df[df["Reference_Name"] == amplicon]
    t = A.build_table(df, cuts, names, window, A.PAIR_TOLERANCE,
                      callable_lo, callable_hi)
    n = int(t["reads"].sum())
    counts = (t[t["event_signature"] != ""]
              .groupby("event_signature")["reads"].sum().to_dict())
    return counts, n


def compare(edit_df, control_dir, config, amplicon, n_edited,
            window, callable_lo, callable_hi):
    """Returns a DataFrame of per-allele enrichment statistics."""
    ctrl, n_ctrl = _signature_counts(control_dir, config, amplicon, window,
                                     callable_lo, callable_hi)
    rows, pvals = [], []
    for _, r in edit_df.iterrows():
        sig = r["event_signature"]
        if not sig:
            continue
        k_e = int(r["reads"])
        k_c = int(ctrl.get(sig, 0))
        p_e = k_e / n_edited
        p_c = k_c / n_ctrl if n_ctrl else 0.0
        if HAVE_SCIPY:
            _, p = fisher_exact([[k_e, n_edited - k_e],
                                 [k_c, n_ctrl - k_c]], alternative="greater")
        else:
            p = float("nan")
        pvals.append(p)
        rows.append({
            "allele_id": r["allele_id"],
            "event_signature": sig,
            "attribution_class": r["attribution_class"],
            "edited_reads": k_e, "edited_pct": round(100 * p_e, 5),
            "control_reads": k_c, "control_pct": round(100 * p_c, 5),
            "control_upper95_pct": round(100 * wilson_upper(k_c, n_ctrl), 5),
            "absolute_excess_pct": round(100 * (p_e - p_c), 5),
            "normalized_adjusted_pct": round(
                100 * max(0.0, (p_e - p_c) / (1 - p_c)) if p_c < 1 else 0.0, 5),
            "fisher_p_onesided": p,
        })
    out = pd.DataFrame(rows)
    if len(out):
        out["fdr_bh"] = bh_adjust(pvals)
        out["above_control_background"] = (
            out["edited_pct"] > out["control_upper95_pct"])
        out["evidence"] = out.apply(
            lambda r: "enriched (FDR<0.05, above background)"
            if (r["fdr_bh"] < 0.05 and r["above_control_background"])
            else "not distinguishable from background", axis=1)
        out = out.sort_values("edited_reads", ascending=False)
    out.attrs["n_edited"] = n_edited
    out.attrs["n_control"] = n_ctrl
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edited-dir", required=True)
    ap.add_argument("--control-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default="03_summary")
    ap.add_argument("--amplicon", default="Reference")
    ap.add_argument("--window", type=int, default=5)
    a = ap.parse_args()

    if not HAVE_SCIPY:
        print("WARNING: scipy not installed; p-values will be NaN. "
              "pip install scipy")

    import analyze_alleles as A
    cfg = A.load_cfg(a.config)
    p = cfg.get("params", {})
    amp = cfg["_amp"]
    lo = p.get("exclude_bp_from_left", 15) + 1
    hi = len(amp) - p.get("exclude_bp_from_right", 15)

    cuts = [g["cut_site_1based"] for g in cfg["guides"]]
    names = [g["name"] for g in cfg["guides"]]
    df = A.load_alleles(a.edited_dir)
    if "Reference_Name" in df.columns:
        df = df[df["Reference_Name"] == a.amplicon]
    t = A.build_table(df, cuts, names, a.window, A.PAIR_TOLERANCE, lo, hi)
    n_e = int(t["reads"].sum())
    edit = (t.groupby(["allele_id", "event_signature", "attribution_class"])
             .agg(reads=("reads", "sum")).reset_index())

    out = compare(edit, a.control_dir, a.config, a.amplicon, n_e, a.window, lo, hi)
    os.makedirs(a.outdir, exist_ok=True)
    dst = os.path.join(a.outdir, "control_comparison.tsv")
    out.to_csv(dst, sep="\t", index=False)
    print("edited N = %d | control N = %d" % (n_e, out.attrs["n_control"]))
    print(out.head(20).to_string(index=False))
    print("\nwrote %s" % dst)
    print("NOTE: without UMIs these are screening statistics, not independent "
          "molecular observations. Use biological replicates for inference.")


if __name__ == "__main__":
    main()
