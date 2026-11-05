#!/usr/bin/env python3
"""
Put CRISPResso2's Pass-A numbers next to ont_guide_panel.py's, per sample.

Framework §16.2 calls the CRISPResso2 run a cross-check, not the primary result
for ONT data. This script makes the check readable: two independent engines,
one number each, and an explicit concordance verdict.

The two are NOT expected to agree exactly, and the reasons are enumerated in the
output so a disagreement can be diagnosed rather than argued about:

  * CRISPResso2 aligns globally with Needleman-Wunsch; ont_guide_panel.py chains
    unique k-mer anchors. On reads with a large deletion the two place the
    junction by different criteria.
  * CRISPResso2 counts a read as modified if ANY indel intersects the
    quantification window (a SPAN rule); the panel analyser requires a
    BREAKPOINT inside the window (an ENDPOINT rule, framework §6.2). SPAN is
    the looser of the two, so CRISPResso2 should read equal or higher.
  * ont_guide_panel.py drops indels below 3 bp and 1-2 bp indels adjacent to
    homopolymer runs (§16.3); CRISPResso2 counts every indel in the window.
    This pushes CRISPResso2 higher again on ONT data.
  * Denominators differ: CRISPResso2 uses reads aligned to the reference;
    the panel analyser uses reads with >=5 unique k-mer anchors.

Expect CRISPResso2 >= panel, with the gap widening on the noisiest libraries.
A LOWER CRISPResso2 number, or a gap above ~15 points on a clean library, is
worth investigating — the usual cause is large-deletion reads failing to align.

Usage:
  python3 compare_crispresso.py --config config.yaml \
      --crispresso-dir 02_crispresso --summary-dir 03_summary
"""
import argparse
import os
import zipfile

import pandas as pd
import yaml


def read_quant(d):
    """Read CRISPResso2's quantification summary from a run directory."""
    for name in ("CRISPResso_quantification_of_editing_frequency.txt",
                 "Quantification_of_editing_frequency.txt"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return pd.read_csv(p, sep="\t")
    z = os.path.join(d, "CRISPResso2_info.zip")
    if os.path.exists(z):
        with zipfile.ZipFile(z) as zf:
            for n in zf.namelist():
                if "quantification_of_editing_frequency" in n.lower():
                    return pd.read_csv(zf.open(n), sep="\t")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--crispresso-dir", default="02_crispresso")
    ap.add_argument("--summary-dir", default="03_summary")
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config))
    panel = pd.read_csv(os.path.join(a.summary_dir, "guide_efficiency.tsv"),
                        sep="\t").set_index("sample")

    rows = []
    for s in cfg["samples"]:
        d = os.path.join(a.crispresso_dir, "CRISPResso_on_%s_indel" % s["name"])
        q = read_quant(d)
        if q is None:
            rows.append({"sample": s["name"], "guide": s["guide"],
                         "crispresso_status": "MISSING (%s)" % d})
            continue
        r = q[q["Amplicon"] == cfg["amplicon"]["name"]] if "Amplicon" in q else q
        r = r.iloc[0]
        aligned = float(r.get("Reads_aligned", float("nan")))
        modified = float(r.get("Modified%", float("nan")))
        unmod = float(r.get("Unmodified%", float("nan")))
        pr = panel.loc[s["name"]]
        gap = modified - pr["core_indel_pct"]
        rows.append({
            "sample": s["name"], "guide": s["guide"],
            "crispresso_reads_aligned": int(aligned) if aligned == aligned else None,
            "crispresso_modified_pct": round(modified, 3),
            "crispresso_unmodified_pct": round(unmod, 3),
            "panel_analyzable_reads": int(pr["analyzable_reads"]),
            "panel_core_indel_pct": round(float(pr["core_indel_pct"]), 3),
            "panel_any_indel_pct": round(float(pr["any_indel_pct"]), 3),
            "difference_pts": round(gap, 3),
            "direction_expected": "yes" if gap >= -1.0 else "NO - investigate",
            "crispresso_status": "ok",
        })

    df = pd.DataFrame(rows)
    out = os.path.join(a.summary_dir, "crispresso_concordance.tsv")
    df.to_csv(out, sep="\t", index=False)

    print("=" * 78)
    print("CRISPResso2 cross-check  —  %s" % cfg.get("project"))
    print("=" * 78)
    print(df.to_string(index=False))
    print()
    if "difference_pts" in df and df["difference_pts"].notna().any():
        ok = df[df["direction_expected"] == "yes"]
        print("CRISPResso2 >= panel in %d/%d samples (expected: SPAN rule and no "
              "ONT indel filtering both push CRISPResso2 higher)."
              % (len(ok), len(df)))
        print("Median gap: %.1f points." % df["difference_pts"].median())
        big = df[df["difference_pts"].abs() > 15]
        if len(big):
            print("\nGaps above 15 points — inspect the alignment of large-deletion "
                  "reads in these runs:")
            print(big[["sample", "guide", "crispresso_modified_pct",
                       "panel_core_indel_pct", "difference_pts"]].to_string(index=False))
    print("\nRanking check: the cross-check passes if the ORDER of guides is the "
          "same under both engines. Absolute values are not expected to match.")
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
