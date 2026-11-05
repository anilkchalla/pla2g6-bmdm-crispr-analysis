#!/usr/bin/env python3
"""
Read-accounting ledger and QC gate for a CRISPResso2 run.

Reconciles input reads -> preprocessed -> aligned -> quantified, and applies
the gate thresholds from ANALYSIS_FRAMEWORK.md section 2.

Usage:
  python3 qc_gate.py --crispresso-dir 02_crispresso/CRISPResso_on_<sample>
"""
import argparse
import json
import os
import sys

GATES = {
    "preprocessing_retention_pct": 80.0,   # G2
    "alignment_rate_pct": 90.0,            # G3: aligned / preprocessed
    "min_aligned_reads": 10000,            # G4
}


def read_mapping_stats(d):
    p = os.path.join(d, "CRISPResso_mapping_statistics.txt")
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        vals = fh.readline().rstrip("\n").split("\t")
    return {h: int(v) for h, v in zip(hdr, vals)}


def read_quant(d):
    p = os.path.join(d, "CRISPResso_quantification_of_editing_frequency.txt")
    rows = []
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def read_fastp(d):
    p = os.path.join(d, "fastp_report.json")
    if not os.path.exists(p):
        return None
    j = json.load(open(p))
    return {
        "before": j.get("summary", {}).get("before_filtering", {}),
        "after": j.get("summary", {}).get("after_filtering", {}),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crispresso-dir", required=True)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any gate fails")
    a = ap.parse_args()
    d = a.crispresso_dir

    ms = read_mapping_stats(d)
    quant = read_quant(d)
    fp = read_fastp(d)

    n_in = ms["READS IN INPUTS"]
    n_pre = ms["READS AFTER PREPROCESSING"]
    n_aln = ms["READS ALIGNED"]

    print("=" * 70)
    print("READ ACCOUNTING LEDGER")
    print("=" * 70)
    print("%-32s %12s %10s %10s" % ("stage", "reads", "% of in", "% of prev"))
    prev = n_in
    for label, n in [("reads in input", n_in),
                     ("after preprocessing (fastp)", n_pre),
                     ("aligned to a reference", n_aln)]:
        print("%-32s %12d %9.2f%% %9.2f%%"
              % (label, n, 100 * n / n_in, 100 * n / prev if prev else 0))
        prev = n
    print()
    print("%-32s %12d" % ("lost in preprocessing", n_in - n_pre))
    print("%-32s %12d" % ("preprocessed but unaligned", n_pre - n_aln))
    if fp and fp["before"]:
        tr = fp["before"].get("total_reads")
        if tr:
            print("%-32s %12d" % ("fastp saw (pairs x2)", tr))

    print()
    print("=" * 70)
    print("PER-AMPLICON QUANTIFICATION")
    print("=" * 70)
    print("%-14s %10s %12s %12s %10s" %
          ("amplicon", "reads", "unmod %", "mod %", "% of aln"))
    total_all = None
    for r in quant:
        total_all = int(r["Reads_aligned_all_amplicons"])
        n = int(r["Reads_aligned"])
        print("%-14s %10d %11.4f%% %11.4f%% %9.4f%%"
              % (r["Amplicon"], n, float(r["Unmodified%"]),
                 float(r["Modified%"]), 100 * n / total_all))

    print()
    print("=" * 70)
    print("QC GATES")
    print("=" * 70)
    results = []
    retention = 100 * n_pre / n_in if n_in else 0
    aln_rate = 100 * n_aln / n_pre if n_pre else 0
    for name, val, thr, op in [
        ("G2 preprocessing retention", retention, GATES["preprocessing_retention_pct"], ">="),
        ("G3 alignment rate", aln_rate, GATES["alignment_rate_pct"], ">="),
        ("G4 aligned read depth", n_aln, GATES["min_aligned_reads"], ">="),
    ]:
        ok = val >= thr
        results.append(ok)
        unit = "%" if "rate" in name or "retention" in name else ""
        print("  [%s] %-30s %10.2f%s  (need %s %g%s)"
              % ("PASS" if ok else "FAIL", name, val, unit, op, thr, unit))

    print()
    if all(results):
        print("STATUS: all gates passed")
    else:
        print("STATUS: %d gate(s) FAILED -- interpret results with caution"
              % results.count(False))
        if a.strict:
            sys.exit(1)


if __name__ == "__main__":
    main()
