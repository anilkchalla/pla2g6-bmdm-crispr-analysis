#!/usr/bin/env python3
"""Stage 3 - reclassify vendor per-allele tables across BOTH cut sites.

Why this is needed: vendor pipelines anchor each sample to a single NHEJ group and
score editing inside one narrow quantification window. An allele edited only at the
*other* cut site does not touch that window and is not counted, so the headline
edit_rate is scoped to one site and the two samples of a dual-guide experiment are
not reporting the same quantity.

This script re-parses the CIGAR of every reported allele, maps each indel event to
its reference interval, and assigns it to cut 1, cut 2, or both. It also classifies
each allele by reading frame over the CDS.

A caveat the output cannot repair: because the vendor clusters alleles around one
cut, single-site alleles at the other cut may be absent from the table entirely,
producing mirror-image zeros between samples. Per-guide single-site efficiencies
should not be read off these tables; use the anchor recount instead.

Usage:
    python3 y452b3_reclassify_alleles.py --edit-analysis data/Y452B3_edit-analysis
"""

import argparse
import glob
import os
import re

import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import CDS, CUT1, CUT2

# Vendor "clustered core" start for each sample id, read from the order summary.
CORE_START_COL = "core_start"


def parse_cigar(cigar, ref_start):
    """Return [(op, ref_start, ref_end, length)] for I/D events, 1-based inclusive."""
    pos = ref_start
    events = []
    for n, op in re.findall(r"(\d+)([MIDN])", cigar):
        n = int(n)
        if op == "M":
            pos += n
        elif op in "DN":
            events.append(("D", pos, pos + n - 1, n))
            pos += n
        elif op == "I":
            events.append(("I", pos, pos - 1, n))
    return events


def touches(events, cut, slack=1):
    """True if any indel disrupts the blunt cut between `cut` and `cut`+1."""
    for op, s, e, n in events:
        if op == "D" and not (e < cut - slack or s > cut + slack):
            return True
        if op == "I" and abs(s - (cut + 1)) <= 3:
            return True
    return False


def frame_effect(events):
    """Net coding-sequence length change, restricted to the CDS."""
    net = 0
    for op, s, e, n in events:
        if op == "D":
            overlap = max(0, min(e, CDS[1]) - max(s, CDS[0]) + 1)
            net -= overlap
        elif CDS[0] <= s <= CDS[1]:
            net += n
    return net


def read_core_starts(edit_dir):
    """Map sample id -> clustered-core start from the order summary."""
    summary = glob.glob(os.path.join(edit_dir, "*order_summary.tsv"))
    if not summary:
        raise FileNotFoundError("no *order_summary.tsv in " + edit_dir)
    cores = {}
    with open(summary[0]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_id, idx_core = header.index("osid"), header.index(CORE_START_COL)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            cores[p[idx_id]] = int(p[idx_core])
    return cores


def analyse(path, core_start, excision_tol=2):
    """`excision_tol` is the end-processing slack allowed around each cut site
    when deciding whether one contiguous deletion is the programmed excision.
    Observed deletion boundaries cluster at 74 and 210 with 71-75 and 208-211
    shoulders, and the excision fraction plateaus from tol=2 upward, so 2 bp is
    the smallest value that captures the real spread rather than a tuned one."""
    counts = {"unedited": 0, "cut1_only": 0, "cut2_only": 0, "dual": 0, "low_support": 0}
    frame = {"frameshift": 0, "in_frame": 0}
    excision_sizes = {}
    total = 0

    for line in open(path):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "Molecule_ID":
            continue
        n_reads = int(parts[6])
        cigar = parts[13].strip()
        total += n_reads
        if "very low support" in line:
            counts["low_support"] += n_reads
            continue

        events = parse_cigar(cigar, core_start)
        if not events:
            counts["unedited"] += n_reads
            continue

        # The excision test is applied to every allele, not only to those the
        # touches() test calls "dual". A deletion starting at 75 with the same
        # right boundary is still the programmed excision with 1 bp of extra
        # resection, but falls outside touches()'s 1 bp slack at cut 1.
        for op, s, e, n in events:
            if op == "D" and s <= CUT1 + 1 + excision_tol and e >= CUT2 - excision_tol:
                excision_sizes[n] = excision_sizes.get(n, 0) + n_reads
                break

        h1, h2 = touches(events, CUT1), touches(events, CUT2)
        if h1 and h2:
            counts["dual"] += n_reads
        elif h1:
            counts["cut1_only"] += n_reads
        elif h2:
            counts["cut2_only"] += n_reads
        else:
            counts["unedited"] += n_reads

        net = frame_effect(events)
        if net == 0:
            pass
        elif net % 3 == 0:
            frame["in_frame"] += n_reads
        else:
            frame["frameshift"] += n_reads

    edited = counts["cut1_only"] + counts["cut2_only"] + counts["dual"]
    excision = sum(v for k, v in excision_sizes.items())
    return {
        "total_reads": total,
        "counts": counts,
        "edited_pct": round(100 * edited / total, 2),
        "excision_pct": round(100 * excision / total, 2),
        "frameshift_pct": round(100 * frame["frameshift"] / total, 2),
        "in_frame_pct": round(100 * frame["in_frame"] / total, 2),
        "excision_sizes": dict(sorted(excision_sizes.items(), key=lambda x: -x[1])[:8]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edit-analysis", required=True,
                    help="directory holding *_editing_results.tsv and *order_summary.tsv")
    ap.add_argument("--excision-tolerance", type=int, default=2,
                    help="end-processing slack (bp) allowed around each cut site (default 2)")
    args = ap.parse_args()

    cores = read_core_starts(args.edit_analysis)
    files = sorted(glob.glob(os.path.join(args.edit_analysis, "*_editing_results.tsv")))
    if not files:
        raise SystemExit("no *_editing_results.tsv found")

    print(f"{'sample':<14}{'reads':>8}{'edited':>9}{'excision':>10}"
          f"{'frameshift':>12}{'in-frame':>10}")
    print("-" * 63)
    for path in files:
        base = os.path.basename(path)
        sid = base.replace("_editing_results.tsv", "")
        core = cores.get(sid)
        if core is None:
            print(f"  ! no core_start for {sid}, skipping")
            continue
        r = analyse(path, core, args.excision_tolerance)
        print(f"{sid:<14}{r['total_reads']:>8}{r['edited_pct']:>8}%"
              f"{r['excision_pct']:>9}%{r['frameshift_pct']:>11}%{r['in_frame_pct']:>9}%")

    print("\nNote: cut1_only / cut2_only counts are NOT comparable between samples; "
          "see the module docstring.")


if __name__ == "__main__":
    main()
