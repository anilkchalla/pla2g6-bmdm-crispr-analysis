#!/usr/bin/env python3
"""Stage 2 - alignment-free allele calling from raw reads.

Three anchors are used: one left of cut 1, one between the two cuts, one right of
cut 2. The distance between consecutive anchors gives the net indel at each cut
site independently. When the middle anchor is absent and the outer span has
shortened by at least 50 bp, the read carries the programmed excision.

Being alignment-free matters here: a global aligner with the CRISPResso2 default
--needleman_wunsch_gap_extend -2 penalises the 137 bp deletion by about -294 and
diverts those reads out of the output, so the intended product of the experiment
disappears. An anchor method pays nothing for deletion length.

Usage:
    python3 y452b3_anchor_recount.py --fastq-dir data/fastq \
        --ref refs/Mm_Pla2g6_F2R2-amplicon.fa [--json out.json]
"""

import argparse
import collections
import json
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import (CUT1, CUT2, EXCISION_BP, SAMPLES, fastq_path,
                           fuzzy_find, load_reference, orient, read_fastq)

# anchors: L left of cut1, M between the cuts, R right of cut2
L_SPAN = (21, 44)
M_SPAN = (130, 153)
R_SPAN = (240, 263)
MIN_EXCISION_BP = 50
THRESHOLDS = (1, 2, 3, 5, 10)


def analyse(ref, path):
    aL = ref[L_SPAN[0] - 1:L_SPAN[1]]
    aM = ref[M_SPAN[0] - 1:M_SPAN[1]]
    aR = ref[R_SPAN[0] - 1:R_SPAN[1]]
    gap_LM = M_SPAN[0] - L_SPAN[1] - 1
    gap_MR = R_SPAN[0] - M_SPAN[1] - 1
    gap_LR = R_SPAN[0] - L_SPAN[1] - 1

    callable_reads = 0
    excision = 0
    excision_sizes = collections.Counter()
    worst_indel = collections.Counter()
    per_site = {"cut1_only": 0, "cut2_only": 0, "both": 0, "unedited": 0}
    dropped = {"no_left": 0, "no_right": 0, "ambiguous_middle": 0}

    for seq, _q in read_fastq(path):
        read = orient(seq, aL)
        if read is None:
            dropped["no_left"] += 1
            continue
        iL = fuzzy_find(read, aL)
        endL = iL + len(aL)
        iR = fuzzy_find(read[endL:], aR)
        if iR < 0:
            dropped["no_right"] += 1
            continue
        absR = endL + iR
        iM = fuzzy_find(read[endL:absR], aM)

        if iM < 0:
            net = (absR - endL) - gap_LR
            if net <= -MIN_EXCISION_BP:
                callable_reads += 1
                excision += 1
                excision_sizes[-net] += 1
            else:
                # middle anchor failed to match but the span is near full length:
                # an anchor-match failure, not an excision. Never assigned.
                dropped["ambiguous_middle"] += 1
            continue

        absM = endL + iM
        endM = absM + len(aM)
        n1 = (absM - endL) - gap_LM
        n2 = (absR - endM) - gap_MR
        callable_reads += 1
        worst_indel[max(abs(n1), abs(n2))] += 1
        if n1 and n2:
            per_site["both"] += 1
        elif n1:
            per_site["cut1_only"] += 1
        elif n2:
            per_site["cut2_only"] += 1
        else:
            per_site["unedited"] += 1

    if callable_reads == 0:
        raise RuntimeError(f"no callable reads in {path}")

    curve = {}
    for t in THRESHOLDS:
        n = sum(v for k, v in worst_indel.items() if k >= t) + excision
        curve[t] = round(100 * n / callable_reads, 2)

    return {
        "callable_reads": callable_reads,
        "excision_reads": excision,
        "excision_pct": round(100 * excision / callable_reads, 2),
        "excision_sizes": dict(excision_sizes.most_common()),
        "permissive_curve_pct": curve,
        "per_site": per_site,
        "dropped": dropped,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--json")
    args = ap.parse_args()

    ref = load_reference(args.ref)
    results = {}

    hdr = f"{'sample':<20}{'callable':>9}{'excision':>9}{'exc%':>8}" + \
          "".join(f"{'>=' + str(t) + 'bp':>9}" for t in THRESHOLDS)
    print(hdr)
    print("-" * len(hdr))
    for stem, label, _cond, _strain, _nuc in SAMPLES:
        r = analyse(ref, fastq_path(args.fastq_dir, stem))
        results[label] = r
        row = f"{label:<20}{r['callable_reads']:>9}{r['excision_reads']:>9}{r['excision_pct']:>7}%"
        row += "".join(f"{r['permissive_curve_pct'][t]:>8}%" for t in THRESHOLDS)
        print(row)

    print(f"\nblunt-blunt prediction: {EXCISION_BP} bp deleted "
          f"(positions {CUT1 + 1}-{CUT2}), {298 - EXCISION_BP} bp product")
    for label, r in results.items():
        top = list(r["excision_sizes"].items())[:5]
        if top:
            print(f"  {label:<20} modal excision sizes (bp:reads): {top}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
