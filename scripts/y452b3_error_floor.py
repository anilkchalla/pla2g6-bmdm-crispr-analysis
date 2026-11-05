#!/usr/bin/env python3
"""Stage 4 - background estimation: matched control vs within-amplicon surrogate.

Two ways to estimate the sequencing error floor, scored with identical anchor
geometry so they are directly comparable:

  1. the matched Cas9-only control, at the ACTUAL cut sites  -> the correct answer
  2. no-cut "sham" windows elsewhere in the same amplicon    -> the fallback when
     no control exists

Running both on the same control reads shows the fallback overestimates cut-site
background by two- to threefold, and that it varies sevenfold depending on which
windows are chosen (9.53% at position 30 vs 1.20% at position 265 in the same
298 bp amplicon). It is a conservative bound, not a measurement.

Note that this geometry deliberately scores only small indels local to each site:
a read carrying the programmed excision loses its inner anchors and is excluded,
so these numbers are not editing efficiencies. Use y452b3_anchor_recount.py for
those.

Usage:
    python3 y452b3_error_floor.py --fastq-dir data/fastq \
        --ref refs/Mm_Pla2g6_F2R2-amplicon.fa
"""

import argparse
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import (ANCHOR_LEN, ANCHOR_OFFSET, CUT1, CUT2, NO_CUT_WINDOWS,
                           SAMPLES, anchors_for_site, fastq_path, fuzzy_find,
                           load_reference, orient, read_fastq)


def score_site(ref, reads, site):
    """Apparent indel rate at one site, from pre-loaded reads."""
    left, right, gap = anchors_for_site(ref, site, ANCHOR_LEN, ANCHOR_OFFSET)
    n = any_indel = ge3 = 0
    for seq in reads:
        read = orient(seq, left)
        if read is None:
            continue
        i = fuzzy_find(read, left)
        end = i + ANCHOR_LEN
        j = fuzzy_find(read[end:], right)
        if j < 0:
            continue
        net = j - gap
        n += 1
        if net != 0:
            any_indel += 1
        if abs(net) >= 3:
            ge3 += 1
    if n == 0:
        return n, float("nan"), float("nan")
    return n, 100 * any_indel / n, 100 * ge3 / n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--ref", required=True)
    args = ap.parse_args()

    ref = load_reference(args.ref)
    sites = [(f"no-cut window {w}", w) for w in NO_CUT_WINDOWS] + \
            [(f"CUT 1 ({CUT1}|{CUT1 + 1})", CUT1), (f"CUT 2 ({CUT2}|{CUT2 + 1})", CUT2)]

    table = {}
    for stem, label, cond, _strain, _nuc in SAMPLES:
        reads = [s for s, _q in read_fastq(fastq_path(args.fastq_dir, stem))]
        table[label] = {name: score_site(ref, reads, site) for name, site in sites}

    labels = list(table)
    print(f"{'site':<22}" + "".join(f"{lab:>26}" for lab in labels))
    print(f"{'':<22}" + "".join(f"{'n    any%   >=3bp%':>26}" for _ in labels))
    print("-" * (22 + 26 * len(labels)))
    for name, _site in sites:
        row = f"{name:<22}"
        for lab in labels:
            n, a, g = table[lab][name]
            row += f"{n:>9}{a:>8.2f}{g:>9.2f}"
        print(row)

    print()
    for lab in labels:
        sham = [table[lab][f"no-cut window {w}"][1] for w in NO_CUT_WINDOWS]
        sham3 = [table[lab][f"no-cut window {w}"][2] for w in NO_CUT_WINDOWS]
        c1 = table[lab][f"CUT 1 ({CUT1}|{CUT1 + 1})"][1]
        c2 = table[lab][f"CUT 2 ({CUT2}|{CUT2 + 1})"][1]
        mean = sum(sham) / len(sham)
        spread = max(sham) / min(sham) if min(sham) else float("nan")
        print(f"{lab:<22} sham mean {mean:5.2f}% (>=3bp {sum(sham3)/len(sham3):4.2f}%), "
              f"spread {spread:.1f}x | real cut sites {c1:.2f}% / {c2:.2f}%")


if __name__ == "__main__":
    main()
