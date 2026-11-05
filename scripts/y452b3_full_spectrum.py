#!/usr/bin/env python3
"""SO-C main-text analyses: allele spectrum, reading frame, clustering, rarefaction.

Produces every number the main text needs from the control-matched dataset, with
the Cas9-only control scored identically at every step so each figure carries its
own background.

Outputs a JSON consumed by make_figure4.py / make_figure5.py / make_figure6.py.

Usage:
    python3 y452b3_full_spectrum.py --fastq-dir data/fastq \
        --ref refs/Mm_Pla2g6_F2R2-amplicon.fa --json so_c_spectrum.json
"""

import argparse
import collections
import json
import random
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pla2g6_common import (ANCHOR_LEN, ANCHOR_OFFSET, CDS, CUT1, CUT2, EXCISION_BP,
                           NO_CUT_WINDOWS, SAMPLES, anchors_for_site, fastq_path,
                           fuzzy_find, load_reference, orient, read_fastq)

L_SPAN, M_SPAN, R_SPAN = (21, 44), (130, 153), (240, 263)
MIN_EXCISION_BP = 50
RAREFACTION_DEPTHS = [100, 200, 400, 700, 1000, 1400, 1800, 2200]
RAREFACTION_REPS = 40


def classify_reads(ref, path):
    """Per-read calls: net indel at each cut site, or the excision class.

    Returns a list of records: ("excision", size) or ("local", n1, n2).
    """
    aL = ref[L_SPAN[0] - 1:L_SPAN[1]]
    aM = ref[M_SPAN[0] - 1:M_SPAN[1]]
    aR = ref[R_SPAN[0] - 1:R_SPAN[1]]
    gap_LM = M_SPAN[0] - L_SPAN[1] - 1
    gap_MR = R_SPAN[0] - M_SPAN[1] - 1
    gap_LR = R_SPAN[0] - L_SPAN[1] - 1

    records = []
    for seq, _q in read_fastq(path):
        read = orient(seq, aL)
        if read is None:
            continue
        iL = fuzzy_find(read, aL)
        endL = iL + len(aL)
        iR = fuzzy_find(read[endL:], aR)
        if iR < 0:
            continue
        absR = endL + iR
        iM = fuzzy_find(read[endL:absR], aM)
        if iM < 0:
            net = (absR - endL) - gap_LR
            if net <= -MIN_EXCISION_BP:
                records.append(("excision", -net))
            continue
        absM = endL + iM
        endM = absM + len(aM)
        records.append(("local", (absM - endL) - gap_LM, (absR - endM) - gap_MR))
    return records


def local_spectrum(ref, path, site):
    """Indel-size histogram at one site using tight local anchors.

    Used for the per-site spectra and for the matched-control background trace.
    """
    left, right, gap = anchors_for_site(ref, site, ANCHOR_LEN, ANCHOR_OFFSET)
    hist = collections.Counter()
    n = 0
    for seq, _q in read_fastq(path):
        read = orient(seq, left)
        if read is None:
            continue
        i = fuzzy_find(read, left)
        j = fuzzy_find(read[i + ANCHOR_LEN:], right)
        if j < 0:
            continue
        hist[j - gap] += 1
        n += 1
    return n, hist


def classes_from_records(records, threshold=1):
    """`threshold` is the minimum |indel| in bp counted as an edit at a site.

    At threshold 1 the Cas9-only control reads ~23% edited, which is nanopore
    error, so every figure using threshold 1 must show the control alongside.
    Threshold 3 drops the control to ~4%.
    """
    c = {"unedited": 0, "cut1_only": 0, "cut2_only": 0, "both_independent": 0,
         "excision": 0}
    for r in records:
        if r[0] == "excision":
            c["excision"] += 1
            continue
        _, n1, n2 = r
        n1 = n1 if abs(n1) >= threshold else 0
        n2 = n2 if abs(n2) >= threshold else 0
        if n1 and n2:
            c["both_independent"] += 1
        elif n1:
            c["cut1_only"] += 1
        elif n2:
            c["cut2_only"] += 1
        else:
            c["unedited"] += 1
    return c


def frame_composition(records):
    """Frameshift vs in-frame among edited reads, by net CDS length change."""
    fs = inf = 0
    for r in records:
        if r[0] == "excision":
            if r[1] % 3:
                fs += 1
            else:
                inf += 1
            continue
        _, n1, n2 = r
        net = n1 + n2
        if net == 0:
            continue
        if net % 3:
            fs += 1
        else:
            inf += 1
    return {"frameshift": fs, "in_frame": inf}


def clustering(records):
    """Joint vs marginal cutting, with excision counted as cut at both sites."""
    n = len(records)
    s1 = s2 = joint = 0
    for r in records:
        if r[0] == "excision":
            s1 += 1
            s2 += 1
            joint += 1
            continue
        _, n1, n2 = r
        if n1:
            s1 += 1
        if n2:
            s2 += 1
        if n1 and n2:
            joint += 1
    p1, p2, pj = s1 / n, s2 / n, joint / n
    return {"n": n, "p_site1": round(100 * p1, 2), "p_site2": round(100 * p2, 2),
            "p_joint": round(100 * pj, 2),
            "p_independent": round(100 * p1 * p2, 3),
            "enrichment": round(pj / (p1 * p2), 2) if p1 * p2 else None,
            "p_site2_given_site1": round(100 * pj / p1, 1) if p1 else None}


def rarefaction(records, depths, reps, seed=0):
    rng = random.Random(seed)
    out = {}
    for d in depths:
        if d > len(records):
            continue
        acc = collections.defaultdict(list)
        distinct = []
        for _ in range(reps):
            sub = rng.sample(records, d)
            c = classes_from_records(sub)
            for k, v in c.items():
                acc[k].append(100 * v / d)
            sizes = {r[1] for r in sub if r[0] == "excision"}
            sizes |= {r[1] for r in sub if r[0] == "local" and r[1]}
            distinct.append(len(sizes))
        out[d] = {k: {"mean": round(float(np.mean(v)), 2),
                      "sd": round(float(np.std(v)), 2)} for k, v in acc.items()}
        out[d]["distinct_alleles"] = {"mean": round(float(np.mean(distinct)), 1),
                                      "sd": round(float(np.std(distinct)), 1)}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--json", default="so_c_spectrum.json")
    args = ap.parse_args()

    ref = load_reference(args.ref)
    out = {}
    for stem, label, cond, strain, nuc in SAMPLES:
        path = fastq_path(args.fastq_dir, stem)
        records = classify_reads(ref, path)
        n1, h1 = local_spectrum(ref, path, CUT1)
        n2, h2 = local_spectrum(ref, path, CUT2)
        entry = {
            "label": label, "condition": cond, "strain": strain, "nucleofection": nuc,
            "scorable": len(records),
            "classes": classes_from_records(records),
            "classes_min3bp": classes_from_records(records, threshold=3),
            "frame": frame_composition(records),
            "clustering": clustering(records),
            "site1_spectrum": {"n": n1, "hist": {str(k): v for k, v in sorted(h1.items())}},
            "site2_spectrum": {"n": n2, "hist": {str(k): v for k, v in sorted(h2.items())}},
            "excision_sizes": {str(r[1]): 0 for r in records if r[0] == "excision"},
        }
        sizes = collections.Counter(r[1] for r in records if r[0] == "excision")
        entry["excision_sizes"] = {str(k): v for k, v in sorted(sizes.items())}
        if cond == "dual_guide":
            entry["rarefaction"] = rarefaction(records, RAREFACTION_DEPTHS,
                                               RAREFACTION_REPS)
        out[label] = entry

    with open(args.json, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"{'sample':<20}{'scorable':>9}{'unedit':>8}{'cut1':>7}{'cut2':>7}"
          f"{'both':>7}{'excis':>7}{'FS':>7}{'IF':>7}{'P(s2|s1)':>10}")
    for label, e in out.items():
        c, f, cl = e["classes"], e["frame"], e["clustering"]
        n = e["scorable"]
        print(f"{label:<20}{n:>9}{100*c['unedited']/n:>7.1f}%{100*c['cut1_only']/n:>6.1f}%"
              f"{100*c['cut2_only']/n:>6.1f}%{100*c['both_independent']/n:>6.1f}%"
              f"{100*c['excision']/n:>6.1f}%{100*f['frameshift']/n:>6.1f}%"
              f"{100*f['in_frame']/n:>6.1f}%{cl['p_site2_given_site1']:>9}%")
    print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
