"""Shared helpers for the Pla2g6 dual-guide amplicon analysis.

All coordinates are 1-based inclusive on the 298 bp amplicon.
Cas9 cuts between spacer positions 17 and 18 (3 bp from the PAM):
  guide 1 cuts between amplicon positions 73 and 74
  guide 2 cuts between amplicon positions 210 and 211
A precise dual-cut product removes positions 74-210 (137 bp) -> 161 bp.
"""

import os
import numpy as np

CUT1 = 73      # blunt cut lies between CUT1 and CUT1+1
CUT2 = 210
CDS = (51, 260)
EXCISION_BP = CUT2 - CUT1          # 137
PRODUCT_BP = 298 - EXCISION_BP     # 161
NO_CUT_WINDOWS = (30, 140, 170, 265)
ANCHOR_LEN = 20
ANCHOR_OFFSET = 5

_TR = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_TR)[::-1]


def load_reference(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    if len(seq) != 298:
        raise ValueError(f"expected a 298 bp amplicon, got {len(seq)} bp from {path}")
    return seq


def read_fastq(path, max_len=700, min_mean_q=12):
    """Yield (sequence, mean_quality) for reads passing basic filters.

    Reads longer than `max_len` are concatemers/chimeras and are dropped.
    """
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip().upper()
            fh.readline()
            qual = fh.readline().strip()
            if len(seq) > max_len:
                continue
            q = np.frombuffer(qual.encode(), dtype=np.uint8) - 33
            if q.size == 0 or q.mean() < min_mean_q:
                continue
            yield seq, float(q.mean())


def fuzzy_find(read, anchor, max_mismatch=2):
    """Index of `anchor` in `read`, allowing up to `max_mismatch` substitutions.

    Exact match is tried first; the Hamming scan is only reached when that fails.
    Returns -1 if no position is within the mismatch budget.
    """
    i = read.find(anchor)
    if i >= 0:
        return i
    n, m = len(read), len(anchor)
    if n < m:
        return -1
    r = np.frombuffer(read.encode(), dtype=np.uint8)
    a = np.frombuffer(anchor.encode(), dtype=np.uint8)
    mm = (np.lib.stride_tricks.sliding_window_view(r, m) != a).sum(1)
    j = int(mm.argmin())
    return j if mm[j] <= max_mismatch else -1


def anchors_for_site(ref, cut, anchor_len=ANCHOR_LEN, offset=ANCHOR_OFFSET):
    """Two anchors placed `offset` bp outside a cut site.

    Returns (left_anchor, right_anchor, expected_inner_gap).
    `cut` is the 1-based position immediately 5' of the blunt cut.
    """
    left = ref[cut - offset - anchor_len:cut - offset]
    right = ref[cut + offset:cut + offset + anchor_len]
    return left, right, 2 * offset


def orient(read, probe_anchor):
    """Return the read in reference orientation, or None if neither strand matches."""
    for candidate in (read, revcomp(read)):
        if fuzzy_find(candidate, probe_anchor) >= 0:
            return candidate
    return None


SAMPLES = [
    ("Y452B3_1_Cas9_Ctrl_1", "Cas9 ctrl C57BL/6", "cas9_only", "C57BL/6", 5),
    ("Y452B3_2_Pla2g6_1",    "Pla2g6 C57BL/6",    "dual_guide", "C57BL/6", 5),
    ("Y452B3_3_Cas9_Ctrl_2", "Cas9 ctrl NOD",     "cas9_only", "NOD", 6),
    ("Y452B3_4_Pla2g6_2",    "Pla2g6 NOD",        "dual_guide", "NOD", 6),
]


def fastq_path(fastq_dir, stem):
    for ext in (".fastq", ".fq", ".fastq.gz"):
        p = os.path.join(fastq_dir, stem + ext)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"no FASTQ for {stem} in {fastq_dir}")
