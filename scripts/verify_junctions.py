#!/usr/bin/env python3
"""
Independent verification of called alleles - framework section 16.6.4 extended.

The allele caller uses k-mer anchor chaining. This script re-counts the SAME
alleles by a completely different route: reconstruct the predicted product
sequence for each called allele, take a 30 bp string centred on its junction,
and count exact occurrences of that string (and its reverse complement) in the
raw FASTQ. No alignment, no chaining, no shared code path.

If the two counts agree, the allele is real and the chaining coordinates are
right. If the junction string is absent from the reads, the call is an
alignment artifact.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ont_guide_panel import rc, read_fastq
import yaml

FLANK = 15


def product(ref, sig):
    """Rebuild the predicted allele sequence from a signature string."""
    ev = []
    for part in sig.split("|"):
        kind = part[0]
        body = part[1:].split(":")
        sz = int(body[0])
        if kind == "D":
            s, e = (int(x) for x in body[1].split("-"))
            ev.append(("D", s, e, ""))
        else:
            s = int(body[1].split("_")[0])
            ev.append(("I", s, s, body[2] if len(body) > 2 else "N" * sz))
    ev.sort(key=lambda x: x[1])
    out, pos, juncs = [], 1, []
    for kind, s, e, seq in ev:
        out.append(ref[pos - 1:s - 1] if kind == "D" else ref[pos - 1:s])
        juncs.append(sum(len(x) for x in out))
        if kind == "I":
            out.append(seq)
            juncs.append(sum(len(x) for x in out))
            pos = s + 1
        else:
            pos = e + 1
    out.append(ref[pos - 1:])
    return "".join(out), juncs


def main(cfg_path, data_json, out_tsv):
    cfg = yaml.safe_load(open(cfg_path))
    ref = "".join(cfg["amplicon"]["sequence"].upper().split())
    data = json.load(open(data_json))
    root = os.path.dirname(os.path.abspath(cfg_path))
    for _s in cfg["samples"]:
        if not os.path.isabs(_s["fastq"]):
            _s["fastq"] = os.path.join(root, _s["fastq"])

    rows = []
    for s in cfg["samples"]:
        reads = [seq for seq, q in read_fastq(s["fastq"])]
        blob = "\n".join(reads)
        blob_rc = "\n".join(rc(r) for r in reads)
        for al in data["samples"][s["name"]]["top_alleles"]:
            prod, juncs = product(ref, al["sig"])
            j = juncs[0]
            probe = prod[max(0, j - FLANK):j + FLANK]
            if len(probe) < 2 * FLANK:
                continue
            n = blob.count(probe) + blob_rc.count(probe)
            wt_hits = ref.count(probe)
            rows.append({
                "sample": s["name"], "guide": s["guide"],
                "allele": al["sig"], "chaining_reads": al["reads"],
                "junction_probe_30bp": probe,
                "probe_reads_exact": n,
                "probe_in_wt_reference": wt_hits,
                "ratio_probe_over_chaining": round(n / al["reads"], 3),
                "predicted_product_len": al["predicted_len"],
                "observed_median_len": al["median_len"],
                "len_delta": al["median_len"] - al["predicted_len"],
            })
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(df.to_string(index=False))
    ok = df[df["probe_in_wt_reference"] == 0]
    print("\nalleles checked            : %d" % len(df))
    print("junction found in raw reads: %d/%d"
          % ((ok["probe_reads_exact"] > 0).sum(), len(ok)))
    print("median probe/chaining ratio: %.2f" % ok["ratio_probe_over_chaining"].median())
    print("median |length delta|      : %.1f bp" % ok["len_delta"].abs().median())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
