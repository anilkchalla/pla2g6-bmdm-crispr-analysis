#!/usr/bin/env python3
"""
Stage 0 pre-flight guide validator.

CRISPResso2 finds a spacer by string match on either strand but does NOT verify
that a PAM is present, and computes the cut site assuming the orientation in
which it matched. A spacer entered in the wrong orientation therefore yields a
silently misplaced quantification window.

This script checks, for each spacer:
  - forward and reverse-complement match position in the amplicon
  - candidate PAM on each strand
  - PAM validity for the chosen nuclease
  - predicted blunt cleavage position
  - for >=2 guides, the predicted excision product

Usage:
  python3 check_guides.py --amplicon amplicon.fa \
      --guides CGTCAACACCCTCAGTAGCG,GAATGTCTCCAATCGCACCT \
      [--nuclease SpCas9]
  python3 check_guides.py --config config.yaml
"""
import argparse
import re
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# name -> (pam_regex, pam_side, spacer_len, cut_offset_from_3prime_end_of_spacer)
# cut_offset follows CRISPResso2's --quantification_window_center convention:
# position relative to the 3' end of the spacer (negative = inside the spacer).
NUCLEASES = {
    "SpCas9":  (r"[ACGT]GG",      "3prime", 20, -3),
    "SpCas9-NAG": (r"[ACGT]AG",   "3prime", 20, -3),
    "SpG":     (r"[ACGT][ACGT]G", "3prime", 20, -3),
    "SpRY":    (r"[ACGT][ACGT][ACGT]", "3prime", 20, -3),
    "SaCas9":  (r"[ACGT][ACGT]G[AG][AG]T", "3prime", 21, -3),
    "Cas12a":  (r"TTT[ACG]",      "5prime", 23, 1),
}


def rc(s):
    return s.translate(COMP)[::-1]


def find_all(hay, needle):
    return [m.start() for m in re.finditer("(?=%s)" % re.escape(needle), hay)]


def read_fasta(path):
    seq, name = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line)
    return name, "".join(seq).upper()


def evaluate(amp, spacer, nuclease):
    pam_re, pam_side, _slen, cut_off = NUCLEASES[nuclease]
    pam_len = len(re.findall(r"\[[^\]]*\]|[ACGT]", pam_re))
    results = []

    for strand, hits in (("+", find_all(amp, spacer)), ("-", find_all(amp, rc(spacer)))):
        for start in hits:
            end = start + len(spacer)  # exclusive, 0-based, top-strand coords
            # `cut` = 0-based slice boundary of the nick == number of top-strand
            # bases lying 5' of it. The nick sits between 1-based `cut` and `cut+1`.
            if strand == "+":
                pam = (amp[end:end + pam_len] if pam_side == "3prime"
                       else amp[start - pam_len:start] if start >= pam_len else "")
                cut = end + cut_off
            else:
                pam = (rc(amp[start - pam_len:start]) if pam_side == "3prime"
                       else rc(amp[end:end + pam_len]))
                if pam_side == "3prime" and start < pam_len:
                    pam = ""
                cut = start - cut_off
            ok = bool(pam) and len(pam) == pam_len and re.fullmatch(pam_re, pam) is not None
            results.append({
                "strand": strand,
                "start0": start, "end0": end,
                "start1": start + 1, "end1": end,
                "pam": pam or "(off amplicon)",
                "pam_ok": ok,
                # nick sits between 1-based positions `cut` and `cut+1`;
                # equivalently at 0-based slice boundary `cut`
                "cut_between_1based": (cut, cut + 1),
                "cut0": cut,
            })
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amplicon", help="FASTA file or raw sequence")
    ap.add_argument("--guides", help="comma-separated spacer sequences, no PAM")
    ap.add_argument("--nuclease", default="SpCas9", choices=sorted(NUCLEASES))
    ap.add_argument("--config", help="YAML config with amplicon/guides/nuclease")
    args = ap.parse_args()

    if args.config:
        try:
            import yaml
        except ImportError:
            sys.exit("pyyaml required for --config: pip install pyyaml")
        cfg = yaml.safe_load(open(args.config))
        amp = "".join(cfg["amplicon"]["sequence"].upper().split())
        amp_name = cfg["amplicon"].get("name", "amplicon")
        guides = [g["spacer"].upper() for g in cfg["guides"]]
        gnames = [g.get("name", "g%d" % (i + 1)) for i, g in enumerate(cfg["guides"])]
        nuclease = cfg.get("nuclease", "SpCas9")
    else:
        if not args.amplicon or not args.guides:
            sys.exit("need --amplicon and --guides (or --config)")
        if args.amplicon.endswith((".fa", ".fasta", ".fna")):
            amp_name, amp = read_fasta(args.amplicon)
        else:
            amp_name, amp = "amplicon", args.amplicon.upper()
        guides = [g.strip().upper() for g in args.guides.split(",")]
        gnames = ["g%d" % (i + 1) for i in range(len(guides))]
        nuclease = args.nuclease

    print("=" * 74)
    print("GUIDE PRE-FLIGHT VALIDATION")
    print("=" * 74)
    print("Amplicon : %s (%d bp)" % (amp_name, len(amp)))
    print("Nuclease : %s" % nuclease)
    print()

    valid_cuts = []
    problems = []

    for name, spacer in zip(gnames, guides):
        print("-" * 74)
        print("%s : %s (%d nt)" % (name, spacer, len(spacer)))
        hits = evaluate(amp, spacer, nuclease)
        if not hits:
            print("  !! NOT FOUND in amplicon on either strand")
            problems.append("%s: not found in amplicon" % name)
            continue
        good = [h for h in hits if h["pam_ok"]]
        for h in hits:
            flag = "OK " if h["pam_ok"] else "BAD"
            print("  [%s] strand %s  pos %d-%d (1-based)  PAM=%s  cut between %d/%d"
                  % (flag, h["strand"], h["start1"], h["end1"], h["pam"],
                     *h["cut_between_1based"]))
        if not good:
            print("  !! NO VALID PAM for the spacer as entered")
            # Is the reverse complement a valid protospacer at the same locus?
            alt = rc(spacer)
            alt_hits = [h for h in evaluate(amp, alt, nuclease) if h["pam_ok"]]
            if alt_hits:
                a = alt_hits[0]
                print()
                print("  ** ORIENTATION ERROR DETECTED")
                print("     The reverse complement IS a valid protospacer:")
                print("       spacer : %s" % alt)
                print("       strand %s  pos %d-%d (1-based)  PAM=%s  cut between %d/%d"
                      % (a["strand"], a["start1"], a["end1"], a["pam"],
                         *a["cut_between_1based"]))
                wrong = hits[0]["cut0"]
                print("     CRISPResso2, given the spacer as entered, would centre the")
                print("     quantification window at 1-based %d instead of %d "
                      "-- an error of %d bp."
                      % (wrong + 1, a["cut0"] + 1, abs(a["cut0"] - wrong)))
                print("     -> ACTION: enter this spacer as %s" % alt)
                problems.append("%s: WRONG ORIENTATION, use %s" % (name, alt))
                valid_cuts.append((name + "(corrected)", a["cut0"]))
            else:
                print("     -> neither orientation has a valid %s PAM at this locus."
                      % nuclease)
                print("        Check the amplicon, the spacer, or the nuclease "
                      "(SpG/SpRY/Cas12a?).")
                print("     -> reverse complement is: %s" % alt)
                problems.append("%s: no valid PAM in either orientation" % name)
        else:
            h = good[0]
            if h["strand"] == "-":
                print("     (bottom-strand protospacer; CRISPResso2 detects the")
                print("      orientation and places the cut correctly)")
            valid_cuts.append((name, h["cut0"]))
        print()

    if len(valid_cuts) >= 2:
        print("=" * 74)
        print("PREDICTED DUAL-GUIDE EXCISION PRODUCT")
        print("=" * 74)
        vc = sorted(valid_cuts, key=lambda x: x[1])
        for (n1, c1), (n2, c2) in zip(vc, vc[1:]):
            prod = amp[:c1] + amp[c2:]
            print("%s (cut after 1-based %d) x %s (cut after 1-based %d)"
                  % (n1, c1, n2, c2))
            print("  deleted segment : 1-based %d-%d  (%d bp)"
                  % (c1 + 1, c2, c2 - c1))
            print("  product length  : %d bp  (WT %d bp)" % (len(prod), len(amp)))
            print("  product sequence:")
            for i in range(0, len(prod), 70):
                print("    " + prod[i:i + 70])
            print()

    print("=" * 74)
    if problems:
        print("STATUS: %d issue(s) require attention" % len(problems))
        for p in problems:
            print("  - " + p)
        sys.exit(1)
    print("STATUS: all guides validated")


if __name__ == "__main__":
    main()
