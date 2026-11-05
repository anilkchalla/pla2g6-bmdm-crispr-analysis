#!/usr/bin/env bash
# Config-driven CRISPResso2 runner.
# Usage: bash scripts/run_crispresso.sh projects/<project>/config.yaml
set -euo pipefail

CFG="${1:?usage: run_crispresso.sh <config.yaml>}"
PROJ="$(cd "$(dirname "$CFG")" && pwd)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# --- parse config (requires python3 + pyyaml) --------------------------------
eval "$(python3 - "$CFG" <<'PY'
import sys, yaml, shlex
c = yaml.safe_load(open(sys.argv[1]))
amp = c["amplicon"]
refs = [("Reference", "".join(amp["sequence"].split()))]
for e in c.get("extra_amplicons", []):
    refs.append((e["name"], "".join(e["sequence"].split())))
guides = c["guides"]
p = c.get("params", {})
def emit(k, v): print("%s=%s" % (k, shlex.quote(str(v))))
emit("SAMPLE",   c["sample"])
emit("R1",       c["reads"]["r1"])
emit("R2",       c["reads"].get("r2", ""))
emit("AMPSEQ",   ",".join(s for _, s in refs))
emit("AMPNAME",  ",".join(n for n, _ in refs))
emit("GUIDESEQ", ",".join(g["spacer"] for g in guides))
emit("GUIDENAME",",".join(g["name"] for g in guides))
emit("CUTSITES", ",".join(str(g["cut_site_1based"]) for g in guides))
emit("W",        p.get("quantification_window_size", 5))
emit("WC",       p.get("quantification_window_center", -3))
emit("MINQ",     p.get("min_average_read_quality", 30))
emit("MINBPQ",   p.get("min_single_bp_quality", 20))
emit("GAPEXT",   p.get("needleman_wunsch_gap_extend", 0))
emit("EXCL_L",   p.get("exclude_bp_from_left", 15))
emit("EXCL_R",   p.get("exclude_bp_from_right", 15))
emit("PLOTW",    p.get("plot_window_size", 30))
emit("MINFREQ",  p.get("min_frequency_alleles_around_cut_to_plot", 0.1))
emit("THREADS",  p.get("threads", 4))
PY
)"

cd "$PROJ"
mkdir -p 01_qc 02_crispresso 03_summary 04_report

# --- Stage 0: pre-flight guide validation ------------------------------------
echo ">>> Stage 0: guide validation"
python3 "$ROOT/scripts/check_guides.py" --config "$CFG" \
  | tee 01_qc/guide_validation.txt || {
    echo "!! Guide validation reported issues. Review 01_qc/guide_validation.txt."
    echo "!! Set force: true in config to proceed anyway."
    if [ "$(python3 -c "import yaml,sys;print(yaml.safe_load(open('$CFG')).get('force',False))")" != "True" ]; then
      exit 1
    fi
  }

# --- Stage 1: raw QC ----------------------------------------------------------
echo ">>> Stage 1: FastQC"
if command -v fastqc >/dev/null; then
  fastqc -q -o 01_qc "$R1" ${R2:+"$R2"}
else
  echo "    fastqc not found, skipping"
fi

# --- Stage 2-4: CRISPResso2 ---------------------------------------------------
echo ">>> Stage 2-4: CRISPResso2"
CRISPResso \
  -r1 "$R1" ${R2:+-r2 "$R2"} \
  -a  "$AMPSEQ" --amplicon_name "$AMPNAME" \
  -g  "$GUIDESEQ" --guide_name "$GUIDENAME" \
  -w "$W" -wc "$WC" \
  -q "$MINQ" -s "$MINBPQ" \
  --needleman_wunsch_gap_extend "$GAPEXT" \
  --exclude_bp_from_left "$EXCL_L" --exclude_bp_from_right "$EXCL_R" \
  --write_detailed_allele_table \
  --plot_window_size "$PLOTW" \
  --min_frequency_alleles_around_cut_to_plot "$MINFREQ" \
  -n "$SAMPLE" -o 02_crispresso -p "$THREADS"

OUT="02_crispresso/CRISPResso_on_${SAMPLE}"

# --- Stage 5: QC gate + summary -----------------------------------------------
echo ">>> Stage 5: QC gate"
python3 "$ROOT/scripts/qc_gate.py" --crispresso-dir "$OUT" \
  | tee 01_qc/read_accounting_ledger.txt

echo ">>> Stage 5: allele summary"
python3 "$ROOT/scripts/summarize_alleles.py" \
  --crispresso-dir "$OUT" \
  --cut-sites "$CUTSITES" --guide-names "$GUIDENAME" \
  --window "$W" --outdir 03_summary

echo
echo "Done. Outputs:"
echo "  01_qc/       guide validation, FastQC, read ledger"
echo "  02_crispresso/ CRISPResso2 report ($OUT)"
echo "  03_summary/  allele tables, guide efficiency, figures"
