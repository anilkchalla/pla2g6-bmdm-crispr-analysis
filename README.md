# Pla2g6 dual-guide CRISPR analysis in primary mouse BMDMs

Analysis code and reference files for:

> Challa AK, Nair R, Arumugam G, Huynh TYL, Folmar JR, Ahamed A, Chowdarapu S, Ramanadham S.
> *Silent failure modes in CRISPR-Cas9 gene targeting of primary mouse macrophages: guide RNA
> chemistry, heteroduplex screening and amplicon-sequencing analysis.* Preprint, 2026.

This repository contains everything needed to reproduce the reported numbers from raw reads.
Raw sequencing data are deposited separately (see [Data availability](#data-availability)).

---

## What the analysis does

The experiment places two Cas9 cut sites 137 bp apart in a 298 bp amplicon spanning the exon 2
region of mouse *Pla2g6*, so that cutting at both sites excises the intervening fragment. The
analysis has to answer three questions that a single "percent edited" cannot:

1. How much editing occurred at each site, and at both together?
2. How much of the apparent editing is sequencing or PCR artifact?
3. How many edited alleles are actually null alleles?

**All coordinates in this repository are 1-based and inclusive**, referring to the 298 bp amplicon
in `refs/Mm_Pla2g6_F2R2-amplicon.fa`.

### Target geometry

Cas9 cleaves between positions 17 and 18 of the spacer, i.e. 3 bp from the PAM.

| | Spacer | Strand | Amplicon coords | PAM | Blunt cut |
|---|---|---|---|---|---|
| Guide 1 | `CGCTACTGAGGGTGTTG`\|`ACG` | minus | 71–90 | `AGG` (plus 68–70 = `CCT`) | 73\|74 |
| Guide 2 | `GAATGTCTCCAATCGCA`\|`CCT` | plus | 194–213 | `GGG` (plus 214–216) | 210\|211 |

A precise dual-cut product removes positions **74–210** (137 bp), leaving a **161 bp** product.
Positions 74 and 211 are both `C`, so the junction carries 1 bp of microhomology and the gap can be
left- or right-aligned with identical sequence; we report the left-aligned form because it
coincides with the cut sites. Aligners that right-align will report deletion of 75–211 for the
same allele.

### Method, in order

**Stage 0 — guide validation (`scripts/check_guides.py`).**
Locates each spacer in both orientations and *requires an NGG PAM at the inferred position* before
any quantification window is defined. This is the one step that must never be skipped. CRISPResso2
matches a spacer by string search on either strand but does not verify the PAM, so a spacer entered
in the wrong orientation yields a silently misplaced window — in the Dataset A run this displaced guide
1's window by 14 bp with no warning. Vendor pipelines that *infer* the spacer from the data can
land one base away, which changes the predicted product by one base. Supply spacers explicitly and
check them.

**Stage 1 — read accounting (`scripts/qc_gate.py`).**
Every read is accounted for: total, passing quality, matching the primary amplicon, spanning the
region of interest, scorable. Reads that drop out at each step are reported rather than silently
excluded, because dropout is usually non-random — reads carrying an anchor error are enriched for
indels.

**Stage 2 — anchor-based allele calling (`scripts/y452b3_anchor_recount.py`).**
Two 20 bp anchors are placed 5 bp outside each cut site; indel size is the net change in the
distance between them. This is alignment-free, so a 137 bp deletion costs nothing to detect — in
contrast to global alignment, where the default `--needleman_wunsch_gap_extend -2` penalises that
deletion by about −294 and diverts such reads out of the output entirely. A read is called as
carrying the excision when the inner anchor is absent and the outer span has shortened by at least
50 bp. Reads with a failed anchor match are reported as a separate non-scorable class, never
assigned.

**Stage 3 — allele reclassification (`scripts/y452b3_reclassify_alleles.py`).**
Vendor per-allele tables are anchored on a single quantification window, so an allele edited only
at the *other* cut site is not counted. This script re-parses the CIGAR of every allele and assigns
each indel event to cut 1, cut 2, or both by its reference interval, recovering total editing
across both sites. It also classifies each allele by reading frame over the CDS to separate
frameshift from in-frame outcomes.

**Stage 4 — background estimation (`scripts/y452b3_error_floor.py`).**
Two estimates, which do not agree:

- *Matched Cas9-only control* at the actual cut sites. This is the correct denominator.
- *Within-amplicon no-cut windows* at positions 30, 140, 170, 265, scored with identical anchor
  geometry. This is the fallback when no control exists.

Running both on the same control reads shows the fallback **overestimates cut-site background by
two- to threefold** and varies **sevenfold** with window placement (9.53% at position 30 vs 1.20%
at position 265 in the same amplicon). Treat it as a conservative bound, not a measurement.

**Stage 5 — substitution audit (`scripts/snp_audit.py`).**
Substitutions are never counted as editing. Omitting this inflated the Dataset A estimate by 10.8%
relative.

**Stage 6 — figures (`scripts/make_figure7.py`, `scripts/make_figures.py`).**

### Why thresholds matter

Scored permissively — counting any net length change across the cut-site interval — a Cas9-only
control reads **22–23% "edited"**. Requiring indels of at least 5 bp drops it below 1%; requiring
breakpoints at a cut site drops it to 0.07–0.28%. The dual-guide samples stay above 44% at every
threshold. Permissive calling is not the conservative choice; it manufactures signal, and the size
of the effect is invisible without a control.

---

## Datasets

Guide chemistry and sequencing platform differ together between the three datasets, so each is named
by both rather than by date. No two differ in only one variable; none of the comparisons between them
is controlled.

| Short form | Dataset | Guide chemistry | Platform | Matched control | Samples |
|---|---|---|---|---|---|
| **Dataset A** | IVT–Illumina | in vitro transcribed | Illumina MiSeq | none | 1 |
| **Dataset B** | synthetic–ONT | chemically synthesised | ONT | none | 1 |
| **Dataset C** | synthetic–ONT, controlled | chemically synthesised | ONT | yes, per nucleofection | 4 |

The scripts prefixed `y452b3_` operate on **Dataset C** (vendor order Y452B3).

## Headline results reproduced by this code

**Dataset C**: two independent nucleofections in two mouse backgrounds, each with a
matched Cas9-only control:

| Sample | Analysable reads | Total edited | Excision | Frameshift | In-frame |
|---|---|---|---|---|---|
| Cas9 control, Nuc#5 (C57BL/6) | 1,417 | 0.28% | 0.00% | 0.28% | 0% |
| Pla2g6 RNP, Nuc#5 (C57BL/6) | 1,752 | 54.7% | 33.1% | 34.8% | 20.0% |
| Cas9 control, Nuc#6 (NOD) | 2,966 | 0.07% | 0.00% | 0.07% | 0% |
| Pla2g6 RNP, Nuc#6 (NOD) | 2,380 | 44.5% | 26.2% | 31.5% | 12.9% |

Across 4,383 Cas9-only control reads, **zero** carried a deletion spanning both cut sites.

*Control values are the vendor breakpoint-anchored `edit_rate`. `y452b3_reclassify_alleles.py`
reports 0.00% rather than 0.07% for the NOD control because it sets aside two reads flagged as
very-low-support; the difference is two reads and immaterial to any conclusion.*

The excision allele is identical between backgrounds at base resolution: a 137 bp deletion of
positions 74–210, junction `...TCGGACGCCTCGTC / CTGGGACTGTGTCC...`, giving the 161 bp product seen
as a discrete band by HMA.

---

## Quick start

```bash
conda env create -f environment.yml
conda activate pla2g6-analysis

# 1. validate guides BEFORE anything else
python3 scripts/check_guides.py --config configs/pla2g6_dual_guide.yaml

# 2. anchor-based recount from raw FASTQ
python3 scripts/y452b3_anchor_recount.py --fastq-dir data/fastq --ref refs/Mm_Pla2g6_F2R2-amplicon.fa

# 3. reclassify vendor allele tables across both cut sites
python3 scripts/y452b3_reclassify_alleles.py --edit-analysis data/Y452B3_edit-analysis

# 4. background: matched control vs within-amplicon no-cut windows
python3 scripts/y452b3_error_floor.py --fastq-dir data/fastq --ref refs/Mm_Pla2g6_F2R2-amplicon.fa

# 5. figure
python3 scripts/make_figure7.py --fastq-dir data/fastq --out Figure7_matched_control.png
```

`data/` is empty in this repository; see below.

---

## Repository layout

```
├── README.md                             this file — the method
├── LICENSE                               MIT
├── environment.yml                       conda environment
├── CITATION.cff                          how to cite
├── .gitignore                            excludes FASTQ, PDFs, __pycache__
├── refs/
│   └── Mm_Pla2g6_F2R2-amplicon.fa        298 bp reference amplicon
├── configs/
│   └── pla2g6_dual_guide.yaml            spacers, PAMs, cut sites, CDS bounds
├── scripts/
│   ├── check_guides.py                   Stage 0 — guide/PAM/cut-site validator
│   ├── qc_gate.py                        Stage 1 — read-accounting ledger
│   ├── y452b3_anchor_recount.py          Stage 2 — alignment-free allele calling
│   ├── y452b3_reclassify_alleles.py      Stage 3 — CIGAR reclassification, both cut sites
│   ├── y452b3_error_floor.py             Stage 4 — matched control vs no-cut windows
│   ├── snp_audit.py                      Stage 5 — substitution audit
│   ├── positional_scan.py                length-unbiased sham-cut-site scan
│   ├── control_compare.py                Fisher + BH vs matched control
│   ├── analyze_alleles.py                allele attribution, evidence tiers, Wilson CIs
│   ├── ont_amplicon.py                   long-read analyzer
│   ├── verify_junctions.py               junction verification against assembly
│   ├── compare_crispresso.py             CRISPResso2 parameter comparison
│   ├── make_figures.py                   Figures 4–6
│   ├── make_figure7.py                   Figure 7 — matched control
│   └── run_crispresso.sh                 config-driven CRISPResso2 runner
└── data/                                 (empty — see Data availability)
```

---

## Data availability

Raw reads are **not** in this repository. Deposit and link them here:

| Dataset | Platform | Run | Accession |
|---|---|---|---|
| Dataset A | Illumina MiSeq, 2×251 bp | — | PRJNA1554390, run `[SRR to be added]` |
| Dataset B | ONT (Plasmidsaurus) | 97K3PK | PRJNA1554390, run `[SRR to be added]` |
| Dataset C | ONT (Plasmidsaurus) | Y452B3 | PRJNA1554390, runs `[SRR to be added]` |

To reproduce, download the reads into `data/fastq/` and the vendor edit-analysis tables into
`data/Y452B3_edit-analysis/`.

---

## A note on the Dataset A excision figure

The Dataset A library gives 0.019–0.132% excision while HMA on the same PCR product shows a clear
161 bp band, and the Dataset B data give 26–33%. Two mechanisms depress the Illumina number and
both act in the same direction: default global-alignment gap penalties discard the 137 bp deletion,
and size-selective bead cleanup during library preparation depletes a 161 bp product that HMA and
nanopore sequencing, which read the raw PCR product, both retain. The Dataset A value is a floor and
should not be quoted as an editing efficiency or compared directly with the Dataset B figures.

---

## License

MIT — see [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).

## Contact

Anil K. Challa — akchalla@uab.edu
