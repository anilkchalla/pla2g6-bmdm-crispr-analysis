# Audit — Pla2g6 dual-guide run (CRISPResso2 v2.3.3, July 2024 data)

**Revision 2 — 6 August 2026.** Reanalysed under framework v2.0.
**Revision 1 contained an error in the central biological conclusion. See §3.**

All figures recomputed from the run's own output files. Coordinates are 1-based inclusive.
No matched untreated control exists for this sample, so **every editing figure below is putative** (framework §1).

---

## 1. Headline result

| Metric | Value | 95% CI |
|---|---|---|
| Reads in input | 570,586 | |
| Analyzable reads (`N`) | 401,563 (70.4% of input) | |
| Exact-reference WT | 81.54% | 81.42–81.66% |
| SNP-only (not editing) | 12.03% | 11.93–12.13% |
| **Raw core indel efficiency (±5 bp)** | **1.78%** | 1.74–1.82% |
| Guide 1 marginal, endpoint-anchored | 0.41% | 0.39–0.43% |
| Guide 2 marginal, endpoint-anchored | 0.70% | 0.67–0.73% |
| Clean 137 bp excision | 0.019% | |
| Distinct edit alleles | 622 | |
| Edit alleles meeting main threshold (≥402 reads) | **0** | |

**Both guides cut. Both are weak.** Neither reaches 1%. The paired excision at 0.019% is what two sub-1% guides predict. The limiting factor is editing efficiency at this locus overall, not one defective guide.

---

## 2. Finding 1 — Guide 1's spacer was entered in the wrong orientation

The spacer used was `CGTCAACACCCTCAGTAGCG`. It matches the **top** strand at 71–90, but the three nucleotides 3′ of it are **TCA** — not an NGG PAM. Its reverse complement, `CGCTACTGAGGGTGTTGACG`, is a valid protospacer on the **bottom** strand with PAM **AGG** (top-strand `CCT` at 68–70).

CRISPResso2 matched the top strand and recorded cut points at 0-based `[86, 209]`. The true guide 1 cut is at 0-based 73 — a **14 bp error** in window placement.

**Independent confirmation.** The predicted blunt excision product using corrected cut sites (delete 74–210, 137 bp) is 161 bp and **byte-identical** to the sequence supplied as `--expected_hdr_amplicon_seq`. The design was right; only the entered orientation was wrong.

```
        68  71                    90  93
top  ...CCT CGTCAACACCCTCAGTAGCG TCA...
        ^^^ ^^^^^^^^^^^^^^^^^^^^ ^^^
        |   spacer as entered    no PAM
        PAM (AGG on bottom strand)
              ▲ cut 73/74  ← true          ▲ cut 87/88  ← as run
```

CRISPResso2 never checks for a PAM, so this passed silently. `scripts/check_guides.py` now blocks it.

---

## 3. CORRECTION — Guide 1 is active, contrary to revision 1

**Revision 1 of this audit concluded guide 1 was inactive and recommended redesigning it. That was wrong.**

Revision 1 ranked guides using CRISPResso2's `Effect_vector_deletion.txt`, which counts **every deleted base** in every read. A single 53 bp deletion contributes to 53 positions. Guide 2 produces many large deletions and so dominated the per-base signal; guide 1 produces mostly small discrete indels and appeared to sit at background (1.38× vs 6.34×). The background estimate was itself inflated by the same large-deletion shoulders.

The length-unbiased test (`positional_scan.py`) counts each read **once**, at its indel breakpoint, at every candidate nick position across the amplicon. Non-guide positions act as sham cut sites; their median is the empirical background.

| Position | Anchored reads | % of `N` | × amplicon median | Rank of 257 |
|---|---:|---:|---:|---:|
| Sham background (median) | 312 | 0.0777% | 1.0× | — |
| **Guide 1 (74)** | **1,663** | **0.4141%** | **5.3×** | **19** |
| **Guide 2 (211)** | **2,806** | **0.6988%** | **9.0×** | **8** |

The guide 1 profile is a discrete peak, not a plateau of background: 1.1–1.3× at positions 60–62, rising to 5.3–5.8× across 68–80, falling to 1.3–1.5× by 84–94. Of the top 25 positions amplicon-wide, 14 are within ±15 bp of guide 2 and **8 of the remaining 11 are in the guide 1 region**. See `03_summary/fig5_positional_scan.png`.

Corroborating allele-level evidence: the single most abundant edit allele in the entire sample is **`D3:72-74`, 351 reads (0.087%)** — a 3 bp deletion sitting exactly on the guide 1 nick.

**Revised recommendation.** Guide 1 does not need redesigning on this evidence. Both guides cut at sub-1% efficiency, which points at delivery, RNP dose, exposure time, or locus accessibility rather than guide sequence. Correcting the orientation and re-running will place guide 1's window correctly, but the efficiency will remain low.

---

## 4. Finding 2 — Substitutions were counted as editing

The run did not use `--ignore_substitutions`.

| | Reads | % of `N` |
|---|---:|---:|
| CRISPResso2 `Modified%` as reported | 6,050 | 1.5066% |
| of which **substitution-only** | 591 | 0.1472% |
| indel-containing (correct numerator) | 5,459 | 1.3594% |

Editing was overstated by **10.8% relative**. The `1.78%` headline in §1 is higher than 1.36% only because it uses a ±5 bp window instead of CRISPResso2's ±1 bp; at ±1 bp the framework reproduces 1.685% span-rule / 1.36% CRISPResso-equivalent.

Window sensitivity (span rule, substitutions excluded):

| Core window | Reads | % of `N` | g1 | g2 |
|---|---:|---:|---:|---:|
| ±1 bp | 6,768 | 1.685% | 0.581% | 1.201% |
| ±5 bp | 7,148 | 1.780% | 0.631% | 1.257% |
| ±20 bp | 7,345 | 1.829% | 0.654% | 1.285% |

Editing is not strongly window-dependent, which is itself reassuring.

## 5. Finding 3 — The excision product was declared as an HDR donor

Passed via `--expected_hdr_amplicon_seq`. 1,381 reads were assigned to the HDR reference; only **76 were unmodified**. The reported `HDR Modified% = 94.5%` is an artifact — HDR reads are scored over a window not anchored to a cut site, so any substitution anywhere marks a read modified. The interpretable figure is 76 reads = **0.019%**.

Correct declaration: a second entry in `-a` with `--amplicon_name "Reference,Excision"`. Applied in `config.yaml`.

Consistency check: if the two cuts occurred independently, excision ≈ 0.41% × 0.70% ≈ 0.003%. Observed 0.019% is ~6× higher, consistent with correlated cutting in the same cells at the same RNP dose. The numbers hang together.

## 6. Finding 4 — Alignment parameters suppress large deletions

`--needleman_wunsch_gap_extend` was left at the default −2. A 137 bp deletion then costs ≈ −20 + (137 × −2) = −294, pushing true excision reads below the alignment threshold. This biases the excision estimate **downward** by an unknown amount, so 0.019% is a floor. Set to `0`.

## 7. Finding 5 — QC failures

| Gate | Value | Threshold | Result |
|---|---:|---:|---|
| Preprocessing retention | 77.4% | ≥80% | **investigate** |
| Alignment rate | 91.2% | ≥90% | pass |
| Analyzable depth | 401,563 | ≥10,000 | pass |

128,774 reads (22.6%) lost before alignment. From `CRISPResso_RUNNING_LOG.txt`:

- `Your FASTQ may be invalid, please check the tail` — R2 truncated or corrupted at the tail.
- `different read numbers of the 1789 pack / Read1 256, Read2 194 / Ignore the unmatched reads` — R1 and R2 not properly paired; fastp saw 458,178 pairs against 570,586 R1 records.

Re-download, verify checksums, re-pair with `seqkit pair`.

## 8. Finding 6 — "Unmodified" is not wild-type

| Class | Distinct | Reads | % of `N` |
|---|---:|---:|---:|
| Exact-reference WT | 1 | 327,452 | **81.54%** |
| SNP-only | 12,771 | 48,315 | 12.03% |
| Unmodified, other (N-containing) | 491 | 17,734 | 4.42% |
| CRISPR-associated indel | 622 | 7,148 | 1.78% |
| Distal/background indel | — | 914 | 0.23% |

The headline `Unmodified% = 98.49` must not be read as 98.49% perfect wild type. Only 81.5% of reads match the reference exactly. The 12% SNP-only fraction is sequencing error across a 268 bp callable interval, not polymorphism — the substitution audit found **0 positions** meeting the ≥1% and ≥20-read reporting threshold.

## 9. Finding 7 — No mock control

With editing at ~1.8% and no untreated sample, every figure here is **putative** under framework §1. The positional scan is suggestive but is not a substitute. A mock control sequenced on the same run and analysed identically is the single highest-value addition to this experiment.

## 10. Finding 8 — Editing is diffuse, not dominated by one allele

622 distinct edit alleles; **none** reaches the main-report threshold of ≥10 reads *and* ≥0.1% of `N` (≥402 reads). The largest is `D3:72-74` at 351 reads (0.087%). 324 alleles meet the ≥3-read appendix threshold.

At this efficiency the repair outcome is a smear rather than a dominant product, which matters for downstream cloning: no single allele is common enough to expect it reliably in a small set of picked clones.

---

## 11. Recommended actions

| # | Action | Effort | Changes the answer? |
|---|---|---|---|
| 1 | Sequence and analyse a matched mock control | medium | **Yes — converts putative to established** |
| 2 | Fix the FASTQ files (re-download, checksum, re-pair) | low | Yes — recovers up to 23% of reads |
| 3 | Re-run both passes from `config.yaml` (corrected orientation, Excision as reference amplicon, `--ignore_substitutions`, `gap_extend 0`) | low | Yes — correct windows and numerator |
| 4 | Improve delivery / RNP dose / exposure before changing guides | high | **Yes — this is the actual blocker** |
| 5 | Confirm excision by ddPCR or long-read sequencing | high | Confirmatory; current 0.019% is a floor |
| 6 | Biological replicates | high | Required for any inferential claim |

Note the reordering from revision 1, which put "redesign guide 1" at rank 4. That recommendation is withdrawn.

---

## 12. Provenance

Source run: `../../CRISPResso_output/CRISPResso_on_Pla2g6_dual_guide/` (CRISPResso2 v2.3.3)

```
CRISPResso -r1 ..._R1_001.fastq -r2 ..._R2_001.fastq \
  -a <298bp WT> -g CGTCAACACCCTCAGTAGCG,GAATGTCTCCAATCGCACCT \
  -e <161bp excision> --expected_hdr_amplicon_seq <161bp excision> \
  -n Pla2g6_dual_guide -o CRISPResso_output --plot_window_size 30
```

**Parser validation.** The framework's alignment parser reproduces CRISPResso2's own indel count to within 1 read of 5,459 (99.98% concordance) when restricted to the same ±1 bp quantification window. Differences at wider windows are scope, not disagreement.

Reproduce:

```bash
cd scripts
python3 check_guides.py    --config ../projects/Pla2g6_dual_guide_2024-07-15/config.yaml
python3 qc_gate.py         --crispresso-dir ../CRISPResso_output/CRISPResso_on_Pla2g6_dual_guide
python3 analyze_alleles.py --crispresso-dir ../CRISPResso_output/CRISPResso_on_Pla2g6_dual_guide \
                           --config ../projects/Pla2g6_dual_guide_2024-07-15/config.yaml --outdir ../projects/.../03_summary
python3 positional_scan.py --crispresso-dir ... --config ... --outdir ...
python3 snp_audit.py       --crispresso-dir ... --config ... --outdir ...
```
