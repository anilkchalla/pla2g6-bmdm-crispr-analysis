# data/

Raw reads are not stored in this repository. Download them here before running the pipeline.

```
data/
├── fastq/                    # raw reads, one file per sample
└── Y452B3_edit-analysis/     # vendor per-allele tables for the SO-C dataset
```

## Where to get the reads

All reads are deposited in the NCBI Sequence Read Archive under BioProject
**PRJNA1554390**. Download the following into `data/fastq/`:

| Dataset | File | Reads |
|---|---|---:|
| II | `Pla2g6.15.July.2024_S1_L001_R1_001.fastq` | 570,586 |
| II | `Pla2g6.15.July.2024_S1_L001_R2_001.fastq` | 570,586 |
| SO | `97K3PK_1_BMDM-Pla2g6.fastq` | 5,462 |
| SO-C | `Y452B3_1_Cas9_Ctrl_1.fastq` | 3,595 |
| SO-C | `Y452B3_2_Pla2g6_1.fastq` | 3,779 |
| SO-C | `Y452B3_3_Cas9_Ctrl_2.fastq` | 4,871 |
| SO-C | `Y452B3_4_Pla2g6_2.fastq` | 3,903 |

The scripts expect these exact filenames; `scripts/pla2g6_common.py` resolves them by sample stem.

The vendor per-allele tables consumed by `y452b3_reclassify_alleles.py` are provided as
supplementary material with the preprint.
