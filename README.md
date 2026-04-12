<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white" />
  <img src="https://img.shields.io/badge/PyTorch-2.x-EE4C2C?logo=pytorch&logoColor=white" />
  <img src="https://img.shields.io/badge/🤗_DNABERT--2-117M-FFD21E" />
  <img src="https://img.shields.io/badge/ESMFold-REST_API-6366f1" />
  <img src="https://img.shields.io/badge/mapDamage2-PMD_Profiling-16a34a" />
  <img src="https://img.shields.io/badge/W%26B-Experiment_Tracking-FFBE00?logo=weightsandbiases&logoColor=black" />
  <img src="https://img.shields.io/badge/License-MIT-22c55e" />
  <img src="https://img.shields.io/badge/Status-Manuscript_Prep-f97316" />
</p>

<h1 align="center">Phoenix-LM</h1>
<p align="center"><b>Damage-Aware Ancient DNA Infilling with Protein Stability Validation</b></p>
<p align="center">Angshuman Chakravertty · B.Tech CSE (Data Science), NMIMS Hyderabad</p>

---

## The Finding

Standard masked language model fine-tuning on ancient DNA is not just suboptimal — **at the innermost terminal positions (peak PMD zone), it performs worse than random chance.**

| Method | T_END = 3 nt | T_END = 5 nt | T_END = 10 nt |
|--------|-------------|-------------|--------------|
| Zero-shot DNABERT-2 | 15.5% | 19.2% | 22.8% |
| Random (chance) | 27.7% | 25.7% | 25.0% |
| MLM fine-tuning | 20.5% ⚠️ | 24.0% | 26.5% |
| **DAM (proposed)** | **30.8%** ✓ | **31.2%** ✓ | **32.2%** ✓ |

DAM is the only method that consistently exceeds random. All six terminal zone widths: **p < 0.001** (paired t-test, n = 626 windows across 7 genes).

---

## What This Is

Ancient DNA carries a forensic signature: cytosines deaminate to thymine preferentially at 5′ fragment ends; guanines oxidise at 3′ ends. This post-mortem damage (PMD) is non-random, position-dependent, and empirically characterised by mapDamage2. Every DNA language model trained before this work ignores this structure entirely — they apply uniform random masking during training, building no specialised reconstruction ability at the precise positions where ancient sequence information is most degraded.

**Phoenix-LM makes one principled change:** replace uniform masking with **Damage-Aware Masking (DAM)** — a custom PyTorch `DataCollator` that samples per-position masking probabilities directly from real mapDamage2 PMD profiles. The masking gradient mirrors the biological damage gradient. Fine-tuning on this objective produces a model that has seen disproportionately more training signal at the positions that matter most in ancient sequence reconstruction.

The reconstructed coding sequences are then validated at the protein level via ESMFold, establishing that neither model introduces fold-disrupting substitutions — reconstruction is biologically viable.

---

## Figure 1 — Empirical PMD Damage Profile (Two Woolly Mammoth Specimens)

![Figure 1 — PMD Damage Profile](results/figures/fig1_damage_profile.png)

C→T deamination peaks at position 1 from the 5′ end (~0.35–0.41% per specimen) and decays exponentially to background by position ~10. G→A at the 3′ end mirrors this pattern. ERR852028 (~44,800 BP) shows consistently higher terminal damage than ERR855944 (~4,300 BP), confirming the age-damage relationship. Both specimens converge to a shared background of 0.001 by position 10. The mean of these two profiles is what DAM's masking probabilities are sampled from.

---

## Figure 2 — Training Curves: MLM vs DAM

![Figure 2 — Training Curves](results/figures/fig2_training_curves.png)

DAM achieves **13% lower validation loss** on C/G-only masked positions (3.2736 vs 3.7568) — the model learns the damage grammar more efficiently when the masking distribution matches the damage distribution. Both models converge at epoch 17 (checkpoint-595).

---

## Figure 3 — Four-Baseline Terminal Zone Sensitivity

![Figure 3 — Terminal Zone Sensitivity](results/figures/fig3_damage_scaling.png)

The DAM advantage is **localized to the terminal damage zone and decays inward**, matching the shape of the PMD profile. Δ(DAM−MLM) ranges from +10.4 pp (T_END=3) to +4.2 pp (T_END=15), then stabilises. At every terminal zone width, the 95% bootstrap CI is entirely positive. Zero-shot DNABERT-2 is consistently the worst performer — fine-tuning is necessary; but only damage-aware fine-tuning actually works at the peak-damage zone.

### Figure 3b — Innermost Terminal Zone Spotlight (T_END = 3)

![Figure 3b — Headline](results/figures/fig3b_headline.png)

At the innermost 3 nt — where PMD probability peaks — MLM fine-tuning (20.5%) falls **below random chance** (27.7%). DAM (30.8%) is the only method above random. This is the central claim of the paper in one plot.

---

## Three-Part Validation

### 1 — Nucleotide Recovery (Phase 3)

Evaluation dataset: 626 windows across 7 genes — TRPV3, KCNK9, HBB (training genes, validation split) + TRPA1, UCP1, ADRB3, FASN (held-out genes, novel evaluation). Four baselines evaluated at T_END ∈ {3, 5, 10, 15, 20, 25}.

**At background positions (authentic PMD rates, n=43 windows):** MLM 65.7% vs DAM 65.1%, p = 0.892. No difference — the models are equivalent in the background zone, confirming the advantage is terminal-specific, not a generic artefact.

**At terminal positions:** DAM > MLM at all six T_END values, p < 0.001. The effect generalizes fully to four new genes not present in training.

### 2 — Per-Position Decomposition (Phase 3 ext)

![Per-Position Analysis](results/figures/fig_per_position.png)

Fine-grained per-nucleotide analysis (d = 1..25 from each end, 5′-C and 3′-G strands separated). At d = 1 from the 5′ end — the single nucleotide with highest PMD probability — DAM reconstructs 5/16 cytosines correctly vs 1/16 for MLM (5× improvement). The gradient decays with distance, matching the mapDamage2 C→T profile.

### 3 — Protein Structural Validation (Phase 4)

ESMFold REST API used to fold all reconstructed protein sequences for TRPV3, KCNK9, and HBB. TM-score and Cα-RMSD computed via pure-Python Kabsch superposition (Zhang & Skolnick 2004).

| Gene | Method | pLDDT (ref) | pLDDT (recon) | TM-score | Cα-RMSD |
|------|--------|-------------|---------------|----------|---------|
| TRPV3† | MLM | 77.37 | 77.44 | 0.9808 | 1.092 Å |
| TRPV3† | DAM | 77.37 | 77.44 | 0.9808 | 1.092 Å |
| KCNK9 | MLM | 67.94 | 67.02 | 0.9520 | 1.732 Å |
| KCNK9 | DAM | 67.94 | 67.02 | 0.9520 | 1.732 Å |
| HBB | MLM | 79.70 | 80.78 | 0.9728 | 0.815 Å |
| HBB | DAM | 79.70 | 81.34 | 0.9714 | 0.840 Å |

† TRPV3 full-length (791 aa) truncated to N-terminal ankyrin repeat domain (aa 1–400; ESMFold API limit).

All TM-scores > 0.95. All substitutions are conservative (R↔K, E↔D, S↔T type). Neither model introduces fold-disrupting mutations — reconstruction is biologically viable regardless of which model is used.

---

## Pipeline

```
Woolly Mammoth aDNA (SRA: ERP008929)          Asian Elephant genome (GCF_024166365.1)
   ERR855944 · ERR852028                            elephant_ref.fa
          │                                               │
     BWA-aln (aDNA params)                         gene_coords.csv
          │                                    (TRPV3, KCNK9, HBB, TRPA1, UCP1, ADRB3, FASN)
     mammoth.bam                                          │
          │                                         build_dataset.py
     mapDamage2                                           │
          │                                        344 windows (train/val split)
   damage_profile.npy                                     │
          │                                               │
          └──────────────────────────────────────────┐   │
                                                     ▼   ▼
                                        ┌─────────────────────────────┐
                                        │    DNABERT-2 Fine-tuning    │
                                        │  ┌─────────┐  ┌─────────┐  │
                                        │  │   MLM   │  │   DAM   │  │
                                        │  │ uniform │  │ PMD-    │  │
                                        │  │ masking │  │ weighted│  │
                                        │  └─────────┘  └─────────┘  │
                                        └─────────────────────────────┘
                                                     │
                        ┌────────────────────────────┼────────────────────────┐
                        ▼                            ▼                        ▼
               Nucleotide recovery           T_END sensitivity          ESMFold (REST API)
               at authentic PMD              analysis + baselines        pLDDT · TM-score
               (Table 1 + 2)                 (Figure 3, 626 windows)     (Table 3)
```

---

## Key Implementation Details

**DamageAwareDataCollator** (`masking/collator_dam.py`): The core contribution. Computes a per-position, per-base-type masking probability matrix from the mapDamage2 output. The matrix is rescaled so C/G positions average 15% masking (equal total density to MLM), preserving the damage gradient as the only difference between the two training objectives. A/T tokens are never masked — biologically correct, since deamination is C/G-specific.

**BPE cascade avoidance** (`evaluation/simulate_damage.py`): DNABERT-2 uses byte-pair encoding. Substituting one nucleotide changes token boundaries across a large surrounding context (tested: 2 substitutions → 31 changed token positions). The evaluation uses `offset_mapping` from the original tokenization to locate which token spans each damaged nucleotide, then masks only that token in `gt_ids`. The damaged sequence is never re-tokenized.

**Pure-Python structural comparison** (`protein/tmalign_compare.py`): TM-score (Zhang & Skolnick 2004) and Cα-RMSD computed via Kabsch superposition — no TMalign binary dependency. Validated against published TM-score values for known structure pairs.

---

## Species

| Role | Species | Accession |
|------|---------|-----------|
| Ancient DNA | *Mammuthus primigenius* (Woolly Mammoth) | NCBI SRA · ERP008929 (Palkopoulou et al., 2015) |
| Reference | *Elephas maximus* (Asian Elephant) | NCBI RefSeq · GCF_024166365.1 (EanMak 1.0) |

**Training genes:** TRPV3 (cold sensing), KCNK9 (temperature-gated ion channel), HBB (hemoglobin β-subunit)  
**Evaluation-only genes:** TRPA1, UCP1, ADRB3, FASN — all cold-adaptation relevant, none seen during training

---

## Repository Structure

```
Phoenix-LM/
├── damage/
│   ├── parse_profiles.py          # misincorporation.txt → damage_profile.npy
│   ├── visualize_damage.py        # Figure 1
│   └── gene_coords.csv            # Genomic coordinates for all 7 genes
│
├── masking/
│   ├── collator_mlm.py            # Baseline: uniform random masking
│   └── collator_dam.py            # Core contribution: Damage-Aware Masking
│
├── training/
│   ├── build_dataset.py           # 2000bp sliding window dataset (344 windows)
│   ├── train.py                   # HuggingFace Trainer fine-tuning
│   ├── config_mlm.yaml
│   └── config_dam.yaml
│
├── evaluation/
│   ├── simulate_damage.py         # PMD-rate damage simulation (BPE-safe)
│   ├── evaluate_reconstruction.py # Background-site recovery + BLEU-4
│   ├── evaluate_terminal.py       # Terminal-zone recovery (T_END=10)
│   ├── evaluate_scaling.py        # T_END sweep + 4-baseline comparison
│   ├── evaluate_per_position.py   # Per-position d=1..25 decomposition
│   └── expand_validation.py       # Fetch extra genes + build held-out eval set
│
├── protein/
│   ├── fetch_cds.py               # NCBI RefSeq mRNA CDS via Entrez elink
│   ├── damage_reconstruct_cds.py  # CDS damage + reconstruction + translation
│   ├── run_esmfold.py             # ESMFold REST API → 9 PDB files
│   └── tmalign_compare.py         # Kabsch + TM-score (pure Python)
│
├── results/figures/               # All paper figures (PDF + PNG)
├── Documentation/                 # Phase-by-phase research logs (08 files)
└── Dataset/                       # Gitignored: BAM, reference FASTA
```

---

## Reproduce

```bash
# Environment
conda create -n phoenix python=3.10 && conda activate phoenix
pip install transformers datasets biopython wandb sacrebleu torch scipy matplotlib
conda install -c bioconda samtools sra-tools mapdamage2

# Download mammoth reads (SRA)
prefetch ERR855944 ERR852028
fastq-dump --outdir Dataset/mammoth/ ERR855944 ERR852028

# Download Asian Elephant reference
datasets download genome accession GCF_024166365.1 --include genome \
  --filename Dataset/elephant_reference/ncbi_dataset.zip

# Phase 1 — align + damage profile
bwa index Dataset/elephant_reference/elephant_ref.fa
bwa aln -l 16500 -n 0.01 -t 4 Dataset/elephant_reference/elephant_ref.fa \
  Dataset/mammoth/ERR855944.fastq | bwa samse ... | samtools sort -o mammoth_ERR855944.bam
mapDamage2 -i Dataset/mammoth/mammoth_ERR855944.bam \
  -r Dataset/elephant_reference/elephant_ref.fa \
  --folder Dataset/results/results_mammoth_ERR855944/
python3 damage/parse_profiles.py

# Phase 2 — fine-tune both models (~11–14h each on Apple M1 Pro)
python3 training/build_dataset.py
python3 training/train.py --config training/config_mlm.yaml
python3 training/train.py --config training/config_dam.yaml

# Phase 3 — evaluation
python3 evaluation/simulate_damage.py
python3 evaluation/evaluate_reconstruction.py      # Table 1
python3 evaluation/evaluate_terminal.py            # Table 2
NCBI_EMAIL=you@email.com python3 evaluation/expand_validation.py  # held-out genes
python3 evaluation/evaluate_scaling.py             # Figure 3

# Phase 3 ext — per-position decomposition
python3 evaluation/evaluate_per_position.py

# Phase 4 — protein structural validation
NCBI_EMAIL=you@email.com python3 protein/fetch_cds.py
python3 protein/damage_reconstruct_cds.py
python3 protein/run_esmfold.py
python3 protein/tmalign_compare.py                 # Table 3
```

---

## Citation

```bibtex
@misc{chakravertty2026phoenix,
  author = {Chakravertty, Angshuman},
  title  = {Phoenix-LM: Damage-Aware Ancient DNA Infilling with Protein Stability Validation},
  year   = {2026},
  note   = {Manuscript in preparation}
}
```

---

## License

Code: MIT · Data: subject to NCBI SRA terms of use (accessions ERP008929, GCF_024166365.1)
