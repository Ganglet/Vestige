# Phoenix-LM

**Damage-Aware Ancient DNA Infilling with Protein Stability Validation**

Angshuman Chakravertty · B.Tech CSE (Data Science), NMIMS Hyderabad  
Target venue: *Bioinformatics* (Oxford) · *Briefings in Bioinformatics* · ISMB/ECCB Workshop

---

## What This Is

Ancient DNA (aDNA) degrades in a precise, non-random pattern. Cytosines deaminate to thymine preferentially at fragment ends, guanines oxidize at the 3′ end, and hydrolysis shatters reads to 30–70 bp. Existing DNA language models (DNABERT-2, Nucleotide Transformer) were pre-trained with uniform random masking — they have no concept of this damage structure.

Phoenix-LM makes one core methodological change: replace random masking with **Damage-Aware Masking (DAM)** — a position-weighted masking scheduler whose probabilities are sampled directly from real mapDamage2 damage profiles. The hypothesis is that a model fine-tuned this way will reconstruct damaged ancient sequences more accurately and produce more structurally plausible proteins.

The reconstructed coding sequences are validated via ESMFold (pLDDT, TM-score vs. Asian Elephant reference). A lightweight biosecurity classifier screens all outputs against regulated sequences.

---

## Research Contributions

| # | Contribution | What Gets Built | Metric |
|---|-------------|-----------------|--------|
| 1 | **Damage-Aware Masking (DAM)** | Custom PyTorch `DataCollator` sampling from mapDamage2 profiles | Nucleotide recovery rate; BLEU-4 vs. MLM baseline |
| 2 | **Protein Stability Validation** | ESMFold pipeline: reconstructed CDS → 3D structure → pLDDT + TM-score | Mean pLDDT; TM-score vs. Asian Elephant reference protein |
| 3 | **Biosecurity Safety Layer** | 1D CNN classifier trained on NCBI Select Agent + benign sequences | AUC; false positive rate on mammoth cold-adaptation genes |

---

## Species & Target Genes

| Role | Species | Source |
|------|---------|--------|
| Ancient DNA | *Mammuthus primigenius* (Woolly Mammoth) | NCBI SRA · PRJEB44109 (van der Valk et al., 2021, *Nature*) |
| Reference genome | *Elephas maximus* (Asian Elephant) | NCBI RefSeq · GCF_024166365.1 (EanMak 1.0) |

**Target coding sequences:** TRPV3 (cold sensing), KCNK9 (temperature-gated ion channel), HBB (hemoglobin β-subunit — oxygen affinity at cold temperatures). These are the mammoth cold-adaptation genes with known functional significance.

---

## Tech Stack

| Layer | Tool |
|-------|------|
| Base model | DNABERT-2 (`zhihan1996/DNABERT-2-117M`) |
| Fine-tuning | HuggingFace Transformers + Trainer API |
| Damage profiling | mapDamage2 |
| Read alignment | BWA-aln (aDNA-tuned parameters) |
| Genome I/O | Biopython |
| Protein folding | ESMFold (`facebook/esmfold_v1`) |
| Structure similarity | TMalign |
| Experiment tracking | Weights & Biases |
| Safety classifier | PyTorch 1D CNN |

---

## Repository Structure

```
Phoenix-LM/
├── Dataset/
│   ├── mammoth/              # Raw reads, BAM, SRA files (gitignored — re-download from NCBI)
│   ├── elephant_reference/   # Reference FASTA + BWA indexes (gitignored — re-download from NCBI)
│   └── results/              # mapDamage2 output tables and plots (tracked)
│
├── damage/
│   ├── parse_profiles.py     # misincorporation.txt → damage_profile.npy
│   └── visualize_damage.py   # C→T / G→A frequency plots (Figure 1)
│
├── masking/
│   ├── collator_mlm.py       # Baseline: uniform random masking
│   └── collator_dam.py       # Core contribution: Damage-Aware Masking
│
├── training/
│   ├── config_mlm.yaml       # Baseline fine-tuning config
│   ├── config_dam.yaml       # DAM fine-tuning config
│   └── train.py              # HuggingFace Trainer fine-tuning script
│
├── evaluation/
│   ├── simulate_damage.py    # Apply damage profiles to modern CDS (eval set construction)
│   ├── evaluate_reconstruction.py  # Nucleotide recovery rate + BLEU-4
│   └── protein_validation.py       # ESMFold + TMalign pipeline
│
├── biosecurity/
│   ├── build_dataset.py      # Select Agent + benign sequence dataset
│   ├── train_classifier.py   # 1D CNN training
│   └── scan_reconstructions.py     # Screen all model outputs
│
├── results/
│   ├── figures/              # All paper figures
│   └── tables/               # CSV exports of all metric tables
│
├── Documentation/            # Step-by-step research documentation (journal-grade)
│
└── paper/
    └── phoenix_lm_draft.tex  # LaTeX manuscript
```

---

## Data Availability

Large genomic files are excluded from this repository (see `.gitignore`). To reproduce the dataset:

**Mammoth aDNA reads:**
```bash
# SRA Toolkit
prefetch ERR855944 ERR852028
fastq-dump --outdir Dataset/mammoth/ ERR855944

# Subsample to 5M reads (seed fixed for reproducibility)
seqtk sample -s 42 Dataset/mammoth/ERR855944.fastq 5000000 > Dataset/mammoth/ERR855944_5M.fastq
```

**Asian Elephant reference genome:**
```bash
datasets download genome accession GCF_024166365.1 --include genome,gff3 \
  --filename Dataset/elephant_reference/ncbi_dataset.zip
unzip Dataset/elephant_reference/ncbi_dataset.zip -d Dataset/elephant_reference/
```

**Alignment:**
```bash
bwa index Dataset/elephant_reference/elephant_ref.fa
bwa aln -l 16500 -n 0.01 -t 4 Dataset/elephant_reference/elephant_ref.fa \
  Dataset/mammoth/ERR855944_5M.fastq > Dataset/mammoth/mammoth_ERR855944.sai
bwa samse Dataset/elephant_reference/elephant_ref.fa \
  Dataset/mammoth/mammoth_ERR855944.sai Dataset/mammoth/ERR855944_5M.fastq \
  | samtools sort -@ 4 -o Dataset/mammoth/mammoth_ERR855944.bam
samtools index Dataset/mammoth/mammoth_ERR855944.bam
```

---

## Setup

```bash
conda create -n phoenix python=3.10
conda activate phoenix
pip install transformers datasets biopython mapdamage2 wandb sacrebleu torch
conda install -c bioconda samtools sra-tools
```

---

## Documentation

Step-by-step research logs are in `Documentation/`. Each file covers one pipeline stage with full command records, parameter justifications, and intermediate results — structured for journal methods reproducibility.

---

## Citation

> Chakravertty, A. (2026). Phoenix-LM: Damage-Aware Ancient DNA Infilling with Protein Stability Validation. *Manuscript in preparation.*

---

## License

Code: MIT  
Data: subject to NCBI SRA terms of use. See individual dataset accessions.
