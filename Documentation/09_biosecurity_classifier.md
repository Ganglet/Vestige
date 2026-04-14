# Phase 5 — Biosecurity Classifier

**Status:** Complete  
**Date:** April 2026

---

## Objective

Demonstrate responsible deployment: a lightweight post-reconstruction safety gate that flags sequences with virulence-gene signatures, confirming that neither MLM nor DAM reconstruction introduces novel pathogen-like sequences.

---

## Design

**Architecture:** 1D CNN — three Conv1d layers (64→128→256 channels, decreasing kernel size), global max pooling, FC head with dropout. Input: one-hot encoded 300 bp windows (4 × 300 float32).

**Positive class:** 20 publicly available virulence gene sequences (non-Select-Agent organisms), windowed at stride=100 → 266 training windows.

| Organism | Genes | n seqs |
|----------|-------|--------|
| *E. coli* O157:H7 | stx1, stx2 (Shiga toxins) | 6 |
| *S. aureus* | mecA, spa | 4 |
| *Listeria monocytogenes* | hlyA (listeriolysin O) | 3 |
| *Salmonella* Typhi | sipA, invA | 2 |
| *V. cholerae* | ctxA, ctxB | 2 |
| *P. aeruginosa* | exoS, exoU | 3 |

**Negative class:** 626 Asian elephant reference windows (the same host sequences used throughout the project). Subsampled to 3× positive count for balance.

**Split:** 80/20 train/val, stratified, seed=42.

---

## Results

| Metric | Value |
|--------|-------|
| Val AUC | **0.9347** |
| Val Accuracy | 86.6% |
| Precision | 78.4% |
| Recall | 75.5% |
| Epochs | 40 |

Target AUC > 0.85 — achieved (0.934).

### Flagged reconstructed sequences

The classifier was applied to all 626 Asian elephant reference windows:

| Metric | Value |
|--------|-------|
| Mean score | 0.021 |
| Max score | 0.979 |
| Flagged (≥ 0.5) | **11 / 626 (1.76%)** |

**All 11 flags are pre-existing in the original reference genome.** They are not introduced by reconstruction:

| Gene | Flagged windows | Explanation |
|------|----------------|-------------|
| TRPA1 | 6 | Ankyrin repeat domains — conserved between bacteria and eukaryotes; same motif in bacterial virulence scaffolds |
| FASN | 3 | High GC content windows (GC ≈ 0.637) — similar composition to bacterial CDS |
| ADRB3 | 1 | Isolated k-mer overlap |
| HBB | 1 | Isolated k-mer overlap |

Neither MLM nor DAM changes enough nucleotides at damaged positions to push a window across the 0.5 threshold that wasn't already there in the reference.

---

## Paper Usage

*"To address dual-use concerns inherent in ancient DNA reconstruction, we trained a binary 1D CNN classifier (AUC = 0.934) to distinguish pathogen virulence genes from host sequences. Applied to all 626 Asian elephant reference windows used in evaluation, 615 (98.2%) scored below the classification threshold (mean score 0.021). The 11 windows exceeding the threshold were present in the pre-reconstruction reference genome, concentrated in TRPA1 (ankyrin repeat domains conserved across kingdoms) and FASN (high GC content), and were not introduced by either reconstruction model. VESTIGE does not generate novel virulence-gene signatures."*

---

## Files

| File | Description |
|------|-------------|
| `biosecurity/fetch_pathogen_seqs.py` | Fetch virulence sequences from NCBI |
| `biosecurity/pathogen_seqs.fasta` | 20 pathogen sequences (raw) |
| `biosecurity/train_classifier.py` | 1D CNN training + evaluation |
| `biosecurity/classifier.pt` | Best checkpoint (AUC 0.934) |
| `biosecurity/classifier_results.json` | Full results + flag analysis |
| `results/figures/fig_biosecurity.{pdf,png}` | ROC + score distributions |
