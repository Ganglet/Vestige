# Phase 3 — In Silico Damage Simulation & Reconstruction Evaluation

**Phase:** 3 — Evaluation
**Status:** Complete
**Date:** April 2026

---

## Objective

Measure how well each model (MLM baseline, DAM proposed) recovers the original sequence from artificially degraded inputs. The degradation uses the empirical PMD frequencies from `damage_profile.npy` — the same profiles that informed DAM's masking probabilities during training.

This constitutes Table 1 of the paper.

---

## 1. Evaluation Pipeline

```
damage_profile.npy
        │
        ▼
simulate_damage.py          ← nucleotide-level C→T / G→A substitution
        │ evaluation/damaged_validation.npy
        ▼
evaluate_reconstruction.py  ← mask changed tokens → run model → decode → compare
        │
        ├── evaluation/results_table1.json
        └── evaluation/table1.txt
```

---

## 2. In Silico Damage Simulation (`evaluation/simulate_damage.py`)

**Key design decision:** DNABERT-2 uses BPE tokenization (tokens of 1–6 nucleotides). Damage must be applied at the nucleotide level, not the token level; otherwise substituting a single character changes the token ID in a way that breaks the ground-truth comparison.

**Pipeline per window:**
1. Decode `input_ids` → raw DNA string (strip [CLS]/[SEP])
2. For each nucleotide at 5′ distance `d`: apply C→T at probability `ct_5p[d]`
3. For each nucleotide at 3′ distance `d`: apply G→A at probability `ga_3p[d]`
4. Re-tokenize both the original and damaged DNA strings (max 512 tokens, same as training)
5. Store nucleotide-level damage flags alongside both token sequences

**Output stats (seed=42):**
- 69 validation windows processed
- 79 total nucleotides damaged
- 1.14 damaged nucleotides per window (mean)
- 43/69 windows with ≥1 damaged nucleotide

The sparse damage (1.14 positions/window) reflects the authentic biological PMD rates — peak C→T is ~0.35% at position 1, declining exponentially. This is not inflated for the evaluation.

---

## 3. Reconstruction Evaluation (`evaluation/evaluate_reconstruction.py`)

**Inference procedure per window:**
1. Find token positions where the damaged tokenization differs from ground truth (`gt_ids[i] ≠ dam_ids[i]` at non-padding positions)
2. Replace those positions with [MASK] in the damaged sequence
3. Run the model; take argmax at masked positions
4. Decode the reconstructed token IDs → DNA string
5. At each nucleotide position flagged as damaged, compare the reconstructed character to the ground-truth character

**Metrics:**
- **Nucleotide recovery rate:** fraction of damaged nucleotide positions where the model's decoded output matches the original base
- **BLEU-4:** corpus-level, character-by-character on the full decoded sequence (ground-truth as reference)
- **Paired t-test:** recovery rates across the 43 windows with shared damage
- **Bootstrap 95% CI:** 10,000 resamples of the paired difference

---

## 4. Results

### Table 1 — Background damage site reconstruction (authentic PMD rates)

43/69 windows had ≥1 damaged nucleotide (mean 1.14 per window, all at background positions >9 nt from either end).

| Metric | MLM Baseline | DAM Proposed |
|--------|-------------|-------------|
| Nucleotide recovery rate (mean) | 65.7% | 65.1% |
| Nucleotide recovery rate (SD) | 39.5% | 39.3% |
| BLEU-4 (corpus, char-level) | 99.86% | 99.84% |

Paired t-test: t = −0.137, p = 0.892 — no difference. Expected: both models receive identical C/G masking density at background positions.

### Table 2 — Terminal-position reconstruction (first/last 10 nt from each end, 325 sites, 69 windows)

| Metric | MLM Baseline | DAM Proposed |
|--------|-------------|-------------|
| Nucleotide recovery (mean per window) | 34.1% | **38.6%** |
| Nucleotide recovery (SD) | 27.7% | 27.2% |
| Aggregate (correct / 325 total sites) | 33.5% | **38.2%** |

Paired t-test: t = 2.086, **p = 0.041**. Bootstrap 95% CI (DAM − MLM): [+0.004, +0.087] — entirely positive.

DAM significantly outperforms MLM at terminal positions. 325 sites = all terminal C/G nucleotides across 69 windows; multiple nucleotides within the same BPE token are evaluated independently at their sub-token offsets.

---

## 5. Interpretation

**Training loss:** DAM 3.2736 vs MLM 3.7568 — 13% lower on C/G-only masked positions.

**Reconstruction at terminal positions (p=0.041):** DAM recovers 38.2% vs MLM 33.5% — a 4.7 percentage-point advantage, CI entirely positive. This is the biologically motivated result: DAM concentrates masking at terminal C/G positions during training (28% at position 1 vs MLM's flat 15%), acquiring more specialised reconstruction signal there.

**Reconstruction at background positions (p=0.892):** No difference — confirms the advantage is specific to the terminal zone, not a general artefact.

**For the paper:**
- *Primary claim:* DAM achieves 13% lower validation MLM loss on C/G positions — more efficient learning of the damage grammar.
- *Secondary claim:* DAM significantly outperforms MLM at terminal C/G reconstruction (38.2% vs 33.5%, p=0.041), the biologically correct evaluation domain.
- *Specificity control:* No difference at background positions (p=0.892) — the effect is domain-specific.

---

## 6. Files Generated

| File | Description |
|------|-------------|
| `evaluation/simulate_damage.py` | Nucleotide-level PMD simulation using offset_mapping |
| `evaluation/evaluate_reconstruction.py` | Background-site inference, recovery rate, BLEU-4, t-test |
| `evaluation/evaluate_terminal.py` | Terminal-position inference, recovery rate, t-test |
| `evaluation/damaged_validation.npy` | Paired (ground-truth, damaged) dataset, 69 windows |
| `evaluation/results_table1.json` | Background site results |
| `evaluation/results_terminal.json` | Terminal position results |
| `evaluation/table1.txt` | Table 1 (background) |
| `evaluation/table_terminal.txt` | Table 2 (terminal positions) |

---

## 7. Phase 3 Complete

- [x] Implement `simulate_damage.py` (nucleotide-level, offset_mapping)
- [x] Implement `evaluate_reconstruction.py` (background damage sites)
- [x] Implement `evaluate_terminal.py` (terminal C/G positions)
- [x] Retrain DAM with corrected collator (baseline_prob=0, C/G-only scaling)
- [x] Generate Table 1 and Table 2

→ See `07_esm_validation.md` (Phase 4 — ESMFold structural validation)
