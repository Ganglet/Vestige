# Fine-Tuning Runs — MLM Baseline & DAM Proposed

**Phase:** 2 — DAM Implementation & Fine-Tuning
**Status:** MLM baseline complete; DAM proposed pending
**Date:** April 2026

---

## Objective

Run two controlled fine-tuning experiments on DNABERT-2 using the TRPV3/KCNK9/HBB dataset. The only variable between the two runs is the data collator — everything else (model, dataset, hyperparameters, random seed) is identical. This is the ablation study that constitutes the core experimental claim of the paper.

---

## 1. Experimental Setup

**Model:** `zhihan1996/DNABERT-2-117M` — 117M parameter BERT-style encoder, pre-trained on multi-species genomic sequences

**Dataset:** `training/dataset/` — 344 windows (275 train / 69 val) tiled from TRPV3, KCNK9, HBB CDS sequences (see `04_dataset_construction.md`)

**Hardware:** Apple M1 Pro (MPS backend) — ~29s/step at batch size 8; actual wall-clock time per run: **11 hours 24 minutes**

### Hyperparameters (identical for both runs)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Epochs | 20 | Sufficient for visible plateau on a 275-example dataset |
| Batch size | 8 | Stable on M1 Pro 16GB unified memory |
| Learning rate | 2×10⁻⁵ | Standard BERT fine-tuning range |
| LR schedule | Cosine with 10% warmup | Prevents early instability |
| Weight decay | 0.01 | Light regularization |
| `save_total_limit` | 3 | Avoids ~9 GB disk usage from 20 checkpoints |
| Seed | 42 | Fixed for reproducibility |

### The only difference

| | Run 1 | Run 2 |
|--|-------|-------|
| **Config** | `config_mlm.yaml` | `config_dam.yaml` |
| **Collator** | `MLMCollator` — 15% uniform random masking | `DamageAwareDataCollator` — position/base-conditioned masking, scaled to 15% mean |
| **W&B run** | `mlm-baseline` | `dam-proposed` |

---

## 2. MLM Baseline Results

**Status:** Complete
**Checkpoint:** `training/checkpoints/mlm-baseline/checkpoint-595` (best)
**Final eval loss (epoch 20):** 3.8350
**Best eval loss:** 3.7568 (epoch 17)

### Validation loss per epoch

| Epoch | Eval Loss | Epoch | Eval Loss |
|-------|-----------|-------|-----------|
| 1 | 5.3914 | 11 | 4.1006 |
| 2 | 5.1830 | 12 | 4.0240 |
| 3 | 5.0984 | 13 | 3.9608 |
| 4 | 4.9974 | 14 | 3.8005 |
| 5 | 4.8330 | 15 | 3.8815 |
| 6 | 4.6699 | 16 | 3.8683 |
| 7 | 4.5055 | 17 | **3.7568** |
| 8 | 4.4364 | 18 | 3.8777 |
| 9 | 4.2395 | 19 | 3.9107 |
| 10 | 4.1701 | 20 | 3.8350 |

### Convergence assessment

- Loss drops steeply from 5.39 → 4.02 over epochs 1–12 (~1.37 reduction)
- Plateau begins at epoch 14 with oscillation of ±0.04 through epoch 20
- Loss delta over last 3 epochs: **0.043** — model has converged
- No sign of overfitting: eval loss tracks training loss throughout; no divergence

The 20-epoch run is well-justified — the plateau is clearly visible and the curve will be defensible to reviewers.

**Figure 2:** `results/figures/fig2_training_curves.pdf` — training and validation loss curves (DAM curve will be added after the second run completes)

---

## 3. DAM Proposed — Pending

Run the DAM fine-tuning after the MLM baseline:

```bash
python3 training/train.py training/config_dam.yaml
```

Results will be added to this document and the training curves figure will be regenerated with both curves overlaid.

---

## 4. Dependency Issues Encountered

See `problems_and_decisions.md` P12 for the full account. Summary:
- `tf-keras` had to be installed separately (`pip install tf-keras`)
- `accelerate` was not auto-installed by transformers (`pip install accelerate`)
- `use_mps_device` argument removed from `TrainingArguments` in recent transformers — Trainer detects MPS automatically
- `save_safetensors=False` required due to DNABERT-2's tied embedding weights (standard BERT architecture)
- All fixes are in `requirements.txt` and `train.py`

---

## 5. Remaining Phase 2 Steps

- [x] Implement collators — `03_dam_collator.md`
- [x] Build dataset — `04_dataset_construction.md`
- [x] Write `train.py`, `config_mlm.yaml`, `config_dam.yaml`
- [x] Run 1 — MLM baseline (20 epochs, best eval loss 3.7568)
- [ ] Run 2 — DAM proposed
- [ ] Regenerate Figure 2 with both curves
- **DELIVERABLE:** Two model checkpoints + Figure 2 → Phase 3 can begin

→ See `06_evaluation.md`
