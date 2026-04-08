# Fine-Tuning Dataset Construction

**Phase:** 2 — DAM Implementation & Fine-Tuning
**Status:** Complete
**Date:** April 2026

---

## Objective

Extract CDS sequences for TRPV3, KCNK9, and HBB from the Asian Elephant reference genome, tile into overlapping windows compatible with DNABERT-2's 512-token context limit, and produce a train/validation split for both fine-tuning runs (MLM baseline and DAM proposed).

**Script:** `training/build_dataset.py`
**Output:** `training/dataset/` (HuggingFace `DatasetDict`), `training/dataset_stats.txt`

---

## 1. Sequence Extraction

Coordinates from `damage/gene_coords.csv` (produced in Phase 1 via NCBI Entrez). Sequences extracted from `Dataset/elephant_reference/elephant_ref.fa` using Biopython `SeqIO.index()` (1-based inclusive slicing). Minus-strand genes (HBB) are reverse-complemented before tiling.

**Soft-masking issue:** `elephant_ref.fa` uses lowercase letters for RepeatMasker soft-masked repeats. DNABERT-2's BPE tokenizer has no lowercase vocabulary — lowercase input triggers character-level fallback and approximately doubles the token count per window (tested: 2000 bp uppercase → ~400 tokens; 2000 bp with 31% lowercase → ~910 tokens). All sequences uppercased before tokenization.

---

## 2. Tiling Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `WINDOW_BP` | 2,000 | Tokenizes to ~400 tokens across all three genes — within the 512-token limit with headroom |
| `STRIDE_BP` | 200 | 10× overlap; adjacent windows share ~1,800 bp (~360 tokens) of context |
| `MAX_TOKENS` | 512 | DNABERT-2 positional embedding limit |

Tiling is done in bp rather than token space because `DamageAwareDataCollator` uses in-window token positions as proxies for fragment-terminus distance. Tiling in bp keeps those positions biologically interpretable.

---

## 3. Dataset Statistics

```
Gene       CDS bp  Windows Filtered     Kept
----------------------------------------------
TRPV3      46,944      226        0      226
KCNK9      21,014       97        0       97
HBB         5,821       21        0       21
----------------------------------------------
TOTAL                                    344

Train windows : 275  (80%)
Val   windows :  69  (20%)
```

Zero windows exceeded 512 tokens after uppercasing. All three genes are represented in both splits (verified by checking unique gene labels post-shuffle).

---

## 4. Dataset Format

Saved as a HuggingFace `DatasetDict` to `training/dataset/`. Each example has four fields:

| Field | Type | Description |
|-------|------|-------------|
| `gene` | string | Source gene label (`TRPV3`, `KCNK9`, `HBB`) |
| `input_ids` | list[int] | Token IDs, padded to 512 |
| `attention_mask` | list[int] | 1 for real tokens, 0 for padding |
| `token_type_ids` | list[int] | All zeros (single-sequence input) |

**Load:**
```python
from datasets import load_from_disk
dataset = load_from_disk("training/dataset")
train_ds = dataset["train"]       # 275 examples
val_ds   = dataset["validation"]  # 69 examples
```

**Pass to Trainer:**
```python
from training.build_dataset import main  # or load directly
# Both collators receive the same dataset — only data_collator differs
trainer = Trainer(
    model=model,
    train_dataset=dataset["train"],
    eval_dataset=dataset["validation"],
    data_collator=dam_collator,   # or mlm_collator for the baseline run
    ...
)
```

---

## 5. Remaining Phase 2 Steps

- [x] Implement `DamageAwareDataCollator` — see `03_dam_collator.md`
- [x] Construct fine-tuning dataset — this document
- [ ] Write `training/train.py` — HuggingFace Trainer fine-tuning script
- [ ] Write `training/config_mlm.yaml` and `training/config_dam.yaml`
- [ ] Fine-tune DNABERT-2: Run 1 (MLM baseline) → W&B run `mlm-baseline`
- [ ] Fine-tune DNABERT-2: Run 2 (DAM proposed) → W&B run `dam-proposed`
- **DELIVERABLE:** Two model checkpoints + W&B training curves + validation loss comparison

→ See `05_finetuning.md`
