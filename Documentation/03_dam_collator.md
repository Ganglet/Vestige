# Damage-Aware Masking — DataCollator Implementation

**Phase:** 2 — DAM Implementation & Fine-Tuning
**Status:** In Progress — collator implemented and unit-tested; dataset construction next
**Date:** April 2026

---

## Objective

Implement `DamageAwareDataCollator` as a subclass of HuggingFace `DataCollatorForLanguageModeling`, overriding the masking method to assign per-position, per-base masking probabilities derived from the mapDamage2 profiles produced in Phase 1. This is the primary methodological contribution of VESTIGE.

---

## 1. Files

| File | Role |
|------|------|
| `masking/collator_dam.py` | `DamageAwareDataCollator` — primary contribution |
| `masking/collator_mlm.py` | Baseline — re-exports `DataCollatorForLanguageModeling` as `MLMCollator` |

---

## 2. Masking Logic

The parent class `DataCollatorForLanguageModeling.torch_mask_tokens()` fills a `(batch, seq_len)` probability matrix uniformly with `mlm_probability` (default 0.15). `DamageAwareDataCollator` overrides this with a position- and base-conditioned matrix:

```
C at token position p from 5′ end  →  p_mask = ct_5p[ min(p, 69) ]
G at token position p from 3′ end  →  p_mask = ga_3p[ min(seq_len−1−p, 69) ]
A / T                               →  p_mask = baseline_prob  (default 0.03)
```

`ct_5p` and `ga_3p` are the mean damage frequency arrays from `damage_profile.npy` (70 positions each, both specimens averaged).

After building the raw probability matrix, it is uniformly scaled so E[masking rate] = `scale_to` (default 0.15), preserving the relative damage gradient while matching the MLM baseline training density. This makes the ablation a controlled comparison — the only variable between DAM and MLM is the masking distribution, not the total proportion of masked tokens.

The 80% → `[MASK]`, 10% → random token, 10% → unchanged schedule is kept identical to the base class.

---

## 3. Implementation Details

### `_build_base_lookup()`

Called once at initialisation. Iterates over the full tokenizer vocabulary and maps each `token_id` to an integer encoding `{0:other, 1:A, 2:C, 3:G, 4:T}`. BPE prefix markers (`##`, `▁`) are stripped before reading the first character. The result is a `vocab_size` uint8 tensor stored as `self._base_lookup`.

This avoids any per-token string decoding at training time.

### `_prob_matrix()`

Fully vectorized — no Python loops over positions or batch elements. The computation:

```python
pos     = torch.arange(seq_len)
pos_5p  = pos.clamp(max=69)
pos_3p  = (seq_len - 1 - pos).clamp(max=69)

base_ids = self._base_lookup[inputs]          # (batch, seq_len)
prob     = torch.full(..., self.baseline_prob)
prob     = torch.where(is_C, ct[pos_5p], prob)
prob     = torch.where(is_G, ga[pos_3p], prob)

# Optional scaling to target masking rate
if self.scale_to:
    prob = (prob * scale_to / prob.mean()).clamp(max=1.0)
```

### `scale_to` parameter

| `scale_to` | Behaviour |
|------------|-----------|
| `0.15` (default) | Scales probability matrix so E[masking] ≈ 15%. Relative damage pattern preserved. Use for the ablation study. |
| `None` | Raw damage frequencies used as-is (sparse masking ~1–2%). Biologically literal; useful for analysis but not for the ablation. |

The raw damage frequencies from our specimens peak at 0.35–0.41% at position 1 — far below the 15% MLM rate. Without scaling, DAM would see ~10× fewer masked tokens per batch, making a direct loss comparison invalid. Scaling corrects this.

---

## 4. Unit Test

**Claim:** for an all-cytosine input, masking density is significantly elevated at terminal positions relative to central positions.

**Method:** chi-squared test of independence on a 2×2 contingency table (terminal vs. central) × (masked vs. not-masked). 10,000 trials, `scale_to=None` (raw damage probabilities), synthetic exponentially decaying profile `ct_5p[i] = 0.05 × e^{−0.3i} + 0.001`.

Groups:
- Terminal: positions 1–5 (5 positions, highest expected deamination)
- Central: positions 21–50 (30 positions, expected background)

```
Terminal masking rate (pos  1– 5):  0.03098
Central  masking rate (pos 21–50):  0.00110
Ratio terminal / central:           28.1×
χ²(1) = 7154.8,  p = 0.00e+00
PASS: masking is significantly concentrated at terminal positions (p < 0.001) ✓
```

Assertion: `p < 0.001` and `terminal_rate > central_rate`. Run anytime:
```bash
python masking/collator_dam.py
```

---

## 5. Usage

```python
from transformers import AutoTokenizer
from masking.collator_dam import DamageAwareDataCollator
from masking.collator_mlm import MLMCollator

tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)

# DAM — proposed
dam_collator = DamageAwareDataCollator(
    tokenizer=tokenizer,
    damage_profile_path="damage/damage_profile.npy",
    baseline_prob=0.03,
    scale_to=0.15,
)

# MLM — baseline (identical HuggingFace collator, no subclassing needed)
mlm_collator = MLMCollator(tokenizer=tokenizer, mlm=True, mlm_probability=0.15)
```

Both collators are passed directly to `HuggingFace Trainer` via the `data_collator` argument.

---

## 6. Remaining Phase 2 Steps

- [ ] Construct fine-tuning dataset: extract CDS sequences for TRPV3, KCNK9, HBB from `elephant_ref.fa` using `gene_coords.csv`; tile into 512-token windows with 50-token stride; 80/20 train/val split
- [ ] Write `training/train.py` — HuggingFace Trainer fine-tuning script
- [ ] Write `training/config_mlm.yaml` and `training/config_dam.yaml`
- [ ] Fine-tune DNABERT-2: Run 1 (MLM baseline) → W&B run `mlm-baseline`
- [ ] Fine-tune DNABERT-2: Run 2 (DAM proposed) → W&B run `dam-proposed`
- **DELIVERABLE:** Two model checkpoints + W&B training curves + validation loss comparison figure

→ See `04_dataset_construction.md`
