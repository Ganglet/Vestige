"""
Baseline collator — standard HuggingFace DataCollatorForLanguageModeling.
15% uniform random masking, no damage conditioning. Used as the ablation
control against DamageAwareDataCollator.

Usage:
    from masking.collator_mlm import MLMCollator
    collator = MLMCollator(tokenizer=tokenizer, mlm=True, mlm_probability=0.15)
"""
from transformers import DataCollatorForLanguageModeling as MLMCollator  # noqa: F401
