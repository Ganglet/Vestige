"""
Fine-tuning script for Phoenix-LM Phase 2 ablation.

Runs either the MLM baseline or DAM proposed fine-tuning run depending
on the config file passed. Both runs use identical hyperparameters —
only the data_collator differs.

Usage (from project root):
    python training/train.py training/config_mlm.yaml
    python training/train.py training/config_dam.yaml

Outputs (per run):
    training/checkpoints/<run_name>/          model checkpoints
    training/checkpoints/<run_name>/results.json  final eval metrics
"""
import json
import sys
import types
from pathlib import Path

import yaml

# triton is CUDA-only — unavailable on macOS. Both transformers and DNABERT-2
# check for it at import time. The stub must have __spec__ set or
# importlib.util.find_spec() raises ValueError.
import importlib.util, importlib.machinery
if importlib.util.find_spec("triton") is None:
    _stub = types.ModuleType("triton")
    _stub.__spec__ = importlib.machinery.ModuleSpec("triton", loader=None)
    _stub.__version__ = "0.0.0"
    sys.modules["triton"] = _stub

import torch
from datasets import load_from_disk
from transformers import (
    AutoModelForMaskedLM,
    AutoTokenizer,
    Trainer,
    TrainingArguments,
)

# Local imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from masking.collator_dam import DamageAwareDataCollator
from masking.collator_mlm import MLMCollator


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def build_collator(cfg: dict, tokenizer):
    if cfg["collator"] == "mlm":
        return MLMCollator(tokenizer=tokenizer, mlm=True, mlm_probability=0.15)
    elif cfg["collator"] == "dam":
        return DamageAwareDataCollator(
            tokenizer=tokenizer,
            damage_profile_path="damage/damage_profile.npy",
            baseline_prob=0.0,   # A/T never masked — only C/G positions trained
            scale_to=0.15,
        )
    raise ValueError(f"Unknown collator: {cfg['collator']}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python training/train.py <config.yaml>")
        sys.exit(1)

    cfg = load_config(sys.argv[1])
    print(f"\n=== Run: {cfg['run_name']} | Collator: {cfg['collator']} ===\n")

    # Device
    if torch.backends.mps.is_available():
        device = "mps"
    elif torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"
    print(f"Device: {device}")

    # Tokenizer + model
    tokenizer = AutoTokenizer.from_pretrained(cfg["model_name"], trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(cfg["model_name"], trust_remote_code=True)

    # Dataset
    dataset = load_from_disk(cfg["dataset_path"])
    dataset = dataset.remove_columns(["gene"])   # Trainer expects only tensor columns

    # Collator
    collator = build_collator(cfg, tokenizer)

    # Training arguments
    training_args = TrainingArguments(
        output_dir=cfg["output_dir"],
        run_name=cfg["run_name"],
        num_train_epochs=cfg["num_train_epochs"],
        per_device_train_batch_size=cfg["per_device_train_batch_size"],
        per_device_eval_batch_size=cfg["per_device_eval_batch_size"],
        learning_rate=cfg["learning_rate"],
        warmup_ratio=cfg["warmup_ratio"],
        weight_decay=cfg["weight_decay"],
        lr_scheduler_type=cfg["lr_scheduler_type"],
        eval_strategy=cfg["eval_strategy"],
        save_strategy=cfg["save_strategy"],
        load_best_model_at_end=cfg["load_best_model_at_end"],
        metric_for_best_model=cfg["metric_for_best_model"],
        seed=cfg["seed"],
        dataloader_num_workers=cfg["dataloader_num_workers"],
        logging_steps=cfg["logging_steps"],
        report_to=cfg["report_to"],
        save_safetensors=False,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=dataset["train"],
        eval_dataset=dataset["validation"],
        data_collator=collator,
    )

    trainer.train()

    # Save best checkpoint and eval results
    trainer.save_model(cfg["output_dir"])

    metrics = trainer.evaluate()
    results_path = Path(cfg["output_dir"]) / "results.json"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with open(results_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nFinal eval loss: {metrics['eval_loss']:.4f}")
    print(f"Checkpoint saved → {cfg['output_dir']}")
    print(f"Results   saved → {results_path}")


if __name__ == "__main__":
    main()
