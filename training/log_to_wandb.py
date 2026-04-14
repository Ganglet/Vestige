"""
Retroactively logs both fine-tuning runs to W&B from trainer_state.json files.

Produces the same dashboard as live logging would have — per-step train loss,
per-epoch eval loss, final hyperparameters, best metric, run config.

Run from project root AFTER `wandb login`:
    python training/log_to_wandb.py
"""
import json
from pathlib import Path
import wandb

WANDB_PROJECT = "vestige"

RUNS = [
    # MLM baseline already logged (run ID: vfdklmzb) — uncomment to re-log if needed
    # {
    #     "run_name":    "mlm-baseline",
    #     "state_path":  "training/checkpoints/mlm-baseline/checkpoint-700/trainer_state.json",
    #     "config_path": "training/config_mlm.yaml",
    #     "collator":    "DataCollatorForLanguageModeling (uniform 15%)",
    # },
    {
        # Corrected DAM run: baseline_prob=0, C/G-only scaling, retrained 2026-04-11
        "run_name":    "dam-proposed-corrected",
        "state_path":  "training/checkpoints/dam-proposed/checkpoint-700/trainer_state.json",
        "config_path": "training/config_dam.yaml",
        "collator":    "DamageAwareDataCollator (baseline_prob=0, C/G-only scaling, scale_to=0.15)",
    },
]

# Hyperparameters common to both runs
COMMON_CONFIG = {
    "model":                  "zhihan1996/DNABERT-2-117M",
    "num_train_epochs":       20,
    "per_device_train_batch_size": 8,
    "learning_rate":          2e-5,
    "warmup_ratio":           0.1,
    "weight_decay":           0.01,
    "lr_scheduler_type":      "cosine",
    "seed":                   42,
    "dataset":                "TRPV3/KCNK9/HBB CDS — Elephas maximus EanMak 1.0",
    "train_windows":          275,
    "val_windows":            69,
    "window_bp":              2000,
    "stride_bp":              200,
    "hardware":               "Apple M1 Pro (MPS)",
    "wall_clock_hours":       14.53,
    "baseline_prob":          0.0,
    "scale_to":               0.15,
    "note":                   "Corrected collator — original run (oun1npv4) had baseline_prob=0.03 bug",
}

import yaml

for run_meta in RUNS:
    with open(run_meta["state_path"]) as f:
        state = json.load(f)

    with open(run_meta["config_path"]) as f:
        yaml_cfg = yaml.safe_load(f)

    config = {**COMMON_CONFIG, "collator": run_meta["collator"], **yaml_cfg}

    run = wandb.init(
        project=WANDB_PROJECT,
        name=run_meta["run_name"],
        config=config,
        reinit=True,
    )

    # Log every entry in log_history — train loss at each logging step,
    # eval loss at each epoch boundary
    for entry in state["log_history"]:
        log_dict = {}
        step = entry.get("step")
        if "loss" in entry:
            log_dict["train/loss"]     = entry["loss"]
            log_dict["train/lr"]       = entry.get("learning_rate", 0)
            log_dict["train/grad_norm"] = entry.get("grad_norm", 0)
            log_dict["epoch"]          = entry.get("epoch")
        if "eval_loss" in entry:
            log_dict["eval/loss"]              = entry["eval_loss"]
            log_dict["eval/runtime"]           = entry.get("eval_runtime")
            log_dict["eval/samples_per_second"] = entry.get("eval_samples_per_second")
            log_dict["epoch"]                  = entry.get("epoch")
        if log_dict:
            wandb.log(log_dict, step=step)

    # Summary metrics
    wandb.run.summary["best_eval_loss"]       = state["best_metric"]
    wandb.run.summary["best_checkpoint"]      = state["best_model_checkpoint"]
    wandb.run.summary["final_eval_loss"]      = [
        e["eval_loss"] for e in state["log_history"] if "eval_loss" in e
    ][-1]
    wandb.run.summary["total_steps"]          = state["global_step"]
    wandb.run.summary["total_flos"]           = state.get("total_flos", 0)

    run.finish()
    print(f"Logged {run_meta['run_name']} → wandb.ai/{wandb.api.default_entity}/{WANDB_PROJECT}")
