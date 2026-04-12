"""
Log all evaluation results to W&B — Phoenix-LM.

Reads from existing JSON files (no model inference needed).
Creates a single W&B run with:
  - Phase 2: training loss comparison (Table)
  - Phase 3: background + terminal reconstruction (Tables 1 & 2)
  - Phase 3 ext: T_END sensitivity sweep (line chart — 4 baselines, Figure 3)
  - Phase 3 ext: per-position d=1..25 decomposition (line chart)
  - Phase 4: ESMFold pLDDT + TM-score (Table 3)
  - Figures: uploads all PNGs from results/figures/

Run from project root:
    python3 results/log_eval_to_wandb.py
"""
import json
import math
import wandb
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
TABLE1_JSON       = Path("evaluation/results_table1.json")
TERMINAL_JSON     = Path("evaluation/results_terminal.json")
SCALING_JSON      = Path("evaluation/scaling_results.json")
PER_POS_JSON      = Path("evaluation/per_position_results.json")
TABLE3_JSON       = Path("protein/table3.json")
FIGURES_DIR       = Path("results/figures")
# ───────────────────────────────────────────────────────────────────────────────


def safe(v):
    """Replace NaN/None with None (W&B accepts None, not float NaN)."""
    if v is None:
        return None
    try:
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return v


def pct(v):
    """Format as percentage string."""
    return f"{v:.1%}" if v is not None else "—"


def main():
    wandb.init(
        project="phoenix-lm",
        name="evaluation-results",
        notes=(
            "Full evaluation results for Phoenix-LM: "
            "DAM vs MLM fine-tuning on ancient DNA reconstruction. "
            "Includes T_END sensitivity sweep (626 windows, 7 genes, 4 baselines), "
            "per-position decomposition, and ESMFold structural validation."
        ),
        tags=["aDNA", "DNABERT-2", "DAM", "evaluation", "ancient-dna"],
        config={
            "base_model":       "zhihan1996/DNABERT-2-117M",
            "training_genes":   ["TRPV3", "KCNK9", "HBB"],
            "eval_genes":       ["TRPV3", "KCNK9", "HBB", "TRPA1", "UCP1", "ADRB3", "FASN"],
            "ancient_specimen": "Mammuthus primigenius (ERP008929)",
            "reference":        "Elephas maximus (GCF_024166365.1)",
            "window_bp":        2000,
            "stride_bp":        200,
            "max_tokens":       512,
            "val_fraction":     0.20,
            "n_train_windows":  275,
            "n_val_windows":    69,
            "n_extra_windows":  557,
            "t_end_values":     [3, 5, 10, 15, 20, 25],
            "seed":             42,
        },
    )

    # Line-chart metrics need independent x-axes
    wandb.define_metric("t_end_sweep/*", step_metric="t_end_sweep/t_end")
    wandb.define_metric("per_position/*", step_metric="per_position/d")

    # ── Phase 2 — training loss comparison (Table) ────────────────────────────
    p2_table = wandb.Table(columns=["method", "val_loss", "improvement_pct", "best_checkpoint"])
    p2_table.add_data("MLM baseline", 3.7568, 0.0,  595)
    p2_table.add_data("DAM proposed", 3.2736, 13.0, 595)
    wandb.log({"phase2/training_loss_comparison": p2_table})
    print("✓ Phase 2 logged")

    # ── Phase 3 — background reconstruction (Table 1) ─────────────────────────
    t1 = json.loads(TABLE1_JSON.read_text())
    ps = t1.get("paired_stats", {})

    bg_table = wandb.Table(columns=[
        "method", "recovery_mean", "bleu4", "valid_windows", "p_value", "mean_diff_DAM_minus_MLM"
    ])
    bg_table.add_data(
        "MLM baseline",
        safe(t1["mlm_baseline"]["mean_recovery"]),
        safe(t1["mlm_baseline"]["bleu4"]),
        t1["mlm_baseline"]["valid_windows"],
        safe(ps.get("p_value")),
        None,
    )
    bg_table.add_data(
        "DAM proposed",
        safe(t1["dam_proposed"]["mean_recovery"]),
        safe(t1["dam_proposed"]["bleu4"]),
        t1["dam_proposed"]["valid_windows"],
        safe(ps.get("p_value")),
        safe(ps.get("mean_diff_dam_minus_mlm")),
    )
    wandb.log({"phase3_background/table1_reconstruction": bg_table})
    print("✓ Phase 3 background logged")

    # ── Phase 3 — terminal reconstruction (Table 2) ───────────────────────────
    term = json.loads(TERMINAL_JSON.read_text())
    tps  = term.get("paired_stats", {})
    ci   = tps.get("bootstrap_ci_95") or [None, None]

    term_table = wandb.Table(columns=[
        "method", "recovery_mean", "aggregate_recovery", "total_sites",
        "valid_windows", "t_end_nt", "p_value", "ci_lo", "ci_hi", "mean_diff"
    ])
    term_table.add_data(
        "MLM baseline",
        safe(term["mlm_baseline"]["mean_recovery"]),
        safe(term["mlm_baseline"]["aggregate_recovery"]),
        term["mlm_baseline"]["total_sites"],
        term["mlm_baseline"]["valid_windows"],
        term.get("terminal_window_nt", 10),
        safe(tps.get("p_value")),
        safe(ci[0]), safe(ci[1]),
        None,
    )
    term_table.add_data(
        "DAM proposed",
        safe(term["dam_proposed"]["mean_recovery"]),
        safe(term["dam_proposed"]["aggregate_recovery"]),
        term["dam_proposed"]["total_sites"],
        term["dam_proposed"]["valid_windows"],
        term.get("terminal_window_nt", 10),
        safe(tps.get("p_value")),
        safe(ci[0]), safe(ci[1]),
        safe(tps.get("mean_diff_dam_minus_mlm")),
    )
    wandb.log({"phase3_terminal/table2_terminal_reconstruction": term_table})
    print("✓ Phase 3 terminal logged")

    # ── Phase 3 ext — T_END sweep (Figure 3) ──────────────────────────────────
    # Line chart: x = t_end, y = recovery rate for each of 4 methods
    scaling = json.loads(SCALING_JSON.read_text())

    sweep_table = wandb.Table(columns=[
        "t_end", "n_sites", "n_windows",
        "random", "zeroshot", "mlm", "dam", "delta_pp", "p_two", "p_one"
    ])

    for t_str, r in sorted(scaling.items(), key=lambda x: int(x[0])):
        t_end = int(t_str)
        sweep_table.add_data(
            t_end,
            r["n_sites"], r["n_windows"],
            safe(r.get("rnd_mean")),
            safe(r.get("zs_mean")),
            safe(r["mlm_mean"]),
            safe(r["dam_mean"]),
            safe(r["delta"]),
            safe(r["p_value_two"]),
            safe(r["p_value_one"]),
        )
        wandb.log({
            "t_end_sweep/t_end":    t_end,
            "t_end_sweep/random":   safe(r.get("rnd_mean")),
            "t_end_sweep/zeroshot": safe(r.get("zs_mean")),
            "t_end_sweep/mlm":      safe(r["mlm_mean"]),
            "t_end_sweep/dam":      safe(r["dam_mean"]),
            "t_end_sweep/delta_pp": safe(r["delta"]),
            "t_end_sweep/p_two":    safe(r["p_value_two"]),
            "t_end_sweep/p_one":    safe(r["p_value_one"]),
            "t_end_sweep/n_sites":  r["n_sites"],
        })

    wandb.log({"t_end_sweep/results_table": sweep_table})
    print("✓ T_END sweep logged")

    # ── Phase 3 ext — per-position decomposition ──────────────────────────────
    # Line chart: x = d (distance from end), y = recovery rate 5'/3' strands
    if PER_POS_JSON.exists():
        per_pos = json.loads(PER_POS_JSON.read_text())

        pp_table = wandb.Table(columns=["d", "strand", "n", "mlm", "dam", "delta_pp", "p_value"])

        for d_str, strands in sorted(per_pos.items(), key=lambda x: int(x[0])):
            d = int(d_str)
            for strand, s in strands.items():
                pp_table.add_data(
                    d, strand,
                    s["n"],
                    safe(s["mlm"]), safe(s["dam"]),
                    safe(s["delta"]), safe(s["p_value"]),
                )
            wandb.log({
                "per_position/d":        d,
                "per_position/5p_mlm":   safe(strands["5p"]["mlm"]),
                "per_position/5p_dam":   safe(strands["5p"]["dam"]),
                "per_position/5p_delta": safe(strands["5p"]["delta"]),
                "per_position/3p_mlm":   safe(strands["3p"]["mlm"]),
                "per_position/3p_dam":   safe(strands["3p"]["dam"]),
                "per_position/3p_delta": safe(strands["3p"]["delta"]),
            })

        wandb.log({"per_position/results_table": pp_table})
        print("✓ Per-position results logged")

    # ── Phase 4 — ESMFold structural validation (Table 3) ─────────────────────
    t3 = json.loads(TABLE3_JSON.read_text())

    t3_table = wandb.Table(columns=[
        "gene", "method", "plddt_ref", "plddt_recon", "tm_score", "rmsd_ca_angstrom", "aligned_aa"
    ])

    for gene, methods in t3.items():
        for method, r in methods.items():
            if not r:
                continue
            t3_table.add_data(
                gene, method,
                safe(r.get("plddt_ref")),
                safe(r.get("plddt_method")),
                safe(r.get("tm_score")),
                safe(r.get("rmsd_ca")),
                r.get("aligned_aa"),
            )

    wandb.log({"phase4/structural_validation_table3": t3_table})
    print("✓ Phase 4 ESMFold results logged")

    # ── Figures ───────────────────────────────────────────────────────────────
    fig_labels = {
        "fig1_damage_profile.png":  "Figure 1 — PMD Damage Profile",
        "fig2_training_curves.png": "Figure 2 — Training Curves (MLM vs DAM)",
        "fig3_damage_scaling.png":  "Figure 3 — T_END Sensitivity (4 baselines, 626 windows)",
        "fig3b_headline.png":       "Figure 3b — T_END=3 Spotlight (MLM < Random)",
        "fig_per_position.png":     "Figure S1 — Per-Position Decomposition (d=1..25)",
    }
    images = []
    for fname, caption in fig_labels.items():
        fpath = FIGURES_DIR / fname
        if fpath.exists():
            images.append(wandb.Image(str(fpath), caption=caption))
        else:
            print(f"  (skipping {fname} — not found)")

    if images:
        wandb.log({"figures": images})
        print(f"✓ {len(images)} figures uploaded")

    run_url = wandb.run.url if wandb.run else ""
    wandb.finish()
    print(f"\nDone. View run at: {run_url}")


if __name__ == "__main__":
    main()
