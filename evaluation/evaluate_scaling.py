"""
Terminal-zone sensitivity analysis — Phoenix-LM Phase 3 (Figure 3).

Runs the terminal-position reconstruction evaluation at T_END ∈ {3, 5, 10, 15, 20, 25} nt.
T_END defines how many nucleotides from each window end are treated as "terminal" —
the zone where mapDamage2 shows exponentially decaying C→T and G→A deamination.

DAM masking probability peaks at position 1 (~28%) and decays toward the window
interior. If DAM's advantage is genuine, it should be strongest at small T_END
(highest-damage zone) and diminish as T_END grows to include background-like positions.

Method: mask ALL terminal C/G positions in each validation window (no stochastic
damage applied), matching the approach in evaluate_terminal.py. This gives well-
powered evaluation (n = hundreds of sites) at every T_END.

Outputs:
    evaluation/scaling_results.json
    results/figures/fig3_damage_scaling.{pdf,png}

Run from project root:
    python3 evaluation/evaluate_scaling.py
"""
import sys, types, importlib.util, importlib.machinery

if importlib.util.find_spec("triton") is None:
    _stub = types.ModuleType("triton")
    _stub.__spec__ = importlib.machinery.ModuleSpec("triton", loader=None)
    _stub.__version__ = "0.0.0"
    sys.modules["triton"] = _stub

import json
import numpy as np
import torch
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from transformers import AutoModelForMaskedLM, AutoTokenizer

DAMAGED_VAL    = "evaluation/damaged_validation.npy"
TOKENIZER_NAME = "zhihan1996/DNABERT-2-117M"
CHECKPOINTS = {
    "MLM baseline": "training/checkpoints/mlm-baseline/checkpoint-595",
    "DAM proposed":  "training/checkpoints/dam-proposed/checkpoint-595",
}
OUT_JSON = Path("evaluation/scaling_results.json")
OUT_FIG  = Path("results/figures/fig3_damage_scaling")

T_ENDS = [3, 5, 10, 15, 20, 25]
N_BOOT = 10_000
SEED   = 42
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"


def load_model(ckpt, device):
    model = AutoModelForMaskedLM.from_pretrained(
        TOKENIZER_NAME, trust_remote_code=True, dtype=torch.float32,
    )
    state = torch.load(Path(ckpt) / "pytorch_model.bin",
                       map_location="cpu", weights_only=True)
    model.load_state_dict(state, strict=False)
    return model.eval().to(device)


def find_terminal(dna: str, offsets: list, t_end: int):
    n = len(dna)
    pairs, seen, tok_positions = [], set(), []
    for nuc_pos, base in enumerate(dna):
        if not ((base == "C" and nuc_pos < t_end) or
                (base == "G" and nuc_pos >= n - t_end)):
            continue
        for tok_idx, (s, e) in enumerate(offsets):
            if s == e:
                continue
            if s <= nuc_pos < e:
                pairs.append((nuc_pos, tok_idx))
                if tok_idx not in seen:
                    tok_positions.append(tok_idx)
                    seen.add(tok_idx)
                break
    return pairs, tok_positions


@torch.no_grad()
def predict_at(model, tokenizer, gt_ids, tok_positions, attn_mask, device):
    mask_id = tokenizer.convert_tokens_to_ids(tokenizer.mask_token)
    masked  = gt_ids[:]
    for tp in tok_positions:
        masked[tp] = mask_id
    logits = model(
        input_ids=torch.tensor([masked], dtype=torch.long, device=device),
        attention_mask=torch.tensor([attn_mask], dtype=torch.long, device=device),
    ).logits
    log_probs = torch.log_softmax(logits[0], dim=-1)  # (seq_len, vocab)
    preds    = {tp: int(logits[0, tp].argmax()) for tp in tok_positions}
    lp_correct = {tp: float(log_probs[tp, gt_ids[tp]]) for tp in tok_positions}
    return preds, lp_correct


def nuc_correct(dna, nuc_pos, tok_pos, pred_id, offsets, tokenizer) -> bool:
    tok = tokenizer.convert_ids_to_tokens([pred_id])[0]
    if not tok:
        return False
    s   = tok.lstrip("#").lstrip("▁").upper()
    off = nuc_pos - offsets[tok_pos][0]
    return 0 <= off < len(s) and s[off] == dna[nuc_pos]


def bootstrap_ci(diffs, n, rng):
    means = np.array([rng.choice(diffs, size=len(diffs), replace=True).mean()
                      for _ in range(n)])
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def run_t_end(t_end, windows, mlm_model, dam_model, tokenizer, rng):
    mlm_rates, dam_rates = [], []
    mlm_lp,   dam_lp    = [], []   # mean log-prob of correct token per window
    total_sites = 0

    for w in windows:
        offsets = [tuple(o) for o in w["offsets"]]
        pairs, tok_positions = find_terminal(w["dna_gt"], offsets, t_end)
        if not pairs:
            continue
        total_sites += len(pairs)

        for model, rates, lp_list in [
            (mlm_model, mlm_rates, mlm_lp),
            (dam_model, dam_rates, dam_lp),
        ]:
            preds, lp_correct = predict_at(model, tokenizer, w["gt_ids"],
                                           tok_positions, w["attention_mask"], DEVICE)
            correct = sum(
                nuc_correct(w["dna_gt"], np_, tp, preds[tp], offsets, tokenizer)
                for np_, tp in pairs if tp in preds
            )
            rates.append(correct / len(pairs))
            # mean log-prob of correct token across all terminal tok_positions in window
            lp_list.append(float(np.mean([lp_correct[tp] for tp in tok_positions])))

    r_mlm = np.array(mlm_rates)
    r_dam = np.array(dam_rates)
    lp_mlm = np.array(mlm_lp)
    lp_dam = np.array(dam_lp)
    n_w   = len(r_mlm)

    t_stat = p_val = p_one = ci_lo = ci_hi = float("nan")
    lp_t = lp_p = lp_p_one = float("nan")
    if n_w >= 2:
        # Recovery rate — one-sided (alternative: DAM > MLM)
        t_stat, p_two = stats.ttest_rel(r_dam, r_mlm)
        p_val  = float(p_two)
        p_one  = float(p_two / 2) if t_stat > 0 else float(1 - p_two / 2)
        diffs  = r_dam - r_mlm
        ci_lo, ci_hi = bootstrap_ci(diffs, N_BOOT, rng)
        # Log-probability — one-sided
        lp_t, lp_p_two = stats.ttest_rel(lp_dam, lp_mlm)
        lp_p     = float(lp_p_two)
        lp_p_one = float(lp_p_two / 2) if lp_t > 0 else float(1 - lp_p_two / 2)

    return {
        "t_end":          t_end,
        "n_sites":        total_sites,
        "n_windows":      n_w,
        "mlm_mean":       float(r_mlm.mean())  if n_w else float("nan"),
        "dam_mean":       float(r_dam.mean())  if n_w else float("nan"),
        "delta":          float((r_dam - r_mlm).mean()) if n_w else float("nan"),
        "t_stat":         float(t_stat),
        "p_value_two":    p_val,
        "p_value_one":    p_one,
        "ci_95":          [ci_lo, ci_hi],
        "mlm_logprob":    float(lp_mlm.mean()) if n_w else float("nan"),
        "dam_logprob":    float(lp_dam.mean()) if n_w else float("nan"),
        "delta_logprob":  float((lp_dam - lp_mlm).mean()) if n_w else float("nan"),
        "lp_t_stat":      float(lp_t),
        "lp_p_two":       lp_p,
        "lp_p_one":       lp_p_one,
    }


def main():
    rng       = np.random.default_rng(SEED)
    tokenizer = AutoTokenizer.from_pretrained(TOKENIZER_NAME, trust_remote_code=True)
    windows   = np.load(DAMAGED_VAL, allow_pickle=True).tolist()

    print("Loading MLM baseline...")
    mlm_model = load_model(CHECKPOINTS["MLM baseline"], DEVICE)
    print("Loading DAM proposed...")
    dam_model = load_model(CHECKPOINTS["DAM proposed"], DEVICE)

    print(f"\n{'T_END':>6}  {'Sites':>7}  {'Win':>5}  "
          f"{'MLM':>8}  {'DAM':>8}  {'Δ':>8}  {'p(2-sided)':>11}  {'p(1-sided)':>11}")
    print("-" * 76)

    all_results = {}
    for t_end in T_ENDS:
        r = run_t_end(t_end, windows, mlm_model, dam_model, tokenizer, rng)
        all_results[t_end] = r

        def ps(v): return f"{v:.3f}" if not np.isnan(v) else "  N/A"
        def sig(v): return " *" if not np.isnan(v) and v < 0.05 else "  "

        print(f"{t_end:>5}nt  {r['n_sites']:>7}  {r['n_windows']:>5}  "
              f"{r['mlm_mean']:>8.4f}  {r['dam_mean']:>8.4f}  {r['delta']:>+8.4f}  "
              f"{ps(r['p_value_two'])}{sig(r['p_value_two'])}  "
              f"{ps(r['p_value_one'])}{sig(r['p_value_one'])}")

    OUT_JSON.write_text(json.dumps({str(k): v for k, v in all_results.items()}, indent=2))

    # ── Figure 3 ────────────────────────────────────────────────────────────────
    t_vals  = [r["t_end"]        for r in all_results.values()]
    mlm_v   = [r["mlm_mean"]    for r in all_results.values()]
    dam_v   = [r["dam_mean"]    for r in all_results.values()]
    deltas  = [r["delta"]       for r in all_results.values()]
    p_vals  = [r["p_value_one"] for r in all_results.values()]   # one-sided
    n_sites = [r["n_sites"]     for r in all_results.values()]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1.plot(t_vals, mlm_v, "o-", color="#2980b9", lw=2, ms=7, label="MLM baseline")
    ax1.plot(t_vals, dam_v, "s-", color="#c0392b", lw=2, ms=7, label="DAM proposed")
    for tv, p, d, m in zip(t_vals, p_vals, dam_v, mlm_v):
        if not np.isnan(p) and p < 0.05:
            ax1.annotate("*", xy=(tv, max(d, m) + 0.006), ha="center",
                         fontsize=13, color="#27ae60", fontweight="bold")
    ax1.set_xlabel("Terminal zone width (nt from each end)", fontsize=11)
    ax1.set_ylabel("Recovery rate (mean per window)", fontsize=11)
    ax1.set_title("Recovery rate vs terminal zone width", fontsize=12, fontweight="bold")
    ax1.set_xticks(t_vals)
    ax1.legend(fontsize=9, frameon=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    colors = ["#27ae60" if d > 0 else "#e74c3c" for d in deltas]
    ax2.bar(t_vals, deltas, color=colors, alpha=0.75, width=1.6)
    ax2.axhline(0, color="black", lw=0.8, ls="--")
    for tv, p, d in zip(t_vals, p_vals, deltas):
        if not np.isnan(p) and p < 0.05:
            ax2.annotate("*", xy=(tv, d + 0.003 if d >= 0 else d - 0.006),
                         ha="center", fontsize=13, color="#27ae60", fontweight="bold")
    ax2b = ax2.twinx()
    ax2b.plot(t_vals, n_sites, "D--", color="gray", lw=1.2, ms=5, label="Sites (n)")
    ax2b.set_ylabel("Sites evaluated (n)", fontsize=10, color="gray")
    ax2b.tick_params(axis="y", labelcolor="gray")
    ax2.set_xlabel("Terminal zone width (nt from each end)", fontsize=11)
    ax2.set_ylabel("DAM − MLM recovery (Δ)", fontsize=11)
    ax2.set_title("DAM advantage vs terminal zone width", fontsize=12, fontweight="bold")
    ax2.set_xticks(t_vals)
    ax2.spines["top"].set_visible(False)
    ax2b.spines["top"].set_visible(False)

    fig.suptitle(
        "Terminal C/G reconstruction · 69 validation windows · TRPV3 / KCNK9 / HBB"
        " · * p < 0.05 (paired t-test)",
        fontsize=9, color="gray", y=1.01,
    )
    plt.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=300, bbox_inches="tight")
    print(f"\nFigure 3 → {OUT_FIG}.{{pdf,png}}")
    print(f"JSON     → {OUT_JSON}")


if __name__ == "__main__":
    main()
