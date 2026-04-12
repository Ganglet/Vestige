"""
Phase 4 — Step 2: Damage CDS sequences, reconstruct with MLM and DAM,
translate to protein.

For each gene CDS (from protein/sequences/{gene}_cds.fa):
  1. Apply PMD-profile damage at DAMAGE_SCALE × authentic rates (default 10×).
     Amplification ensures enough damaged positions for a meaningful comparison.
  2. Tile the CDS in 2000 bp windows (stride 200 bp) — same as training.
  3. For each damaged nucleotide position p, identify the "best" window (one
     where p is nearest the 1000 bp centre) and record that window/token.
  4. Run each model once per needed window; collect predictions.
  5. Stitch: replace damaged positions in the damaged CDS string with the
     model's predicted nucleotide.
  6. Translate the stitched CDS → protein; save to protein/sequences/.

Outputs:
    protein/sequences/{gene}_mlm_protein.fa
    protein/sequences/{gene}_dam_protein.fa
    protein/damage_reconstruct_log.txt

Run from project root:
    python protein/damage_reconstruct_cds.py
"""
import sys, types, importlib.util, importlib.machinery

if importlib.util.find_spec("triton") is None:
    _stub = types.ModuleType("triton")
    _stub.__spec__ = importlib.machinery.ModuleSpec("triton", loader=None)
    _stub.__version__ = "0.0.0"
    sys.modules["triton"] = _stub

import numpy as np
import torch
from pathlib import Path
from Bio.Seq import Seq
from transformers import AutoModelForMaskedLM, AutoTokenizer

TOKENIZER_NAME  = "zhihan1996/DNABERT-2-117M"
PROFILE_PATH    = "damage/damage_profile.npy"
SEQ_DIR         = Path("protein/sequences")
OUT_DIR         = Path("protein/sequences")
LOG_PATH        = Path("protein/damage_reconstruct_log.txt")

WINDOW_BP       = 2000
STRIDE_BP       = 200
MAX_TOKENS      = 512
SEED            = 42

# Per-gene damage scale — multiply authentic PMD rates by this factor.
# HBB CDS is only 489 nt; at 10× seed=42 produced 0 damaged positions.
DAMAGE_SCALE_PER_GENE = {
    "TRPV3": 10,
    "KCNK9": 10,
    "HBB":   30,   # amplified further — short CDS + unlucky seed
}

CHECKPOINTS = {
    "mlm": "training/checkpoints/mlm-baseline/checkpoint-595",
    "dam": "training/checkpoints/dam-proposed/checkpoint-595",
}

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
GENES  = ["TRPV3", "KCNK9", "HBB"]


# ── Model loading ──────────────────────────────────────────────────────────────

def load_model(ckpt_path: str, device: str):
    model = AutoModelForMaskedLM.from_pretrained(
        TOKENIZER_NAME, trust_remote_code=True, dtype=torch.float32,
    )
    state = torch.load(
        Path(ckpt_path) / "pytorch_model.bin",
        map_location="cpu", weights_only=True,
    )
    missing, unexpected = model.load_state_dict(state, strict=False)
    if missing or unexpected:
        raise RuntimeError(f"Weight mismatch — {missing}, {unexpected}")
    model.eval().to(device)
    return model


# ── Damage application ─────────────────────────────────────────────────────────

def apply_damage(cds: str, ct_5p: np.ndarray, ga_3p: np.ndarray,
                 scale: float, rng) -> tuple[str, list[int]]:
    """Apply scaled PMD damage; return (damaged_cds, list_of_damaged_nuc_positions)."""
    seq = list(cds)
    n   = len(seq)
    damaged = []
    for i, base in enumerate(seq):
        if base == "C":
            prob = min(ct_5p[min(i, len(ct_5p) - 1)] * scale, 1.0)
            if rng.random() < prob:
                seq[i] = "T"
                damaged.append(i)
        elif base == "G":
            dist = n - 1 - i
            prob = min(ga_3p[min(dist, len(ga_3p) - 1)] * scale, 1.0)
            if rng.random() < prob:
                seq[i] = "A"
                damaged.append(i)
    return "".join(seq), damaged


# ── Window tiling ──────────────────────────────────────────────────────────────

def make_window_starts(cds_len: int, window: int, stride: int) -> list[int]:
    starts = list(range(0, max(1, cds_len - window + 1), stride))
    # Ensure last window covers the tail
    if cds_len > window and (cds_len - window) % stride != 0:
        starts.append(cds_len - window)
    return starts


def best_window_for_pos(pos: int, window_starts: list[int], window: int) -> int | None:
    """Return the index into window_starts for which pos is nearest centre."""
    centre = window // 2
    best_idx = None
    best_dist = float("inf")
    for idx, start in enumerate(window_starts):
        end = start + window
        if start <= pos < end:
            dist = abs(pos - start - centre)
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
    return best_idx


# ── Token-position lookup ──────────────────────────────────────────────────────

def nuc_to_tok(local_nuc_pos: int, offsets: list) -> int | None:
    for tok_idx, (s, e) in enumerate(offsets):
        if s == e:
            continue
        if s <= local_nuc_pos < e:
            return tok_idx
    return None


# ── Inference ─────────────────────────────────────────────────────────────────

@torch.no_grad()
def predict_nuc(model, tokenizer, window_seq: str, local_nuc_pos: int, device: str) -> str | None:
    """
    Mask the token containing local_nuc_pos in window_seq, run the model,
    decode the predicted nucleotide at the sub-token offset. Returns None if fails.
    """
    enc = tokenizer(
        window_seq, return_offsets_mapping=True,
        padding="max_length", truncation=True, max_length=MAX_TOKENS,
    )
    offsets  = enc["offset_mapping"]
    gt_ids   = enc["input_ids"]
    attn     = enc["attention_mask"]

    tok_pos = nuc_to_tok(local_nuc_pos, offsets)
    if tok_pos is None:
        return None

    mask_id  = tokenizer.convert_tokens_to_ids(tokenizer.mask_token)
    masked   = gt_ids[:]
    masked[tok_pos] = mask_id

    inp    = torch.tensor([masked],   dtype=torch.long, device=device)
    am     = torch.tensor([attn],     dtype=torch.long, device=device)
    logits = model(input_ids=inp, attention_mask=am).logits  # (1, seq, vocab)

    pred_id    = int(logits[0, tok_pos].argmax())
    pred_token = tokenizer.convert_ids_to_tokens([pred_id])[0]
    if pred_token is None:
        return None

    pred_nuc_str = pred_token.lstrip("#").lstrip("▁").upper()
    nuc_offset   = local_nuc_pos - offsets[tok_pos][0]

    if nuc_offset < 0 or nuc_offset >= len(pred_nuc_str):
        return None
    return pred_nuc_str[nuc_offset]


# ── Translation ────────────────────────────────────────────────────────────────

def translate_cds(cds: str) -> str:
    """Translate CDS (ATG...stop), strip trailing *, handle partial codons."""
    # Trim to complete codons
    trimmed = cds[:len(cds) - len(cds) % 3]
    return str(Seq(trimmed).translate(to_stop=True))


# ── Main ───────────────────────────────────────────────────────────────────────

def run():
    rng     = np.random.default_rng(SEED)
    profile = np.load(PROFILE_PATH, allow_pickle=True).item()
    ct_5p   = profile["ct_5p"]
    ga_3p   = profile["ga_3p"]

    tokenizer = AutoTokenizer.from_pretrained(TOKENIZER_NAME, trust_remote_code=True)
    log_lines = []

    for gene in GENES:
        cds_path = SEQ_DIR / f"{gene}_cds.fa"
        if not cds_path.exists():
            print(f"SKIP {gene} — {cds_path} not found. Run fetch_cds.py first.")
            continue

        # Read CDS (skip header)
        lines = cds_path.read_text().strip().splitlines()
        cds_ref = "".join(l for l in lines if not l.startswith(">")).upper()
        print(f"\n── {gene} ({len(cds_ref)} nt) ─────────────────────────────────────")

        # Apply damage
        scale = DAMAGE_SCALE_PER_GENE.get(gene, 10)
        cds_damaged, damaged_positions = apply_damage(cds_ref, ct_5p, ga_3p, scale, rng)
        print(f"  Damaged positions ({scale}×): {len(damaged_positions)}")
        log_lines.append(f"{gene}: {len(damaged_positions)} damaged positions at {scale}× PMD scale")

        if not damaged_positions:
            print(f"  No positions damaged — consider increasing DAMAGE_SCALE_PER_GENE for {gene}.")
            continue

        window_starts = make_window_starts(len(cds_ref), WINDOW_BP, STRIDE_BP)

        # Assign each damaged position to its best window
        pos_to_win = {}  # global_nuc_pos → window_start_index
        for pos in damaged_positions:
            best = best_window_for_pos(pos, window_starts, WINDOW_BP)
            if best is not None:
                pos_to_win[pos] = best

        # Group by window
        win_to_positions: dict[int, list[int]] = {}
        for pos, win_idx in pos_to_win.items():
            win_to_positions.setdefault(win_idx, []).append(pos)

        for model_key, ckpt in CHECKPOINTS.items():
            print(f"\n  [{model_key}] Loading {ckpt}...")
            model = load_model(ckpt, DEVICE)

            reconstructed = list(cds_damaged)

            for win_idx, positions in win_to_positions.items():
                win_start = window_starts[win_idx]
                win_end   = min(win_start + WINDOW_BP, len(cds_damaged))
                window_seq = cds_damaged[win_start:win_end]

                for global_pos in positions:
                    local_pos = global_pos - win_start
                    pred = predict_nuc(model, tokenizer, window_seq, local_pos, DEVICE)
                    if pred is not None:
                        reconstructed[global_pos] = pred

            reconstructed_cds = "".join(reconstructed)

            # Count how many positions were correctly restored
            n_correct = sum(
                1 for p in damaged_positions
                if reconstructed_cds[p] == cds_ref[p]
            )
            pct = n_correct / len(damaged_positions) * 100
            print(f"  [{model_key}] Reconstructed {n_correct}/{len(damaged_positions)} damaged positions correctly ({pct:.1f}%)")
            log_lines.append(f"  {model_key}: {n_correct}/{len(damaged_positions)} correct ({pct:.1f}%)")

            # Translate
            protein = translate_cds(reconstructed_cds)
            print(f"  [{model_key}] Protein length: {len(protein)} aa")

            out_path = OUT_DIR / f"{gene}_{model_key}_protein.fa"
            out_path.write_text(f">{gene}_{model_key}_reconstructed\n{protein}\n")
            print(f"  [{model_key}] Saved → {out_path}")

            del model
            if DEVICE == "mps":
                torch.mps.empty_cache()

    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    LOG_PATH.write_text("\n".join(log_lines) + "\n")
    print(f"\nLog → {LOG_PATH}")


if __name__ == "__main__":
    run()
