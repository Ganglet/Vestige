"""
VESTIGE — End-to-end pipeline schematic for CIBM submission.
Run from project root:
    python3 results/plot_pipeline.py
Outputs: results/figures/fig_pipeline.{pdf,png}
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

OUT_DIR = Path("results/figures")

# ── Palette ──────────────────────────────────────────────────────────────────
C_DATA   = "#1e3a5f"
C_DATA_L = "#d6e4f0"
C_PROC   = "#2d6a4f"
C_PROC_L = "#d8f3dc"
C_CORE   = "#7b2d8b"
C_CORE_L = "#f3e8ff"
C_BASE   = "#b5451b"
C_BASE_L = "#fde8d8"
C_EVAL   = "#1a5276"
C_EVAL_L = "#d6eaf8"
C_OUT    = "#145a32"
C_OUT_L  = "#d5f5e3"
C_ARROW  = "#4a4a4a"
FONT     = "DejaVu Sans"

# ── Global size tokens ────────────────────────────────────────────────────────
FS_L   = 11.0   # box label
FS_S   =  9.0   # box sublabel  (FS_L - 2)
FS_SEC =  9.5   # section badge
FS_NOTE=  8.5   # italic notes
LW_BOX = 1.9    # default box edge
LW_ARR = 1.9    # default arrow
MS     = 11     # arrowhead mutation_scale


def box(ax, x, y, w, h, label, sublabel=None,
        fc="#f0f0f0", ec="#333333", fontsize=FS_L,
        bold=True, radius=0.025, lw=LW_BOX, label_color="black"):
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=f"round,pad=0.01,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3
    )
    ax.add_patch(patch)
    fw = "bold" if bold else "normal"
    if sublabel:
        ax.text(x, y + h * 0.15, label, ha="center", va="center",
                fontsize=fontsize, fontweight=fw, color=label_color,
                fontfamily=FONT, zorder=4)
        ax.text(x, y - h * 0.20, sublabel, ha="center", va="center",
                fontsize=fontsize - 2, color="#444444",
                fontfamily=FONT, zorder=4, style="italic")
    else:
        ax.text(x, y, label, ha="center", va="center",
                fontsize=fontsize, fontweight=fw, color=label_color,
                fontfamily=FONT, zorder=4)
    return patch


def arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=LW_ARR, style="-|>",
          shrink=3, connectionstyle="arc3,rad=0.0"):
    ax.annotate("",
        xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle=style, color=color, lw=lw,
            mutation_scale=MS,
            shrinkA=shrink, shrinkB=shrink,
            connectionstyle=connectionstyle,
        ),
        zorder=5
    )


def section_label(ax, x, y, text, color):
    ax.text(x, y, text, ha="left", va="center",
            fontsize=FS_SEC, color=color, fontweight="bold",
            fontfamily=FONT,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color,
                      lw=1.2, alpha=0.9),
            zorder=6)


def divider(ax, y, xmin=0.01, xmax=0.99, color="#cccccc", ls="--"):
    ax.axhline(y, xmin=xmin, xmax=xmax, color=color, lw=1.0, ls=ls, zorder=1)


def make_pipeline():
    fig, ax = plt.subplots(figsize=(18, 12.5))   # larger canvas for 100% readability
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10.1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 1 — DATA INPUTS  (y=9.3, h=0.82)
    # box top=9.71, bottom=8.89
    # ═══════════════════════════════════════════════════════════════════════
    section_label(ax, 0.15, 9.65, "① DATA", C_DATA)

    box(ax, 2.5, 9.3, 3.2, 0.82,
        "Woolly Mammoth aDNA",
        "NCBI SRA · ERP008929\nERR855944 (4,300 BP)  ·  ERR852028 (44,800 BP)",
        fc=C_DATA_L, ec=C_DATA, label_color=C_DATA)

    box(ax, 7.5, 9.3, 3.2, 0.82,
        "Asian Elephant Reference",
        "GCF_024166365.1 (EanMak 1.0)\ngene_coords.csv  (NCBI Entrez lookup)",
        fc=C_DATA_L, ec=C_DATA, label_color=C_DATA)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 2 — PREPROCESSING  (y=8.1, h=0.70)
    # box top=8.45, bottom=7.75
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 8.72, color="#bbbbbb")
    section_label(ax, 0.15, 8.52, "② PREPROCESSING", C_PROC)

    box(ax, 2.5,  8.1, 2.4, 0.70, "BWA-aln",
        "aDNA params\n-l 16500  -n 0.01",
        fc=C_PROC_L, ec=C_PROC, label_color=C_PROC)

    box(ax, 5.5,  8.1, 2.4, 0.70, "mapDamage2",
        "misincorporation.txt\n→ damage_profile.npy",
        fc=C_PROC_L, ec=C_PROC, label_color=C_PROC)

    box(ax, 8.8,  8.1, 2.6, 0.70, "Gene Extraction",
        "NCBI Entrez · 7 genes · 2000 bp windows\nstride 200 bp · 626 eval windows",
        fc=C_PROC_L, ec=C_PROC, label_color=C_PROC)

    ax.text(2.5, 7.67, "seqtk subsample  ·  5M reads/specimen  ·  seed 42",
            ha="center", va="center", fontsize=FS_NOTE, color=C_PROC,
            fontfamily=FONT, style="italic")

    # Row 1 → Row 2 arrows  (data bottom=8.89 → preproc top=8.45)
    arrow(ax, 2.5, 8.89, 2.5,  8.45)
    arrow(ax, 3.7, 8.1,  4.3,  8.1)
    arrow(ax, 7.5, 8.89, 8.8,  8.45)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 3 — DATASET + MASKING  (y=6.85, h=0.85)
    # box top=7.275, bottom=6.425
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 7.55, color="#bbbbbb")
    section_label(ax, 0.15, 7.32, "③ DATASET & MASKING", C_CORE)

    box(ax, 3.0, 6.85, 3.4, 0.85, "Dataset",
        "Train 275 · Val 69 · Held-out 557\n80/20 split · 3 training + 4 held-out genes",
        fc="#f8f4ff", ec=C_CORE, label_color=C_CORE)

    box(ax, 7.0, 6.85, 2.6, 0.85, "MLM Collator  (baseline)",
        "Uniform 15% masking\nAll positions equal probability",
        fc=C_BASE_L, ec=C_BASE, label_color=C_BASE)

    box(ax, 10.4, 6.85, 2.9, 0.85, "DAM Collator  ★  (proposed)",
        "PMD-weighted masking\nC/G avg 15% · A/T never masked",
        fc=C_CORE_L, ec=C_CORE, label_color=C_CORE, lw=2.5)

    # ── Preprocessing → Masking: L-shaped orthogonal routes ──────────────
    # Gap between preproc bottom (7.75) and masking top (7.275):
    #   y=7.62  Gene Extraction → Dataset
    #   y=7.44  mapDamage2 → DAM Collator

    def lpath(xs, ys, color, lw=1.7):
        ax.plot(xs, ys, color=color, lw=lw, zorder=5,
                solid_capstyle="round", solid_joinstyle="round")

    def ltip(x1, y1, x2, y2, color, lw=1.7):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                   lw=lw, mutation_scale=MS,
                                   shrinkA=0, shrinkB=2), zorder=6)

    lpath([8.8, 8.8, 3.0], [7.75, 7.62, 7.62], C_PROC)
    ltip(3.0, 7.62, 3.0, 7.275, C_PROC)

    lpath([5.5, 5.5, 10.4], [7.75, 7.44, 7.44], C_PROC)
    ltip(10.4, 7.44, 10.4, 7.275, C_PROC)

    # Dataset → MLM (direct horizontal; MLM left edge = 7.0 - 1.3 = 5.7)
    arrow(ax, 4.7, 6.85, 5.7, 6.85)

    # Dashed "same dataset" connector between collators (MLM right=8.3, DAM left=8.95)
    ax.plot([8.3, 8.95], [6.85, 6.85], color=C_ARROW, lw=1.3,
            ls="dashed", zorder=4, solid_capstyle="round")
    ax.text(8.625, 6.95, "same dataset", ha="center", va="bottom",
            fontsize=8.0, color=C_ARROW, fontfamily=FONT, style="italic", zorder=5)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 4 — FINE-TUNING  (y=5.65, h=0.80)
    # box top=6.05, bottom=5.25
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 6.25, color="#bbbbbb")
    section_label(ax, 0.15, 6.02, "④ FINE-TUNING", C_BASE)

    box(ax, 7.0,  5.65, 2.6, 0.80, "DNABERT-2  (MLM)",
        "117M params · HuggingFace Trainer\nepoch 17 · loss 3.757",
        fc=C_BASE_L, ec=C_BASE, label_color=C_BASE)

    box(ax, 10.4, 5.65, 2.9, 0.80, "DNABERT-2  (DAM)  ★",
        "117M params · HuggingFace Trainer\nepoch 17 · loss 3.274  (↓13%)",
        fc=C_CORE_L, ec=C_CORE, label_color=C_CORE, lw=2.5)

    # Masking bottom (6.425) → Fine-tuning top (6.05)
    arrow(ax, 7.0,  6.425, 7.0,  6.05)
    arrow(ax, 10.4, 6.425, 10.4, 6.05)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 5 — EVALUATION  (y=4.3, h=0.82)
    # box top=4.71, bottom=3.89
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 5.08, color="#bbbbbb")
    section_label(ax, 0.15, 4.74, "⑤ EVALUATION", C_EVAL)

    primary_eval = [
        (1.5,  4.3, "Nucleotide\nRecovery",  "Background\n+ Terminal"),
        (4.0,  4.3, "Terminal Zone\nSweep",  "T_END ∈ {3,5,10,15,20,25}\n6 widths · 626 windows"),
        (12.0, 4.3, "In silico\nDamage Sim", "BPE-safe\noffset mapping"),
    ]
    for ex, ey, el, esl in primary_eval:
        box(ax, ex, ey, 2.2, 0.82, el, esl,
            fc=C_EVAL_L, ec=C_EVAL, label_color=C_EVAL)

    C_SUPP_L, C_SUPP = "#f0f0f0", "#999999"
    supp_eval = [
        (6.5, 4.3, "Per-position\nAnalysis",  "d = 1…25 · 5′-C / 3′-G\n(Fig. S1)"),
        (9.0, 4.3, "Damage\nIntensity Sweep", "5% – 40% peak rates\n(Fig. S2)"),
    ]
    for ex, ey, el, esl in supp_eval:
        patch = FancyBboxPatch(
            (ex - 1.1, ey - 0.41), 2.2, 0.82,
            boxstyle="round,pad=0.01,rounding_size=0.025",
            facecolor=C_SUPP_L, edgecolor=C_SUPP,
            linewidth=1.2, linestyle="dashed", zorder=3
        )
        ax.add_patch(patch)
        ax.text(ex, ey + 0.82 * 0.12, el, ha="center", va="center",
                fontsize=10, fontweight="bold", color=C_SUPP,
                fontfamily=FONT, zorder=4)
        ax.text(ex, ey - 0.82 * 0.22, esl, ha="center", va="center",
                fontsize=8.5, color=C_SUPP, fontfamily=FONT,
                zorder=4, style="italic")

    all_eval_x = [e[0] for e in primary_eval] + [e[0] for e in supp_eval]

    # ── Bus routing: models → horizontal bus → eval boxes ────────────────
    BUS_Y  = 4.98
    BUS_X0 = 1.5
    BUS_X1 = 12.0
    BUS_MID= (7.0 + 10.4) / 2   # 8.7

    # Fine-tuning bottoms (5.25) drop to bus
    ax.annotate("", xy=(7.0,  BUS_Y), xytext=(7.0,  5.25),
                arrowprops=dict(arrowstyle="-", color=C_BASE, lw=1.7), zorder=5)
    ax.annotate("", xy=(10.4, BUS_Y), xytext=(10.4, 5.25),
                arrowprops=dict(arrowstyle="-", color=C_CORE, lw=1.7), zorder=5)

    ax.plot([BUS_X0, BUS_MID], [BUS_Y, BUS_Y],
            color=C_BASE, lw=1.7, zorder=5, solid_capstyle="round")
    ax.plot([BUS_MID, BUS_X1], [BUS_Y, BUS_Y],
            color=C_CORE, lw=1.7, zorder=5, solid_capstyle="round")

    for ex in all_eval_x:
        col = C_BASE if ex <= BUS_MID else C_CORE
        ax.annotate("", xy=(ex, 4.71), xytext=(ex, BUS_Y),
                    arrowprops=dict(arrowstyle="-|>", color=col,
                                   lw=1.5, mutation_scale=MS), zorder=5)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 6 — VALIDATION  (y=2.9, h=0.90)
    # box top=3.35, bottom=2.45
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 3.66, color="#bbbbbb")
    section_label(ax, 0.15, 3.33, "⑥ VALIDATION", C_OUT)

    box(ax, 4.5,  2.9, 3.8, 0.90,
        "Protein Structural Validation",
        "ESMFold REST API · Kabsch superposition\nTM-score > 0.95 · Cα-RMSD  (TRPV3, KCNK9, HBB)",
        fc=C_OUT_L, ec=C_OUT, label_color=C_OUT)

    box(ax, 10.5, 2.9, 3.8, 0.90,
        "Biosecurity Classifier",
        "1D CNN · AUC 0.934\n10 virulence gene classes · 98.2% windows cleared",
        fc=C_OUT_L, ec=C_OUT, label_color=C_OUT)

    # Eval bottom (3.89) → Validation top (3.35)
    arrow(ax, 4.0, 3.89, 4.5,  3.35, color=C_OUT)
    arrow(ax, 9.0, 3.89, 10.5, 3.35, color=C_OUT)

    # ═══════════════════════════════════════════════════════════════════════
    # ROW 7 — KEY RESULT
    # ═══════════════════════════════════════════════════════════════════════
    divider(ax, 2.27, color="#bbbbbb")

    result_patch = FancyBboxPatch(
        (1.5, 1.18), 13.0, 0.95,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor="#fffbeb", edgecolor="#b45309", linewidth=2.2, zorder=3
    )
    ax.add_patch(result_patch)

    ax.text(8.0, 1.87,
            "Central Result",
            ha="center", va="center", fontsize=11.5, fontweight="bold",
            color="#92400e", fontfamily=FONT, zorder=4)

    ax.text(8.0, 1.48,
            "At peak PMD positions (T_END = 3 nt):  "
            "Zero-shot 15.5%  ·  MLM 20.5%  <  Random chance 27.7%  ·  "
            "DAM 30.8%   (p < 0.001,  n = 626 windows,  7 genes)",
            ha="center", va="center", fontsize=10.5, color="#1a1a2e",
            fontfamily=FONT, zorder=4)

    # Validation bottom (2.45) → Result top (2.13)
    arrow(ax, 4.5,  2.45, 5.5,  2.13, color=C_OUT)
    arrow(ax, 10.5, 2.45, 10.5, 2.13, color=C_OUT)

    # ── W&B badge ────────────────────────────────────────────────────────
    ax.text(15.85, 0.18, "W&B · Full reproducibility · seed 42",
            ha="right", va="bottom", fontsize=7.5, color="#888888",
            fontfamily=FONT, style="italic")

    # ── Legend ───────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(fc=C_DATA_L, ec=C_DATA, label="Data input"),
        mpatches.Patch(fc=C_PROC_L, ec=C_PROC, label="Preprocessing"),
        mpatches.Patch(fc=C_BASE_L, ec=C_BASE, label="MLM baseline"),
        mpatches.Patch(fc=C_CORE_L, ec=C_CORE, label="DAM (proposed)  ★"),
        mpatches.Patch(fc=C_EVAL_L, ec=C_EVAL, label="Evaluation"),
        mpatches.Patch(fc=C_OUT_L,  ec=C_OUT,  label="Validation / Output"),
    ]
    ax.legend(handles=legend_patches, loc="lower left",
              bbox_to_anchor=(0.0, 0.0),
              fontsize=9.0, frameon=True, framealpha=0.95,
              edgecolor="#cccccc", ncol=3,
              handlelength=1.4, handleheight=1.0,
              columnspacing=1.2)

    plt.tight_layout(pad=0.3)

    out = OUT_DIR / "fig_pipeline"
    fig.savefig(str(out) + ".pdf", bbox_inches="tight", dpi=300)
    fig.savefig(str(out) + ".png", dpi=300, bbox_inches="tight")
    print(f"Pipeline figure → {out}.{{pdf,png}}")
    plt.close(fig)


if __name__ == "__main__":
    make_pipeline()
