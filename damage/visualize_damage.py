"""
Figure 1: C→T (5′) and G→A (3′) damage frequency profiles from mapDamage2.
Overlays both specimens with mean shown as solid line.

Run from project root: python damage/visualize_damage.py
Requires damage/damage_profile.npy (run parse_profiles.py first).
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

p = np.load("damage/damage_profile.npy", allow_pickle=True).item()
pos = p["positions"]

SAMPLES = {"ERR855944": "#c0392b", "ERR852028": "#e74c3c"}
MEAN_COLORS = {"ct": "#922b21", "ga": "#1a5276"}
BLUE_SAMPLES = {"ERR855944": "#2980b9", "ERR852028": "#5dade2"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

for name, color in SAMPLES.items():
    ax1.plot(pos, p[f"ct_5p_{name}"], color=color, linewidth=1, alpha=0.5, label=name)
for name, color in BLUE_SAMPLES.items():
    ax2.plot(pos, p[f"ga_3p_{name}"], color=color, linewidth=1, alpha=0.5, label=name)

ax1.plot(pos, p["ct_5p"], color=MEAN_COLORS["ct"], linewidth=2, label="mean")
ax1.fill_between(pos, p["ct_5p"], alpha=0.10, color=MEAN_COLORS["ct"])
ax2.plot(pos, p["ga_3p"], color=MEAN_COLORS["ga"], linewidth=2, label="mean")
ax2.fill_between(pos, p["ga_3p"], alpha=0.10, color=MEAN_COLORS["ga"])

for ax, title, xlabel in [
    (ax1, "C→T  (5′ end)", "Position from 5′ end (bp)"),
    (ax2, "G→A  (3′ end)", "Position from 3′ end (bp)"),
]:
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Substitution frequency", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(pos[0], pos[-1])
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, frameon=False)

fig.suptitle(
    "Palkopoulou et al. 2015  ·  ERR855944 + ERR852028  ·  mapDamage2 v2.1.0",
    fontsize=9, color="gray", y=1.01,
)
plt.tight_layout()

out = Path("results/figures")
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "fig1_damage_profile.pdf", bbox_inches="tight")
fig.savefig(out / "fig1_damage_profile.png", dpi=300, bbox_inches="tight")
print("Saved results/figures/fig1_damage_profile.{pdf,png}")
