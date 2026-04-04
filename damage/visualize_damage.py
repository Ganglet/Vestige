"""
Figure 1: C→T (5′) and G→A (3′) damage frequency profiles from mapDamage2.

Run from project root: python damage/visualize_damage.py
Requires damage/damage_profile.npy (run parse_profiles.py first).
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

profile = np.load("damage/damage_profile.npy", allow_pickle=True).item()
ct_5p = profile["ct_5p"]
ga_3p = profile["ga_3p"]
pos = profile["positions"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

RED = "#c0392b"
BLUE = "#2980b9"

for ax, freq, color, title, xlabel in [
    (ax1, ct_5p, RED,  "C→T  (5′ end)", "Position from 5′ end (bp)"),
    (ax2, ga_3p, BLUE, "G→A  (3′ end)", "Position from 3′ end (bp)"),
]:
    ax.plot(pos, freq, color=color, linewidth=1.5)
    ax.fill_between(pos, freq, alpha=0.12, color=color)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Substitution frequency", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(pos[0], pos[-1])
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle(
    "ERR855944  ·  Palkopoulou et al. 2015  ·  mapDamage2 v2.1.0",
    fontsize=9, color="gray", y=1.01
)
plt.tight_layout()

out = Path("results/figures")
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "fig1_damage_profile.pdf", bbox_inches="tight")
fig.savefig(out / "fig1_damage_profile.png", dpi=300, bbox_inches="tight")
print("Saved results/figures/fig1_damage_profile.{pdf,png}")
