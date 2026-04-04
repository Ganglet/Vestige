"""
Parse mapDamage2 misincorporation.txt → damage_profile.npy

Aggregates C→T (5') and G→A (3') substitution frequencies
across all chromosomes and strands, grouped by read-end position.

Run from project root: python damage/parse_profiles.py
"""
import numpy as np
import pandas as pd
from pathlib import Path

MISINCORP = Path("Dataset/results/results_mammoth_ERR855944/misincorporation.txt")
OUT = Path("damage/damage_profile.npy")

df = pd.read_csv(MISINCORP, sep="\t", comment="#")
# mapDamage2 uses ">" in column names — pandas reads them as-is
ct_col = "C>T"
ga_col = "G>A"

# Aggregate across all chromosomes and strands
five_p = df[df["End"] == "5p"].groupby("Pos")[[ct_col, "Total"]].sum()
three_p = df[df["End"] == "3p"].groupby("Pos")[[ga_col, "Total"]].sum()

ct_freq = (five_p[ct_col] / five_p["Total"]).values
ga_freq = (three_p[ga_col] / three_p["Total"]).values
positions = five_p.index.values  # 1..70

profile = {"ct_5p": ct_freq, "ga_3p": ga_freq, "positions": positions}
np.save(OUT, profile, allow_pickle=True)

print(f"Saved {OUT}  —  {len(positions)} positions")
print(f"  C→T pos 1 (5′): {ct_freq[0]:.4f}")
print(f"  G→A pos 1 (3′): {ga_freq[0]:.4f}")
