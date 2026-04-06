"""
Parse mapDamage2 misincorporation.txt → damage_profile.npy

Parses both specimens (ERR855944, ERR852028), stores individual profiles
and their mean. The DAM collator uses ct_5p / ga_3p (mean).

Run from project root: python damage/parse_profiles.py
"""
import numpy as np
import pandas as pd
from pathlib import Path

SAMPLES = {
    "ERR855944": Path("Dataset/results/results_mammoth_ERR855944/misincorporation.txt"),
    "ERR852028": Path("Dataset/results/results_mammoth_ERR852028/misincorporation.txt"),
}
OUT = Path("damage/damage_profile.npy")

CT, GA = "C>T", "G>A"


def parse(path: Path) -> tuple:
    df = pd.read_csv(path, sep="\t", comment="#")
    five_p  = df[df["End"] == "5p"].groupby("Pos")[[CT, "Total"]].sum()
    three_p = df[df["End"] == "3p"].groupby("Pos")[[GA, "Total"]].sum()
    return (
        (five_p[CT]  / five_p["Total"]).values,
        (three_p[GA] / three_p["Total"]).values,
        five_p.index.values,
    )


profiles = {}
for name, path in SAMPLES.items():
    ct, ga, pos = parse(path)
    profiles[f"ct_5p_{name}"] = ct
    profiles[f"ga_3p_{name}"] = ga
    print(f"{name}  C→T pos 1: {ct[0]:.4f}  G→A pos 1: {ga[0]:.4f}")

names = list(SAMPLES.keys())
profiles["ct_5p"]    = np.mean([profiles[f"ct_5p_{n}"] for n in names], axis=0)
profiles["ga_3p"]    = np.mean([profiles[f"ga_3p_{n}"] for n in names], axis=0)
profiles["positions"] = pos

np.save(OUT, profiles, allow_pickle=True)
print(f"\nSaved {OUT}  —  {len(pos)} positions")
print(f"  mean C→T pos 1: {profiles['ct_5p'][0]:.4f}")
print(f"  mean G→A pos 1: {profiles['ga_3p'][0]:.4f}")
