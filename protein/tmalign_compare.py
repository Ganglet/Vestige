"""
Phase 4 — Step 4: Structural comparison of reconstructed vs reference proteins.

Computes TM-score and Cα-RMSD using the standard formulae on the ESMFold PDB files.
No external binary required — implemented in pure Python / NumPy.

TM-score formula (Zhang & Skolnick 2004):
    TM-score = (1/L_ref) * Σ_i [ 1 / (1 + (d_i / d0)^2) ]
    d0 = 1.24 * (L_ref - 15)^(1/3) - 1.8   (clamped to ≥0.5)

The query structure is superposed onto the reference via the Kabsch algorithm
(optimal rigid-body alignment minimising RMSD) before TM-score is calculated.

Outputs:
    protein/table3.txt      — Table 3 for the paper
    protein/table3.json     — machine-readable

Run from project root:
    python3 protein/tmalign_compare.py
"""
import json
import numpy as np
from pathlib import Path

STRUCT_DIR = Path("protein/structures")
PLDDT_JSON = STRUCT_DIR / "plddt_scores.json"
OUT_JSON   = Path("protein/table3.json")
OUT_TXT    = Path("protein/table3.txt")

GENES   = ["TRPV3", "KCNK9", "HBB"]
METHODS = ["mlm", "dam"]


# ── PDB parsing ────────────────────────────────────────────────────────────────

def read_ca_coords(pdb_path: Path) -> np.ndarray:
    """Return (N, 3) array of Cα coordinates from an ATOM record PDB file."""
    coords = []
    seen   = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        resnum = int(line[22:26].strip())
        if resnum in seen:
            continue
        seen.add(resnum)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])
    return np.array(coords, dtype=np.float64)


# ── Kabsch superposition ───────────────────────────────────────────────────────

def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> tuple[float, np.ndarray]:
    """
    Kabsch algorithm: optimally rotate P onto Q.
    Returns (RMSD after superposition, rotated P).
    Both P and Q must have the same shape (N, 3).
    """
    # Centre
    p_c = P - P.mean(axis=0)
    q_c = Q - Q.mean(axis=0)

    # Covariance matrix
    H   = p_c.T @ q_c
    U, S, Vt = np.linalg.svd(H)

    # Ensure proper rotation (handle reflection)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])

    R  = Vt.T @ D @ U.T
    P_rot = p_c @ R.T + Q.mean(axis=0)

    rmsd = float(np.sqrt(((P_rot - Q) ** 2).sum(axis=1).mean()))
    return rmsd, P_rot


# ── TM-score ──────────────────────────────────────────────────────────────────

def tm_score(P_superposed: np.ndarray, Q: np.ndarray, L_ref: int) -> float:
    """
    TM-score of query (P_superposed) vs reference (Q), normalised by L_ref.
    P_superposed and Q must already be aligned (same residue ordering).
    """
    d0 = max(1.24 * (L_ref - 15) ** (1.0 / 3.0) - 1.8, 0.5)
    n  = min(len(P_superposed), len(Q))
    d  = np.sqrt(((P_superposed[:n] - Q[:n]) ** 2).sum(axis=1))
    return float((1.0 / L_ref) * (1.0 / (1.0 + (d / d0) ** 2)).sum())


# ── Main ───────────────────────────────────────────────────────────────────────

def run():
    if not PLDDT_JSON.exists():
        print(f"ERROR: {PLDDT_JSON} not found. Run run_esmfold.py first.")
        return

    plddt   = json.loads(PLDDT_JSON.read_text())
    results = {}

    for gene in GENES:
        ref_pdb = STRUCT_DIR / f"{gene}_reference.pdb"
        if not ref_pdb.exists():
            print(f"SKIP {gene} — reference PDB not found.")
            continue

        Q       = read_ca_coords(ref_pdb)
        L_ref   = len(Q)
        results[gene] = {}

        for method in METHODS:
            query_pdb = STRUCT_DIR / f"{gene}_{method}.pdb"
            label_q   = f"{gene}_{method}"
            label_ref = f"{gene}_reference"

            if not query_pdb.exists():
                print(f"SKIP {gene}/{method} — PDB not found.")
                results[gene][method] = {}
                continue

            P = read_ca_coords(query_pdb)

            # Align to same length
            n = min(len(P), len(Q))
            P_trim, Q_trim = P[:n], Q[:n]

            rmsd, P_rot = kabsch_rmsd(P_trim, Q_trim)
            tms         = tm_score(P_rot, Q_trim, L_ref)

            ref_plddt  = plddt.get(label_ref, {}).get("mean_plddt", float("nan"))
            meth_plddt = plddt.get(label_q,   {}).get("mean_plddt", float("nan"))

            results[gene][method] = {
                "plddt_ref":    round(ref_plddt,  2),
                "plddt_method": round(meth_plddt, 2),
                "tm_score":     round(tms,  4),
                "rmsd_ca":      round(rmsd, 3),
                "aligned_aa":   n,
                "L_ref":        L_ref,
            }
            print(f"{gene:6s} {method}  pLDDT={meth_plddt:5.2f}  "
                  f"TM={tms:.4f}  RMSD={rmsd:.3f} Å  (n={n}/{L_ref})")

    # ── JSON ────────────────────────────────────────────────────────────────────
    OUT_JSON.write_text(json.dumps(results, indent=2))

    # ── Table 3 ─────────────────────────────────────────────────────────────────
    lines = [
        "\nTable 3 — ESMFold structural validation (Phase 4)",
        "=" * 76,
        f"{'Gene':<8}  {'Method':<5}  {'pLDDT(ref)':>10}  {'pLDDT(recon)':>12}  "
        f"{'TM-score':>9}  {'Cα-RMSD':>8}  {'n':>5}",
        "-" * 76,
    ]
    for gene in GENES:
        if gene not in results:
            continue
        for method in METHODS:
            r = results[gene].get(method, {})
            if not r:
                continue
            lines.append(
                f"{gene:<8}  {method:<5}  {r.get('plddt_ref', float('nan')):>10.2f}  "
                f"{r.get('plddt_method', float('nan')):>12.2f}  "
                f"{r.get('tm_score', float('nan')):>9.4f}  "
                f"{r.get('rmsd_ca', float('nan')):>8.3f}  "
                f"{r.get('aligned_aa', 0):>5}"
            )
        lines.append("")

    lines += [
        "=" * 76,
        "",
        "Notes:",
        "  TM-score: Zhang & Skolnick (2004) formula, normalised by reference length",
        "  TM-score > 0.5 = same fold topology; > 0.9 = near-identical structure",
        "  Cα-RMSD after Kabsch superposition",
        "  TRPV3 folded on N-terminal domain (aa 1–400); full-length 791 aa exceeds API limit",
        f"  Damage applied at 10–30× authentic PMD rates (amplified for testability)",
        "  pLDDT in 0–100 scale (ESMFold Cα B-factor × 100)",
    ]

    table = "\n".join(lines) + "\n"
    print(table)
    OUT_TXT.write_text(table)
    print(f"Table 3 → {OUT_TXT}")
    print(f"JSON    → {OUT_JSON}")


if __name__ == "__main__":
    run()
