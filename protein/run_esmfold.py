"""
Phase 4 — Step 3: Fold protein sequences with ESMFold via the public REST API.

Uses the ESMFold API endpoint (no model download required):
    POST https://api.esmatlas.com/foldSequence/v1/pdb/
    Body: raw protein sequence string
    Response: PDB format text

Folds 9 proteins:
    {gene}_reference  — from protein/sequences/{gene}_protein.fa
    {gene}_mlm        — from protein/sequences/{gene}_mlm_protein.fa
    {gene}_dam        — from protein/sequences/{gene}_dam_protein.fa

Outputs:
    protein/structures/{gene}_{source}.pdb
    protein/structures/plddt_scores.json  — mean pLDDT per structure

Run from project root:
    python3 protein/run_esmfold.py
"""
import json
import ssl
import time
import urllib.request
import urllib.error
from pathlib import Path

# Bypass SSL cert issues on macOS Python.org installs
ssl._create_default_https_context = ssl._create_unverified_context

SEQ_DIR    = Path("protein/sequences")
STRUCT_DIR = Path("protein/structures")
PLDDT_JSON = STRUCT_DIR / "plddt_scores.json"

STRUCT_DIR.mkdir(parents=True, exist_ok=True)

ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
RETRY_WAIT  = 5    # seconds between retries
MAX_RETRIES = 3

GENES   = ["TRPV3", "KCNK9", "HBB"]
SOURCES = {
    "reference": "{gene}_protein.fa",
    "mlm":       "{gene}_mlm_protein.fa",
    "dam":       "{gene}_dam_protein.fa",
}

# TRPV3 is 791 aa — exceeds ESMFold API limit (~600 aa).
# Truncate to N-terminal ankyrin repeat domain (aa 1–400).
# The 13 damaged nucleotide positions in the TRPV3 CDS map to
# aa 1–8 (5′ damage) and aa 784–791 (3′ damage); the N-terminal
# domain captures the 5′-damaged residues.
TRUNCATE = {"TRPV3": 400}


def read_fasta(path: Path) -> str | None:
    if not path.exists():
        return None
    lines = path.read_text().strip().splitlines()
    return "".join(l for l in lines if not l.startswith(">"))


def call_esmfold_api(sequence: str) -> str | None:
    """POST sequence to ESMFold API; return PDB string or None on failure."""
    data = sequence.encode("utf-8")
    req  = urllib.request.Request(
        ESMFOLD_API,
        data=data,
        method="POST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} — {e.reason} (attempt {attempt}/{MAX_RETRIES})")
        except Exception as e:
            print(f"    Error: {e} (attempt {attempt}/{MAX_RETRIES})")
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_WAIT)
    return None


def extract_plddt(pdb_str: str) -> list[float]:
    """
    Extract per-residue pLDDT from B-factor column using Cα atoms.
    ESMFold API returns B-factors in 0–1 scale; multiply by 100 for standard reporting.
    Falls back to first ATOM per residue if Cα not found.
    """
    ca_values = {}
    fallback  = {}
    for line in pdb_str.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        resnum    = int(line[22:26].strip())
        bfactor   = float(line[60:66].strip()) * 100.0   # convert to 0–100 scale
        if atom_name == "CA":
            ca_values[resnum] = bfactor
        if resnum not in fallback:
            fallback[resnum] = bfactor
    # Prefer Cα; fall back to first atom if Cα absent for a residue
    combined = {**fallback, **ca_values}
    return [combined[k] for k in sorted(combined)]


def run():
    plddt_scores = {}

    for gene in GENES:
        for source, template in SOURCES.items():
            fname   = template.format(gene=gene)
            fa_path = SEQ_DIR / fname
            seq     = read_fasta(fa_path)
            label   = f"{gene}_{source}"

            if seq is None:
                print(f"SKIP {label} — {fa_path} not found.")
                continue

            # Truncate oversized sequences to their N-terminal domain
            trunc = TRUNCATE.get(gene)
            if trunc and len(seq) > trunc:
                seq = seq[:trunc]
                print(f"  Truncating {label} to first {trunc} aa (N-terminal domain)")

            pdb_out = STRUCT_DIR / f"{label}.pdb"
            if pdb_out.exists():
                print(f"SKIP {label} — PDB already exists.")
                pdb_str = pdb_out.read_text()
                plddt   = extract_plddt(pdb_str)
                mean_p  = sum(plddt) / len(plddt) if plddt else float("nan")
                plddt_scores[label] = {
                    "mean_plddt":    round(mean_p, 2),
                    "residues":      len(plddt),
                    "plddt_per_res": [round(v, 2) for v in plddt],
                }
                print(f"  {label}: cached  mean pLDDT={mean_p:.2f}")
                continue

            print(f"Folding {label} ({len(seq)} aa) via ESMFold API...")
            pdb_str = call_esmfold_api(seq)

            if pdb_str is None:
                print(f"  FAILED to fold {label}.")
                continue

            pdb_out.write_text(pdb_str)
            plddt   = extract_plddt(pdb_str)
            mean_p  = sum(plddt) / len(plddt) if plddt else float("nan")
            plddt_scores[label] = {
                "mean_plddt":    round(mean_p, 2),
                "residues":      len(plddt),
                "plddt_per_res": [round(v, 2) for v in plddt],
            }
            print(f"  {label}: mean pLDDT={mean_p:.2f}  ({len(plddt)} residues)  → {pdb_out}")
            time.sleep(1)   # be polite to the API

    PLDDT_JSON.write_text(json.dumps(plddt_scores, indent=2))
    print(f"\npLDDT scores → {PLDDT_JSON}")

    # Quick summary
    print("\n── pLDDT summary ─────────────────────────────────────────────")
    for label, data in plddt_scores.items():
        print(f"  {label:<25}  {data['mean_plddt']:>6.2f}  ({data['residues']} aa)")


if __name__ == "__main__":
    run()
