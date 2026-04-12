# Phase 4 — CDS Translation & ESMFold Structural Validation

**Phase:** 4 — Protein Structural Validation
**Status:** In progress
**Date:** April 2026

---

## Objective

Provide a second, modality-orthogonal line of evidence: does DAM-reconstructed DNA produce more structurally plausible proteins than MLM-reconstructed DNA? Sequence recovery (Phase 3) measures nucleotide accuracy; structural validation measures whether the downstream protein retains its fold.

Three genes: TRPV3, KCNK9, HBB — same as fine-tuning.

---

## Why spliced mRNA, not the genomic windows

The fine-tuning dataset used genomic intervals from `gene_coords.csv` (full chromosomal spans including introns). Those cannot be translated directly. Phase 4 fetches the **spliced mRNA CDS** from NCBI RefSeq for each gene in *Elephas maximus* (GCF_024166365.1 / EanMak 1.0). The models generalise: DNABERT-2 encodes DNA agnostically of splicing, and the exonic regions are a subset of what the fine-tuning data covered.

---

## Pipeline

```
protein/fetch_cds.py           NCBI RefSeq mRNA → {gene}_cds.fa + {gene}_protein.fa
        │
        ▼
protein/damage_reconstruct_cds.py
        │  Apply PMD damage at 10× authentic rates (amplified for testability)
        │  Tile CDS → 2000 bp windows → run MLM / DAM → stitch → translate
        │
        ├── {gene}_mlm_protein.fa
        └── {gene}_dam_protein.fa
        │
        ▼
protein/run_esmfold.py         ESMFold (facebook/esmfold_v1) → 9 PDB files + pLDDT
        │
        ▼
protein/tmalign_compare.py     TMalign (reconstructed vs reference) → Table 3
```

---

## 1. CDS Fetching (`protein/fetch_cds.py`)

Searches NCBI nucleotide database for RefSeq mRNA records (NM_* or XM_*) for each gene in *Elephas maximus* using Biopython Entrez. Extracts the CDS feature from the GenBank record and translates using `/translation` qualifier where available.

**Output:** `protein/sequences/{gene}_cds.fa`, `protein/sequences/{gene}_protein.fa`

---

## 2. Damage, Reconstruction, Translation (`protein/damage_reconstruct_cds.py`)

**Damage:** Authentic PMD profiles from `damage/damage_profile.npy`, scaled by `DAMAGE_SCALE=10`. At 10× authentic rates, peak C→T probability at position 1 is ~3.5% (vs 0.35% authentic). This ensures ≥10 damaged positions per gene for a meaningful comparison.

**Reconstruction:** For each damaged nucleotide position `p`:
1. Identify the "best" window — the 2000 bp window in which `p` is nearest the centre (position 1000). This minimises edge effects.
2. Mask the BPE token spanning `p` in that window.
3. Run the model; decode the predicted nucleotide at `p`'s sub-token offset.

**Stitching:** Replace damaged positions in the damaged CDS string with the model's predicted nucleotide. Undamaged positions are unchanged (identical for both models — same source CDS).

**Translation:** `Bio.Seq.translate(to_stop=True)` on the stitched CDS. Trimmed to complete codons.

---

## 3. ESMFold (`protein/run_esmfold.py`)

Model: `facebook/esmfold_v1` (~2.7 GB). Runs on CPU (MPS lacks some required operations). Folds 9 protein sequences (3 genes × reference + MLM + DAM). pLDDT extracted from B-factor column of ATOM records.

| Gene | Reference length | Expected fold time (CPU, M1 Pro) |
|------|-----------------|----------------------------------|
| HBB  | ~146 aa         | ~5 min                           |
| KCNK9| ~335 aa         | ~10 min                          |
| TRPV3| ~760 aa         | ~20 min                          |

Total: ~35 min estimated.

---

## 4. Structural Comparison (`protein/tmalign_compare.py`)

**TM-align** (Zhang lab binary): pairwise alignment of reconstructed vs reference PDB. TM-score normalised by reference structure length.

**Metrics:**
- Mean pLDDT — ESMFold confidence per structure (from B-factor column)
- TM-score — structural similarity to reference (>0.5 = same fold family)
- RMSD — Cα RMSD in Ångström
- Aligned length

---

## 5. Results

**Status:** Complete.

### Table 3 — ESMFold structural validation

| Gene | Method | pLDDT (ref) | pLDDT (recon) | TM-score | Cα-RMSD (Å) | n |
|------|--------|-------------|---------------|----------|-------------|---|
| TRPV3† | MLM | 77.37 | 77.44 | 0.9808 | 1.092 | 400 |
| TRPV3† | DAM | 77.37 | 77.44 | 0.9808 | 1.092 | 400 |
| KCNK9 | MLM | 67.94 | 67.02 | 0.9520 | 1.732 | 347 |
| KCNK9 | DAM | 67.94 | 67.02 | 0.9520 | 1.732 | 347 |
| HBB | MLM | 79.70 | 80.78 | 0.9728 | 0.815 | 162 |
| HBB | DAM | 79.70 | 81.34 | 0.9714 | 0.840 | 162 |

† TRPV3 full-length (791 aa) exceeds ESMFold API limit (~600 aa); folded on N-terminal ankyrin repeat domain (aa 1–400). The 5′-damaged residues (aa 1–8) are within this domain.

**Amino acid changes (reconstructed vs reference):**

| Gene | MLM aa Δ | DAM aa Δ |
|------|----------|----------|
| TRPV3 | 7 / 791 | 10 / 791 |
| KCNK9 | 3 / 347 | 3 / 347 (identical proteins) |
| HBB   | 3 / 162 | 4 / 162 |

All substitutions are conservative (R↔K, E↔D, S↔T, G↔D type).

### Interpretation

All TM-scores > 0.95 — all reconstructed sequences maintain the reference fold topology. MLM and DAM are structurally indistinguishable at this level of sequence perturbation (3–10 conservative aa changes in 162–791 aa proteins). pLDDT within 2 points of reference for all structures.

**For the paper:** Phase 4 is a fold integrity check, not a ranking. The claim is: *"Both reconstructed sequences maintain the reference protein topology (TM-score > 0.95 for all three genes), confirming that neither model introduces fold-disrupting substitutions under simulated ancient DNA damage."* The structural comparison does not discriminate MLM from DAM — that distinction comes from Phase 3 (p=0.041, terminal position recovery). Phase 4 validates that reconstruction is biologically viable.

---

## 6. Run order

```bash
cd Phoenix_Project

# Step 1 — fetch mRNA CDS from NCBI (requires internet, ~1 min)
NCBI_EMAIL=your@email.com python3 protein/fetch_cds.py

# Step 2 — damage + reconstruct + translate (~5-10 min on MPS)
python3 protein/damage_reconstruct_cds.py

# Step 3 — ESMFold (~35 min on CPU)
python3 protein/run_esmfold.py

# Step 4 — TM-align comparison (requires TMalign binary on PATH)
python3 protein/tmalign_compare.py
```

**TMalign install:**
```bash
conda install -c bioconda tmalign
# or: download binary from https://zhanggroup.org/TM-align/
```

---

## 7. Phase 4 Checklist

- [x] Fetch mRNA CDS for TRPV3, KCNK9, HBB — `protein/fetch_cds.py`
- [x] Damage + reconstruct + translate — `protein/damage_reconstruct_cds.py`
- [x] ESMFold structural prediction — `protein/run_esmfold.py`
- [x] TM-align comparison + Table 3 — `protein/tmalign_compare.py`
- [x] Fill in Table 3 above with actual results
- [x] Update `problems_and_decisions.md` with any issues (P17–P21)

→ See `08_biosecurity.md` (Phase 5 — Biosecurity Classifier + Manuscript)
