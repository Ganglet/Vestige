# Damage Profile Parsing & Gene Coordinate Extraction

**Phase:** 1 — Data Acquisition & Validation
**Status:** Complete
**Date:** April 2026

---

## Objective

Parse the mapDamage2 `misincorporation.txt` into numpy arrays for use by the DAM collator, produce Figure 1 (C→T / G→A frequency by position), and extract chromosomal coordinates for the three target coding sequences (TRPV3, KCNK9, HBB) from the EanMak 1.0 annotation.

---

## 1. Damage Profile Parsing

**Script:** `damage/parse_profiles.py`
**Input:** `Dataset/results/results_mammoth_ERR855944/misincorporation.txt`
**Output:** `damage/damage_profile.npy`

The misincorporation table has 8,960 data rows — 4,480 for each terminus (5p / 3p), covering 16 chromosomes × 2 strands × 70 positions per end. Frequencies are aggregated across all chromosomes and strands by grouping on position:

```
freq_C→T[pos] = Σ(C>T counts) / Σ(Total coverage)   [5p rows]
freq_G→A[pos] = Σ(G>A counts) / Σ(Total coverage)   [3p rows]
```

`damage_profile.npy` is a pickled dict with three keys:

| Key | Shape | Contents |
|-----|-------|----------|
| `ct_5p` | (70,) | C→T frequency, positions 1–70 from 5′ end |
| `ga_3p` | (70,) | G→A frequency, positions 1–70 from 3′ end |
| `positions` | (70,) | Position index array (1..70) |

**Load:**
```python
profile = np.load("damage/damage_profile.npy", allow_pickle=True).item()
```

---

## 2. Damage Profile Values

| Position | C→T (5′) | G→A (3′) |
|----------|----------|----------|
| 1 | 0.0029 | 0.0031 |
| 2 | 0.0020 | 0.0020 |
| 3 | 0.0016 | 0.0016 |
| 4 | 0.0013 | 0.0014 |
| 5 | 0.0012 | 0.0012 |
| 20–70 (mean) | 0.0010 | 0.0010 |

End-damage elevates substitution frequency ~3× above background at position 1 and decays to baseline by position 10. The signal is real but moderate — consistent with permafrost-preserved material from a relatively young specimen.

---

## 3. Figure 1 — Damage Frequency Plot

**Script:** `damage/visualize_damage.py`
**Input:** `damage/damage_profile.npy`
**Output:** `results/figures/fig1_damage_profile.pdf` + `.png`

Two-panel figure: C→T (red, left) and G→A (blue, right) vs. read-end position. Replaces the failed R-generated `Fragmisincorporation_plot.pdf` from mapDamage2 (R ≥ 4.x incompatibility — see Phase 1 doc).

---

## 4. Gene Coordinate Extraction

**Script:** `damage/extract_genes.py`
**Method:** NCBI Entrez (`Bio.Entrez`) — avoids downloading the ~3 GB GFF3 annotation
**Output:** `damage/gene_coords.csv`

The NCBI datasets CLI at `/usr/local/bin/datasets` was non-functional (XML stub installed, not the actual binary), so the GFF-based approach was replaced with direct Entrez queries against the NCBI Gene database.

### Coordinates (EanMak 1.0 / GCF_024166365.1, 1-based inclusive)

| Gene | NCBI Symbol | Chromosome | Start | End | Strand | Length |
|------|-------------|------------|-------|-----|--------|--------|
| TRPV3 | TRPV3 | NC_064837.1 (chr 19) | 17,626,780 | 17,673,723 | + | 46,944 bp |
| KCNK9 | KCNK9 | NC_064833.1 (chr 15) | 4,939,078 | 4,960,091 | + | 21,014 bp |
| HBB | LOC126080006 | NC_064825.1 (chr 7) | 57,243,805 | 57,249,625 | − | 5,821 bp |

### Note on HBB annotation

EanMak 1.0 does not carry a gene annotated as `HBB`. The locus is `LOC126080006`, described as "hemoglobin subunit beta-1/2-like." This is expected: in proboscideans the ancestral β-globin and δ-globin genes fused into a single HBB/HBD hybrid locus (Opazo et al. 2009; Signore et al. 2019). The Asian Elephant genome carries this fused locus rather than separate HBB and HBD entries, so no HGNC-standard HBB symbol is assigned. For the purposes of this project, `LOC126080006` is the correct HBB ortholog and the naming difference is noted wherever the gene is referenced.

---

## 5. Phase 1 Deliverables — Status

| Deliverable | File | Status |
|-------------|------|--------|
| Sorted, indexed BAM | `Dataset/mammoth/mammoth_ERR855944.bam` | ✓ Done |
| mapDamage2 misincorporation table | `Dataset/results/results_mammoth_ERR855944/misincorporation.txt` | ✓ Done |
| Damage profile arrays | `damage/damage_profile.npy` | ✓ Done |
| Figure 1 | `results/figures/fig1_damage_profile.pdf` | ✓ Done |
| Gene coordinates | `damage/gene_coords.csv` | ✓ Done |

**Phase 2 (DAM collator implementation) can begin.**

→ See `03_dam_collator.md`
