# Damage Profile Parsing & Gene Coordinate Extraction

**Phase:** 1 — Data Acquisition & Validation
**Status:** Complete
**Date:** April 2026

---

## Objective

Parse mapDamage2 `misincorporation.txt` outputs from both specimens into numpy arrays for use by the DAM collator, produce Figure 1 (C→T / G→A frequency by position), and extract chromosomal coordinates for the three target coding sequences (TRPV3, KCNK9, HBB) from the EanMak 1.0 annotation.

---

## 1. Alignment — ERR852028

ERR852028 was aligned using identical parameters to ERR855944 (BWA-aln, `-l 16500 -n 0.01 -t 4`, subsampled to 5M reads, seed 42). Both specimens are from Palkopoulou et al. 2015.

**Alignment statistics:**

| Accession | Reads | Mapped | Mapping rate |
|-----------|-------|--------|--------------|
| ERR855944 | 5,000,000 | 4,956,238 | 99.12% |
| ERR852028 | 5,000,000 | 4,966,693 | 99.33% |

Both > 99% — consistent with the known ~6 Myr divergence from the Asian Elephant reference.

---

## 2. Damage Profile Parsing

**Script:** `damage/parse_profiles.py`
**Inputs:** `Dataset/results/results_mammoth_ERR855944/misincorporation.txt`, `Dataset/results/results_mammoth_ERR852028/misincorporation.txt`
**Output:** `damage/damage_profile.npy`

Per-specimen frequencies are aggregated across all chromosomes and strands by grouping on position:

```
freq_C→T[pos] = Σ(C>T counts) / Σ(Total coverage)   [5p rows]
freq_G→A[pos] = Σ(G>A counts) / Σ(Total coverage)   [3p rows]
```

The output file stores individual profiles for both specimens plus their mean. The DAM collator reads `ct_5p` and `ga_3p` (the mean arrays).

`damage_profile.npy` is a pickled dict:

| Key | Shape | Contents |
|-----|-------|----------|
| `ct_5p` | (70,) | C→T frequency — mean of both specimens (used by DAM) |
| `ga_3p` | (70,) | G→A frequency — mean of both specimens (used by DAM) |
| `ct_5p_ERR855944` | (70,) | C→T, ERR855944 individual |
| `ga_3p_ERR855944` | (70,) | G→A, ERR855944 individual |
| `ct_5p_ERR852028` | (70,) | C→T, ERR852028 individual |
| `ga_3p_ERR852028` | (70,) | G→A, ERR852028 individual |
| `positions` | (70,) | Position index array (1..70) |

**Load:**
```python
profile = np.load("damage/damage_profile.npy", allow_pickle=True).item()
```

---

## 3. Damage Profile Values

| Position | ERR855944 C→T | ERR852028 C→T | ERR855944 G→A | ERR852028 G→A |
|----------|--------------|--------------|--------------|--------------|
| 1 | 0.0029 | 0.0041 | 0.0031 | 0.0042 |
| 2 | 0.0020 | 0.0026 | 0.0020 | 0.0028 |
| 3 | 0.0016 | 0.0021 | 0.0016 | 0.0022 |
| 4 | 0.0013 | 0.0018 | 0.0014 | 0.0019 |
| 5 | 0.0012 | 0.0016 | 0.0012 | 0.0016 |
| 20–70 (mean) | 0.0010 | 0.0010 | 0.0010 | 0.0010 |

End-damage elevates substitution frequency ~3–4× above background at position 1 and decays to a shared baseline of 0.0010 by position ~10. ERR852028 shows slightly higher terminal damage, but both specimens converge to identical background — the damage patterns are consistent across both specimens and follow the canonical aDNA deamination signature.

---

## 4. Figure 1 — Damage Frequency Plot

**Script:** `damage/visualize_damage.py`
**Input:** `damage/damage_profile.npy`
**Output:** `results/figures/fig1_damage_profile.pdf` + `.png`

Two-panel figure: C→T (red, left) and G→A (blue, right) vs. read-end position. Both specimens are plotted as individual lines with the mean overlaid as a solid darker line. Replaces the failed R-generated `Fragmisincorporation_plot.pdf` from mapDamage2 (R ≥ 4.x incompatibility — see Phase 1 doc).

---

## 5. Gene Coordinate Extraction

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

## 6. Phase 1 Deliverables — Status

| Deliverable | File | Status |
|-------------|------|--------|
| Sorted, indexed BAM (ERR855944) | `Dataset/mammoth/mammoth_ERR855944.bam` | ✓ Done |
| Sorted, indexed BAM (ERR852028) | `Dataset/mammoth/mammoth_ERR852028.bam` | ✓ Done |
| mapDamage2 tables | `Dataset/results/results_mammoth_ERR855944/` + `ERR852028/` | ✓ Done |
| Damage profile arrays (both specimens + mean) | `damage/damage_profile.npy` | ✓ Done |
| Figure 1 | `results/figures/fig1_damage_profile.pdf` | ✓ Done |
| Gene coordinates | `damage/gene_coords.csv` | ✓ Done |

**Phase 2 (DAM collator implementation) can begin.**

→ See `03_dam_collator.md`
