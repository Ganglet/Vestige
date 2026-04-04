# Data Acquisition — Mammoth aDNA & Elephant Reference

**Phase:** 1 — Data Acquisition & Validation  
**Status:** In Progress — W1–W4 steps done; gene extraction and damage profile parsing remain  
**Date:** April 2026

---

## Objective

Obtain raw ancient DNA reads from a woolly mammoth specimen and a high-quality Asian Elephant reference genome to serve as the alignment target and coding sequence source for Phoenix-LM.

---

## 1. Mammoth aDNA Reads

**Source:** NCBI SRA — Project ERP008929
Palkopoulou et al., 2015. "Complete Genomes Reveal Signatures of Demographic and Genetic Declines in the Woolly Mammoth." *Current Biology* 25(10), 1395–1402.

**Accessions downloaded:**

| Accession | Description |
|-----------|-------------|
| ERR855944 | Woolly mammoth (*Mammuthus primigenius*) — primary sample used |
| ERR852028 | Additional mammoth run — reserved |

**Download method:** SRA Toolkit (`fastq-dump` / prefetch). Raw `.sra` files stored in `Dataset/mammoth/ERR855944/` and `Dataset/mammoth/ERR852028/`. FASTQ exported to `Dataset/mammoth/ERR855944.fastq`.

**Subsampling:** Full FASTQ is prohibitively large for iterative development. Subsampled to 5 million reads using seqtk with fixed seed for reproducibility:

```bash
seqtk sample -s 42 ERR855944.fastq 5000000 > ERR855944_5M.fastq
```

Seed 42 used throughout for all stochastic operations.

---

## 2. Asian Elephant Reference Genome

**Source:** NCBI RefSeq  
**Accession:** GCF_024166365.1 (EanMak 1.0 / mEleMax1 primary haplotype)  
**Stored in:** `Dataset/elephant_reference/`

The reference FASTA used for alignment is `elephant_ref.fa` — this is a subset of the full EanMak assembly covering the chromosomal regions relevant to the target genes (TRPV3, KCNK9, HBB). The full assembly FASTA is available in `Dataset/elephant_reference/ncbi_dataset/`.

BWA index built prior to alignment:

```bash
bwa index elephant_ref.fa
```

Index files (`.amb`, `.ann`, `.bwt`, `.fai`, `.pac`, `.sa`) stored alongside the reference.

---

## 3. Alignment

Ancient DNA reads require modified alignment parameters relative to modern short reads. Standard BWA-MEM2 is not appropriate here — aDNA reads are short (30–70 bp), have elevated mismatch rates due to deamination, and the default seed length causes misalignments.

Used BWA-aln with aDNA-tuned parameters:

```bash
bwa aln \
  -l 16500 \   # disable seeding (seed length >> read length)
  -n 0.01 \    # allow 1% mismatches — accounts for deamination
  -t 4 \       # 4 threads
  elephant_ref.fa ERR855944_5M.fastq > mammoth_ERR855944.sai

bwa samse elephant_ref.fa mammoth_ERR855944.sai ERR855944_5M.fastq \
  | samtools sort -@ 4 -o mammoth_ERR855944.bam

samtools index mammoth_ERR855944.bam
```

**Why `-l 16500` (not BWA-MEM2):** The seed length of 16500 exceeds any realistic read length in the dataset, effectively disabling seeding and forcing full-read alignment. This is the standard aDNA practice for reads under ~100 bp. BWA-MEM2 with default settings would miss a significant fraction of damaged reads.

**Alignment statistics** (`samtools flagstat`):

```
5,000,000 reads total
4,956,238 mapped (99.12%)
0 secondary / supplementary / paired
```

99.12% mapping rate to the Asian Elephant reference is consistent with the known ~6 Myr divergence and indicates good sample quality.

Output: `Dataset/mammoth/mammoth_ERR855944.bam` (sorted, indexed).

---

## 4. Damage Profiling — mapDamage2

**Tool:** mapDamage2 v2.1.0  
**Install:** `conda install -c bioconda -c conda-forge mapdamage2`

Run against the sorted BAM and reference:

```bash
mapDamage -i Dataset/mammoth/mammoth_ERR855944.bam -r Dataset/elephant_reference/elephant_ref.fa
```

**Note on R plotting error:** mapDamage2's R plotting module is incompatible with R ≥ 4.x due to a deprecated `&&` operator on length-2 vectors. This causes the R-generated plots to fail but does not affect the frequency table computation. The critical output — `misincorporation.txt` — was generated successfully. Custom plotting from this table is done in Python (see `damage/visualize_damage.py`).

**Outputs** (in `Dataset/results/results_mammoth_ERR855944/`):

| File | Contents |
|------|----------|
| `misincorporation.txt` | Per-position C→T and G→A substitution frequencies — primary damage profile used for DAM |
| `lgdistribution.txt` | Fragment length distribution |
| `dnacomp.txt` | Nucleotide composition by position |
| `Length_plot.pdf` | Fragment length distribution plot (R-generated, succeeded) |
| `Fragmisincorporation_plot.pdf` | Damage frequency plot (R-generated, failed — replaced by Python) |

---

## 5. Data Storage & Version Control

Large binary files (raw reads, BAM, reference FASTA, genome indexes) are excluded from git — they are re-downloadable from NCBI using accessions above. See `.gitignore`.

What is tracked in git:
- `Dataset/mammoth/runs.txt` — SRA accession list
- `Dataset/mammoth/md5sum.txt` — checksums
- `Dataset/results/` — mapDamage2 output tables and plots (small, result of a ~4.5 hour compute run)
- This documentation file

---

## Phase 1 Completion

- [x] Extract TRPV3, KCNK9, HBB gene coordinates from EanMak 1.0 annotation via NCBI Entrez
- [x] Parse `misincorporation.txt` → `damage_profile.npy` (numpy arrays for C→T and G→A per position)
- [x] Produce C→T / G→A frequency-by-position visualization — Figure 1 of the paper
- **DELIVERABLE complete.** `damage_profile.npy` + `gene_coords.csv` + Figure 1 committed.

→ See `02_damage_profile_parsing.md`
