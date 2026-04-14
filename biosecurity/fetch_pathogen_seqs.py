"""
Fetch pathogen virulence gene sequences for biosecurity classifier training.

Uses non-Select-Agent pathogens with unrestricted public sequences.
Fetches by direct accession number — more reliable than text queries.

Positive class (virulence genes):
  - E. coli O157:H7    stx1/stx2   (Shiga toxins — standard biosecurity ref)
  - S. aureus          mecA / spa  (MRSA resistance + surface protein)
  - Listeria mono.     hlyA        (listeriolysin O, pore-forming toxin)
  - Salmonella Typhi   sipA / invA (invasion machinery)
  - V. cholerae        ctxA/ctxB   (cholera toxin)
  - Pseudomonas aerug. exoS/exoU   (type III secretion effectors)

Output: biosecurity/pathogen_seqs.fasta

Run:
    NCBI_EMAIL=your@email.com python3 biosecurity/fetch_pathogen_seqs.py
"""
import os, sys, time
from pathlib import Path
from Bio import Entrez, SeqIO

EMAIL = os.environ.get("NCBI_EMAIL", "")
if not EMAIL:
    sys.exit("Set NCBI_EMAIL environment variable before running.")

Entrez.email = EMAIL
Entrez.tool  = "VESTIGE-biosecurity"

OUT_FASTA = Path("biosecurity/pathogen_seqs.fasta")
MIN_LEN   = 200
MAX_LEN   = 6000

# Direct accession numbers — reliable, unrestricted, peer-reviewed sequences
ACCESSIONS = {
    # E. coli O157:H7 — Shiga toxins (widely used biosecurity benchmark)
    "Ecoli_stx1":   ["X07865",  "AF461168", "AF500189"],
    "Ecoli_stx2":   ["AY286000","AF043627", "EU078440"],
    # Staphylococcus aureus — MRSA virulence
    "Saureus_mecA": ["Y00688",  "AB505628", "AF033191"],
    "Saureus_spa":  ["X52143",  "M18104"],
    # Listeria monocytogenes — hlyA (listeriolysin O)
    "Listeria_hlyA":["X15127",  "M28372",   "AY878649"],
    # Salmonella — invasion genes
    "Salmonella_sipA": ["AF043295", "AE006468"],
    "Salmonella_invA": ["M90846",   "L04718"],
    # Vibrio cholerae — cholera toxin (already in public RefSeq)
    "Vcholerae_ctxA":  ["X00171",   "AF325734"],
    "Vcholerae_ctxB":  ["M35586",   "AE003852"],
    # Pseudomonas aeruginosa — type III effectors
    "Pseudo_exoS":  ["AF013406", "U88545"],
    "Pseudo_exoU":  ["AF043312", "U88545"],
}


def fetch_by_accessions(accession_list: list[str], label: str) -> list[tuple[str, str]]:
    seqs = []
    for acc in accession_list:
        try:
            handle  = Entrez.efetch(db="nucleotide", id=acc,
                                    rettype="fasta", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            for rec in records:
                s = str(rec.seq).upper().replace("N", "")
                if MIN_LEN <= len(s) <= MAX_LEN and set(s) <= {"A", "T", "G", "C"}:
                    seqs.append((f"{label}|{acc}", s))
            time.sleep(0.35)   # stay under NCBI 3 req/sec limit
        except Exception as e:
            print(f"    ✗ {acc}: {e}")
    return seqs


def main():
    OUT_FASTA.parent.mkdir(exist_ok=True)

    all_seqs = []
    for label, accessions in ACCESSIONS.items():
        organism, gene = label.split("_", 1)
        print(f"  {organism} / {gene} ...", end=" ", flush=True)
        seqs = fetch_by_accessions(accessions, label)
        print(f"{len(seqs)} sequences")
        all_seqs.extend(seqs)

    with open(OUT_FASTA, "w") as fh:
        for header, seq in all_seqs:
            fh.write(f">{header}\n{seq}\n")

    print(f"\nTotal pathogen sequences: {len(all_seqs)}")
    print(f"Saved → {OUT_FASTA}")
    if len(all_seqs) < 20:
        print("WARNING: fewer than 20 sequences — check network/accessions")


if __name__ == "__main__":
    main()
