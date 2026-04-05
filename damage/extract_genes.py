"""
Extract TRPV3, KCNK9, HBB coding sequence coordinates for Elephas maximus
(EanMak 1.0 / GCF_024166365.1) via NCBI Entrez.

Outputs damage/gene_coords.csv with columns: gene, ncbi_name, chr, start, end, strand

Notes:
  - HBB is unannotated under that symbol in GCF_024166365.1; retrieved as
    LOC126080006 ("hemoglobin subunit beta-1/2-like") via species-level search.
  - Coordinates are 1-based, inclusive (BED-style conversion applied).
  - Minus-strand genes: ChrStart > ChrStop in NCBI summary; we swap to report
    smaller coord as start.

Run from project root: python damage/extract_genes.py
"""
import csv
import os
import sys
import time
from pathlib import Path
from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "your.email@example.com")

# gene label → Entrez query term (all for Elephas maximus, taxid 9783)
QUERIES = {
    "TRPV3": "TRPV3[Gene Name] AND Elephas maximus[Organism]",
    "KCNK9": "KCNK9[Gene Name] AND Elephas maximus[Organism]",
    "HBB":   "hemoglobin subunit beta[Description] AND txid9783[Organism:exp]",
}
OUT = Path("damage/gene_coords.csv")


def fetch_location(label: str, query: str) -> dict | None:
    handle = Entrez.esearch(db="gene", term=query)
    rec = Entrez.read(handle); handle.close()
    if not rec["IdList"]:
        print(f"  {label}: no hits")
        return None

    gene_id = rec["IdList"][0]
    handle = Entrez.esummary(db="gene", id=gene_id)
    summary = Entrez.read(handle); handle.close()

    try:
        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        loc = doc["GenomicInfo"][0]
        raw_start = int(loc["ChrStart"])
        raw_stop  = int(loc["ChrStop"])
        # NCBI: if ChrStart > ChrStop the gene is on the minus strand
        if raw_start <= raw_stop:
            start, end, strand = raw_start + 1, raw_stop, "+"
        else:
            start, end, strand = raw_stop + 1, raw_start, "-"
        return {
            "gene":      label,
            "ncbi_name": doc["Name"],
            "chr":       loc["ChrAccVer"],
            "start":     start,
            "end":       end,
            "strand":    strand,
        }
    except (KeyError, IndexError) as e:
        print(f"  {label}: parse error — {e}")
        return None


found = []
for label, query in QUERIES.items():
    print(f"Querying {label}...")
    result = fetch_location(label, query)
    if result:
        found.append(result)
        print(f"  {result['ncbi_name']}  {result['chr']}:{result['start']}-{result['end']}  ({result['strand']})")
    time.sleep(0.4)

if not found:
    print("No genes retrieved. Check NCBI_EMAIL and network access.")
    sys.exit(1)

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene", "ncbi_name", "chr", "start", "end", "strand"])
    w.writeheader()
    w.writerows(found)

print(f"\nSaved {OUT}")
