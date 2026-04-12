"""
Phase 4 — Step 1: Fetch spliced mRNA CDS from NCBI for TRPV3, KCNK9, HBB
in Elephas maximus, translate to protein.

Strategy:
  1. esearch gene database for each gene name + Elephas maximus → get gene ID
  2. elink gene → refseq mRNA (gene_nuccore_refseqrna)
  3. efetch top mRNA accession in GenBank format → extract CDS feature

Outputs (protein/sequences/):
    {gene}_cds.fa          spliced coding sequence (ATG … stop codon inclusive)
    {gene}_protein.fa      translated amino acid sequence

Run from project root:
    NCBI_EMAIL=your@email.com python3 protein/fetch_cds.py
"""
import os
import ssl
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# macOS Python.org installs lack system CA certs — bypass SSL for NCBI queries
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = os.environ.get("NCBI_EMAIL", "your.email@example.com")

OUT_DIR = Path("protein/sequences")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Gene label → NCBI gene search query (taxid 9783 = Elephas maximus)
GENE_QUERIES = {
    "TRPV3": "TRPV3[Gene Name] AND txid9783[Organism:exp]",
    "KCNK9": "KCNK9[Gene Name] AND txid9783[Organism:exp]",
    "HBB":   "(HBB[Gene Name] OR LOC126080006[Gene Name] OR hemoglobin subunit beta[Gene Name]) AND txid9783[Organism:exp]",
}


def gene_search(query: str) -> list[str]:
    """Search NCBI gene database; return list of gene UIDs."""
    handle = Entrez.esearch(db="gene", term=query, retmax=5)
    rec = Entrez.read(handle)
    handle.close()
    return rec["IdList"]


def elink_gene_to_mrna(gene_uid: str) -> list[str]:
    """Return nucleotide UIDs for RefSeq mRNA linked to this gene."""
    handle = Entrez.elink(
        dbfrom="gene", db="nucleotide",
        id=gene_uid, linkname="gene_nuccore_refseqrna",
    )
    records = Entrez.read(handle)
    handle.close()
    uids = []
    for linkset in records:
        for db_link in linkset.get("LinkSetDb", []):
            for link in db_link.get("Link", []):
                uids.append(link["Id"])
    return uids


def fetch_genbank(uid: str):
    """Fetch a GenBank record by UID; return SeqRecord or None."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text")
        rec = SeqIO.read(handle, "genbank")
        handle.close()
        return rec
    except Exception as e:
        print(f"    efetch failed for {uid}: {e}")
        return None


def extract_cds(record) -> tuple[str, str] | None:
    """Extract first CDS feature → (cds_nt, protein_aa) or None."""
    for feat in record.features:
        if feat.type != "CDS":
            continue
        cds_seq = str(feat.extract(record.seq)).upper()
        if "translation" in feat.qualifiers:
            protein = feat.qualifiers["translation"][0]
        else:
            protein = str(Seq(cds_seq).translate(to_stop=True))
        return cds_seq, protein
    return None


def run():
    for gene, query in GENE_QUERIES.items():
        print(f"\n── {gene} ───────────────────────────────────────────")

        # Step 1: find gene UID
        gene_uids = gene_search(query)
        time.sleep(0.4)
        if not gene_uids:
            print(f"  No gene hits for: {query}")
            continue
        gene_uid = gene_uids[0]
        print(f"  Gene UID: {gene_uid}")

        # Step 2: elink gene → RefSeq mRNA
        mrna_uids = elink_gene_to_mrna(gene_uid)
        time.sleep(0.4)
        if not mrna_uids:
            print(f"  No linked RefSeq mRNA for gene UID {gene_uid}.")
            continue
        print(f"  Linked mRNA UIDs: {mrna_uids[:5]} ({len(mrna_uids)} total)")

        # Step 3: efetch each mRNA until we get a usable CDS
        cds_nt = protein_aa = None
        for uid in mrna_uids[:5]:
            print(f"  Fetching UID {uid}...")
            record = fetch_genbank(uid)
            time.sleep(0.4)
            if record is None:
                continue
            # Skip chromosomal records
            if record.id.startswith("NC_") or record.id.startswith("NW_"):
                print(f"    Skipping genomic record {record.id}")
                continue
            result = extract_cds(record)
            if result is not None:
                cds_nt, protein_aa = result
                print(f"  Accession: {record.id}  |  CDS: {len(cds_nt)} nt  |  Protein: {len(protein_aa)} aa")
                break
            else:
                print(f"    No CDS in {record.id}")

        if cds_nt is None:
            print(f"  WARNING: could not get CDS for {gene}.")
            continue

        (OUT_DIR / f"{gene}_cds.fa").write_text(f">{gene}_cds_Elephas_maximus\n{cds_nt}\n")
        protein_clean = protein_aa.rstrip("*")
        (OUT_DIR / f"{gene}_protein.fa").write_text(f">{gene}_protein_Elephas_maximus\n{protein_clean}\n")
        print(f"  Saved {gene}_cds.fa  +  {gene}_protein.fa  ({len(protein_clean)} aa)")

    print("\nDone. Check protein/sequences/")


if __name__ == "__main__":
    run()
