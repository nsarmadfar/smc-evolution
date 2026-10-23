#!/usr/bin/env python3
"""Download proteomes, reference sequences and outgroups from UniProt.

Proteomes: the full UniProt reference proteome of each ingroup species.
  SMC proteins are found later by sequence similarity (02_classify.py), not
  by annotation, because many SMC proteins in non-model species are
  annotated inconsistently.

References: reviewed human and budding yeast SMC1-6 and RAD50. They are the
  search queries and are used to label the hits.

Outgroups: the reviewed bacterial or archaeal Smc of each outgroup species.

Output (in OUTDIR): references.tsv, outgroups.tsv, proteomes/<code>.faa.gz

Usage: python scripts/01_fetch_uniprot.py config/species.tsv data/raw
"""
import gzip
import io
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

API = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,id,gene_primary,organism_id,organism_name,length,sequence"

REFERENCE_QUERIES = {
    "human": ("9606", ["SMC1A", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6", "RAD50"]),
    "yeast": ("559292", ["SMC1", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6", "RAD50"]),
}


def get(params, retries=4):
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(API, params=params, timeout=600)
            r.raise_for_status()
            return r
        except requests.RequestException as err:
            print(f"  attempt {attempt} failed: {err}", file=sys.stderr)
            time.sleep(5 * attempt)
    raise RuntimeError(f"UniProt query failed: {params['query']}")


def uniprot_tsv(query):
    text = get({"query": query, "format": "tsv", "fields": FIELDS}).text
    if not text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(text), sep="\t")


def family_of(gene):
    gene = str(gene).upper()
    if gene.startswith("RAD50"):
        return "RAD50"
    if gene.startswith("SMC") and len(gene) >= 4 and gene[3].isdigit():
        return "SMC" + gene[3]
    return None


def read_fasta(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def fetch_proteome(taxon_id, path):
    """Download a reference proteome; keep one strain (the largest)."""
    r = get({"query": f"taxonomy_id:{taxon_id} AND keyword:KW-1185",
             "format": "fasta", "compressed": "true"})
    records = list(read_fasta(gzip.decompress(r.content).decode()))
    if not records:
        return 0, None
    ox = [re.search(r"OX=(\d+)", h).group(1) for h, _ in records]
    top = Counter(ox).most_common(1)[0][0]
    with gzip.open(path, "wt") as fh:
        n = 0
        for (header, seq), strain in zip(records, ox):
            if strain == top:
                fh.write(f">{header}\n{seq}\n")
                n += 1
    return n, top


def main(species_file, outdir):
    outdir = Path(outdir)
    (outdir / "proteomes").mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(species_file, sep="\t", dtype=str)

    # 1. References
    refs = []
    for label, (taxon, genes) in REFERENCE_QUERIES.items():
        gene_q = " OR ".join(f"gene_exact:{g}" for g in genes)
        df = uniprot_tsv(f"({gene_q}) AND organism_id:{taxon} AND reviewed:true")
        df["family"] = df["Gene Names (primary)"].map(family_of)
        refs.append(df.dropna(subset=["family"]))
        print(f"references ({label}): {len(df)} sequences")
    pd.concat(refs).to_csv(outdir / "references.tsv", sep="\t", index=False)

    # 2. Proteomes (ingroup) and Smc proteins (outgroup)
    outgroups = []
    for sp in species.itertuples():
        if sp.role == "ingroup":
            path = outdir / "proteomes" / f"{sp.code}.faa.gz"
            if path.exists():
                print(f"{sp.code:5s} {sp.species:32s} already downloaded")
                continue
            n, strain = fetch_proteome(sp.taxon_id, path)
            if n == 0:
                print(f"WARNING: no proteome for {sp.species} ({sp.taxon_id})")
            else:
                print(f"{sp.code:5s} {sp.species:32s} {n:6d} proteins "
                      f"(strain taxon {strain})")
        else:
            df = uniprot_tsv(f"gene_exact:smc AND organism_id:{sp.taxon_id} "
                             "AND reviewed:true")
            if df.empty:
                df = uniprot_tsv(f"gene_exact:smc AND organism_id:{sp.taxon_id}")
            if df.empty:
                print(f"WARNING: no Smc for outgroup {sp.species}")
                continue
            df["code"] = sp.code
            outgroups.append(df)
            print(f"{sp.code:5s} {sp.species:32s} {len(df):6d} Smc (outgroup)")
        time.sleep(1)
    pd.concat(outgroups).to_csv(outdir / "outgroups.tsv", sep="\t", index=False)
    print(f"\nDone. Data in {outdir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
