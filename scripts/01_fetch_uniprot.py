#!/usr/bin/env python3
"""Download SMC candidate proteins and reference sequences from UniProt.

Candidates (ingroup): every protein in the species' UniProt reference proteome
that carries the SMC N-terminal domain (Pfam PF02463, SMC_N). This domain is
shared by SMC1-6, RAD50 and some other ATPases, so the next step classifies them.

Outgroup: the reviewed bacterial or archaeal Smc protein of each outgroup species.

References: reviewed human and budding yeast SMC1-6 and RAD50, used only to
label the candidates.

Usage: python scripts/01_fetch_uniprot.py config/species.tsv data/raw
"""
import io
import sys
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,id,gene_primary,organism_id,organism_name,length,sequence"

REFERENCE_QUERIES = {
    # family -> (organism_id, gene names)
    "human": ("9606", ["SMC1A", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6", "RAD50"]),
    "yeast": ("559292", ["SMC1", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6", "RAD50"]),
}


def uniprot_tsv(query, retries=4):
    """Run a UniProt query and return the result as a data frame."""
    params = {"query": query, "format": "tsv", "fields": FIELDS}
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(API, params=params, timeout=120)
            r.raise_for_status()
            if not r.text.strip():
                return pd.DataFrame()
            return pd.read_csv(io.StringIO(r.text), sep="\t")
        except requests.RequestException as err:
            print(f"  attempt {attempt} failed: {err}", file=sys.stderr)
            time.sleep(5 * attempt)
    raise RuntimeError(f"UniProt query failed: {query}")


def family_of(gene):
    """Map a reference gene name (SMC1A, SMC1, RAD50 ...) to a family label."""
    gene = str(gene).upper()
    if gene.startswith("RAD50"):
        return "RAD50"
    if gene.startswith("SMC") and len(gene) >= 4 and gene[3].isdigit():
        return "SMC" + gene[3]
    return None


def main(species_file, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(species_file, sep="\t", dtype=str)

    # 1. References
    refs = []
    for label, (taxon, genes) in REFERENCE_QUERIES.items():
        gene_q = " OR ".join(f"gene_exact:{g}" for g in genes)
        df = uniprot_tsv(f"({gene_q}) AND organism_id:{taxon} AND reviewed:true")
        df["family"] = df["Gene Names (primary)"].map(family_of)
        refs.append(df.dropna(subset=["family"]))
        print(f"references ({label}): {len(df)} sequences")
    refs = pd.concat(refs)
    refs.to_csv(outdir / "references.tsv", sep="\t", index=False)

    # 2. Candidates and outgroups
    rows = []
    for sp in species.itertuples():
        if sp.role == "ingroup":
            q = f"xref:pfam-PF02463 AND taxonomy_id:{sp.taxon_id} AND keyword:KW-1185"
        else:
            q = f"gene_exact:smc AND organism_id:{sp.taxon_id} AND reviewed:true"
        df = uniprot_tsv(q)
        if df.empty and sp.role == "outgroup":
            # fall back to unreviewed entries
            df = uniprot_tsv(f"gene_exact:smc AND organism_id:{sp.taxon_id}")
        if df.empty:
            print(f"WARNING: no sequences for {sp.species} ({sp.taxon_id})")
            continue
        # A species can have several strains with a reference proteome.
        # Keep only the strain with the most hits so each species counts once.
        top_strain = df["Organism (ID)"].value_counts().idxmax()
        df = df[df["Organism (ID)"] == top_strain].copy()
        df["code"] = sp.code
        df["role"] = sp.role
        rows.append(df)
        print(f"{sp.code:5s} {sp.species:32s} {len(df):3d} sequences "
              f"(strain taxon {top_strain})")
        time.sleep(0.5)  # be polite to the API

    cand = pd.concat(rows)
    cand.to_csv(outdir / "candidates.tsv", sep="\t", index=False)
    print(f"\nWrote {len(refs)} references and {len(cand)} candidates to {outdir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
