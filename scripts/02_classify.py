#!/usr/bin/env python3
"""Find SMC proteins in each proteome by sequence similarity and label them.

1. Search: every protein of every proteome is compared with the human and
   yeast reference SMC1-6 and RAD50 proteins (DIAMOND blastp). A protein is
   a candidate when it aligns to the N-terminal ATPase head of a reference
   (alignment starting within the first HEAD_WINDOW residues). This keeps
   true SMC-like proteins and drops proteins that only share coiled coils.
2. Label: each candidate gets the family of its best hit. It is
   'ambiguous' when the best and second-best families score too similarly,
   and 'fragment' when shorter than min_length.
3. Redundancy: within a species, a sequence >= 95% identical over >= 90% of
   its length to a longer one is 'redundant' (isoforms, duplicate records).

Outputs (in OUTDIR):
  classification.tsv  one row per candidate with its label and scores
  tree_input.faa      SMC1-6 sequences plus prokaryotic outgroups

Usage: python scripts/02_classify.py data/raw config/species.tsv results [min_length] [threads]
"""
import gzip
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

MIN_RATIO = 1.15     # best family bitscore / second-best family bitscore
HEAD_WINDOW = 200    # hit must start within this many residues of the reference
TREE_FAMILIES = ["SMC1", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6"]
COLS = ["query", "subject", "pident", "aln_len", "sstart", "evalue", "bitscore"]


def write_fasta(records, path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def parse_header(header):
    """UniProt FASTA header -> accession, entry name, gene name."""
    first = header.split()[0]
    parts = first.split("|")
    acc, entry = (parts[1], parts[2]) if len(parts) >= 3 else (first, first)
    gene = re.search(r"GN=(\S+)", header)
    return acc, entry, gene.group(1) if gene else ""


def run_diamond(query, refs, out, threads):
    db = out.with_suffix("")
    subprocess.run(["diamond", "makedb", "--quiet", "--in", str(refs),
                    "-d", str(db)], check=True)
    subprocess.run(["diamond", "blastp", "--quiet", "--more-sensitive",
                    "-q", str(query), "-d", str(db), "-o", str(out),
                    "-k", "50", "-e", "1e-5", "-p", str(threads),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "sstart", "evalue", "bitscore"], check=True)
    return pd.read_csv(out, sep="\t", names=COLS)


def main(rawdir, species_file, outdir, min_length=700, threads=4):
    rawdir, outdir = Path(rawdir), Path(outdir)
    work = outdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(species_file, sep="\t", dtype=str)
    organism = dict(zip(species["code"], species["species"]))

    refs = pd.read_csv(rawdir / "references.tsv", sep="\t")
    write_fasta(zip(refs["family"] + "|" + refs["Entry"], refs["Sequence"]),
                work / "references.faa")

    # 1. Search all proteomes
    proteins = []
    for code in species.loc[species["role"] == "ingroup", "code"]:
        path = rawdir / "proteomes" / f"{code}.faa.gz"
        if not path.exists():
            print(f"WARNING: no proteome for {code}")
            continue
        for header, seq in read_fasta(path):
            acc, entry, gene = parse_header(header)
            proteins.append({"name": f"{code}_{acc}", "code": code,
                             "Entry": acc, "Entry Name": entry,
                             "Gene Names (primary)": gene,
                             "Organism": organism[code], "Length": len(seq),
                             "Sequence": seq})
    proteins = pd.DataFrame(proteins)
    write_fasta(zip(proteins["name"], proteins["Sequence"]), work / "proteomes.faa")
    print(f"Searching {len(proteins)} proteins from "
          f"{proteins['code'].nunique()} proteomes")
    hits = run_diamond(work / "proteomes.faa", work / "references.faa",
                       work / "diamond_hits.tsv", threads)
    hits["family"] = hits["subject"].str.split("|").str[0]
    head_hits = set(hits.loc[hits["sstart"] <= HEAD_WINDOW, "query"])
    hits = hits[hits["query"].isin(head_hits)]

    # 2. Label by best family
    best = (hits.groupby(["query", "family"])["bitscore"].max()
                .reset_index().sort_values("bitscore", ascending=False))
    rows = []
    for q, grp in best.groupby("query", sort=False):
        top = grp.iloc[0]
        second = grp.iloc[1]["bitscore"] if len(grp) > 1 else 0.0
        rows.append({"name": q, "best_family": top["family"],
                     "best_bitscore": top["bitscore"],
                     "second_bitscore": second,
                     "ratio": round(top["bitscore"] / second, 2)
                              if second > 0 else float("inf")})
    table = proteins.merge(pd.DataFrame(rows), on="name", how="inner")

    def label(r):
        if r["Length"] < min_length:
            return "fragment"
        if r["ratio"] < MIN_RATIO:
            return "ambiguous"
        return r["best_family"]

    table["label"] = table.apply(label, axis=1)

    # 3. Redundant entries within a species
    write_fasta(zip(table["name"], table["Sequence"]), work / "candidates.faa")
    self_hits = run_diamond(work / "candidates.faa", work / "candidates.faa",
                            work / "self_hits.tsv", threads)
    lengths = dict(zip(table["name"], table["Length"]))
    self_hits = self_hits[(self_hits["query"] != self_hits["subject"]) &
                          (self_hits["pident"] >= 95) &
                          (self_hits["query"].str.split("_").str[0] ==
                           self_hits["subject"].str.split("_").str[0])]
    redundant = set()
    for r in self_hits.sort_values("bitscore", ascending=False).itertuples():
        lq, ls = lengths[r.query], lengths[r.subject]
        shorter, longer = ((r.query, r.subject) if (lq, r.query) < (ls, r.subject)
                           else (r.subject, r.query))
        if r.aln_len >= 0.9 * min(lq, ls) and longer not in redundant:
            redundant.add(shorter)
    table.loc[table["name"].isin(redundant), "label"] = "redundant"

    keep_cols = ["code", "Entry", "Entry Name", "Gene Names (primary)",
                 "Organism", "Length", "best_family", "best_bitscore",
                 "second_bitscore", "ratio", "label"]
    table[keep_cols].sort_values(["code", "label", "best_bitscore"],
                                 ascending=[True, True, False]).to_csv(
        outdir / "classification.tsv", sep="\t", index=False)

    # Sequences for the tree: SMC1-6 from all species + outgroups
    tree = table[table["label"].isin(TREE_FAMILIES)].copy()
    tree["tip"] = tree["label"] + "_" + tree["code"] + "_" + tree["Entry"]
    out = pd.read_csv(rawdir / "outgroups.tsv", sep="\t")
    out["tip"] = "prokSMC_" + out["code"] + "_" + out["Entry"]
    records = list(zip(tree["tip"], tree["Sequence"])) + \
              list(zip(out["tip"], out["Sequence"]))
    write_fasta(records, outdir / "tree_input.faa")

    print(table["label"].value_counts().to_string())
    print(f"\n{len(records)} sequences written to {outdir / 'tree_input.faa'}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3],
         int(sys.argv[4]) if len(sys.argv) > 4 else 700,
         int(sys.argv[5]) if len(sys.argv) > 5 else 4)
