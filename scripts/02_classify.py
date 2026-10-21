#!/usr/bin/env python3
"""Label each candidate protein as SMC1-6, RAD50 or 'other'.

Each candidate is compared with the human and yeast reference proteins
(DIAMOND blastp). It gets the family of its best hit. A label is called
'ambiguous' when the best and second-best families score too similarly,
and short sequences (likely fragments) are removed.

Outputs (in OUTDIR):
  classification.tsv  one row per candidate with its label and scores
  tree_input.faa      SMC1-6 sequences plus prokaryotic outgroups, for the tree

Usage: python scripts/02_classify.py data/raw results [min_length] [threads]
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

MIN_RATIO = 1.15  # best family bitscore / second family bitscore
TREE_FAMILIES = ["SMC1", "SMC2", "SMC3", "SMC4", "SMC5", "SMC6"]


def write_fasta(records, path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def run_diamond(query, refs, out, threads):
    db = out.with_suffix("")
    subprocess.run(["diamond", "makedb", "--quiet", "--in", str(refs),
                    "-d", str(db)], check=True)
    subprocess.run(["diamond", "blastp", "--quiet", "--more-sensitive",
                    "-q", str(query), "-d", str(db), "-o", str(out),
                    "-k", "50", "-e", "1e-5", "-p", str(threads),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "evalue", "bitscore"], check=True)
    cols = ["query", "subject", "pident", "aln_len", "evalue", "bitscore"]
    return pd.read_csv(out, sep="\t", names=cols)


def main(rawdir, outdir, min_length=700, threads=4):
    rawdir, outdir = Path(rawdir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "work"
    work.mkdir(exist_ok=True)

    refs = pd.read_csv(rawdir / "references.tsv", sep="\t")
    cand = pd.read_csv(rawdir / "candidates.tsv", sep="\t")
    cand["name"] = cand["code"] + "_" + cand["Entry"]

    write_fasta(zip(refs["family"] + "|" + refs["Entry"], refs["Sequence"]),
                work / "references.faa")
    ingroup = cand[cand["role"] == "ingroup"]
    write_fasta(zip(ingroup["name"], ingroup["Sequence"]), work / "candidates.faa")

    hits = run_diamond(work / "candidates.faa", work / "references.faa",
                       work / "diamond_hits.tsv", threads)
    hits["family"] = hits["subject"].str.split("|").str[0]

    # best bitscore per candidate and family
    best = (hits.groupby(["query", "family"])["bitscore"].max()
                .reset_index().sort_values("bitscore", ascending=False))
    rows = []
    for q, grp in best.groupby("query", sort=False):
        top = grp.iloc[0]
        second = grp.iloc[1]["bitscore"] if len(grp) > 1 else 0.0
        ratio = top["bitscore"] / second if second > 0 else float("inf")
        rows.append({"name": q, "best_family": top["family"],
                     "best_bitscore": top["bitscore"],
                     "second_bitscore": second, "ratio": round(ratio, 2)})
    calls = pd.DataFrame(rows)

    table = ingroup.merge(calls, on="name", how="left")
    table["best_family"] = table["best_family"].fillna("other")

    def label(r):
        if r["best_family"] == "other":
            return "other"
        if r["Length"] < min_length:
            return "fragment"
        if r["ratio"] < MIN_RATIO:
            return "ambiguous"
        return r["best_family"]

    table["label"] = table.apply(label, axis=1)
    keep_cols = ["code", "Entry", "Entry Name", "Gene Names (primary)",
                 "Organism", "Length", "best_family", "best_bitscore",
                 "second_bitscore", "ratio", "label"]
    table[keep_cols].sort_values(["code", "label"]).to_csv(
        outdir / "classification.tsv", sep="\t", index=False)

    # Sequences for the tree: SMC1-6 from all ingroup species + outgroups
    tree = table[table["label"].isin(TREE_FAMILIES)].copy()
    tree["tip"] = tree["label"] + "_" + tree["code"] + "_" + tree["Entry"]
    out = cand[cand["role"] == "outgroup"].copy()
    out["tip"] = "prokSMC_" + out["code"] + "_" + out["Entry"]
    records = list(zip(tree["tip"], tree["Sequence"])) + \
              list(zip(out["tip"], out["Sequence"]))
    write_fasta(records, outdir / "tree_input.faa")

    print(table["label"].value_counts().to_string())
    print(f"\n{len(records)} sequences written to {outdir / 'tree_input.faa'}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 700,
         int(sys.argv[4]) if len(sys.argv) > 4 else 4)
