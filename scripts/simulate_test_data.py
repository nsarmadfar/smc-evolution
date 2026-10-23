#!/usr/bin/env python3
"""Simulate a test data set where the true answer is known.

An ancestral SMC sequence evolves along a known history:
  prokaryotic SMC  vs  eukaryotic stem
  eukaryotic stem -> 3 complexes -> 6 families (SMC1-6), all before species split
  each family then evolves along the same species tree.
RAD50 and unrelated proteins are added as decoys. To test the pipeline, one
species loses SMC5 and SMC6, one species gains an extra SMC1 copy, and one
SMC4 is a short fragment.

Each species also gets 200 unrelated proteins, and one species gets a
"coiled-coil" decoy that shares only the middle of SMC1 (no ATPase head),
which the search should ignore.

The output has the same format as 01_fetch_uniprot.py (proteomes, references,
outgroups), so the rest of the pipeline runs unchanged. A file with the true
labels is written for checking.

Usage: python scripts/simulate_test_data.py config/species.tsv data/test [seed]
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AA = list("ACDEFGHIKLMNPQRSTVWY")
FREQ = np.array([7.4, 2.5, 5.4, 5.4, 4.7, 7.4, 2.6, 6.8, 5.8, 9.9,
                 2.5, 4.5, 3.9, 3.4, 5.2, 5.7, 5.1, 7.3, 1.3, 3.2])
FREQ = FREQ / FREQ.sum()
LENGTH = 1100

rng = None
site_rates = None


def random_seq(n):
    return rng.choice(len(AA), size=n, p=FREQ)


def evolve(seq, branch_length):
    """Substitute each site with probability 1 - exp(-rate * branch length)."""
    p = 1 - np.exp(-site_rates[: len(seq)] * branch_length)
    hit = rng.random(len(seq)) < p
    new = seq.copy()
    new[hit] = rng.choice(len(AA), size=hit.sum(), p=FREQ)
    return new


def to_str(seq):
    return "".join(AA[i] for i in seq)


# Species tree as nested tuples of (subtree, branch length); leaves are codes.
SPECIES_TREE = [
    ([([("Hsap", .15), ("Scer", .25)], .05), ("Mbre", .2)], .1),
    ("Ddis", .3),
    ([("Atha", .15), ("Crei", .2)], .1),
    ([([("Pinf", .2), ("Tpse", .2)], .1), ([("Pfal", .35), ("Tthe", .3)], .1)], .05),
    ([("Tbru", .3), ("Ngru", .3)], .1),
    ([("Gint", .45), ("Tvag", .4)], .1),
]


def along_species_tree(seq, subtree, out):
    for child, bl in subtree:
        child_seq = evolve(seq, bl)
        if isinstance(child, str):
            out[child] = child_seq
        else:
            along_species_tree(child_seq, child, out)


def main(species_file, outdir, seed=1):
    global rng, site_rates
    rng = np.random.default_rng(seed)
    site_rates = rng.gamma(shape=0.8, scale=1 / 0.8, size=LENGTH)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(species_file, sep="\t", dtype=str).set_index("code")

    root = random_seq(LENGTH)
    prok = evolve(root, 0.4)
    euk_stem = evolve(root, 0.4)
    families = {}
    for complex_, pair in [("cohesin", ("SMC1", "SMC3")),
                           ("condensin", ("SMC2", "SMC4")),
                           ("smc56", ("SMC5", "SMC6"))]:
        c = evolve(euk_stem, 0.3)
        for fam in pair:
            families[fam] = evolve(c, 0.3)
    families["RAD50"] = evolve(random_seq(LENGTH), 0.1)   # unrelated decoy family

    proteomes = {c: [] for c in species.index[species["role"] == "ingroup"]}
    outgroups, truth, refs = [], [], []
    n = 0

    def add(code, seq, fam, role="ingroup", gene=""):
        nonlocal n
        n += 1
        acc = f"SIM{n:05d}"
        sp = species.loc[code]
        if role == "ingroup":
            header = (f"tr|{acc}|{acc}_{code.upper()} simulated protein "
                      f"OS={sp.species} OX={sp.taxon_id}"
                      + (f" GN={gene}" if gene else ""))
            proteomes[code].append((header, to_str(seq)))
        else:
            outgroups.append({"Entry": acc, "Entry Name": f"{acc}_{code.upper()}",
                              "Gene Names (primary)": gene,
                              "Organism (ID)": sp.taxon_id, "Organism": sp.species,
                              "Length": len(seq), "Sequence": to_str(seq),
                              "code": code})
        truth.append({"Entry": acc, "code": code, "true_family": fam})
        return acc

    for fam, anc in families.items():
        leaves = {}
        along_species_tree(anc, SPECIES_TREE, leaves)
        for code, seq in leaves.items():
            if code == "Gint" and fam in ("SMC5", "SMC6"):
                continue                                   # simulated loss
            if code == "Pfal" and fam == "SMC4":
                seq = seq[:400]                            # simulated fragment
            add(code, seq, fam)
            if code == "Tthe" and fam == "SMC1":
                add(code, evolve(seq, 0.15), fam)          # simulated duplication
            if code in ("Hsap", "Scer"):
                gene = fam + ("A" if fam == "SMC1" and code == "Hsap" else "")
                refs.append({"Entry": f"REF_{code}_{fam}", "Entry Name": "",
                             "Gene Names (primary)": gene,
                             "Organism (ID)": species.loc[code].taxon_id,
                             "Organism": species.loc[code].species,
                             "Length": len(seq), "Sequence": to_str(seq),
                             "family": fam})
    # unrelated background proteins
    for code in proteomes:
        for _ in range(200):
            add(code, random_seq(int(rng.integers(150, 1500))), "other")
    # coiled-coil decoy: middle of SMC1 only, no ATPase head
    middle = evolve(families["SMC1"], 0.3)[400:800]
    add("Ddis", np.concatenate([random_seq(500), middle, random_seq(300)]), "other")
    # prokaryotic outgroups
    for code in species.index[species["role"] == "outgroup"]:
        add(code, evolve(prok, 0.3), "prokSMC", role="outgroup", gene="smc")

    (outdir / "proteomes").mkdir(exist_ok=True)
    for code, recs in proteomes.items():
        with gzip.open(outdir / "proteomes" / f"{code}.faa.gz", "wt") as fh:
            for header, seq in recs:
                fh.write(f">{header}\n{seq}\n")
    pd.DataFrame(refs).to_csv(outdir / "references.tsv", sep="\t", index=False)
    pd.DataFrame(outgroups).to_csv(outdir / "outgroups.tsv", sep="\t", index=False)
    pd.DataFrame(truth).to_csv(outdir / "truth.tsv", sep="\t", index=False)
    print(f"Simulated {n} proteins in {len(proteomes)} proteomes, "
          f"{len(refs)} references and {len(outgroups)} outgroups in {outdir}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
