#!/usr/bin/env python3
"""Root the SMC tree on the prokaryotic outgroup, test whether each SMC
family is monophyletic, and draw the tree coloured by family.

Usage: python scripts/03_plot_tree.py results/smc.treefile config/species.tsv results
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

COLOURS = {
    "SMC1": "#1b6ca8", "SMC3": "#6fb1e4",   # cohesin
    "SMC2": "#c0392b", "SMC4": "#ef8a7c",   # condensin
    "SMC5": "#2e8b57", "SMC6": "#8fd19e",   # SMC5/6
    "prokSMC": "#7f7f7f",
}
COMPLEX = {"SMC1": "cohesin", "SMC3": "cohesin", "SMC2": "condensin",
           "SMC4": "condensin", "SMC5": "SMC5/6", "SMC6": "SMC5/6"}


def family(tip_name):
    return tip_name.split("_")[0]


def main(treefile, species_file, outdir):
    outdir = Path(outdir)
    species = pd.read_csv(species_file, sep="\t", dtype=str).set_index("code")
    tree = Phylo.read(treefile, "newick")
    tips = tree.get_terminals()

    # Root on the prokaryotic SMCs
    outgroup = [t for t in tips if family(t.name) == "prokSMC"]
    if outgroup:
        tree.root_with_outgroup(tree.common_ancestor(outgroup))
    tree.ladderize()

    # Monophyly test. For each family, find its "core clade": the largest
    # clade made only of that family. The family is monophyletic when the
    # core clade holds all its members. Members outside it are listed: they
    # are either misclassified or genuinely odd sequences worth a closer look.
    rows = []
    core_clades = {}
    for fam in sorted({family(t.name) for t in tips} - {"prokSMC"}):
        members = [t.name for t in tips if family(t.name) == fam]
        core = None
        for c in tree.find_clades():
            names = [t.name for t in c.get_terminals()]
            if all(family(n) == fam for n in names):
                if core is None or len(names) > len(core.get_terminals()):
                    core = c
        core_clades[fam] = core
        core_names = {t.name for t in core.get_terminals()}
        outside = [m for m in members if m not in core_names]
        # nearest family of each outside member = family of its sister tips
        placed = []
        for m in outside:
            path = tree.get_path(m)
            parent = path[-2] if len(path) > 1 else tree.root
            sisters = [family(t.name) for t in parent.get_terminals()
                       if t.name != m]
            near = max(set(sisters), key=sisters.count) if sisters else "?"
            placed.append(f"{m} (sits with {near})")
        rows.append({
            "family": fam, "complex": COMPLEX.get(fam, ""),
            "n_sequences": len(members),
            "n_species": len({m.split("_")[1] for m in members}),
            "in_core_clade": len(core_names),
            "monophyletic": not outside,
            "ufboot_support": core.confidence if len(core_names) > 1 else None,
            "outside_core": ";".join(placed),
        })
    mono = pd.DataFrame(rows)
    mono.to_csv(outdir / "family_monophyly.tsv", sep="\t", index=False)
    print(mono.drop(columns="outside_core").to_string(index=False))
    for r in mono[mono["outside_core"] != ""].itertuples():
        print(f"  {r.family} outside core clade: {r.outside_core}")

    # Draw
    def label(clade):
        if not clade.is_terminal():
            return None
        fam, code, acc = clade.name.split("_", 2)
        name = species.loc[code, "species"] if code in species.index else code
        return f"{name} ({acc})"

    label_colour = {label(t): COLOURS.get(family(t.name), "black") for t in tips}
    # Show support values only at the base of each family clade
    family_nodes = {id(c) for c in core_clades.values()}
    for c in tree.find_clades():
        if not c.is_terminal() and id(c) not in family_nodes:
            c.confidence = None

    height = max(8, 0.17 * len(tips))
    fig, ax = plt.subplots(figsize=(11, height))
    plt.rcParams["font.size"] = 7
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label,
               label_colors=lambda lab: label_colour.get(lab, "black"),
               branch_labels=lambda c: (f"{c.confidence:.0f}"
                                        if c.confidence is not None else None))
    handles = [plt.Line2D([], [], color=c, lw=6,
                          label=f"{f} ({COMPLEX[f]})" if f in COMPLEX else
                          "prokaryotic SMC (outgroup)")
               for f, c in COLOURS.items()]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False)
    ax.set_title("Maximum-likelihood tree of SMC proteins across eukaryotes",
                 loc="left")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"smc_tree.{ext}", dpi=200)
    print(f"Tree figure written to {outdir / 'smc_tree.png'}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
