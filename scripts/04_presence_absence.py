#!/usr/bin/env python3
"""Count SMC1-6 copies per species and draw a presence/absence heatmap.

Usage: python scripts/04_presence_absence.py results/classification.tsv config/species.tsv results
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

FAMILIES = ["SMC1", "SMC3", "SMC2", "SMC4", "SMC5", "SMC6"]
GROUP_LABELS = ["cohesin", "cohesin", "condensin", "condensin",
                "SMC5/6", "SMC5/6"]


def main(classification, species_file, outdir):
    outdir = Path(outdir)
    calls = pd.read_csv(classification, sep="\t")
    species = pd.read_csv(species_file, sep="\t", dtype=str)
    species = species[species["role"] == "ingroup"]

    counts = (calls[calls["label"].isin(FAMILIES)]
              .groupby(["code", "label"]).size().unstack(fill_value=0)
              .reindex(index=species["code"], columns=FAMILIES, fill_value=0))
    counts.index = [f"{s} ({g})" for s, g in
                    zip(species["species"], species["supergroup"])]
    counts.to_csv(outdir / "copy_numbers.tsv", sep="\t")
    print(counts.to_string())

    # Cells with no full-length copy, but a fragment or ambiguous hit:
    # these may be annotation problems rather than real losses.
    unsure = calls[calls["label"].isin(["fragment", "ambiguous"])]
    unsure = set(zip(unsure["code"], unsure["best_family"]))
    codes = list(species["code"])

    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(counts) + 1.8))
    shown = counts.clip(upper=3)
    ax.imshow(shown.values, cmap="Blues", vmin=-0.4, vmax=3, aspect="auto")
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            n = counts.iat[i, j]
            text = n
            if n == 0:
                text = ("fragment\nonly" if (codes[i], FAMILIES[j]) in unsure
                        else "absent")
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=8 if n else 6.5,
                    color="#b03a2e" if n == 0 else
                          ("white" if n >= 2 else "black"))
    ax.set_xticks(range(len(FAMILIES)))
    ax.set_xticklabels([f"{f}\n{g}" for f, g in zip(FAMILIES, GROUP_LABELS)],
                       fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Copies of each SMC family per species", loc="left", pad=34)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"copy_numbers.{ext}", dpi=200)
    print(f"Heatmap written to {outdir / 'copy_numbers.png'}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
