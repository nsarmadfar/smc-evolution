#!/usr/bin/env python3
"""Compare pipeline output on simulated data with the known truth.

Usage: python scripts/check_test_results.py data/test results/test
Exits with an error if any check fails.
"""
import sys
from pathlib import Path

import pandas as pd


def main(datadir, resultdir):
    truth = pd.read_csv(Path(datadir) / "truth.tsv", sep="\t")
    calls = pd.read_csv(Path(resultdir) / "classification.tsv", sep="\t")
    mono = pd.read_csv(Path(resultdir) / "family_monophyly.tsv", sep="\t")
    counts = pd.read_csv(Path(resultdir) / "copy_numbers.tsv", sep="\t", index_col=0)

    m = truth[~truth["true_family"].isin(["other", "prokSMC"])].merge(
        calls, on=["Entry", "code"], how="left")
    found_other = truth[truth["true_family"] == "other"].merge(
        calls, on=["Entry", "code"])
    fragment = m["expected"] == "fragment"
    wrong = m[m["expected"] != m["label"]]

    checks = {
        "every SMC found and labelled correctly": wrong.empty,
        "no unrelated or coiled-coil-only proteins picked up": found_other.empty,
        "fragment removed (Pfal SMC4)": (m.loc[fragment, "label"] == "fragment").all(),
        "truncated gene model found (Mbre SMC4)":
            (m.loc[(m["code"] == "Mbre") & (m["true_family"] == "SMC4"),
                   "label"] == "SMC4").all(),
        "isoforms removed (Hsap SMC5 = 1)":
            int(counts.filter(like="Homo", axis=0)["SMC5"].iloc[0]) == 1,
        "all six SMC families monophyletic": bool(mono["monophyletic"].all())
                                              and len(mono) == 6,
        "loss detected (Gint SMC5, SMC6 = 0)":
            counts.filter(like="Giardia", axis=0)[["SMC5", "SMC6"]].eq(0).all().all(),
        "duplication detected (Tthe SMC1 = 2)":
            int(counts.filter(like="Tetrahymena", axis=0)["SMC1"].iloc[0]) == 2,
    }
    for name, ok in checks.items():
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not wrong.empty:
        print(wrong[["Entry", "code", "true_family", "expected", "label", "ratio"]])
    if not all(checks.values()):
        sys.exit(1)
    print("\nAll checks passed.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
