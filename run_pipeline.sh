#!/usr/bin/env bash
# Run the full SMC phylogeny pipeline.
#
#   ./run_pipeline.sh          real data from UniProt  -> results/
#   ./run_pipeline.sh --test   simulated data          -> results/test/
#
# Set THREADS to change the number of CPU cores (default 4).
set -euo pipefail

THREADS="${THREADS:-4}"
MIN_LENGTH=700
SPECIES=config/species.tsv

if [[ "${1:-}" == "--test" ]]; then
    MODE=test
    DATA=data/test
    OUT=results/test
    MODEL="LG+G4"          # fixed model keeps the test quick
else
    MODE=real
    DATA=data/raw
    OUT=results
    MODEL="MFP"            # let ModelFinder choose the best model
fi
mkdir -p "$DATA" "$OUT"

echo "== 1. Get sequences ($MODE)"
if [[ "$MODE" == test ]]; then
    python scripts/simulate_test_data.py "$SPECIES" "$DATA"
elif [[ -s "$DATA/references.tsv" && -s "$DATA/outgroups.tsv" && -d "$DATA/proteomes" ]]; then
    echo "Using sequences already in $DATA (delete that folder to download again)"
else
    python scripts/01_fetch_uniprot.py "$SPECIES" "$DATA"
fi

echo "== 2. Find SMC proteins by similarity search and label them"
python scripts/02_classify.py "$DATA" "$SPECIES" "$OUT" "$MIN_LENGTH" "$THREADS"

echo "== 3. Align (MAFFT) and trim (trimAl)"
mafft --auto --thread "$THREADS" --quiet "$OUT/tree_input.faa" > "$OUT/smc.aln.faa"
trimal -in "$OUT/smc.aln.faa" -out "$OUT/smc.trim.faa" -automated1

echo "== 4. Maximum-likelihood tree (IQ-TREE, 1000 ultrafast bootstraps)"
iqtree2 -s "$OUT/smc.trim.faa" -m "$MODEL" -B 1000 -T AUTO --threads-max "$THREADS" \
        --prefix "$OUT/work/smc" -redo -quiet
cp "$OUT/work/smc.treefile" "$OUT/smc.treefile"
grep -m1 "Best-fit model" "$OUT/work/smc.iqtree" || grep -m1 "Model of substitution" "$OUT/work/smc.iqtree"

echo "== 5. Figures and tables"
python scripts/03_plot_tree.py "$OUT/smc.treefile" "$SPECIES" "$OUT"
python scripts/04_presence_absence.py "$OUT/classification.tsv" "$SPECIES" "$OUT"

if [[ "$MODE" == test ]]; then
    echo "== 6. Check results against the simulated truth"
    python scripts/check_test_results.py "$DATA" "$OUT"
fi
echo "Done. Results are in $OUT/"
