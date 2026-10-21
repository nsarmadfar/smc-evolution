# smc-evolution

Evolution of SMC proteins, the motors that fold chromosomes, across the eukaryotic tree of life.

## Why SMC proteins

SMC (Structural Maintenance of Chromosomes) proteins form the core of three complexes that shape how chromosomes are folded in 3D:

| Complex | SMC pair | Main role |
| --- | --- | --- |
| Cohesin | SMC1 + SMC3 | holds sister chromatids together, forms loops and domains in interphase |
| Condensin | SMC2 + SMC4 | compacts chromosomes for cell division |
| SMC5/6 | SMC5 + SMC6 | DNA repair and genome stability |

Bacteria and archaea have a single SMC protein (or a distant relative such as MukB). Eukaryotes have six. Earlier work (for example Cobbe and Heck, 2004, *Mol. Biol. Evol.* 21:332) suggests that these six arose by gene duplications before the last eukaryotic common ancestor. Most of this evidence comes from model organisms such as yeast, animals and plants. I wanted to test it myself with a set of species that covers all major eukaryotic supergroups, including understudied unicellular eukaryotes, and to see whether any lineage has lost or duplicated an SMC family.

## Questions

1. Does each SMC family (SMC1 to SMC6) form one clade that contains all eukaryotic supergroups? This is what we expect if all six families were already present in the last eukaryotic common ancestor.
2. Which species have lost an SMC family, and which have extra copies?

## Approach

```
UniProt reference proteomes (14 eukaryotes, 6 supergroups)
   │  proteins with the SMC_N domain (Pfam PF02463)
   ▼
Classification: DIAMOND blastp against human and yeast SMC1-6 and RAD50
   │  best-hit family, ambiguous and fragment (< 700 aa) calls removed
   ▼
Alignment (MAFFT) → trimming (trimAl -automated1)
   ▼
Maximum-likelihood tree (IQ-TREE 2, ModelFinder, 1000 ultrafast bootstraps)
   │  rooted on prokaryotic SMC from 2 bacteria and 2 archaea
   ▼
Monophyly test per family + copy-number heatmap
```

The species are listed in [`config/species.tsv`](config/species.tsv). They include opisthokonts (human, budding yeast, the choanoflagellate *Monosiga*), amoebozoa (*Dictyostelium*), plants and algae (*Arabidopsis*, *Chlamydomonas*), stramenopiles (*Phytophthora*, the diatom *Thalassiosira*), alveolates (*Plasmodium*, *Tetrahymena*), discobans (*Trypanosoma*, *Naegleria*) and metamonads (*Giardia*, *Trichomonas*).

### Testing on simulated data first

Before using real data, I test the pipeline on simulated sequences where the true answer is known (`scripts/simulate_test_data.py`). An ancestral sequence evolves along a known history: six families arise by duplication before the species split, then each evolves along the same species tree. I added a simulated gene loss (SMC5 and SMC6 in *Giardia*), a duplication (SMC1 in *Tetrahymena*), a fragment, and unrelated decoy proteins. The pipeline passes when it labels every protein correctly, recovers all six families as clades, and finds the loss and the duplication.

## How to run

```bash
conda env create -f environment.yml
conda activate smc-evolution

./run_pipeline.sh --test     # simulated data, checks the results (about 15 min on 2 cores)
./run_pipeline.sh            # real data from UniProt
```

Use `THREADS=8 ./run_pipeline.sh` to use more CPU cores.

## Output

| File | Content |
| --- | --- |
| `results/classification.tsv` | every candidate protein with its label and scores |
| `results/smc_tree.png` | the rooted tree, coloured by family, support values at the base of each family |
| `results/family_monophyly.tsv` | whether each family forms one clade, and its bootstrap support |
| `results/copy_numbers.png` | copies of SMC1-6 and RAD50 in each species |

## Results

<!-- Fill in after running ./run_pipeline.sh on the real data. -->

## Limitations

- Classification by best BLAST hit against human and yeast proteins can fail for very divergent sequences. The tree is the real test: a protein that is labelled wrongly shows up inside another family's clade.
- An absent family can mean a real gene loss, or a gene that is missing from the proteome annotation. A next step would be to search the genome assemblies directly (for example with HMMER or tblastn) before calling a loss.
- One species per lineage gives a broad but thin sample.

## Repository layout

```
config/species.tsv              species, taxon IDs and supergroups
run_pipeline.sh                 runs all steps
scripts/01_fetch_uniprot.py     download candidates, references and outgroups
scripts/02_classify.py          label candidates with DIAMOND
scripts/03_plot_tree.py         root the tree, test monophyly, draw the tree
scripts/04_presence_absence.py  copy-number table and heatmap
scripts/simulate_test_data.py   simulated data with a known answer
scripts/check_test_results.py   compares test output with the truth
```

## Author

Nargess Sarmadfar ([github.com/nsarmadfar](https://github.com/nsarmadfar))
