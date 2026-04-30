# Tripeptide Angle Analysis

Analysis of signed angles between consecutive Cα→sidechain-centroid vectors for XRX tripeptides in helical regions. Angles are grouped by the size of the left neighbor residue.

## Requirements

- Python 3 with biopython, tqdm, pandas, numpy, scipy, matplotlib
- STRIDE (must be in PATH)
- Snakemake

```
pip install biopython tqdm pandas numpy scipy matplotlib snakemake
```

## Running

Put your `.pdb.gz` files in a folder, then:

```
bash scripts/run_all.sh /path/to/pdbs ARG
```

Output goes to `final/`.

## Results

- `final/angles.tsv` — all computed angles
- `final/pdb_list.txt` — PDB IDs that contributed data
- `final/angle_plot.png` — density plot

![angle distribution](final/angle_plot.png)
