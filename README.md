# 3D Neighborhood Test (3DNT)

Code for the 3D neighborhood test (3DNT), a structure-based method for detecting spatial clustering of missense variants in protein structures, as applied to ATP2B2 in:

> Gerges, S. et al. *Genetic and structural evidence links Ca²⁺ dysregulation and ATP2B2 to neuropsychiatric illness.* (Manuscript submitted, 2026.)

3DNT places a 15 Å sphere around each residue carrying at least one variant and uses a one-sided Fisher's exact test to ask whether case-derived variants are concentrated in that neighborhood relative to the rest of the protein, with Bonferroni correction across tested residues.

This repository contains the implementations used for the two structural models in the paper: the AlphaFold3 prediction of canonical ATP2B2, and the cryo-EM structure of human PMCA2z/a (PDB 28JP).

## Repository layout

```
.
├── README.md
├── LICENSE
├── data/
│   ├── atp2b2_case_variants_schema_asd.tsv
│   ├── atp2b2_control_variants_schema_asd.tsv
│   ├── atp2b2_wt_dec2024_model_0.pdb        (AlphaFold3 model)
│   └── pae.scores_atp2b2_wt_dec2024_model_0.tsv
├── python_code/
│   ├── 3DNT.py                              (AlphaFold3 structure analysis)
│   ├── 3DNT_cryoEM.py                       (cryo-EM PMCA2z/a analysis)
│   └── get_tables.py                        (shared variant loader)
└── results/                                 (output of running the scripts)
```

The cryo-EM PDB file (PDB 28JP) is not included in this repository; download it from rcsb.org after publication and place it at `data/pmca_e1ca_model.pdb` to run `3DNT_cryoEM.py`.

## Requirements

- Python 3.10+
- numpy
- pandas
- scipy
- biopython >= 1.80

Install:

```
pip install numpy pandas scipy biopython
```

## Running

### AlphaFold3 analysis

This script runs 3DNT on the AlphaFold3-predicted structure of human ATP2B2 (UniProt Q01814-1). It uses the AlphaFold3 model (`data/atp2b2_wt_dec2024_model_0.pdb`) and the corresponding predicted aligned error matrix (`data/pae.scores_atp2b2_wt_dec2024_model_0.tsv`) together with the case and control variant tables in `data/`.


From the repository root:

```
python python_code/3DNT.py
```

Reads the case/control variant tables and the AlphaFold3 model in `data/`, applies a residue pLDDT > 50 filter and a residue-pair PAE inflation (> 15 Å in both directions), and prints per-residue Bonferroni-significant centers along with the top-ranked neighborhood.

### Cryo-EM analysis

From the repository root, after obtaining the cryo-EM PDB:

```
python python_code/3DNT_cryoEM.py
```

Maps variants from canonical ATP2B2 numbering (UniProt Q01814-1) to the PMCA2z/a isoform (UniProt Q01814-4) and then to the residue numbering of the cryo-EM model, runs 3DNT on the experimental structure, and writes:

- `results/mutation_mapping_PMCA2za.csv`
- `results/3dnt_cryoem_per_residue.csv`
- `results/3dnt_cryoem_union_stats.csv`

CLI arguments (defaults shown):

```
--case      data/atp2b2_case_variants_schema_asd.tsv
--control   data/atp2b2_control_variants_schema_asd.tsv
--pdb       data/pmca_e1ca_model.pdb
--out-dir   results
--radius    15.0
--chain     A
```

## Method

For each residue r in the protein structure:

1. Define a spherical neighborhood of radius 15 Å around r, using minimum atom-atom distances between r and every other residue (standard amino acid atoms only; alternate locations restricted to the primary conformer).
2. For AlphaFold3 models, exclude residues with mean pLDDT ≤ 50 and inflate distances between residue pairs whose predicted aligned error exceeds 15 Å in both directions.
3. Count case and control variants inside vs. outside the neighborhood.
4. Test for case enrichment with a one-sided Fisher's exact test (`alternative='greater'`).
5. Apply Bonferroni correction across the unique residues tested (those carrying at least one variant after filtering).

## Input formats

Variant TSVs must contain at minimum:

- `Mutation` — substitution string in canonical ATP2B2 numbering, e.g. `E457K`, or
- `aa_pos` — integer residue position in canonical ATP2B2 numbering.

When the cryo-EM script runs, canonical positions are remapped via global pairwise alignment to the PMCA2z/a isoform and then to the cryo-EM model's residue numbering. The script includes a hard-coded sanity check that V885 in the canonical numbering corresponds to V840 in the cryo-EM numbering (the anchor reported in the manuscript).

## Citation

If you use this code, please cite the manuscript above. A Zenodo archive of the released version is available at:

> Gerges, S. *3DNT_ATP2B2: 3D neighborhood test for ATP2B2.* Zenodo. https://doi.org/10.5281/zenodo.20444372

## License

MIT (see `LICENSE`).

## Contact

Sherif Gerges — sherif_gerges@g.harvard.edu
