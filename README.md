# 3DNT — ATP2B2

The 3D neighborhood test applied to ATP2B2, as reported in Gerges et al.

> Gerges S, Straarup NC, El-Brolosy MA, et al. Genetic and structural
> evidence links Ca²⁺ dysregulation and ATP2B2 to neuropsychiatric illness.
> bioRxiv 2025.08.25.672202.

## Method

For each residue carrying at least one variant, every residue within 15 Å is
collected into a neighborhood. A one-sided Fisher's exact test compares the
case:control ratio inside that neighborhood with the ratio across the rest of
the protein.

Four details define the test:

1. **Distances are minimum interatomic**, not Cα–Cα, so side-chain contacts
   count.
2. **A residue lies within its own neighborhood.** Variants at the center
   position contribute to the test at that center. This is set explicitly
   (`np.fill_diagonal(..., 0.0)`) and applies identically to the
   unconditional, conditional, AlphaFold3 and cryo-EM analyses.
3. **Residue pairs with poor mutual confidence are down-weighted.** Where the
   predicted aligned error exceeds 15 Å in *both* directions, the pair is
   pushed to 1000 Å so it cannot form a neighborhood. The diagonal is exempt.
   This applies to AlphaFold3 models only.
4. **Only variant-bearing residues are tested as centers.**

For the AlphaFold3 model, residues with mean pLDDT ≤ 50 are excluded. The
cryo-EM structure carries real temperature factors in the B-factor column, so
no confidence filter is applied there.

**Multiple testing.** Neighborhoods overlap, so per-center p-values are not
independent and Bonferroni over tested residues is conservative. The primary
correction is permutation: case and control labels are shuffled 1,000 times
with variant positions held fixed, the scan re-run at every tested residue,
and the minimum p-value per permutation forms the family-wise null. Both
thresholds are reported.

## Install

```bash
git clone https://github.com/sherifgerges/3DNT_ATP2B2.git
cd 3DNT_ATP2B2
pip install numpy pandas scipy biopython
```

## Run

```bash
python python_code/3DNT.py                                 # AlphaFold3 scan
python python_code/3DNT_cryoEM.py                          # cryo-EM scan
python python_code/3DNT_permutation.py                     # permutation, AlphaFold3
python python_code/3DNT_permutation.py --structure cryoem  # permutation, cryo-EM
```

Outputs are written to `results/`, which is not tracked.

## Expected output

**AlphaFold3** — 324 variants in, 251 after the pLDDT filter, 217 residues
tested, Bonferroni 2.304 × 10⁻⁴:

| center | p |
|---|---|
| 885 | 1.016 × 10⁻⁵ |
| 107 | 3.403 × 10⁻⁵ |
| 448 | 3.403 × 10⁻⁵ |
| 457 | 3.743 × 10⁻⁵ |
| 113 | 5.819 × 10⁻⁵ |
| 446 | 1.444 × 10⁻⁴ |
| 913 | 1.444 × 10⁻⁴ |

Permutation: threshold 9.032 × 10⁻⁴, 0/1000 permutations reached the observed
value at residue 885, FWER-corrected P < 1/1001.

**Cryo-EM** — 253 of 324 variants map to the PMCA2z/a isoform, 219 residues
tested, 62 case and 191 control. Residues 107, 412 and 840 share a p-value of
9.175 × 10⁻⁶ (OR 8.99) because their neighborhoods contain the same variants;
412 and 840 are canonical E457 and V885. Permutation threshold 1.067 × 10⁻³,
0/1000.

## Data

| file | contents |
|---|---|
| `data/atp2b2_case_variants_schema_asd.tsv` | 83 case variants, canonical numbering |
| `data/atp2b2_control_variants_schema_asd.tsv` | 241 control variants |
| `data/atp2b2_wt_dec2024_model_0.pdb` | AlphaFold3 model, pLDDT in the B-factor column |
| `data/pae.scores_atp2b2_wt_dec2024_model_0.tsv` | predicted aligned error, 1243 × 1243 |
| `data/pmca_e1ca_model.pdb` | cryo-EM structure, human PMCA2z/a, E1-Ca state |
| `data/mutation_mapping_PMCA2za.csv` | canonical ATP2B2 → PMCA2z/a residue correspondence |

Variants observed in both cases and controls are excluded, as are synonymous
changes; recurrent observations of the same substitution are collapsed to one
row.

**Residue numbering.** Canonical ATP2B2 and PMCA2z/a numbering differ by a
non-constant offset: 0 before the z/a splice deletion, 45 after. E457 is E412
in the structure, V885 is V840, T107 is unchanged. The correspondence is
shipped as `mutation_mapping_PMCA2za.csv` rather than recomputed by sequence
alignment. The cryo-EM script also *writes* a file of that name to `results/`;
the copy in `data/` is the input.

## Limits

- The test identifies regions where case variants concentrate. It does not
  identify individual causal variants; treat the output as prioritization.
- Overlapping neighborhoods mean the number of significant centers exceeds the
  number of distinct clusters. Residues 446 and 448 are two apart; 107, 113 and
  117 sit together; 885 and 913 are both in the Ca²⁺ site.
- The pLDDT and PAE filters exclude disordered regions, so signal there is
  missed rather than reported as absent.
- Variants seen only in controls are not necessarily benign.

## License

See `LICENSE`.