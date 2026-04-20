# 3D Neighborhood Test (ATP2B2)

This code performs a 3D neighborhood enrichment test to identify spatial clustering of missense variants in ATP2B2, comparing case vs control variants.

## Structure

3DNT/
- python_code/3DNT_atp2b2.py : main analysis script  
- data/ : input files  
- results/ : output (optional)

## Requirements

- Python 3
- numpy
- pandas
- scipy
- biopython

Install:
pip install numpy pandas scipy biopython

## Running

From the 3DNT directory:

cd 3DNT  
python python_code/3DNT_atp2b2.py

## Input files (in data/)

- atp2b2_case_variants_schema_asd.tsv
- atp2b2_control_variants_schema_asd.tsv
- atp2b2_wt_dec2024_model_0.pdb
- pae.scores_atp2b2_wt_dec2024_model_0.tsv

Variant files must contain either:
- aa_pos  
or  
- Mutation (e.g. E457K)

## Method

- Defines 15 Å residue-centered neighborhoods  
- Filters residues with pLDDT ≤ 50  
- Uses one-sided Fisher’s exact test  
- Applies Bonferroni correction  

## Output

Results are printed to the terminal.
