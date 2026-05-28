import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.stats import fisher_exact
import os

PAE_PATH = "data/pae.scores_atp2b2_wt_dec2024_model_0.tsv"
PLDDT_CUTOFF = 50
PAE_CUTOFF = 15


def get_schema_asd_from_mark():
    case_df = pd.read_csv("data/atp2b2_case_variants_schema_asd.tsv", sep="\t")
    control_df = pd.read_csv("data/atp2b2_control_variants_schema_asd.tsv", sep="\t")

    case_df["is_case"] = 1
    control_df["is_case"] = 0

    df = pd.concat([case_df, control_df], ignore_index=True)
    return df


# ---------------------------------------------------------------
# Compute pairwise residue distances (atom-level) + inflate by PAE
# ---------------------------------------------------------------
def get_pairwise_distances(pdb_file, pae_path=PAE_PATH, pae_cutoff=PAE_CUTOFF):
    pae_scores = pd.read_csv(pae_path, sep="\t").values

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    residue_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_atoms.append([atom.get_coord() for atom in residue])

    num_residues = len(residue_atoms)
    pairwise_distances = np.full((num_residues, num_residues), np.inf)

    for i in range(num_residues):
        ai = np.array(residue_atoms[i])
        for j in range(i + 1, num_residues):
            aj = np.array(residue_atoms[j])
            d = np.linalg.norm(ai[:, None, :] - aj[None, :, :], axis=-1).min()
            pairwise_distances[i, j] = d
            pairwise_distances[j, i] = d

    # Inflate distances for residue pairs with poor mutual PAE
    bad = (pae_scores > pae_cutoff) & (pae_scores.T > pae_cutoff)
    pairwise_distances[bad] = 1000

    return pairwise_distances


# ---------------------------------------------------------------
# Spatial enrichment scan + return all significant neighborhoods
# ---------------------------------------------------------------
def neighborhood_test(df, pdb_file, radius, print_neighborhoods=False, print_top_neighborhood=True):
    df = df.copy()

    # Ensure aa_pos exists
    if "aa_pos" not in df.columns:
        if "Mutation" not in df.columns:
            raise ValueError("df must contain either 'aa_pos' or 'Mutation' column.")
        df["aa_pos"] = df["Mutation"].str.extract(r"(\d+)").astype(int)

    # Load structure and compute average residue pLDDT from B-factors
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    residue_plddt = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()[1]
                avg_plddt = float(np.mean([atom.bfactor for atom in residue]))
                residue_plddt[residue_id] = avg_plddt

    # Restrict analysis to well-modeled residues
    valid_residues = {res_id for res_id, v in residue_plddt.items() if v > PLDDT_CUTOFF}

    total_variants = len(df)
    df_filtered = df[df["aa_pos"].isin(valid_residues)].copy()
    kept_variants = len(df_filtered)
    dropped_variants = total_variants - kept_variants

    dropped_residues = set(df.loc[~df["aa_pos"].isin(valid_residues), "aa_pos"].unique())
    kept_unique_residues = sorted(df_filtered["aa_pos"].unique())

    print("\n--- Filtering summary ---")
    print(f"Variants before filtering: {total_variants}")
    print(f"Variants after filtering:  {kept_variants}")
    print(f"Variants dropped:         {dropped_variants}")
    print(f"Unique residue positions retained: {len(kept_unique_residues)}")
    print(f"Unique residue positions dropped:  {len(dropped_residues)}")
    print("-------------------------\n")

    # Bonferroni across unique tested residues
    num_residues_tested = len(kept_unique_residues)
    p_cutoff = 0.05 / num_residues_tested if num_residues_tested > 0 else np.nan
    print(f"Number of residues to test: {num_residues_tested}")
    print(f"Bonferroni cutoff is {p_cutoff}")

    pairwise_distances = get_pairwise_distances(pdb_file)

    n_case = int(df_filtered["is_case"].sum())
    n_ctrl = int((1 - df_filtered["is_case"]).sum())

    pval_dict = {}
    significant_centers = []
    significant_nbhd_tables = []

    # Scan each unique residue once
    for r in kept_unique_residues:
        neighborhood = np.where(pairwise_distances[r - 1, :] <= radius)[0] + 1
        in_nbhd = df_filtered["aa_pos"].isin(neighborhood)

        nbhd_case = int(df_filtered.loc[in_nbhd, "is_case"].sum())
        nbhd_ctrl = int((1 - df_filtered.loc[in_nbhd, "is_case"]).sum())

        if nbhd_case + nbhd_ctrl == 0:
            continue

        table = [[nbhd_case, nbhd_ctrl], [n_case - nbhd_case, n_ctrl - nbhd_ctrl]]
        _, p_value = fisher_exact(table, alternative="greater")
        pval_dict[r] = float(p_value)

        if p_value < p_cutoff:
            significant_centers.append(r)

            nbhd_df = df_filtered.loc[in_nbhd].copy()
            nbhd_df["center"] = r
            nbhd_df["center_p_value"] = float(p_value)
            significant_nbhd_tables.append(nbhd_df)

            print(r, p_value)
            if print_neighborhoods:
                print(nbhd_df)

    # Residue-level p-values table
    df_pvals = (
        pd.DataFrame({"aa_pos": list(pval_dict.keys()), "p_value": list(pval_dict.values())})
        .sort_values("p_value", ascending=True)
        .reset_index(drop=True)
    )

    # Merge residue-level p-values back onto all filtered variants
    df_results = df_filtered.merge(df_pvals, on="aa_pos", how="inner")
    df_results = df_results.drop(columns=["ddG"], errors="ignore")

    print(f"Variants included in analysis: {len(df_results)}")
    print(f"Unique residues tested (with p-values): {df_pvals['aa_pos'].nunique()}")

    # Combine all significant neighborhoods into one table
    if significant_nbhd_tables:
        df_sig_nbhd = pd.concat(significant_nbhd_tables, ignore_index=True)
    else:
        df_sig_nbhd = df_filtered.iloc[0:0].copy()
        df_sig_nbhd["center"] = pd.Series(dtype=int)
        df_sig_nbhd["center_p_value"] = pd.Series(dtype=float)

    # Top neighborhood summary
    if df_pvals.empty:
        top_residues = []
        top_p = None
        top_center = None
        df_top_single = df_filtered.iloc[0:0].copy()
    else:
        top_center = int(df_pvals.loc[0, "aa_pos"])
        top_p = float(df_pvals.loc[0, "p_value"])
        top_residues = (np.where(pairwise_distances[top_center - 1, :] <= radius)[0] + 1).tolist()
        df_top_single = df_filtered[df_filtered["aa_pos"].isin(top_residues)].copy()

    print("\n--- Top neighborhood ---")
    print(f"Top center residue: {top_center}")
    print(f"Top neighborhood p-value: {top_p}")
    print(f"Neighborhood size (unique residues): {len(set(top_residues))}")
    print(f"Variants in top neighborhood: {len(df_top_single)}")
    print(f"Number of significant centers: {len(significant_centers)}")

    if print_top_neighborhood:
        cols = ["aa_pos", "is_case"]
        for c in ["Mutation", "Protein_Change", "variant_id", "Sample", "cohort", "source"]:
            if c in df_top_single.columns and c not in cols:
                cols.append(c)
        cols = [c for c in cols if c in df_top_single.columns]

        print("\nVariants in top neighborhood (sorted by aa_pos):")
        if len(df_top_single) == 0:
            print("(none)")
        else:
            print(df_top_single[cols].sort_values("aa_pos").to_string(index=False))

    print("------------------------\n")

    return df_results, df_sig_nbhd, top_residues, top_p


df = load_schema_asd_variants()
pdb_file = "data/atp2b2_wt_dec2024_model_0.pdb"
radius = 15

df_results, df_sig_nbhd, top_residues, top_p = neighborhood_test(
    df, pdb_file, radius, print_neighborhoods=False, print_top_neighborhood=True
)