#!/usr/bin/env python3
"""
Permutation test for 3DNT — family-wise error rate
==================================================

Case and control labels are shuffled while variant positions are held
fixed, so the null preserves the spatial distribution of variants and the
total number of cases and controls. The scan is re-run at every tested
residue on each permutation, and the minimum p-value per permutation forms
the family-wise null.

This is the primary multiple-testing correction reported for 3DNT.
Neighborhoods overlap, so the Bonferroni threshold over tested residues is
conservative; both are printed for comparison.

Usage
-----
  python python_code/3DNT_permutation.py
  python python_code/3DNT_permutation.py --structure cryoem --n-perm 10000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent

PLDDT_CUTOFF = 50.0
PAE_CUTOFF = 15.0
RADIUS = 15.0


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def load_variants(case_path: Path, control_path: Path) -> pd.DataFrame:
    """Case and control tables, concatenated with an is_case label."""
    case_df = pd.read_csv(case_path, sep="\t")
    control_df = pd.read_csv(control_path, sep="\t")
    case_df["is_case"] = 1
    control_df["is_case"] = 0
    df = pd.concat([case_df, control_df], ignore_index=True)

    if "aa_pos" not in df.columns:
        for col in ("Mutation", "HGVSp/c"):
            if col in df.columns:
                df["aa_pos"] = df[col].str.extract(r"(\d+)").astype(int)
                break
        else:
            raise ValueError(
                "Variant tables must contain 'aa_pos', 'Mutation' or 'HGVSp/c'.")
    return df


def load_structure(pdb_path: Path, chain_id: str | None = None):
    """Ordered residues of one chain: (residue_number, coords, mean B-factor)."""
    model = next(PDBParser(QUIET=True).get_structure("s", str(pdb_path)).get_models())
    chains = list(model)
    if chain_id is not None:
        chains = [c for c in chains if c.id == chain_id]
        if not chains:
            raise ValueError(f"chain '{chain_id}' not found in {pdb_path}")
    residues = []
    for residue in chains[0]:
        if residue.get_id()[0] != " ":
            continue
        residues.append((
            residue.get_id()[1],
            np.array([a.get_coord() for a in residue]),
            float(np.mean([a.get_bfactor() for a in residue])),
        ))
    return residues


def load_pae(pae_path: Path, n_expected: int) -> np.ndarray:
    """N x N PAE matrix, tolerant of a header row and a row-name column."""
    rows = [ln.rstrip("\n").split("\t")
            for ln in open(pae_path) if ln.strip()]

    def numeric(tokens):
        try:
            [float(t) for t in tokens]
            return True
        except ValueError:
            return False

    if not numeric(rows[0]):
        rows = rows[1:]
    if len(rows[0]) == n_expected + 1:
        rows = [r[1:] for r in rows]

    pae = np.array(rows, dtype=float)
    if pae.shape != (n_expected, n_expected):
        raise ValueError(f"PAE matrix is {pae.shape}, expected "
                         f"({n_expected}, {n_expected})")
    return pae


# ---------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------
def pairwise_min_distances(residues, pae=None, pae_cutoff=PAE_CUTOFF, far=1000.0):
    """Minimum interatomic distance between every residue pair.

    A residue lies within its own neighborhood: the diagonal is set to zero,
    explicitly and after the PAE step, so the convention is visible and
    identical to the unconditional and cryo-EM analyses.
    """
    n = len(residues)
    d = np.zeros((n, n), dtype=float)
    for i in range(n):
        ai = residues[i][1]
        for j in range(i + 1, n):
            aj = residues[j][1]
            d[i, j] = d[j, i] = np.linalg.norm(
                ai[:, None, :] - aj[None, :, :], axis=-1).min()

    if pae is not None:
        bad = (pae > pae_cutoff) & (pae.T > pae_cutoff)
        np.fill_diagonal(bad, False)
        d[bad] = far

    np.fill_diagonal(d, 0.0)
    return d


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------
def build_membership(df, residues, dist, radius, plddt_cutoff):
    """Which variants fall inside each centre's neighborhood.

    The distance matrix is fixed across permutations, so membership is
    computed once and reused.
    """
    res_nums = [r[0] for r in residues]
    bfac = {r[0]: r[2] for r in residues}
    row_of = {num: i for i, num in enumerate(res_nums)}

    well_modelled = {n for n in res_nums if bfac[n] > plddt_cutoff}
    kept = df[df["aa_pos"].isin(well_modelled)].copy().reset_index(drop=True)
    centers = sorted(kept["aa_pos"].unique().tolist())
    pos = kept["aa_pos"].to_numpy()

    M = np.zeros((len(centers), len(kept)), dtype=bool)
    for k, c in enumerate(centers):
        near = dist[row_of[c], :] <= radius
        nbhd = {res_nums[i] for i in np.where(near)[0]}
        M[k] = np.isin(pos, list(nbhd))
    return M, centers, kept


def scan_labels(M, case_vec, n_case, n_ctrl):
    """One-sided Fisher p at every centre, for one labelling."""
    a = M @ case_vec
    b = M.sum(axis=1) - a
    p = np.ones(M.shape[0])
    for k in range(M.shape[0]):
        if a[k] + b[k] == 0:
            continue
        _, p[k] = fisher_exact(
            [[a[k], b[k]], [n_case - a[k], n_ctrl - b[k]]], alternative="greater")
    return p, a, b


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structure", choices=["alphafold", "cryoem"],
                    default="alphafold")
    ap.add_argument("--case", type=Path,
                    default=ROOT / "data" / "atp2b2_case_variants_schema_asd.tsv")
    ap.add_argument("--control", type=Path,
                    default=ROOT / "data" / "atp2b2_control_variants_schema_asd.tsv")
    ap.add_argument("--mapping", type=Path,
                    default=ROOT / "data" / "mutation_mapping_PMCA2za.csv",
                    help="Used only when --structure cryoem.")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--radius", type=float, default=RADIUS)
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_variants(args.case, args.control)

    if args.structure == "alphafold":
        pdb = ROOT / "data" / "atp2b2_wt_dec2024_model_0.pdb"
        residues = load_structure(pdb)
        pae = load_pae(ROOT / "data" /
                       "pae.scores_atp2b2_wt_dec2024_model_0.tsv", len(residues))
        plddt_cutoff = PLDDT_CUTOFF
    else:
        pdb = ROOT / "data" / "pmca_e1ca_model.pdb"
        residues = load_structure(pdb, chain_id="A")
        pae = None                       # no PAE for an experimental structure
        plddt_cutoff = 0.0               # B-factors are temperature factors
        m = pd.read_csv(args.mapping)
        m = m[m["status"] == "mapped"]
        pos_map = dict(zip(m["original_resi"], pd.to_numeric(m["new_resi"])))
        df["canonical_pos"] = df["aa_pos"]
        df["aa_pos"] = df["canonical_pos"].map(pos_map)
        df = df.dropna(subset=["aa_pos"]).copy()
        df["aa_pos"] = df["aa_pos"].astype(int)

    print(f"structure : {args.structure} ({pdb.name})")
    print(f"residues  : {len(residues)}")

    dist = pairwise_min_distances(residues, pae)
    M, centers, kept = build_membership(df, residues, dist,
                                        args.radius, plddt_cutoff)
    case_vec = kept["is_case"].to_numpy()
    n_case, n_ctrl = int(case_vec.sum()), int((1 - case_vec).sum())

    p_obs, a_obs, b_obs = scan_labels(M, case_vec, n_case, n_ctrl)
    i_top = int(np.argmin(p_obs))
    p_min, top_ctr = float(p_obs[i_top]), centers[i_top]

    print(f"centres   : {len(centers)}   case {n_case} / control {n_ctrl}")
    print(f"observed top centre {top_ctr}, p = {p_min:.3e}\n")

    rng = np.random.default_rng(args.seed)
    perm_min = np.empty(args.n_perm)
    for t in range(args.n_perm):
        perm_min[t] = scan_labels(M, rng.permutation(case_vec),
                                  n_case, n_ctrl)[0].min()
        if (t + 1) % 250 == 0:
            print(f"  {t + 1}/{args.n_perm}")

    thresh = float(np.percentile(perm_min, 5))
    n_at_least = int((perm_min <= p_min).sum())
    fwer_p = (n_at_least + 1) / (args.n_perm + 1)   # never reports exactly 0

    print(f"\n  5th-percentile FWER threshold : {thresh:.3e}")
    print(f"  Bonferroni threshold          : {0.05 / len(centers):.3e}")
    print(f"  permutations reaching p_obs   : {n_at_least}/{args.n_perm}")
    print(f"  FWER-corrected p              : {fwer_p:.3e}"
          f"{f'  (< 1/{args.n_perm + 1})' if n_at_least == 0 else ''}")

    surv = sorted(((centers[k], p_obs[k], int(a_obs[k]), int(b_obs[k]))
                   for k in np.where(p_obs <= thresh)[0]), key=lambda x: x[1])
    print(f"\n  centres passing the permutation threshold: {len(surv)}")
    for c, p, a, b in surv:
        print(f"    {c:>5}   p = {p:.3e}   case_in = {a:>3}   ctrl_in = {b:>3}")

    tag = args.structure
    pd.DataFrame({"center": centers, "p_observed": p_obs,
                  "case_in": a_obs, "ctrl_in": b_obs,
                  "passes_perm_threshold": p_obs <= thresh}
                 ).sort_values("p_observed").to_csv(
        args.out_dir / f"3dnt_permutation_observed_{tag}.csv", index=False)
    pd.DataFrame({"perm_min_p": perm_min}).to_csv(
        args.out_dir / f"3dnt_permutation_min_p_{tag}.csv", index=False)
    print(f"\n  -> {args.out_dir}/3dnt_permutation_observed_{tag}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())