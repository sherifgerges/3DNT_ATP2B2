"""
3D Neighborhood Test (3DNT) for ATP2B2 on the human PMCA2z/a cryo-EM structure.

Workflow
--------
1. Remap missense variants from canonical ATP2B2 numbering (UniProt Q01814-1)
   to the PMCA2z/a isoform (UniProt Q01814-4), then to the residue numbering
   of the cryo-EM model.
2. Compute residue-residue minimum atom-atom distances from the cryo-EM PDB.
3. For each residue carrying at least one variant, define a 15 Angstrom
   spherical neighborhood and test for enrichment of case vs. control variants
   using a one-sided Fisher's exact test (alternative='greater').
4. Apply Bonferroni correction across the unique tested residues. Report
   per-residue statistics, the union of Bonferroni-significant neighborhoods,
   and write all outputs to disk.

Inputs (defaults are relative to the repository root)
-----------------------------------------------------
    data/atp2b2_case_variants_schema_asd.tsv
    data/atp2b2_control_variants_schema_asd.tsv
    data/pmca_e1ca_model.pdb

Outputs (written under --out-dir, default results/)
---------------------------------------------------
    mutation_mapping_PMCA2za.csv
    3dnt_cryoem_per_residue.csv
    3dnt_cryoem_union_stats.csv

Usage
-----
    python python_code/3DNT_atp2b2_cryoem.py \\
        --case data/atp2b2_case_variants_schema_asd.tsv \\
        --control data/atp2b2_control_variants_schema_asd.tsv \\
        --pdb data/pmca_e1ca_model.pdb \\
        --out-dir results \\
        --radius 15 \\
        --chain A
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData
from Bio.PDB import PDBParser, is_aa
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Canonical ATP2B2 (UniProt Q01814-1).
CANONICAL_SEQ = (
    "MGDMTNSDFYSKNQRNESSHGGEFGCTMEELRSLMELRGTEAVVKIKETYGDTEAICRRLKTSPVEGLPGTAPDLEKRKQIFGQNFI"
    "PPKKPKTFLQLVWEALQDVTLIILEIAAIISLGLSFYHPPGEGNEGCATAQGGAEDEGEAEAGWIEGAAILLSVICVVLVTAFNDWS"
    "KEKQFRGLQSRIEQEQKFTVVRAGQVVQIPVAEIVVGDIAQVKYGDLLPADGLFIQGNDLKIDESSLTGESDQVRKSVDKDPMLLSG"
    "THVMEGSGRMLVTAVGVNSQTGIIFTLLGAGGEEEEKKDKKGVKKGDGLQLPAADGAAASNAADSANASLVNGKMQDGNVDASQSKA"
    "KQQDGAAAMEMQPLKSAEGGDADDRKKASMHKKEKSVLQGKLTKLAVQIGKAGLVMSAITVIILVLYFTVDTFVVNKKPWLPECTPV"
    "YVQYFVKFFIIGVTVLVVAVPEGLPLAVTISLAYSVKKMMKDNNLVRHLDACETMGNATAICSDKTGTLTTNRMTVVQAYVGDVHYK"
    "EIPDPSSINTKTMELLINAIAINSAYTTKILPPEKEGALPRQVGNKTECGLLGFVLDLKQDYEPVRSQMPEEKLYKVYTFNSVRKSM"
    "STVIKLPDESFRMYSKGASEIVLKKCCKILNGAGEPRVFRPRDRDEMVKKVIEPMACDGLRTICVAYRDFPSSPEPDWDNENDILNE"
    "LTCICVVGIEDPVRPEVPEAIRKCQRAGITVRMVTGDNINTARAIAIKCGIIHPGEDFLCLEGKEFNRRIRNEKGEIEQERIDKIWP"
    "KLRVLARSSPTDKHTLVKGIIDSTHTEQRQVVAVTGDGTNDGPALKKADVGFAMGIAGTDVAKEASDIILTDDNFSSIVKAVMWGRN"
    "VYDSISKFLQFQLTVNVVAVIVAFTGACITQDSPLKAVQMLWVNLIMDTFASLALATEPPTETLLLRKPYGRNKPLISRTMMKNILG"
    "HAVYQLALIFTLLFVGEKMFQIDSGRNAPLHSPPSEHYTIIFNTFVMMQLFNEINARKIHGERNVFDGIFRNPIFCTIVLGTFAIQI"
    "VIVQFGGKPFSCSPLQLDQWMWCIFIGLGELVWGQVIATIPTSRLKFLKEAGRLTQKEEIPEEELNEDVEEIDHAERELRRGQILWF"
    "RGLNRIQTQIRVVKAFRSSLYEGLEKPESRTSIHNFMAHPEFRIEDSQPHIPLIDDTDLEEDAALKQNSSPPSSLNKNNSAIDSGIN"
    "LTTDTSKSATSSSPGSPIHSLETSL"
)

# PMCA2z/a isoform (UniProt Q01814-4).
ISOFORM_SEQ = (
    "GEFGCTMEELRSLMELRGTEAVVKIKETYGDTEAICRRLKTSPVEGLPGTAPDLEKRKQIFGQNFIPPKKPKTFLQLVWEALQDVTL"
    "IILEIAAIISLGLSFYHPPGEGNEGCATAQGGAEDEGEAEAGWIEGAAILLSVICVVLVTAFNDWSKEKQFRGLQSRIEQEQKFTVV"
    "RAGQVVQIPVAEIVVGDIAQVKYGDLLPADGLFIQGNDLKIDESSLTGESDQVRKSVDKDPMLLSGTHVMEGSGRMLVTAVGVNSQT"
    "GIIFTLLKSVLQGKLTKLAVQIGKAGLVMSAITVIILVLYFTVDTFVVNKKPWLPECTPVYVQYFVKFFIIGVTVLVVAVPEGLPLA"
    "VTISLAYSVKKMMKDNNLVRHLDACETMGNATAICSDKTGTLTTNRMTVVQAYVGDVHYKEIPDPSSINTKTMELLINAIAINSAYT"
    "TKILPPEKEGALPRQVGNKTECGLLGFVLDLKQDYEPVRSQMPEEKLYKVYTFNSVRKSMSTVIKLPDESFRMYSKGASEIVLKKCC"
    "KILNGAGEPRVFRPRDRDEMVKKVIEPMACDGLRTICVAYRDFPSSPEPDWDNENDILNELTCICVVGIEDPVRPEVPEAIRKCQRA"
    "GITVRMVTGDNINTARAIAIKCGIIHPGEDFLCLEGKEFNRRIRNEKGEIEQERIDKIWPKLRVLARSSPTDKHTLVKGIIDSTHTE"
    "QRQVVAVTGDGTNDGPALKKADVGFAMGIAGTDVAKEASDIILTDDNFSSIVKAVMWGRNVYDSISKFLQFQLTVNVVAVIVAFTGA"
    "CITQDSPLKAVQMLWVNLIMDTFASLALATEPPTETLLLRKPYGRNKPLISRTMMKNILGHAVYQLALIFTLLFVGEKMFQIDSGRN"
    "APLHSPPSEHYTIIFNTFVMMQLFNEINARKIHGERNVFDGIFRNPIFCTIVLGTFAIQIVIVQFGGKPFSCSPLQLDQWMWCIFIG"
    "LGELVWGQVIATI"
)

# Anchor: V885 in the canonical sequence is reported as V840 in the PMCA2z/a
# cryo-EM model (see manuscript Figure 5). Used as a sanity check to guarantee
# that the alignment and PDB residue mapping have not drifted.
ANCHOR_CANONICAL_POS = 885
ANCHOR_ISOFORM_POS = 840

# 3-letter -> 1-letter amino acid lookup, plus the set of canonical 1-letter codes.
_AA3_TO_1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}
_AA1_SET = set(IUPACData.protein_letters)


# ---------------------------------------------------------------------------
# Variant I/O
# ---------------------------------------------------------------------------

def load_schema_asd_variants(case_path: Path, control_path: Path) -> pd.DataFrame:
    """Load case and control variant tables and concatenate with an is_case label.

    Each input TSV is expected to contain at least one of:
        - 'aa_pos': integer residue position in canonical ATP2B2 numbering, or
        - 'Mutation': a string such as 'E457K' from which aa_pos can be parsed.

    Returns
    -------
    pandas.DataFrame with columns from the inputs plus an integer 'is_case' column
    (1 for case, 0 for control).
    """
    case_df = pd.read_csv(case_path, sep="\t")
    control_df = pd.read_csv(control_path, sep="\t")
    case_df["is_case"] = 1
    control_df["is_case"] = 0
    df = pd.concat([case_df, control_df], ignore_index=True)

    # Backfill aa_pos from Mutation if needed.
    if "aa_pos" not in df.columns:
        if "Mutation" not in df.columns:
            raise ValueError("Variant tables must contain 'aa_pos' or 'Mutation'.")
        df["aa_pos"] = df["Mutation"].str.extract(r"(\d+)").astype(int)
    return df


# ---------------------------------------------------------------------------
# Mutation parsing
# ---------------------------------------------------------------------------

_SUB_RE = re.compile(r"^([A-Za-z])\s*(\d+)\s*([A-Za-z\*])$")
_DEL_RE = re.compile(r"^([A-Za-z])\s*(\d+)\s*del$", re.IGNORECASE)


def parse_mutation(mut_str: str) -> dict:
    """Parse a mutation string into its components.

    Supports substitutions (E412K), stop gains (G700*), and single-residue
    deletions (A123del). The 3DNT analysis itself uses only substitutions;
    the broader parser is retained so the same routine can be reused upstream.
    """
    s = mut_str.strip()
    m = _SUB_RE.match(s)
    if m:
        ref, pos, alt = m.groups()
        return {
            "ref": ref.upper(),
            "pos": int(pos),
            "alt": alt.upper(),
            "kind": "stop" if alt == "*" else "sub",
            "raw": s,
        }
    m = _DEL_RE.match(s)
    if m:
        ref, pos = m.groups()
        return {
            "ref": ref.upper(),
            "pos": int(pos),
            "alt": "-",
            "kind": "del",
            "raw": s,
        }
    return {"ref": None, "pos": None, "alt": None, "kind": "unknown", "raw": s}


# ---------------------------------------------------------------------------
# Sequence alignment helpers
# ---------------------------------------------------------------------------

def _make_aligner(
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
) -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    return aligner


def align_position_map(seq_a: str, seq_b: str) -> dict[int, int | None]:
    """Globally align seq_a to seq_b and return a 1-indexed position map.

    For every position i in seq_a (1-indexed):
        - returns the corresponding 1-indexed position in seq_b if both residues
          are aligned without a gap, or
        - returns None if position i in seq_a aligns to a gap in seq_b.
    """
    aligner = _make_aligner()
    alignment = aligner.align(seq_a, seq_b)[0]
    aligned_a, aligned_b = str(alignment[0]), str(alignment[1])

    pos_a = pos_b = 0
    mapping: dict[int, int | None] = {}
    for char_a, char_b in zip(aligned_a, aligned_b):
        if char_a != "-":
            pos_a += 1
        if char_b != "-":
            pos_b += 1
        if char_a != "-" and char_b != "-":
            mapping[pos_a] = pos_b
        elif char_a != "-" and char_b == "-":
            mapping[pos_a] = None
    return mapping


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def pdb_sequence_and_index(pdb_path: Path, chain_id: str) -> tuple[str, dict[int, tuple[int, str]]]:
    """Return the PDB sequence (1-letter codes) for a chain plus a mapping from
    the 1-indexed position in that sequence to the (resnum, icode) in the PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", str(pdb_path))
    model = next(structure.get_models())
    chain = model[chain_id]

    chars: list[str] = []
    seq_pos_to_pdb: dict[int, tuple[int, str]] = {}
    seq_pos = 0
    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        aa1 = _AA3_TO_1.get(residue.get_resname().upper(), "X")
        if aa1 not in _AA1_SET:
            aa1 = "X"
        seq_pos += 1
        chars.append(aa1)
        _het, resnum, icode = residue.get_id()
        icode_str = icode.strip() if isinstance(icode, str) else ""
        seq_pos_to_pdb[seq_pos] = (int(resnum), icode_str)
    return "".join(chars), seq_pos_to_pdb


# ---------------------------------------------------------------------------
# Remap variants: canonical -> isoform -> PDB
# ---------------------------------------------------------------------------

def remap_variants_to_pdb(
    df: pd.DataFrame,
    canonical_seq: str,
    isoform_seq: str,
    pdb_path: Path,
    chain_id: str = "A",
) -> pd.DataFrame:
    """Remap variants from canonical ATP2B2 positions to PDB residue numbers."""
    canon_to_iso = align_position_map(canonical_seq, isoform_seq)
    pdb_seq, pdb_seq_pos_to_pdb_resnum = pdb_sequence_and_index(pdb_path, chain_id)
    iso_to_pdb_seq_pos = align_position_map(isoform_seq, pdb_seq)

    # Sanity check: the V885 -> V840 anchor must survive the canonical->isoform map.
    anchor_iso = canon_to_iso.get(ANCHOR_CANONICAL_POS)
    if anchor_iso != ANCHOR_ISOFORM_POS:
        raise RuntimeError(
            f"Alignment anchor check failed: canonical residue "
            f"{ANCHOR_CANONICAL_POS} mapped to isoform {anchor_iso}, "
            f"expected {ANCHOR_ISOFORM_POS}. Inspect sequences before proceeding."
        )

    rows: list[dict] = []
    for _, row in df.iterrows():
        mut_str = row.get("Mutation")
        info = parse_mutation(mut_str) if isinstance(mut_str, str) else {
            "ref": None, "pos": None, "alt": None, "kind": "unknown", "raw": mut_str
        }
        canon_pos = int(row["aa_pos"])
        is_case = int(row["is_case"])

        iso_pos = canon_to_iso.get(canon_pos)
        if iso_pos is None:
            rows.append({
                "raw": info["raw"], "original_resi": canon_pos, "new_resi": None,
                "is_case": is_case, "status": "absent_in_isoform", "remapped_mut": None,
            })
            continue

        pdb_seq_pos = iso_to_pdb_seq_pos.get(iso_pos)
        if pdb_seq_pos is None or pdb_seq_pos not in pdb_seq_pos_to_pdb_resnum:
            rows.append({
                "raw": info["raw"], "original_resi": canon_pos, "new_resi": None,
                "is_case": is_case, "status": "absent_in_structure", "remapped_mut": None,
            })
            continue

        resnum, icode = pdb_seq_pos_to_pdb_resnum[pdb_seq_pos]
        new_resi = f"{resnum}{icode}" if icode else str(resnum)

        remapped_mut: str | None = None
        if info["kind"] in ("sub", "stop"):
            remapped_mut = f"{info['ref']}{new_resi}{info['alt']}"
        elif info["kind"] == "del":
            remapped_mut = f"{info['ref']}{new_resi}del"

        rows.append({
            "raw": info["raw"], "original_resi": canon_pos, "new_resi": new_resi,
            "is_case": is_case, "status": "mapped", "remapped_mut": remapped_mut,
        })

    out = pd.DataFrame(rows)
    n_total = len(out)
    n_mapped = int((out["status"] == "mapped").sum())
    print(f"  Total variants:           {n_total}")
    print(f"  Mapped:                   {n_mapped}")
    print(f"  Absent in isoform:        {int((out['status'] == 'absent_in_isoform').sum())}")
    print(f"  Absent in structure:      {int((out['status'] == 'absent_in_structure').sum())}")
    return out


# ---------------------------------------------------------------------------
# Pairwise minimum atom-atom distances
# ---------------------------------------------------------------------------

def compute_pairwise_distances(pdb_path: Path) -> tuple[np.ndarray, list[tuple[str, int, str]]]:
    """Compute pairwise minimum atom-atom distances across standard residues.

    Standard residues only. Hetero atoms and waters are skipped. Alternate
    locations are restricted to the primary conformer (' ', 'A', '').
    The diagonal is set to 0.0 so each residue is in its own neighborhood.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))

    ids: list[tuple[str, int, str]] = []
    residue_atoms: list[np.ndarray] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if not is_aa(residue, standard=True):
                    continue
                coords = [
                    atom.get_coord()
                    for atom in residue
                    if atom.get_altloc() in (" ", "A", "")
                ]
                if not coords:
                    continue
                ids.append((chain.id, int(residue.id[1]), residue.id[2]))
                residue_atoms.append(np.asarray(coords, dtype=float))

    n = len(residue_atoms)
    distances = np.full((n, n), np.inf, dtype=float)
    for i in range(n):
        ai = residue_atoms[i]
        for j in range(i + 1, n):
            aj = residue_atoms[j]
            diff = ai[:, None, :] - aj[None, :, :]
            d2 = np.sum(diff * diff, axis=2)
            d = float(np.sqrt(np.min(d2)))
            distances[i, j] = d
            distances[j, i] = d
    np.fill_diagonal(distances, 0.0)
    return distances, ids


# ---------------------------------------------------------------------------
# 3D neighborhood test
# ---------------------------------------------------------------------------

def neighborhood_test(
    variants: pd.DataFrame,
    pdb_path: Path,
    radius: float,
    chain: str = "A",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the 3D neighborhood test on a single chain.

    Parameters
    ----------
    variants : DataFrame with integer columns 'aa_pos' (PDB residue number on
               the selected chain) and 'is_case' (1/0).
    pdb_path : path to the PDB file.
    radius   : neighborhood radius in Angstroms.
    chain    : chain ID to scan (single-chain analysis only).

    Returns
    -------
    variants_with_p : input variants merged with per-residue statistics
    per_residue     : one row per residue center tested
    union_stats     : single-row summary of the union of significant neighborhoods
    """
    distances, ids = compute_pairwise_distances(pdb_path)

    # Map (chain, resnum) -> array index, restricted to the requested chain
    # and to residues without insertion codes (the typical scan case).
    pos_to_idx: dict[int, int] = {
        resnum: k
        for k, (c, resnum, icode) in enumerate(ids)
        if c == chain and icode == " "
    }

    df = variants.copy()
    df = df.loc[df["aa_pos"].notna()].copy()
    df["aa_pos"] = df["aa_pos"].astype(int)
    df["is_case"] = df["is_case"].astype(int)
    df = df.loc[df["aa_pos"].isin(pos_to_idx)].copy()
    df["idx"] = df["aa_pos"].map(pos_to_idx)

    n_case = int(df["is_case"].sum())
    n_ctrl = int((1 - df["is_case"]).sum())

    residues_to_scan = sorted(df["aa_pos"].unique().tolist())
    m = len(residues_to_scan)
    p_cutoff = 0.05 / m if m > 0 else float("nan")

    print(f"  Residues scanned:         {m}")
    print(f"  Bonferroni cutoff:        {p_cutoff:.3e}")
    print(f"  Case variants used:       {n_case}")
    print(f"  Control variants used:    {n_ctrl}")

    idxs = df["idx"].to_numpy()
    is_case = df["is_case"].to_numpy()

    results: list[dict] = []
    sig_center_to_nbhd_idx: dict[int, np.ndarray] = {}

    for resnum in residues_to_scan:
        center_idx = pos_to_idx[resnum]
        nbhd_idx = np.where(distances[center_idx, :] <= radius)[0]
        in_nbhd = np.isin(idxs, nbhd_idx)

        nbhd_case = int(is_case[in_nbhd].sum())
        nbhd_ctrl = int((1 - is_case[in_nbhd]).sum())
        if nbhd_case + nbhd_ctrl == 0:
            continue

        contingency = [
            [nbhd_case, nbhd_ctrl],
            [n_case - nbhd_case, n_ctrl - nbhd_ctrl],
        ]
        odds_ratio, p_value = fisher_exact(contingency, alternative="greater")

        results.append({
            "aa_pos": resnum,
            "center_idx": center_idx,
            "nbhd_case": nbhd_case,
            "nbhd_ctrl": nbhd_ctrl,
            "total_case": n_case,
            "total_ctrl": n_ctrl,
            "case_pct": nbhd_case / n_case if n_case > 0 else float("nan"),
            "ctrl_pct": nbhd_ctrl / n_ctrl if n_ctrl > 0 else float("nan"),
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "bonferroni_cutoff": p_cutoff,
            "bonferroni_significant": bool(p_value < p_cutoff),
        })

        if p_value < p_cutoff:
            sig_center_to_nbhd_idx[resnum] = nbhd_idx

    per_residue = (
        pd.DataFrame(results)
        .sort_values("p_value")
        .reset_index(drop=True)
    )

    # Union of Bonferroni-significant neighborhoods.
    if sig_center_to_nbhd_idx:
        union_nbhd_idx = np.unique(np.concatenate(list(sig_center_to_nbhd_idx.values())))
        in_union = np.isin(idxs, union_nbhd_idx)
        union_case = int(is_case[in_union].sum())
        union_ctrl = int((1 - is_case[in_union]).sum())
        union_contingency = [
            [union_case, union_ctrl],
            [n_case - union_case, n_ctrl - union_ctrl],
        ]
        union_or, union_p = fisher_exact(union_contingency, alternative="greater")
        union_stats = {
            "n_significant_centers": len(sig_center_to_nbhd_idx),
            "significant_centers": sorted(sig_center_to_nbhd_idx.keys()),
            "union_n_residues": int(len(union_nbhd_idx)),
            "union_total_resolved_residues": len(ids),
            "union_residue_pct": len(union_nbhd_idx) / len(ids) if ids else float("nan"),
            "union_case": union_case,
            "total_case": n_case,
            "union_case_pct": union_case / n_case if n_case > 0 else float("nan"),
            "union_ctrl": union_ctrl,
            "total_ctrl": n_ctrl,
            "union_ctrl_pct": union_ctrl / n_ctrl if n_ctrl > 0 else float("nan"),
            "union_odds_ratio": union_or,
            "union_p_value": union_p,
        }
    else:
        union_stats = {
            "n_significant_centers": 0,
            "significant_centers": [],
            "union_n_residues": 0,
            "union_total_resolved_residues": len(ids),
            "union_residue_pct": 0.0,
            "union_case": 0,
            "total_case": n_case,
            "union_case_pct": 0.0 if n_case > 0 else float("nan"),
            "union_ctrl": 0,
            "total_ctrl": n_ctrl,
            "union_ctrl_pct": 0.0 if n_ctrl > 0 else float("nan"),
            "union_odds_ratio": float("nan"),
            "union_p_value": float("nan"),
        }

    union_stats_df = pd.DataFrame([union_stats])

    variants_with_p = variants.copy()
    if not per_residue.empty:
        variants_with_p["aa_pos"] = pd.to_numeric(
            variants_with_p["aa_pos"], errors="coerce"
        ).astype("Int64")
        to_merge = per_residue.drop(columns=["center_idx"]).copy()
        to_merge["aa_pos"] = to_merge["aa_pos"].astype("Int64")
        variants_with_p = variants_with_p.merge(to_merge, on="aa_pos", how="left")
    else:
        variants_with_p["p_value"] = float("nan")

    return variants_with_p, per_residue, union_stats_df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _default_repo_root() -> Path:
    # python_code/<script>  -> repo root is the parent of python_code/
    return Path(__file__).resolve().parent.parent


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = _default_repo_root()
    parser = argparse.ArgumentParser(
        description=("Run the 3D neighborhood test on the human PMCA2z/a "
                     "cryo-EM structure."),
    )
    parser.add_argument(
        "--case",
        type=Path,
        default=root / "data" / "atp2b2_case_variants_schema_asd.tsv",
        help="TSV of case variants (canonical ATP2B2 numbering).",
    )
    parser.add_argument(
        "--control",
        type=Path,
        default=root / "data" / "atp2b2_control_variants_schema_asd.tsv",
        help="TSV of control variants (canonical ATP2B2 numbering).",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        default=root / "data" / "pmca_e1ca_model.pdb",
        help="Cryo-EM PDB file (human PMCA2z/a, E1-Ca state).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "results",
        help="Output directory for the mapping CSV and 3DNT result tables.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=15.0,
        help="Neighborhood radius in Angstroms (default: 15).",
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="A",
        help="Chain ID to scan (default: A).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/2] Loading and remapping variants ...")
    variants = load_schema_asd_variants(args.case, args.control)
    mapping = remap_variants_to_pdb(
        variants,
        canonical_seq=CANONICAL_SEQ,
        isoform_seq=ISOFORM_SEQ,
        pdb_path=args.pdb,
        chain_id=args.chain,
    )
    mapping_path = args.out_dir / "mutation_mapping_PMCA2za.csv"
    mapping.to_csv(mapping_path, index=False)
    print(f"  -> {mapping_path}")

    # Keep only successfully mapped variants for the 3DNT scan.
    mapped = mapping.loc[mapping["status"] == "mapped"].copy()
    mapped = mapped.rename(columns={"new_resi": "aa_pos"})
    mapped["aa_pos"] = mapped["aa_pos"].astype(int)
    mapped["is_case"] = mapped["is_case"].astype(int)

    print(f"\n[2/2] Running 3D neighborhood test (radius = {args.radius} A) ...")
    _, per_residue, union_stats = neighborhood_test(
        variants=mapped,
        pdb_path=args.pdb,
        radius=args.radius,
        chain=args.chain,
    )

    per_residue_path = args.out_dir / "3dnt_cryoem_per_residue.csv"
    union_stats_path = args.out_dir / "3dnt_cryoem_union_stats.csv"
    per_residue.to_csv(per_residue_path, index=False)
    union_stats.to_csv(union_stats_path, index=False)
    print(f"  -> {per_residue_path}")
    print(f"  -> {union_stats_path}")

    print("\nTop 10 residues by p-value:")
    cols = ["aa_pos", "nbhd_case", "nbhd_ctrl", "case_pct", "ctrl_pct",
            "odds_ratio", "p_value", "bonferroni_significant"]
    print(per_residue[cols].head(10).to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
