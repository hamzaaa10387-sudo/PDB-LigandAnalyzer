"""Compute drug-likeness descriptors for a ligand and apply Lipinski's rule of five."""
from pathlib import Path
from typing import Dict

from rdkit import Chem
from rdkit.Chem import Descriptors


def load_ligand(pdb_path: str | Path) -> Chem.Mol:
    """Load a ligand from a PDB file.

    Bonds are inferred from atom geometry (``proximityBonding``) because HETATM
    records often lack CONECT bond information. If RDKit still cannot build a
    valid molecule, a ValueError is raised: descriptors are never estimated from
    an incompletely sanitized structure.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"Ligand file not found: {pdb_path}")

    mol = Chem.MolFromPDBFile(
        str(pdb_path),
        sanitize=True,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None:
        raise ValueError(
            f"RDKit could not build a valid molecule from {pdb_path.name}. "
            "The PDB entry may lack suitable bond information. Use a chemically "
            "annotated ligand structure (for example, an SDF from the PDB "
            "chemical component dictionary) instead of estimating descriptors."
        )
    return mol


def evaluate_ligand_properties(pdb_path: str | Path, verbose: bool = True) -> Dict:
    """Return descriptors and Lipinski status for the ligand in `pdb_path`."""
    if verbose:
        print(f"[3/3] Calculating ligand properties with RDKit ({Path(pdb_path).name})...")

    mol = load_ligand(pdb_path)

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # Lipinski's rule of five: each limit is violated when exceeded.
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    compliant = violations <= 1

    if verbose:
        print("\n--- Physicochemical descriptors (Lipinski rule of five) ---")
        print(f"  Molecular weight : {mw:.2f} Da   (limit <= 500)")
        print(f"  LogP             : {logp:.2f}      (limit <= 5)")
        print(f"  H-bond donors    : {hbd}          (limit <= 5)")
        print(f"  H-bond acceptors : {hba}          (limit <= 10)")
        status = "PASSED" if compliant else "FAILED"
        print(f"  Status           : {status} ({violations} violation(s))")
        print(
            "  Note: this is a screening heuristic, not evidence of efficacy, "
            "safety or bioavailability."
        )

    return {
        "mw": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "violations": violations,
        "lipinski_compliant": compliant,
    }
