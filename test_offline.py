"""Offline tests: no network access needed.

Run with:  python tests/test_offline.py     (or: python -m pytest tests)
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402

from chemo_filter import evaluate_ligand_properties  # noqa: E402
from pocket_analyzer import analyze_binding_pocket  # noqa: E402

# Synthetic structure:
#   chain A: ALA 1 close to the ligand (within 5 A)
#   chain B: ALA 1 far away (same residue number, different chain)
#   chain C: water next to the ligand (must be ignored)
#   chain A: GLY 2 far away
TOY_PDB = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      30.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA B   1      30.000   0.000   0.000  1.00  0.00           C
HETATM    4  C1  STI A 500       2.000   0.000   0.000  1.00  0.00           C
HETATM    5  O1  STI A 500       3.000   0.000   0.000  1.00  0.00           O
HETATM    6  O   HOH C 600       1.000   1.000   1.000  1.00  0.00           O
END
"""


def _write(tmp, name, text):
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def test_pocket_is_chain_aware_and_excludes_water():
    """Only chain A's ALA1 is near the ligand; chain B's ALA1 must not be picked up."""
    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = _write(tmp, "toy.pdb", TOY_PDB)
        labels, paths = analyze_binding_pocket(pdb_path, "STI", radius=5.0, output_dir=tmp)
        assert labels == ["A:ALA1"], labels
        for path in paths.values():
            assert os.path.exists(path), path


def test_missing_ligand_raises():
    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = _write(tmp, "toy.pdb", TOY_PDB)
        try:
            analyze_binding_pocket(pdb_path, "XYZ", output_dir=tmp)
        except ValueError:
            return
        raise AssertionError("expected ValueError for a missing ligand")


def _ligand_pdb_from_smiles(tmp, name, smiles):
    """Write a 2D-coordinate PDB (with CONECT records) for a known molecule."""
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    return _write(tmp, name, Chem.MolToPDBBlock(mol))


def test_aspirin_descriptors_and_pass():
    """Aspirin (MW ~180.16) is well inside Lipinski limits."""
    with tempfile.TemporaryDirectory() as tmp:
        path = _ligand_pdb_from_smiles(tmp, "aspirin.pdb", "CC(=O)Oc1ccccc1C(=O)O")
        result = evaluate_ligand_properties(path, verbose=False)
        assert abs(result["mw"] - 180.16) < 0.1, result["mw"]
        assert result["violations"] == 0
        assert result["lipinski_compliant"] is True


def test_violation_count_matches_thresholds():
    """A molecule breaking several rules must be counted and marked as failing."""
    # Large, many-hydroxyl molecule: MW > 500 and HBD > 5.
    smiles = "OCC(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)C(O)CO"
    with tempfile.TemporaryDirectory() as tmp:
        path = _ligand_pdb_from_smiles(tmp, "polyol.pdb", smiles)
        result = evaluate_ligand_properties(path, verbose=False)
        mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=False, proximityBonding=True)
        assert abs(result["mw"] - Descriptors.MolWt(mol)) < 1e-6
        expected = sum([
            result["mw"] > 500,
            result["logp"] > 5,
            result["hbd"] > 5,
            result["hba"] > 10,
        ])
        assert result["violations"] == expected
        assert result["lipinski_compliant"] == (expected <= 1)


if __name__ == "__main__":
    test_pocket_is_chain_aware_and_excludes_water()
    test_missing_ligand_raises()
    test_aspirin_descriptors_and_pass()
    test_violation_count_matches_thresholds()
    print("All offline tests passed.")
