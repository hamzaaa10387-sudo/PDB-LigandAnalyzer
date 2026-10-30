"""PDB-LigandAnalyzer: extract a ligand's binding pocket and screen the ligand.

Usage:
    python main.py                       # default: 1IEP / STI, 5.0 A pocket
    python main.py 1IEP STI              # custom PDB ID and ligand residue name
    python main.py 1IEP STI 6.5          # also set the pocket radius (Angstrom)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from chemo_filter import evaluate_ligand_properties  # noqa: E402
from fetch_pdb import download_pdb_structure  # noqa: E402
from pocket_analyzer import analyze_binding_pocket  # noqa: E402

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
DEFAULT_PDB_ID = "1IEP"
DEFAULT_LIGAND_ID = "STI"
DEFAULT_RADIUS = 5.0


def run_pipeline(pdb_id: str = DEFAULT_PDB_ID,
                 ligand_id: str = DEFAULT_LIGAND_ID,
                 radius: float = DEFAULT_RADIUS) -> int:
    print("=" * 69)
    print("  PDB-LigandAnalyzer: structure-based ligand screening pipeline")
    print("=" * 69)
    print(f"  Target structure : {pdb_id}")
    print(f"  Target ligand    : {ligand_id}")
    print(f"  Pocket radius    : {radius:.1f} A\n")

    try:
        pdb_path = download_pdb_structure(pdb_id, output_dir=DATA_DIR)
        pocket, paths = analyze_binding_pocket(
            pdb_path, ligand_id=ligand_id, radius=radius, output_dir=DATA_DIR
        )
        print("\nBinding-pocket residues:")
        print("  " + ", ".join(pocket))

        evaluate_ligand_properties(paths["ligand"])

        print("\nGenerated files:")
        for label, path in paths.items():
            print(f"  {label:10s}: {path}")
        print("=" * 69)
        return 0
    except Exception as exc:  # report clearly and exit non-zero
        print(f"\n[-] Pipeline failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    args = sys.argv[1:]
    pdb = args[0] if len(args) >= 1 else DEFAULT_PDB_ID
    lig = args[1] if len(args) >= 2 else DEFAULT_LIGAND_ID
    try:
        rad = float(args[2]) if len(args) >= 3 else DEFAULT_RADIUS
    except ValueError:
        print("[-] Radius must be a number, e.g. 5.0", file=sys.stderr)
        sys.exit(1)
    sys.exit(run_pipeline(pdb, lig, rad))
