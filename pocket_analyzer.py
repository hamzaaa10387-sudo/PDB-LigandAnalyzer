"""Locate a bound ligand and extract the surrounding binding pocket."""
import os
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Residue import Residue


class _ResidueSelect(Select):
    """Keep only the residues whose (chain, residue-id) is in `keep`."""

    def __init__(self, keep):
        self.keep = keep

    def accept_residue(self, residue):
        return _key(residue) in self.keep


def _key(residue: Residue) -> Tuple[str, Tuple]:
    chain_id = residue.get_parent().id
    return (chain_id, residue.id)


def analyze_binding_pocket(
    pdb_file: str,
    ligand_id: str,
    radius: float = 5.0,
    output_dir: str = "data",
) -> Tuple[List[str], Dict[str, str]]:
    """Extract the ligand, pocket residues and cleaned receptor.

    Returns a sorted list of pocket residue labels (e.g. "A:LYS271")
    and a dict of output file paths.
    """
    if ligand_id is None or not ligand_id.strip():
        raise ValueError("A ligand residue name is required.")
    ligand_id = ligand_id.strip().upper()

    os.makedirs(output_dir, exist_ok=True)
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]

    ligand_residues = []
    protein_residues = []
    for chain in model:
        for residue in chain:
            name = residue.get_resname().strip()
            if residue.id[0] != " " and name == ligand_id:
                ligand_residues.append(residue)
            elif residue.id[0] == " ":
                # Standard amino-acid / polymer residues only (no waters, no HETATMs).
                protein_residues.append(residue)

    if not ligand_residues:
        raise ValueError(f"Ligand '{ligand_id}' was not found in {pdb_file}.")

    ligand_coords = np.array(
        [atom.get_coord() for res in ligand_residues for atom in res.get_atoms()]
    )
    centroid = ligand_coords.mean(axis=0)
    print(f"[+] Ligand centroid: {np.round(centroid, 2)}")

    pocket_keys = set()
    for residue in protein_residues:
        coords = np.array([atom.get_coord() for atom in residue.get_atoms()])
        if coords.size == 0:
            continue
        # Minimum distance from any atom of the residue to the ligand's atoms.
        dists = np.linalg.norm(
            coords[:, None, :] - ligand_coords[None, :, :], axis=2
        )
        if dists.min() <= radius:
            pocket_keys.add(_key(residue))

    pocket_labels = sorted(
        f"{chain}:{res.get_resname()}{res.id[1]}"
        for chain, res_id in pocket_keys
        for res in [model[chain][res_id]]
    )
    print(f"[✓] {len(pocket_labels)} residues within {radius} Å of the ligand.")

    ligand_keys = {_key(r) for r in ligand_residues}
    paths = {
        "ligand": os.path.join(output_dir, f"{pdb_id.lower()}_ligand_{ligand_id}.pdb"),
        "pocket": os.path.join(output_dir, f"{pdb_id.lower()}_pocket.pdb"),
        "receptor": os.path.join(output_dir, f"{pdb_id.lower()}_receptor_clean.pdb"),
    }

    io = PDBIO()
    io.set_structure(structure)
    io.save(paths["ligand"], _ResidueSelect(ligand_keys))
    io.save(paths["pocket"], _ResidueSelect(pocket_keys))
    protein_keys = {_key(r) for r in protein_residues}
    io.save(paths["receptor"], _ResidueSelect(protein_keys))

    print(f"[✓] Wrote ligand, pocket and receptor files to {output_dir}/")
    return pocket_labels, paths
