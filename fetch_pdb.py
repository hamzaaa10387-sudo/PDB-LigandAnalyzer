"""Download PDB structures from the RCSB archive."""
import os

from Bio.PDB import PDBList


def download_pdb_structure(pdb_id: str, output_dir: str = "data") -> str:
    """Download a PDB entry and return the path to the local .pdb file."""
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError(f"'{pdb_id}' is not a valid 4-character PDB ID.")

    os.makedirs(output_dir, exist_ok=True)
    print(f"[+] Downloading PDB entry {pdb_id}...")

    retriever = PDBList(verbose=False)
    downloaded = retriever.retrieve_pdb_file(
        pdb_id, pdir=output_dir, file_format="pdb"
    )
    if not downloaded or not os.path.exists(downloaded):
        raise FileNotFoundError(f"Could not download PDB entry {pdb_id}.")

    target = os.path.join(output_dir, f"{pdb_id.lower()}.pdb")
    if os.path.abspath(downloaded) != os.path.abspath(target):
        os.replace(downloaded, target)

    print(f"[✓] Saved structure to {target}")
    return target
