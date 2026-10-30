# PDB-LigandAnalyzer

A small, reproducible structure-based pipeline for introductory structural
bioinformatics and cheminformatics.

1. **Download** a structure from the RCSB Protein Data Bank (`src/fetch_pdb.py`).
2. **Extract** the bound ligand, the protein residues within a chosen distance
   of it, and a cleaned receptor file (`src/pocket_analyzer.py`).
3. **Evaluate** the ligand's MW, LogP, H-bond donors and acceptors against
   Lipinski's rule of five with RDKit (`src/chemo_filter.py`).

Default target: PDB entry `1IEP` (ABL kinase) with ligand `STI` (imatinib),
pocket radius 5.0 Å.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Requires Python 3.10+ and an internet connection for the download step.

## Run

```bash
python main.py                     # 1IEP / STI, 5.0 Å
python main.py 1IEP STI            # explicit PDB ID and ligand residue name
python main.py 1IEP STI 6.5        # also set the pocket radius
```

Outputs go to `data/`:

| File | Contents |
| --- | --- |
| `<pdb>.pdb` | downloaded structure |
| `<pdb>_ligand_<LIG>.pdb` | the ligand |
| `<pdb>_pocket.pdb` | protein residues within the radius |
| `<pdb>_receptor_clean.pdb` | standard protein residues (no waters, ligands or ions) |

## Tests

The tests are offline and do not need the RCSB download:

```bash
python tests/test_offline.py       # or: python -m pytest tests
```

They check that:
- pocket residues are selected per chain (a same-numbered residue on another
  chain is not mistaken for the pocket residue) and waters are ignored;
- a missing ligand raises an error;
- aspirin (MW ≈ 180.16) gives the expected descriptors and passes Lipinski;
- violation counts agree with the rule-of-five thresholds.

## How the pocket is defined

A protein residue is in the pocket if any of its atoms lies within the chosen
radius of any ligand atom. Residues are identified by chain and residue number,
so multi-chain entries (homodimers and similar) are handled correctly.

## Limitations

- This is a geometric distance cutoff, not a binding-site prediction method.
- The pipeline does not dock molecules, estimate binding affinity, or predict
  biological activity.
- Descriptors come from the ligand as stored in the PDB file. If RDKit cannot
  build a valid molecule, the program stops with an error rather than guessing.
  Use a chemically annotated ligand (for example, from the PDB chemical
  component dictionary) in that case.
- Lipinski's rule of five is a screening heuristic. Passing it does not establish
  drug-likeness, efficacy, safety or bioavailability.

## Data source

Structures come from the RCSB Protein Data Bank: https://www.rcsb.org/
Follow the database's citation and data-use guidance in research or publications.
