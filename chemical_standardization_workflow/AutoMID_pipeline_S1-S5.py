#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AutoMID pipeline (S1-S5)

Runs SMILES through:
  S1 Canonicalize
  S2 Desalt (keep largest fragment)
  S3 Neutralize
  S4 Canonical tautomer
  S5 Strip isotopes and explicit hydrogens

Local execution:
  1. Place this script and smiles_input.csv in the same folder.
  2. Open Terminal in that folder.
  3. Run:

       caffeinate python3 AutoMID_pipeline_S1-S5.py

Output:
  standardized_smiles_output.xlsx
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


# ---------------------------------------------------------------------
# Local file settings
# ---------------------------------------------------------------------

INPUT_FILE = "smiles_input.csv"
OUTPUT_FILE = "standardized_smiles_output.xlsx"

SMILES_COLUMN = "smiles"
ID_COLUMN = "odo_id"

STAGES = ["S1", "S2", "S3", "S4", "S5"]
ERRORS_ONLY = False


# ---------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------

@dataclass
class Structure:
    odo_id: str
    S0: str
    mol: Optional[Chem.Mol]
    S1: Optional[str] = None
    S2: Optional[str] = None
    S3: Optional[str] = None
    S4: Optional[str] = None
    S5: Optional[str] = None
    error: Optional[str] = None


# ---------------------------------------------------------------------
# Chemistry helpers
# ---------------------------------------------------------------------

def neutralize_atoms(mol: Optional[Chem.Mol]) -> Optional[Chem.Mol]:
    """Neutralize common charged centers where chemically appropriate."""
    if mol is None:
        return None

    pattern = Chem.MolFromSmarts(
        "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
    )

    if pattern is None:
        return mol

    for (idx,) in mol.GetSubstructMatches(pattern):
        atom = mol.GetAtomWithIdx(idx)
        charge = atom.GetFormalCharge()
        hydrogen_count = atom.GetTotalNumHs()

        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(max(hydrogen_count - charge, 0))
        atom.UpdatePropertyCache()

    return mol


class StructurePipeline:
    """Apply the S1-S5 structure-standardization stages."""

    def __init__(self) -> None:
        self.remover = rdMolStandardize.FragmentRemover()
        self.tautomerizer = rdMolStandardize.TautomerEnumerator()
        self.largest_fragment = rdMolStandardize.LargestFragmentChooser()

    def process(self, odo_id: str, smiles: str, stages: List[str]) -> Structure:
        result = Structure(odo_id=odo_id, S0=smiles, mol=None)

        try:
            mol = Chem.MolFromSmiles(str(smiles))

            if mol is None:
                raise ValueError("RDKit failed to parse SMILES")

            result.mol = mol

            # S1 - Canonicalize
            if "S1" in stages:
                result.S1 = Chem.MolToSmiles(result.mol)
                result.mol = Chem.MolFromSmiles(result.S1)

            # S2 - Desalt; if multiple fragments remain, keep the largest
            if "S2" in stages and result.mol is not None:
                desalted = self.remover.remove(result.mol)

                if "." in Chem.MolToSmiles(desalted):
                    desalted = self.largest_fragment.choose(desalted)

                result.mol = desalted
                result.S2 = Chem.MolToSmiles(result.mol)

            # S3 - Neutralize
            if "S3" in stages and result.mol is not None:
                neutralized = neutralize_atoms(result.mol)

                if neutralized is None:
                    raise ValueError("Neutralization failed")

                result.mol = neutralized
                result.S3 = Chem.MolToSmiles(result.mol)

            # S4 - Canonical tautomer
            if "S4" in stages and result.mol is not None:
                tautomer = self.tautomerizer.Canonicalize(result.mol)

                if tautomer is None:
                    raise ValueError("Tautomerization failed")

                result.mol = tautomer
                result.S4 = Chem.MolToSmiles(result.mol)

            # S5 - Remove isotopes and explicit hydrogen atoms
            if "S5" in stages and result.mol is not None:
                for atom in result.mol.GetAtoms():
                    if atom.GetIsotope():
                        atom.SetIsotope(0)

                editable_mol = Chem.RWMol(result.mol)
                hydrogen_indices = [
                    atom.GetIdx()
                    for atom in editable_mol.GetAtoms()
                    if atom.GetAtomicNum() == 1
                ]

                for idx in sorted(hydrogen_indices, reverse=True):
                    editable_mol.RemoveAtom(idx)

                result.mol = editable_mol.GetMol()
                Chem.SanitizeMol(result.mol)
                result.S5 = Chem.MolToSmiles(result.mol)

        except Exception as exc:
            result.error = str(exc)

        return result


# ---------------------------------------------------------------------
# Input/output workflow
# ---------------------------------------------------------------------

def run_pipeline() -> None:
    """Read the local CSV file, standardize structures, and write Excel output."""

    input_path = Path(INPUT_FILE)
    output_path = Path(OUTPUT_FILE)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_path.resolve()}\n"
            f"Place '{INPUT_FILE}' in the same folder as this script."
        )

    dataframe = pd.read_csv(input_path)

    if SMILES_COLUMN not in dataframe.columns:
        raise ValueError(
            f"Missing required SMILES column: '{SMILES_COLUMN}'. "
            f"Available columns: {list(dataframe.columns)}"
        )

    if ID_COLUMN not in dataframe.columns:
        raise ValueError(
            f"Missing required ID column: '{ID_COLUMN}'. "
            f"Available columns: {list(dataframe.columns)}"
        )

    pipeline = StructurePipeline()
    rows = []

    for _, row in dataframe.iterrows():
        odo_id = str(row[ID_COLUMN])
        smiles = str(row[SMILES_COLUMN])

        result = pipeline.process(
            odo_id=odo_id,
            smiles=smiles,
            stages=STAGES,
        )

        rows.append(
            {
                ID_COLUMN: result.odo_id,
                "S0": result.S0,
                "S1": result.S1,
                "S2": result.S2,
                "S3": result.S3,
                "S4": result.S4,
                "S5": result.S5,
                "error": result.error,
            }
        )

    output_dataframe = pd.DataFrame(rows)

    if ERRORS_ONLY:
        output_dataframe = output_dataframe[
            output_dataframe["error"].notna()
            & (output_dataframe["error"] != "")
        ]

    output_dataframe.to_excel(
        output_path,
        index=False,
        engine="openpyxl",
        sheet_name="standardized_structures",
    )

    total_rows = len(output_dataframe)
    error_rows = int(
        (
            output_dataframe["error"].notna()
            & (output_dataframe["error"] != "")
        ).sum()
    )

    print(f"Input file:  {input_path.resolve()}")
    print(f"Output file: {output_path.resolve()}")
    print(f"Rows written: {total_rows}")
    print(f"Rows with errors: {error_rows}")


if __name__ == "__main__":
    run_pipeline()
