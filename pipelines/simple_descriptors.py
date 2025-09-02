#!/usr/bin/env python3
"""
simple_descriptors.py

Compute basic molecular descriptors (InChI, InChIKey, MW, MF) from a CSV
produced by the AutoMID S1–S5 pipeline.

Typical use:
  python pipelines/simple_descriptors.py \
    --in data/example_output.csv \
    --out data/example_descriptors.csv \
    --smiles-col S5 \
    --id-col odo_id

Notes
- Defaults assume you're feeding the pipeline's output (S5 column).
- You can point --smiles-col to S4/S1 or a raw 'smiles' column if needed.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def compute_row(smiles: str):
    """Return (inchi, inchikey, mw, formula) or (None, None, None, None) on failure."""
    try:
        mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None
        if mol is None:
            return (None, None, None, None)
        # InChI / InChIKey (requires RDKit built with InChI support; most conda builds do)
        try:
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
        except Exception:
            # If InChI support is missing, leave these as None but still compute MW/MF
            inchi, inchikey = (None, None)
        mw = Descriptors.MolWt(mol)
        mf = rdMolDescriptors.CalcMolFormula(mol)
        return (inchi, inchikey, mw, mf)
    except Exception:
        return (None, None, None, None)

def main():
    ap = argparse.ArgumentParser(description="Compute InChI, InChIKey, MW, MF from SMILES.")
    ap.add_argument("--in", dest="inp", required=True, help="Input CSV from pipeline (e.g., data/example_output.csv)")
    ap.add_argument("--out", dest="out", required=True, help="Output CSV to write (e.g., data/example_descriptors.csv)")
    ap.add_argument("--smiles-col", dest="scol", default="S5", help="Column containing SMILES (default: S5)")
    ap.add_argument("--id-col", dest="idcol", default="odo_id", help="Identifier column to carry through (default: odo_id)")
    args = ap.parse_args()

    in_path = Path(args.inp)
    if not in_path.exists():
        sys.exit(f"[error] Input file not found: {in_path}")

    try:
        df = pd.read_csv(in_path)
    except Exception as e:
        sys.exit(f"[error] Could not read {in_path}: {e}")

    if args.scol not in df.columns:
        sys.exit(f"[error] SMILES column '{args.scol}' not found. Available columns: {list(df.columns)}")

    # Build result frame
    out_cols = [args.idcol] if args.idcol in df.columns else []
    out = pd.DataFrame()
    if out_cols:
        out[args.idcol] = df[args.idcol]
    out["smiles_used"] = df[args.scol]

    # Compute descriptors row-wise
    desc = df[args.scol].apply(compute_row)
    out["InChI"] = desc.apply(lambda t: t[0])
    out["InChIKey"] = desc.apply(lambda t: t[1])
    out["MW"] = desc.apply(lambda t: t[2])
    out["MF"] = desc.apply(lambda t: t[3])

    # Basic quality flags
    out["parse_ok"] = out["InChI"].notna() | out["MW"].notna() | out["MF"].notna()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"[ok] Wrote descriptors to {out_path}")

if __name__ == "__main__":
    sys.exit(main())

