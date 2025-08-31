#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoMID pipeline (S1–S5)
Runs SMILES through:
  S1 Canonicalize
  S2 Desalt (keep largest fragment)
  S3 Neutralize
  S4 Canonical tautomer
  S5 Strip isotopes & explicit hydrogens

Examples:
  python pipelines/AutoMID_pipeline_S1-S5.py \
    --in data/example_input.csv \
    --out data/example_output.csv \
    --smiles-col smiles \
    --id-col odo_id

  # Using a YAML config for defaults (merged with CLI flags)
  python pipelines/AutoMID_pipeline_S1-S5.py \
    --in data/example_input.csv \
    --out data/example_output.csv \
    --config config/config.yaml
"""

import sys
import argparse
import pathlib
from dataclasses import dataclass
from typing import Optional, List, Dict, Any

import yaml
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


# ---------- Data container ----------
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


# ---------- Chemistry helpers ----------
def neutralize_atoms(mol: Optional[Chem.Mol]) -> Optional[Chem.Mol]:
    """Neutralize common charged centers."""
    if mol is None:
        return None
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    if pattern is None:
        return mol
    for (idx,) in mol.GetSubstructMatches(pattern):
        atom = mol.GetAtomWithIdx(idx)
        charge = atom.GetFormalCharge()
        hcount = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(max(hcount - charge, 0))
        atom.UpdatePropertyCache()
    return mol


class StructurePipeline:
    def __init__(self) -> None:
        self.remover = rdMolStandardize.FragmentRemover()
        self.tautomerizer = rdMolStandardize.TautomerEnumerator()
        self.largestFragment = rdMolStandardize.LargestFragmentChooser()

    def process(self, odo_id: str, smiles: str, stages: List[str]) -> Structure:
        s = Structure(odo_id=odo_id, S0=smiles, mol=None)
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                raise ValueError("RDKit failed to parse SMILES")
            s.mol = mol

            # S1 - Canonicalize
            if "S1" in stages:
                s.S1 = Chem.MolToSmiles(s.mol)
                s.mol = Chem.MolFromSmiles(s.S1)

            # S2 - Desalt; if multiple frags remain, keep largest
            if "S2" in stages and s.mol is not None:
                m = self.remover.remove(s.mol)
                if "." in Chem.MolToSmiles(m):
                    m = self.largestFragment.choose(s.mol)
                s.mol = m
                s.S2 = Chem.MolToSmiles(s.mol)

            # S3 - Neutralize
            if "S3" in stages and s.mol is not None:
                m = neutralize_atoms(s.mol)
                if m is None:
                    s.error = "Neutralization failed"
                else:
                    s.mol = m
                    s.S3 = Chem.MolToSmiles(s.mol)

            # S4 - Canonical tautomer
            if "S4" in stages and s.mol is not None:
                m = self.tautomerizer.Canonicalize(s.mol)
                if m is None:
                    s.error = "Tautomerization failed"
                else:
                    s.mol = m
                    s.S4 = Chem.MolToSmiles(s.mol)

            # S5 - Remove isotopes & explicit H atoms
            if "S5" in stages and s.mol is not None:
                for atom in s.mol.GetAtoms():
                    if atom.GetIsotope():
                        atom.SetIsotope(0)
                rw = Chem.RWMol(s.mol)
                h_indices = [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1]
                for idx in sorted(h_indices, reverse=True):
                    rw.RemoveAtom(idx)
                s.mol = rw.GetMol()
                s.S5 = Chem.MolToSmiles(s.mol)

        except Exception as e:
            s.error = str(e)
        return s


# ---------- I/O helpers ----------
def run_pipeline(inp: str, outp: str, cfg: Dict[str, Any]) -> None:
    smiles_col = cfg.get("smiles_col", "smiles")
    id_col = cfg.get("id_col", "odo_id")
    stages = [s.strip().upper() for s in cfg.get("stages", ["S1", "S2", "S3", "S4", "S5"]) if s.strip()]
    errors_only = bool(cfg.get("errors_only", False))

    df = pd.read_csv(inp)
    if smiles_col not in df.columns:
        raise SystemExit(f"Missing SMILES column: {smiles_col}")
    if id_col not in df.columns:
        raise SystemExit(f"Missing ID column: {id_col}")

    pipe = StructurePipeline()
    rows = []
    for _, row in df.iterrows():
        odo_id = str(row[id_col])
        smi = str(row[smiles_col])
        res = pipe.process(odo_id, smi, stages)
        rows.append({
            id_col: res.odo_id,
            "S0": res.S0,
            "S1": res.S1,
            "S2": res.S2,
            "S3": res.S3,
            "S4": res.S4,
            "S5": res.S5,
            "error": res.error
        })

    out = pd.DataFrame(rows)
    if errors_only:
        out = out[out["error"].notna() & (out["error"] != "")]
    pathlib.Path(outp).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outp, index=False)


# ---------- CLI ----------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="AutoMID S1–S5 pipeline for SMILES normalization (RDKit)."
    )
    p.add_argument("--in", dest="inp", required=True, help="Input CSV")
    p.add_argument("--out", dest="outp", required=True, help="Output CSV")
    p.add_argument("--config", dest="config", default="config/config.yaml",
                   help="Optional YAML config file (overridden by CLI flags)")
    p.add_argument("--smiles-col", default=None, help="SMILES column name")
    p.add_argument("--id-col", default=None, help="ID column name")
    p.add_argument("--stages", default=None,
                   help="Comma-separated subset of S1,S2,S3,S4,S5 (default: all)")
    p.add_argument("--errors-only", action="store_true",
                   help="Only write rows that had an error")
    return p.parse_args(argv)


def load_config(path: str) -> Dict[str, Any]:
    p = pathlib.Path(path)
    if p.exists():
        return yaml.safe_load(p.read_text()) or {}
    return {}


def main(argv=None) -> int:
    args = parse_args(argv)
    # Load config then merge CLI overrides
    cfg = load_config(args.config)
    if args.smiles_col is not None:
        cfg["smiles_col"] = args.smiles_col
    if args.id_col is not None:
        cfg["id_col"] = args.id_col
    if args.stages is not None:
        cfg["stages"] = [s.strip() for s in args.stages.split(",")]
    if args.errors_only:
        cfg["errors_only"] = True

    run_pipeline(args.inp, args.outp, cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
