#!/usr/bin/env python3
"""
Combined SMILES annotation pipeline for CSV files.

Processing order
----------------
1. Calculate RDKit molecular weight from input SMILES.
2. Calculate RDKit molecular formula from input SMILES.
3. Calculate Standard InChI and Standard InChIKey from input SMILES.
4. Retrieve PubChem CID from input SMILES.
5. Retrieve IUPAC name from the PubChem CID.

Everything except the IUPAC name is derived from `canonical_smiles`.
The IUPAC name is retrieved through the PubChem CID.

Output columns (all other input columns are dropped)
----------------------------------------------------
    canonical_smiles
    standard_inchi
    standard_inchikey
    molecular_formula
    molecular_weight
    pubchem_cid
    iupac_name
    rdkit_status
    inchi_status
    pubchem_cid_status
    iupac_status
    annotation_status

Run
---
Place the input CSV in the same directory as this script, then run:

    python3 compound_annotation_pipeline.py

Edit INPUT_FILE, OUTPUT_FILE, and SMILES_COLUMN in the settings section when needed.

Dependencies
------------
pandas, rdkit, pubchempy, tqdm
"""

from __future__ import annotations

import time
import urllib.error
from http.client import IncompleteRead
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import pubchempy as pcp
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from tqdm import tqdm


# -----------------------------------------------------------------------------
# User-editable settings
# -----------------------------------------------------------------------------
INPUT_FILE = Path("smiles_imput_file.csv")
OUTPUT_FILE = Path("compound_processing_output.csv")
SMILES_COLUMN = "canonical_smiles"

MW_DECIMALS = 3

CHECKPOINT_EVERY = 50
REQUEST_DELAY = 0.2
MAX_RETRIES = 3
RETRY_DELAY = 1.0

# Silence RDKit InChI console chatter; status columns carry the audit trail.
RDLogger.DisableLog("rdApp.*")


DEFAULT_SMILES_COLUMNS = (
    "canonical_smiles",
    "chembl_canonical_smiles",
    "SMILES",
    "smiles",
    "rdkit_parent_smiles",
)

# Final column order written to the output and to every checkpoint.
OUTPUT_COLUMNS = (
    "canonical_smiles",
    "standard_inchi",
    "standard_inchikey",
    "molecular_formula",
    "molecular_weight",
    "pubchem_cid",
    "iupac_name",
    "rdkit_status",
    "inchi_status",
    "pubchem_cid_status",
    "iupac_status",
    "annotation_status",
)


def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Remove surrounding whitespace and common invisible characters."""
    df = df.copy()
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.replace("\ufeff", "", regex=False)
        .str.replace("\u200b", "", regex=False)
        .str.replace("\xa0", "", regex=False)
    )
    return df


def resolve_smiles_column(df: pd.DataFrame, requested: Optional[str]) -> str:
    """Use the requested SMILES column or detect a common alternative."""
    if requested:
        cleaned_requested = (
            requested.strip()
            .replace("\ufeff", "")
            .replace("\u200b", "")
            .replace("\xa0", "")
        )
        if cleaned_requested in df.columns:
            return cleaned_requested
        raise ValueError(
            f"SMILES column '{requested}' was not found. "
            f"Available columns: {df.columns.tolist()}"
        )

    for candidate in DEFAULT_SMILES_COLUMNS:
        if candidate in df.columns:
            return candidate

    raise ValueError(
        "No recognized SMILES column was found. Update SMILES_COLUMN in the settings section. "
        f"Available columns: {df.columns.tolist()}"
    )


def normalize_smiles_value(value: Any) -> Optional[str]:
    """Return a stripped SMILES string or None for missing/blank values."""
    if pd.isna(value):
        return None
    smiles = str(value).strip()
    return smiles if smiles else None


MISSING_RDKIT_RESULT: tuple[
    Optional[float], Optional[str], Optional[str], Optional[str], str, str
] = (None, None, None, None, "missing_smiles", "missing_smiles")


def calculate_rdkit_properties(
    smiles: Optional[str],
) -> tuple[Optional[float], Optional[str], Optional[str], Optional[str], str, str]:
    """Calculate MW, formula, Standard InChI and Standard InChIKey from one parse.

    Returns
    -------
    (molecular_weight, molecular_formula, standard_inchi, standard_inchikey,
     rdkit_status, inchi_status)
    """
    if smiles is None:
        return MISSING_RDKIT_RESULT

    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None, None, None, None, "rdkit_parse_error", "skipped_invalid_structure"

    if mol is None:
        return None, None, None, None, "invalid_smiles", "skipped_invalid_structure"

    try:
        molecular_weight = round(float(Descriptors.MolWt(mol)), MW_DECIMALS)
        molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
        rdkit_status = "rdkit_ok"
    except Exception:
        return None, None, None, None, "rdkit_property_error", "skipped_invalid_structure"

    # Standard InChI: default RDKit options produce the standard string.
    try:
        standard_inchi = Chem.MolToInchi(mol)
    except Exception as exc:
        return (
            molecular_weight,
            molecular_formula,
            None,
            None,
            rdkit_status,
            f"inchi_error:{type(exc).__name__}",
        )

    if not standard_inchi:
        return molecular_weight, molecular_formula, None, None, rdkit_status, "inchi_failed"

    try:
        standard_inchikey = Chem.InchiToInchiKey(standard_inchi)
    except Exception as exc:
        return (
            molecular_weight,
            molecular_formula,
            standard_inchi,
            None,
            rdkit_status,
            f"inchikey_error:{type(exc).__name__}",
        )

    if not standard_inchikey:
        return (
            molecular_weight,
            molecular_formula,
            standard_inchi,
            None,
            rdkit_status,
            "inchikey_failed",
        )

    return (
        molecular_weight,
        molecular_formula,
        standard_inchi,
        standard_inchikey,
        rdkit_status,
        "inchi_ok",
    )


def fetch_pubchem_cid(
    smiles: Optional[str],
    max_retries: int,
    retry_delay: float,
) -> tuple[Optional[int], str]:
    """Retrieve the first PubChem CID matching a SMILES string."""
    if smiles is None:
        return None, "missing_smiles"

    retryable_errors = (
        pcp.PubChemHTTPError,
        urllib.error.URLError,
        TimeoutError,
        IncompleteRead,
        ConnectionError,
    )

    for attempt in range(1, max_retries + 1):
        try:
            compounds = pcp.get_compounds(smiles, namespace="smiles")
            if compounds:
                cid = getattr(compounds[0], "cid", None)
                return (int(cid), "cid_found") if cid is not None else (None, "cid_missing")
            return None, "no_pubchem_match"
        except retryable_errors:
            if attempt < max_retries:
                time.sleep(retry_delay * attempt)
        except Exception as exc:
            return None, f"cid_error:{type(exc).__name__}"

    return None, "cid_retry_exhausted"


def fetch_iupac_name(
    cid: Optional[int],
    max_retries: int,
    retry_delay: float,
) -> tuple[Optional[str], str]:
    """Retrieve an IUPAC name using a PubChem CID."""
    if cid is None or pd.isna(cid):
        return None, "missing_cid"

    cid_int = int(cid)
    retryable_errors = (
        pcp.PubChemHTTPError,
        urllib.error.URLError,
        TimeoutError,
        IncompleteRead,
        ConnectionError,
    )

    for attempt in range(1, max_retries + 1):
        try:
            compound = pcp.Compound.from_cid(cid_int)
            iupac_name = getattr(compound, "iupac_name", None)
            if iupac_name:
                return str(iupac_name), "iupac_found"
            return None, "iupac_missing"
        except retryable_errors:
            if attempt < max_retries:
                time.sleep(retry_delay * attempt)
        except Exception as exc:
            return None, f"iupac_error:{type(exc).__name__}"

    return None, "iupac_retry_exhausted"


def select_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return the final schema, tolerating columns not yet populated."""
    frame = df.copy()
    for column in OUTPUT_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA
    return frame[list(OUTPUT_COLUMNS)]


def format_molecular_weight(value: Any) -> Any:
    """Render molecular weight with a fixed number of decimals, no float artifacts."""
    if value is None or pd.isna(value):
        return pd.NA
    return f"{float(value):.{MW_DECIMALS}f}"


def save_checkpoint(df: pd.DataFrame, path: Path) -> None:
    """Save a restart-friendly checkpoint CSV in the final output schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = select_output_columns(df)
    frame["molecular_weight"] = frame["molecular_weight"].map(format_molecular_weight)
    frame.to_csv(path, index=False)


def build_annotation_status(row: pd.Series) -> str:
    """Create a compact row-level audit status."""
    statuses = [
        str(row.get("rdkit_status", "")),
        str(row.get("inchi_status", "")),
        str(row.get("pubchem_cid_status", "")),
        str(row.get("iupac_status", "")),
    ]
    return ";".join(status for status in statuses if status and status != "nan")


def run_pipeline(
    input_path: Path,
    output_path: Path,
    smiles_column: Optional[str],
    checkpoint_every: int,
    request_delay: float,
    max_retries: int,
    retry_delay: float,
) -> pd.DataFrame:
    """Run the complete annotation workflow."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if input_path.suffix.lower() != ".csv":
        raise ValueError("This combined pipeline expects a .csv input file.")

    df = clean_column_names(pd.read_csv(input_path, encoding="utf-8-sig"))
    smiles_column = resolve_smiles_column(df, smiles_column)
    print(f"Using SMILES column: {smiles_column}")
    print(f"Rows loaded: {len(df):,}")

    normalized_smiles = df[smiles_column].map(normalize_smiles_value)

    # Rebuild on the trimmed schema: the input SMILES are carried through as
    # `canonical_smiles` and every other input column is dropped.
    df = pd.DataFrame({"canonical_smiles": df[smiles_column].values})

    # Steps 1-3: molecular weight, molecular formula, Standard InChI and InChIKey.
    rdkit_cache: dict[str, tuple] = {}
    unique_smiles = [s for s in pd.unique(normalized_smiles) if s is not None]

    for smiles in tqdm(unique_smiles, desc="Steps 1-3: RDKit MW, formula, InChI"):
        rdkit_cache[smiles] = calculate_rdkit_properties(smiles)

    def rdkit_field(position: int):
        return normalized_smiles.map(
            lambda s: rdkit_cache.get(s, MISSING_RDKIT_RESULT)[position]
        ).values

    df["molecular_weight"] = rdkit_field(0)
    df["molecular_formula"] = rdkit_field(1)
    df["standard_inchi"] = rdkit_field(2)
    df["standard_inchikey"] = rdkit_field(3)
    df["rdkit_status"] = rdkit_field(4)
    df["inchi_status"] = rdkit_field(5)

    # Step 4: retrieve PubChem CIDs for unique input SMILES.
    cid_cache: dict[str, tuple[Optional[int], str]] = {}
    checkpoint_path = output_path.with_name(f"{output_path.stem}_checkpoint{output_path.suffix}")

    def apply_cid_columns(default_status: str) -> None:
        df["pubchem_cid"] = normalized_smiles.map(
            lambda s: cid_cache.get(s, (None, default_status))[0]
        ).astype("Int64").values
        df["pubchem_cid_status"] = normalized_smiles.map(
            lambda s: cid_cache.get(s, (None, default_status))[1]
        ).values

    for index, smiles in enumerate(
        tqdm(unique_smiles, desc="Step 4: PubChem CID retrieval"), start=1
    ):
        cid_cache[smiles] = fetch_pubchem_cid(
            smiles=smiles,
            max_retries=max_retries,
            retry_delay=retry_delay,
        )

        if request_delay > 0:
            time.sleep(request_delay)

        if checkpoint_every > 0 and index % checkpoint_every == 0:
            apply_cid_columns("not_processed")
            save_checkpoint(df, checkpoint_path)
            print(f"Checkpoint saved after {index:,} unique SMILES: {checkpoint_path}")

    apply_cid_columns("missing_smiles")

    # Step 5: retrieve IUPAC names from the PubChem CIDs produced in Step 4.
    unique_cids = [int(cid) for cid in df["pubchem_cid"].dropna().unique()]
    iupac_cache: dict[int, tuple[Optional[str], str]] = {}

    def apply_iupac_columns(default_status: str) -> None:
        df["iupac_name"] = df["pubchem_cid"].map(
            lambda cid_value: (
                iupac_cache.get(int(cid_value), (None, default_status))[0]
                if pd.notna(cid_value)
                else None
            )
        )
        df["iupac_status"] = df["pubchem_cid"].map(
            lambda cid_value: (
                iupac_cache.get(int(cid_value), (None, default_status))[1]
                if pd.notna(cid_value)
                else "missing_cid"
            )
        )

    for index, cid in enumerate(
        tqdm(unique_cids, desc="Step 5: IUPAC name retrieval"), start=1
    ):
        iupac_cache[cid] = fetch_iupac_name(
            cid=cid,
            max_retries=max_retries,
            retry_delay=retry_delay,
        )

        if request_delay > 0:
            time.sleep(request_delay)

        if checkpoint_every > 0 and index % checkpoint_every == 0:
            apply_iupac_columns("not_processed")
            save_checkpoint(df, checkpoint_path)
            print(f"Checkpoint saved after {index:,} unique CIDs: {checkpoint_path}")

    apply_iupac_columns("missing_cid")

    df["annotation_status"] = df.apply(build_annotation_status, axis=1)

    df = select_output_columns(df)
    save_checkpoint(df, output_path)

    if checkpoint_path.exists():
        checkpoint_path.unlink()

    print(f"Completed. Output saved to: {output_path}")
    print(f"Valid RDKit structures: {(df['rdkit_status'] == 'rdkit_ok').sum():,}/{len(df):,}")
    print(f"Standard InChI generated: {df['standard_inchi'].notna().sum():,}/{len(df):,}")
    print(f"Standard InChIKey generated: {df['standard_inchikey'].notna().sum():,}/{len(df):,}")
    print(f"PubChem CIDs found: {df['pubchem_cid'].notna().sum():,}/{len(df):,}")
    print(f"IUPAC names found: {df['iupac_name'].notna().sum():,}/{len(df):,}")
    return df


def main() -> None:
    """Run the compound annotation pipeline using the settings above."""
    run_pipeline(
        input_path=INPUT_FILE,
        output_path=OUTPUT_FILE,
        smiles_column=SMILES_COLUMN,
        checkpoint_every=CHECKPOINT_EVERY,
        request_delay=REQUEST_DELAY,
        max_retries=MAX_RETRIES,
        retry_delay=RETRY_DELAY,
    )


if __name__ == "__main__":
    main()
