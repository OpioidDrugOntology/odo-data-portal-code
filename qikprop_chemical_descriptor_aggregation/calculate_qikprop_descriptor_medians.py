from pathlib import Path

import pandas as pd


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
INPUT_FILE = Path("Median_calculation_input.csv")
OUTPUT_FILE = Path("Median_calculation_output.xlsx")
SMILES_COLUMN = "canonical_smiles"

MEDIAN_COLUMNS = [
    "dipole",
    "SASA",
    "FISA",
    "donorHB",
    "accptHB",
    "QPlogPw",
    "QPlogPo/w",
    "QPlogS",
]


def calculate_qikprop_medians() -> None:
    """Calculate median QikProp values for each unique canonical SMILES."""

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE.resolve()}"
        )

    df = pd.read_csv(INPUT_FILE)

    required_columns = [SMILES_COLUMN, *MEDIAN_COLUMNS]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(
            "The input file is missing the following required column(s): "
            + ", ".join(missing_columns)
        )

    # Remove rows that do not contain a usable SMILES value.
    df[SMILES_COLUMN] = df[SMILES_COLUMN].astype("string").str.strip()
    df = df[df[SMILES_COLUMN].notna() & df[SMILES_COLUMN].ne("")].copy()

    # Convert descriptor fields to numeric values. Invalid text is treated as
    # missing and is excluded from the corresponding median calculation.
    for column in MEDIAN_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # Calculate one median row per unique canonical SMILES.
    medians = (
        df.groupby(SMILES_COLUMN, as_index=False, sort=False)[MEDIAN_COLUMNS]
        .median()
    )

    # Record the number of source rows represented by each aggregated row.
    source_counts = (
        df.groupby(SMILES_COLUMN, sort=False)
        .size()
        .rename("n_source_rows")
        .reset_index()
    )

    output = medians.merge(source_counts, on=SMILES_COLUMN, how="left")
    output = output[[SMILES_COLUMN, "n_source_rows", *MEDIAN_COLUMNS]]

    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
        output.to_excel(writer, sheet_name="QikProp medians", index=False)

    print(f"Input rows: {len(df):,}")
    print(f"Unique canonical SMILES: {len(output):,}")
    print(f"Output written to: {OUTPUT_FILE.resolve()}")


if __name__ == "__main__":
    calculate_qikprop_medians()
