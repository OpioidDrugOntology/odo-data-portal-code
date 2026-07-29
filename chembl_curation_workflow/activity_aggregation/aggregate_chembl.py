#!/usr/bin/env python3
"""
Aggregate replicate ChEMBL bioactivity rows by collapsing each group to a single
representative: the row whose value is closest to the group median.

Input:  Activity_data_rows_input.csv   (CSV)
Output: aggregated_data_output.xlsx     (Excel workbook, 3 sheets)

Two median paths:
  - pChEMBL path:        rows carrying a 'pChEMBL Value'  -> median on pChEMBL Value
  - Standard Value path: rows whose 'Standard Units' are in STD_VALUE_UNITS
                         -> median on Standard Value
Rows in neither path are left untouched.

Tie-break for even-count / equidistant groups:
  closest to median  ->  earliest Document Year  ->  original file order.
(A missing Document Year is treated as latest, so a dated row is preferred over
an undated one when they tie; if all tied rows are undated, original file order
decides.)

USAGE
-----
Put this script in the SAME folder as the input .csv, open a terminal in that
folder, and run:

    python aggregate_chembl_csv.py

No arguments or file paths are needed. Set INPUT_FILE below to point at a
specific file; otherwise the script uses INPUT_FILE if present, or falls back to
the only .csv it finds in the current folder.

Requires: pandas, openpyxl   ->   pip install pandas openpyxl
"""

import glob
import os
import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# CONFIG  (edit these if needed; all paths are relative to the current folder)
# ---------------------------------------------------------------------------
INPUT_FILE  = "Activity_data_rows_input.csv"
OUTPUT_FILE = "aggregated_data_output.xlsx"

# 9-column grouping key
KEYS = [
    "Smiles", "Molecular Weight", "Standard Type", "Standard Relation",
    "Standard Units", "Target Name", "Target Type", "Target Organism",
    "Assay ChEMBL ID",
]

# Units routed to the Standard Value median path
STD_VALUE_UNITS = ["%", "mg kg-1", "uM", "mg.kg-1"]
# ---------------------------------------------------------------------------


def resolve_input():
    """Use INPUT_FILE if it exists in the current folder; else the only .csv."""
    if os.path.isfile(INPUT_FILE):
        return INPUT_FILE
    candidates = [f for f in glob.glob("*.csv")]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        sys.exit("No .csv file found in this folder. Put the input file here "
                 "or set INPUT_FILE.")
    sys.exit("Multiple .csv files found; set INPUT_FILE to one of:\n  - "
             + "\n  - ".join(candidates))


def read_csv_robust(path):
    """Read a CSV, falling back to latin-1 if the file isn't UTF-8."""
    try:
        return pd.read_csv(path, encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="latin-1")


def select_keep(sub, valuecol):
    """Return the set of orig_idx values to KEEP for one path: one representative
    per group (row closest to the group median), with the documented tie-break."""
    keep = set()
    gsize = sub.groupby(KEYS, dropna=False)[valuecol].transform("size")
    keep.update(sub.loc[gsize == 1, "orig_idx"].tolist())          # singletons kept
    multi = sub[gsize > 1]
    for _, g in multi.groupby(KEYS, dropna=False, sort=False):
        med = g[valuecol].median()
        gg = g.assign(_d=(g[valuecol] - med).abs(),
                      _yr=g["Document Year"].fillna(np.inf))
        gg = gg.sort_values(["_d", "_yr", "orig_idx"])             # tie-break
        keep.add(int(gg["orig_idx"].iloc[0]))
    return keep


def main():
    src = resolve_input()
    print(f"Reading: {src}")
    df = read_csv_robust(src).reset_index(drop=True)
    df["orig_idx"] = np.arange(len(df))
    n = len(df)

    missing = [c for c in KEYS + ["pChEMBL Value", "Standard Value",
                                  "Standard Units", "Document Year"]
               if c not in df.columns]
    if missing:
        sys.exit("Input is missing required columns: " + ", ".join(missing))

    # Path 1: pChEMBL median
    path1 = df[df["pChEMBL Value"].notna()]
    keep1 = select_keep(path1, "pChEMBL Value")

    # Path 2: Standard Value median for the specified units
    path2 = df[df["Standard Units"].isin(STD_VALUE_UNITS) & df["Standard Value"].notna()]
    keep2 = select_keep(path2, "Standard Value")

    path_idx   = set(path1["orig_idx"]) | set(path2["orig_idx"])
    nonpath    = set(df["orig_idx"]) - path_idx            # untouched rows
    final_keep = nonpath | keep1 | keep2

    kept_mask = df["orig_idx"].isin(final_keep)
    out     = df[kept_mask].drop(columns="orig_idx")
    removed = df[~kept_mask].drop(columns="orig_idx")

    removed_p1 = len(path1) - len(keep1)
    removed_p2 = len(path2) - len(keep2)

    summary = pd.DataFrame({
        "Metric": [
            "Original rows",
            "Rows removed by aggregation (total)",
            "  - pChEMBL Value median path",
            "  - Standard Value median path (%, mg kg-1, uM, mg.kg-1)",
            "Rows remaining",
            "Grouping key",
            "Representative rule",
            "Tie-break (even-count / equidistant)",
            "Untouched rows (neither path)",
        ],
        "Value": [
            n, len(removed), removed_p1, removed_p2, len(out),
            " | ".join(KEYS),
            "Row with value closest to group median "
            "(pChEMBL Value for pChEMBL path; Standard Value for the listed units)",
            "Earliest Document Year, then original file order",
            n - len(path1) - len(path2),
        ],
    })

    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as xw:
        out.to_excel(xw, sheet_name="Aggregated", index=False)
        removed.to_excel(xw, sheet_name="Removed rows", index=False)
        summary.to_excel(xw, sheet_name="Aggregation summary", index=False)

    print(f"Original rows:                 {n}")
    print(f"Removed - pChEMBL path:        {removed_p1}")
    print(f"Removed - Standard Value path: {removed_p2}")
    print(f"Removed - total:               {len(removed)}")
    print(f"Rows remaining:                {len(out)}")
    print(f"Written: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
