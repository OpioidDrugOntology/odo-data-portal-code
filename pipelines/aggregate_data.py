import pandas as pd
import numpy as np
import yaml
import argparse

def unique_concat(series, sep=","):
    """Concatenate unique values from a series into one string."""
    vals = [str(x) for x in series.dropna().astype(str) if str(x).strip() != ""]
    seen, out = set(), []
    for v in vals:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return sep.join(out)

def aggregate(input_file, output_file, config_file):
    # Load config
    with open(config_file, "r") as f:
        cfg = yaml.safe_load(f)

    key_cols = cfg["key_columns"]
    value_col = cfg["value_column"]
    numeric_agg = cfg.get("numeric_agg", "median")
    sep = cfg.get("concat_sep", ",")

    # Load dataset
    df = pd.read_excel(input_file) if input_file.endswith(".xlsx") else pd.read_csv(input_file)

    rows = []
    for _, g in df.groupby(key_cols, dropna=False):
        row = {}
        # Numeric aggregation for value_column
        vc = pd.to_numeric(g[value_col], errors="coerce")
        if vc.notna().any():
            if numeric_agg == "median":
                agg_val = float(np.nanmedian(vc.values))
            elif numeric_agg == "mean":
                agg_val = float(np.nanmean(vc.values))
            else:
                raise ValueError(f"Unsupported numeric_agg: {numeric_agg}")
        else:
            agg_val = np.nan

        for c in df.columns:
            if c in key_cols:
                row[c] = g[c].iloc[0]  # keep key column value
            elif c == value_col:
                row[c] = agg_val
            else:
                row[c] = unique_concat(g[c], sep=sep)
        rows.append(row)

    out_df = pd.DataFrame(rows, columns=df.columns)

    # Save
    if output_file.endswith(".xlsx"):
        out_df.to_excel(output_file, index=False)
    else:
        out_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate redundant assay rows.")
    parser.add_argument("--infile", required=True, help="Path to input CSV/XLSX file")
    parser.add_argument("--outfile", required=True, help="Path to save aggregated file")
    parser.add_argument("--config", default="config/aggregate.yaml", help="Config YAML with key/value settings")
    args = parser.parse_args()

    aggregate(args.infile, args.outfile, args.config)
