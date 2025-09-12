import pandas as pd
import numpy as np
import sys

# Usage:
# python scripts/add_endpoint_value_nm.py input.csv output.csv
inp = sys.argv[1]
out = sys.argv[2]

UNIT_TO_NM = {
    "pm": 1e-3,   # pM -> nM
    "nm": 1.0,    # nM -> nM
    "um": 1e3,    # uM/µM/μM -> nM
    "µm": 1e3,
    "μm": 1e3,
    "mm": 1e6,    # mM -> nM
    "m":  1e9,    # M  -> nM
}

def to_nm(val, unit):
    try:
        v = float(val)
    except Exception:
        return np.nan
    if pd.isna(unit):
        return np.nan
    u = str(unit).strip().lower().replace(" ", "")
    return v * UNIT_TO_NM.get(u, np.nan)

# Read (CSV or XLSX)
if inp.lower().endswith(".xlsx"):
    df = pd.read_excel(inp)
else:
    df = pd.read_csv(inp)

# Require these columns
required = ["endpoint", "endpoint_qualifier", "unit_of_measurement", "chembl_p-value_concentration"]
for c in required:
    if c not in df.columns:
        raise SystemExit(f"Missing required column: {c}")

# Compute normalized values
df["endpoint_value_nM"] = [
    to_nm(v, u) for v, u in zip(df["chembl_p-value_concentration"], df["unit_of_measurement"])
]

# Reorder columns: place endpoint_value_nM after endpoint_qualifier
cols = list(df.columns)
if "endpoint_value_nM" in cols:
    cols.remove("endpoint_value_nM")
    insert_pos = cols.index("endpoint_qualifier") + 1
    cols.insert(insert_pos, "endpoint_value_nM")
    df = df[cols]

# Save
if out.lower().endswith(".xlsx"):
    df.to_excel(out, index=False)
else:
    df.to_csv(out, index=False)

print(f"Wrote {out} with endpoint_value_nM inserted after endpoint_qualifier")

