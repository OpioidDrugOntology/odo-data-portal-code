import pandas as pd

# Load the demo input file
df = pd.read_csv("data/aggregate_demo_input.csv")

# Convert p-values to endpoint values (nM) if present
if "chembl_p-value_concentration" not in df.columns:
    raise ValueError("Input file must contain 'chembl_p-value_concentration' column.")

# Calculate endpoint value in nM = 10^(-pValue) * 1e9
df["endpoint_value_nM"] = 10 ** (-df["chembl_p-value_concentration"]) * 1e9

# Save back to the same CSV (overwrites it)
df.to_csv("data/aggregate_demo_input.csv", index=False)

print("Updated aggregate_demo_input.csv with endpoint_value_nM column.")

