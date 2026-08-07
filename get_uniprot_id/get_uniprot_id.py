import pandas as pd
import requests
import time

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

input_file = "Get Uniprot ID.xlsx"
output_file = "Get Uniprot ID_with_UniProt.xlsx"

# ------------------------------------------------------------
# Read input file
# ------------------------------------------------------------

df = pd.read_excel(input_file)

# ------------------------------------------------------------
# UniProt lookup function
# ------------------------------------------------------------

def get_uniprot_id(protein_name, taxonomy_id):
    """
    Retrieve the UniProt accession using protein name
    and NCBI taxonomy ID.

    Reviewed Swiss-Prot entries are searched first.
    """

    if pd.isna(protein_name) or pd.isna(taxonomy_id):
        return None

    protein_name = str(protein_name).strip()
    taxonomy_id = int(taxonomy_id)

    url = "https://rest.uniprot.org/uniprotkb/search"

    # Search reviewed UniProtKB/Swiss-Prot records first
    query = (
        f'protein_name:"{protein_name}" '
        f'AND taxonomy_id:{taxonomy_id} '
        f'AND reviewed:true'
    )

    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,protein_name,organism_name",
        "size": 5
    }

    try:
        response = requests.get(
            url,
            params=params,
            timeout=30
        )

        response.raise_for_status()
        data = response.json()

        results = data.get("results", [])

        if results:
            # Return the first matching reviewed accession
            return results[0]["primaryAccession"]

        # ----------------------------------------------------
        # If no reviewed entry is found, search all UniProtKB
        # ----------------------------------------------------

        query = (
            f'protein_name:"{protein_name}" '
            f'AND taxonomy_id:{taxonomy_id}'
        )

        params["query"] = query

        response = requests.get(
            url,
            params=params,
            timeout=30
        )

        response.raise_for_status()
        data = response.json()

        results = data.get("results", [])

        if results:
            return results[0]["primaryAccession"]

        return None

    except requests.RequestException as e:
        print(
            f"UniProt query failed for "
            f"{protein_name} ({taxonomy_id}): {e}"
        )
        return None


# ------------------------------------------------------------
# Retrieve UniProt IDs
# ------------------------------------------------------------

uniprot_ids = []

for index, row in df.iterrows():

    protein_name = row["protein_name"]
    taxonomy_id = row["ncbi_target_taxonomy_id"]

    print(
        f"Searching: {protein_name} | "
        f"Taxonomy ID: {taxonomy_id}"
    )

    uniprot_id = get_uniprot_id(
        protein_name,
        taxonomy_id
    )

    uniprot_ids.append(uniprot_id)

    print(f"  -> {uniprot_id}")

    # Small delay to avoid sending requests too rapidly
    time.sleep(0.2)


# ------------------------------------------------------------
# Add results to dataframe
# ------------------------------------------------------------

df["uniprot_id"] = uniprot_ids


# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------

df.to_excel(
    output_file,
    index=False
)

print("\nCompleted.")
print(f"Output saved to: {output_file}")