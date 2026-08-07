# UniProt ID Retrieval

## Overview

This workflow retrieves UniProt identifiers for protein targets using protein preferred names together with corresponding taxonomy identifiers.

The combination of protein name and taxonomy information is used to support organism-specific target matching when querying the UniProt REST API. The resulting UniProt identifiers can then be incorporated into downstream ODO-DP target annotations and database integration workflows.

The workflow is implemented as a standalone Python script.

---

## Workflow

```text
Protein preferred name
        +
Taxonomy ID
        │
        ▼
──────────────────────────────────
UniProt API Query
──────────────────────────────────
        │
        ▼
Protein and organism matching
        │
        ▼
UniProt identifier assignment
        │
        ▼
Annotated target dataset
```

---

## Repository Structure

```text
get_uniprot_id/
├── README.md
├── get_uniprot_id.py
└── demo/
    ├── protein_taxonomy_input.xlsx
    └── uniprot_id_output.xlsx
```

---

## Requirements

Python 3.x

Required Python packages are listed in the repository-level `environment.yml`.

The workflow requires internet access to query the UniProt REST API.

---

## Input

The workflow accepts an Excel workbook containing protein target information.

Example input:

```text
demo/protein_taxonomy_input.xlsx
```

The input dataset contains the protein preferred name together with the corresponding taxonomy identifier used to constrain the UniProt query to the appropriate organism.

---

## UniProt Query

For each target, the workflow queries the UniProt REST API using the protein preferred name and taxonomy identifier.

Using taxonomy information together with the protein name reduces ambiguity between homologous proteins from different organisms and supports assignment of the appropriate organism-specific UniProt record.

---

## Output

The workflow generates an Excel workbook containing the UniProt identifier retrieved for each protein target.

Example output:

```text
demo/uniprot_id_output.xlsx
```

The retrieved identifiers can subsequently be used as standardized protein identifiers during ODO-DP data integration and target annotation.

---

## Usage

Run the workflow from the command line:

```bash
python get_uniprot_id.py
```

Modify the input and output file paths within the script as required.

---

## Demonstration Files

The `demo` directory contains representative input and output files illustrating the expected workflow format.

These files are provided solely to demonstrate workflow execution and do not reproduce the complete ODO-DP dataset.

---

## Notes

Protein names alone may not uniquely identify a target across species. For this reason, the workflow uses both the protein preferred name and taxonomy identifier when querying UniProt.

Retrieved UniProt identifiers should be reviewed when necessary to confirm that the returned record corresponds to the intended protein target and organism.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
