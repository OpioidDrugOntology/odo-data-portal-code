# ChEMBL API Reference Retrieval

## Overview

This workflow retrieves publication metadata associated with ChEMBL document identifiers using the ChEMBL web services API. The workflow processes a list of ChEMBL document identifiers, queries the ChEMBL API for available publication metadata, and compiles the results into a structured output dataset suitable for downstream curation within the Opioid Drug Ontology Data Portal (ODO-DP).

The workflow is implemented as a standalone Python script.

---

## Workflow

```text
ChEMBL document identifiers
          │
          ▼
──────────────────────────────────
ChEMBL API Reference Retrieval
──────────────────────────────────
          │
          ▼
Retrieve publication metadata
using the ChEMBL API
          │
          ▼
Reference dataset
          │
          ├──► Excel output
          │
          └──► Provenance metadata
```

---

## Input

The workflow accepts a CSV file containing ChEMBL document identifiers.

Example input:

```text
demo/chembl_doc_id_input.csv
```

---

## Output

The workflow generates an Excel workbook containing the retrieved publication metadata.

Example output:

```text
demo/chembl_api_reference_cascade_output.xlsx
```

---

## Provenance Metadata

A provenance metadata file is generated alongside the Excel output to document information associated with the workflow execution.

Example:

```text
demo/chembl_api_reference_cascade_output.xlsx.provenance.json
```

---

## Repository Structure

```text
reference_retrieval/
│
├── README.md
├── chembl_reference_cascade.py
└── demo/
    ├── chembl_doc_id_input.csv
    ├── chembl_api_reference_cascade_output.xlsx
    └── chembl_api_reference_cascade_output.xlsx.provenance.json
```

---

## Requirements

The workflow requires Python 3.x together with the packages listed in the repository `environment.yml`.

---

## Usage

Run the workflow from the command line:

```bash
python chembl_reference_cascade.py
```

Modify the input and output filenames within the script as required.

---

## Demonstration Files

Representative input and output files are included in the `demo` directory to illustrate the expected file formats and workflow results.

---

## Notes

This workflow provides the publication metadata retrieval component of the ChEMBL curation pipeline for ODO-DP. The retrieved metadata support downstream data curation, annotation, and integration workflows.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
