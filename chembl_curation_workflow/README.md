# ChEMBL API Reference Cascade

## Overview

This workflow retrieves publication reference information associated with ChEMBL document identifiers using the ChEMBL web services API. A cascading retrieval strategy is used to obtain available publication metadata and compile the results into a structured output dataset suitable for downstream processing within the Opioid Drug Ontology Data Portal (ODO-DP).

The workflow is implemented as a standalone Python script.

---

## Workflow

```text
ChEMBL document identifiers
          │
          ▼
──────────────────────────────────
ChEMBL API Reference Cascade
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
chembl_doc_id_input.csv
```

---

## Output

The workflow generates an Excel file containing the retrieved publication reference information.

Example output:

```text
chembl_api_reference_cascade_output.xlsx
```

---

## Provenance Metadata

A provenance file is generated alongside the Excel output to record metadata associated with the generated dataset.

Example:

```text
chembl_api_reference_cascade_output.xlsx.provenance.json
```

---

## Repository Structure

```text
reference_retrieval/
│
├── README.md
├── chembl_reference_cascade.py
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

Input and output filenames can be modified within the script as required.

---

## Demonstration Files

Representative input and output files are included to illustrate the expected file formats and workflow results.

---

## Notes

This workflow forms the reference retrieval component of the ChEMBL curation pipeline for ODO-DP. The retrieved publication metadata are subsequently used during downstream data curation and integration workflows.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
