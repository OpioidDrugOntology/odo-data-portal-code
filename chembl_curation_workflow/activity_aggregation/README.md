# ChEMBL Activity Aggregation

## Overview

This workflow aggregates ChEMBL bioactivity records to generate representative activity values suitable for integration into the Opioid Drug Ontology Data Portal (ODO-DP). Records corresponding to the same experimental context are grouped, non-numeric information is consolidated, and representative activity values are calculated according to predefined aggregation rules.

The workflow is implemented as a standalone Python script.

---

## Workflow

```text
ChEMBL activity records
          │
          ▼
──────────────────────────────────
ChEMBL Activity Aggregation
──────────────────────────────────
          │
          ▼
Group related activity records
          │
          ▼
Aggregate activity values
and associated metadata
          │
          ▼
Aggregated activity dataset
```

---

## Input

The workflow accepts a CSV file containing ChEMBL activity records.

Example input:

```text
demo/Activity_data_rows_input.csv
```

---

## Output

The workflow generates an Excel workbook containing the aggregated activity records.

Example output:

```text
demo/aggregated_data_output.xlsx
```

---

## Aggregation

The workflow groups related activity records and applies predefined aggregation rules to generate representative activity values while preserving associated experimental metadata where appropriate. The resulting dataset provides a single representative record for each aggregated activity group suitable for downstream integration.

---

## Repository Structure

```text
activity_aggregation/
│
├── README.md
├── aggregate_chembl.py
└── demo/
    ├── Activity_data_rows_input.csv
    └── aggregated_data_output.xlsx
```

---

## Requirements

The workflow requires Python 3.x together with the packages listed in the repository `environment.yml`.

---

## Usage

Run the workflow from the command line:

```bash
python aggregate_chembl.py
```

Modify the input and output filenames within the script as required.

---

## Demonstration Files

Representative input and output files are included in the `demo` directory to illustrate the expected file formats and workflow results.

---

## Notes

This workflow provides the activity aggregation component of the ChEMBL curation pipeline for ODO-DP. The aggregated activity records support downstream data curation, annotation, and integration workflows.

Detailed aggregation rules and the rationale for representative record selection are described in the accompanying Zenodo workflow documentation.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
