# Compound Processing Workflow

## Overview

This workflow annotates chemical compounds using canonical SMILES strings as input. Compound identifiers and associated chemical information are retrieved from public chemical resources and compiled into a structured dataset suitable for integration into the Opioid Drug Ontology Data Portal (ODO-DP).

The workflow is implemented as a standalone Python script.

---

## Workflow

```text
Canonical SMILES
          │
          ▼
──────────────────────────────────
Compound Annotation Pipeline
──────────────────────────────────
          │
          ▼
Retrieve compound
annotations
          │
          ▼
Annotated compound dataset
```

---

## Input

The workflow accepts a CSV file containing canonical SMILES strings.

Example input:

```text
demo/smiles_input_file.csv
```

---

## Output

The workflow generates an excel file containing the annotated compound information.

Example output:

```text
demo/compound_processing_output.xlsx
```

---

## Compound Annotation

The workflow retrieves compound annotations associated with each canonical SMILES string and compiles the available information into a structured output dataset for downstream curation and integration.

---

## Repository Structure

```text
compound_processing_workflow/
│
├── README.md
├── compound_annotation_pipeline.py
└── demo/
    ├── smiles_input_file.csv
    └── compound_processioutput.xlsxcsv
---

## Requirements

The workflow requires Python 3.x together with the packages listed in the repository `environment.yml`.

---

## Usage

Run the workflow from the command line:

```bash
python compound_annotation_pipeline.py
```

Modify the input and output filenames within the script as required.

---

## Demonstration Files

Representative input and output files are included in the `demo` directory to illustrate the expected file formats and workflow results.

---

## Notes

This workflow provides the compound annotation component of the ODO-DP data processing pipeline. The annotated compound information supports downstream data curation, annotation, and database integration workflows.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
