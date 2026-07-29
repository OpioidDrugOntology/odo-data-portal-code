# Publication Table Extraction

## Overview

This workflow extracts tabular data from scientific publication images using PaddleX Table Recognition.

The extracted tables are converted into structured Microsoft Excel (`.xlsx`) and HTML (`.html`) formats for manual review and downstream data curation within the Opioid Drug Ontology Data Portal (ODO-DP).

This workflow was developed to improve the efficiency and reproducibility of publication data extraction while preserving table structure for subsequent quality control.

---

## Workflow

```text
Publication table image (.png)
           │
           ▼
──────────────────────────────
PaddleX Table Recognition
──────────────────────────────
           │
           ▼
Structured table extraction
           │
           ├──────────────► Excel (.xlsx)
           │
           └──────────────► HTML (.html)
```

---

## Repository Structure

```text
table_extraction/
│
├── extract_table.py
├── README.md
└── demo/
    ├── demo_table_image.png
    ├── demo_table_output.xlsx
    └── demo_table_output.html
```

---

## Requirements

Python 3.x

Required Python packages include:

- paddlepaddle
- paddlex
- openpyxl
- pandas

Additional dependencies are provided in the repository `environment.yml`.

---

## Input

Input consists of a publication table image.

Supported image formats include:

- PNG
- JPG
- JPEG

The highest extraction accuracy is generally achieved using high-resolution images.

Example input:

```text
demo/demo_table_image.png
```

---

## Output

The workflow generates two output files.

### Excel Output

A structured Microsoft Excel workbook containing the extracted table.

Example:

```text
demo/demo_table_output.xlsx
```

### HTML Output

An HTML representation of the extracted table preserving the detected table layout.

Example:

```text
demo/demo_table_output.html
```

---

## Usage

Run the extraction workflow:

```bash
python extract_table.py
```

If command-line arguments are supported, consult the script header for additional options.

---

## Demonstration Files

The `demo` directory contains representative example files illustrating the expected workflow inputs and outputs.

These files are provided solely to demonstrate workflow execution and do not represent the complete ODO-DP dataset.

---

## Quality Control

Automatic table recognition should always be followed by manual review.

Users should verify:

- table boundaries
- merged cells
- column headers
- row alignment
- numeric values
- special characters
- missing or incorrectly recognized cells

Corrections should be made before downstream data processing.

---

## Notes

This workflow extracts publication tables only.

Additional manual curation and assay description standardization are performed in subsequent workflows within the Manual Curation Workflow.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
