# Manual Curation Workflow

## Overview

The Manual Curation Workflow supports the extraction, standardization, and normalization of publication-derived bioactivity data during construction of the Opioid Drug Ontology Data Portal (ODO-DP).

The workflow converts experimental data extracted from peer-reviewed scientific publications into standardized, machine-readable formats suitable for downstream database integration while preserving essential assay metadata and experimental context.

The workflow consists of three sequential processing steps:

1. Publication Table Extraction
2. Assay Description Standardization
3. ASCII Normalization

Each processing step is implemented as an independent Python workflow and can be executed separately.

---

## Workflow

```text
Publication table image
         │
         ▼
──────────────────────────────────
1. Table Extraction
──────────────────────────────────
         │
         ▼
Extracted publication tables
         │
         ▼
──────────────────────────────────
2. Assay Description
   Standardization
──────────────────────────────────
         │
         ▼
Standardized assay descriptions
         │
         ▼
──────────────────────────────────
3. ASCII Normalization
──────────────────────────────────
         │
         ▼
ASCII-compatible assay descriptions
```

---

## Workflow Components

### 1. Table Extraction

Publication tables contained within scientific articles are converted into structured tabular formats using PaddleX Table Recognition.

Outputs include Microsoft Excel (`.xlsx`) and HTML (`.html`) representations of the extracted tables suitable for manual review and downstream processing.

Directory:

```text
table_extraction/
```

---

### 2. Assay Description Standardization

Raw experimental protocols are transformed into standardized publication-ready assay descriptions using predefined assay templates developed for ODO-DP.

Three assay templates are implemented:

- In vitro affinity assays
- In vitro functional assays
- In vivo functional/behavioral assays

The standardized descriptions preserve critical experimental information while providing consistent formatting across the curated dataset.

Directory:

```text
assay_description_standardization/
```

---

### 3. ASCII Normalization

Standardized assay descriptions are converted to ASCII-compatible text by replacing Unicode characters with equivalent ASCII representations.

This step improves interoperability across software platforms, databases, and downstream analysis tools while preserving the original scientific meaning.

Directory:

```text
ascii_normalization/
```

---

## Repository Structure

```text
manual_curation_workflow/
│
├── README.md
│
├── table_extraction/
│   ├── extract_table.py
│   ├── demo/
│   └── README.md
│
├── assay_description_standardization/
│   ├── ...
│   └── README.md
│
└── ascii_normalization/
    ├── ...
    └── README.md
```

---

## Requirements

The workflows require Python 3.x together with the packages listed in the repository `environment.yml`.

Workflow-specific software requirements and dependencies are described within each workflow directory.

---

## Demonstration Data

Each workflow directory contains representative demonstration files illustrating the expected input and output formats.

These demonstration datasets are provided solely to illustrate workflow execution and are not intended to reproduce the complete ODO-DP dataset.

---

## Execution

Each workflow is executed independently.

Refer to the `README.md` within each workflow directory for:

- Software requirements
- Workflow-specific dependencies
- Required input files
- Generated output files
- Execution instructions
- Demonstration datasets

---

## Notes

The Manual Curation Workflow represents the publication-derived component of the ODO-DP data curation pipeline.

Although these workflows automate multiple processing steps, extracted publication tables and generated assay descriptions should be reviewed before incorporation into downstream database workflows to ensure data quality, consistency, and accurate interpretation of the original publication.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
