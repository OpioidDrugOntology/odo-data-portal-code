# Assay Description Standardization

## Overview

This workflow generates standardized assay descriptions from publication-derived experimental protocols using predefined assay templates developed for the Opioid Drug Ontology Data Portal (ODO-DP).

The standardized descriptions preserve essential experimental information while providing a consistent, publication-ready format suitable for downstream database integration.

To avoid reproducing copyrighted assay descriptions verbatim, standardized descriptions are generated from structured assay templates that retain critical assay concepts and assay conditions.

---

## Workflow

```text
Publication-derived
assay description
          │
          ▼
──────────────────────────────────
Assay Concept Identification
──────────────────────────────────
          │
          ▼
Automatic Template Classification
          │
          ├──► Template A
          │     In vitro affinity binding assays
          │
          ├──► Template B
          │     In vitro functional assays
          │
          └──► Template C
                In vivo functional assays
          │
          ▼
──────────────────────────────────
Template-Based Standardization
──────────────────────────────────
          │
          ▼
Standardized assay description
```

---

## Assay Templates

Three assay templates are implemented.

### Template A

**In vitro affinity assays**

Used for receptor binding and affinity measurements.

---

### Template B

**In vitro functional assays**

Used for functional cellular assays measuring receptor activation or signaling.

---

### Template C

**In vivo functional/behavioral assays**

Used for animal pharmacology and behavioral experiments.

---

## Repository Structure

```text
assay_description_standardization/
│
├── assay_description_standardization.py
├── README.md
└── demo/
    ├── Assay_standardization_input_dataset.csv
    ├── Assay_standardization_output_dataset.xlsx
    ├── Assay_standardization_output_dataset_debug.xlsx
    └── system_prompt.txt
```

---

## Requirements

Python 3.x

Required Python packages are listed in the repository `environment.yml`.

### Anthropic API

This workflow uses Anthropic's Claude API, accessed through the Anthropic Python SDK, to generate standardized assay descriptions.

Before running the workflow, obtain an Anthropic API key and configure it as an environment variable.

**macOS / Linux**

```bash
export ANTHROPIC_API_KEY="your_api_key_here"
```

**Windows (Command Prompt)**

```cmd
set ANTHROPIC_API_KEY=your_api_key_here
```

**Windows (PowerShell)**

```powershell
$env:ANTHROPIC_API_KEY="your_api_key_here"
```

The API key is not included with this repository and must be obtained separately from Anthropic.

---

## Input

Input consists of publication-derived assay descriptions.

Each assay description is analyzed for standardized assay concepts and automatically classified into one of three assay templates:

- **Template A:** In vitro affinity binding assays
- **Template B:** In vitro functional assays
- **Template C:** In vivo functional assays

The selected template is then used to generate a standardized assay description that preserves the critical experimental information reported in the original publication while providing a consistent, standardized format suitable for downstream data integration and analysis.

The assay templates were developed for ODO-DP using a structured, template-based approach designed to generate consistent assay descriptions while preserving the key experimental information reported in the original publication. The template design follows the same general philosophy used by ChEMBL for standardized assay descriptions but was developed independently for ODO-DP.

Example input:

```text
demo/Assay_standardization_input_dataset.csv
```

---

## Output

The workflow generates standardized publication-ready assay descriptions using the appropriate assay template.

Each input assay receives one standardized description.

---

## Usage

After configuring the Anthropic API key, run the workflow:

```bash
python assay_description_standardization.py
```

If command-line arguments are supported, consult the script header for additional options.

---

## Demonstration Files

The `demo` directory contains representative example input and output files illustrating workflow execution.

These files are provided solely to demonstrate workflow functionality and do not reproduce the complete ODO-DP dataset.

---

## Notes

The standardized descriptions are designed to preserve the essential scientific content of the original experimental protocols while providing consistent formatting across the curated dataset.

This workflow is followed by the ASCII Normalization workflow, which converts standardized descriptions into ASCII-compatible text.

---

## Citation

If you use this workflow, please cite the accompanying ODO Data Descriptor publication.

*(Reference will be added following publication.)*
