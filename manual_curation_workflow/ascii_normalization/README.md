# ASCII Normalization Workflow

## Overview

This workflow converts publication-ready assay descriptions into ASCII-compliant text suitable for publication, data repositories, and downstream computational analyses. Non-ASCII characters commonly encountered in biomedical literature—including Greek letters, superscripts, subscripts, degree symbols, and other special characters—are replaced with standardized ASCII equivalents while preserving the scientific meaning of the original text.

The workflow applies context-dependent replacement rules to ensure scientifically accurate normalization. For example, the Greek letter μ is converted to **mu** when referring to the opioid receptor (e.g., mu opioid receptor) but to **u** when used as the SI prefix in measurement units (e.g., μM → uM, μg → ug, μL → uL).

---

## Workflow

**Input assay descriptions**

↓

**Detect non-ASCII characters**

↓

**Apply context-dependent normalization rules**

↓

**Generate ASCII-compliant assay descriptions**

---

## Repository Structure

```text
ascii_normalization/
├── README.md
├── ascii_normalize.py
└── demo/
    ├── ascii_normalization_input_dataset.csv
    └── ASCII_compliant_assay_descriptions_output.xlsx
```

---

## Requirements

- Python 3.10 or later
- pandas
- openpyxl

---

## Demonstration Data

The **demo** directory contains representative files for illustrating the workflow.

**Input**

- `ascii_normalization_input_dataset.csv`

**Output**

- `ASCII_compliant_assay_descriptions_output.xlsx`

---

## Running the Workflow

Place the input dataset in the working directory and execute:

```bash
python ascii_normalize.py
```

The script generates an Excel workbook containing the ASCII-normalized assay descriptions.

---

## Output

The output workbook contains assay descriptions in which non-ASCII characters have been replaced with standardized ASCII representations while preserving the scientific content and terminology required for publication and downstream analyses.

---

## Notes

Normalization rules are context dependent rather than simple character substitutions. Examples include:

| Original | ASCII Output |
|----------|--------------|
| μ opioid receptor | mu opioid receptor |
| μM | uM |
| μg | ug |
| μL | uL |
| α | alpha |
| β | beta |
| γ | gamma |
| δ | delta |
| κ | kappa |
| ± | +/- |
| °C | deg C |
| Ca²⁺ | Ca2+ |
| [³⁵S] | [35S] |

These rules were developed to preserve critical scientific information while ensuring compatibility with software systems and publication platforms that require ASCII-only text.
