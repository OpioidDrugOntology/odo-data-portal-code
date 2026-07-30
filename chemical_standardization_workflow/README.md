# Chemical Structure Standardization Workflow

## Overview

This workflow implements a modified version of the AutoMID chemical structure standardization pipeline using exclusively open-source RDKit functionality. The workflow standardizes molecular structures to generate chemically normalized parent compounds suitable for cheminformatics analyses, virtual screening, and downstream annotation.

The implementation is derived from the AutoMID pipeline, which itself is based on the ChEMBL molecule standardization framework. Unlike the original ChEMBL workflow, this implementation relies entirely on RDKit and does not require proprietary ChEMBL structure-checking components, providing a transparent, portable, and reproducible preprocessing workflow.

The workflow was used during construction of the Opioid Drug Ontology Data Portal (ODO-DP) dataset.

---

## Workflow

Input SMILES are processed through five sequential standardization stages.

### S1 – Canonicalization

Canonical SMILES are generated to remove syntactic variation and standardize atom ordering.

### S2 – Desalting

RDKit's `FragmentRemover` removes inorganic salts and small fragments. When multiple fragments remain, `LargestFragmentChooser` selects the largest organic component.

### S3 – Charge Neutralization

SMARTS-based transformation rules are applied to neutralize formal charges where chemically appropriate, resolving ion-paired species and zwitterions while preserving chemically valid structures.

### S4 – Tautomer Canonicalization

RDKit's `TautomerEnumerator` converts each molecule to its representative canonical tautomer, reducing redundancy arising from alternative tautomeric forms.

### S5 – Hydrogen Removal

Explicit hydrogen atoms are removed to generate standardized parent structures compatible with downstream cheminformatics and virtual screening workflows.

For each compound, the workflow records the following processing stages:

| Stage | Description |
|-------|-------------|
| **S0** | Original SMILES |
| **S1** | Canonicalized |
| **S2** | Desalted |
| **S3** | Neutralized |
| **S4** | Tautomer-standardized |
| **S5** | Final standardized parent SMILES |

---

## Repository Structure

```text
chemical_standardization_workflow/
├── README.md
├── AutoMID_pipeline_S1-S5.py
└── demo/
    ├── smiles_input.csv
    └── standardized_smiles_output.xlsx
```

---

## Requirements

- Python 3.10 or later
- RDKit
- pandas
- openpyxl

All required dependencies can also be installed using the repository-level `environment.yml`.

---

## Demonstration Data

The `demo` directory contains a small example dataset illustrating execution of the workflow.

### Input

**smiles_input.csv**

Representative SMILES strings used for chemical structure standardization.

### Output

**standardized_smiles_output.xlsx**

Excel workbook containing the standardized structures generated during each processing stage (S0–S5).

---

## Running the Workflow

Run the workflow from the command line:

```bash
python AutoMID_pipeline_S1-S5.py
```

Modify the input and output file paths within the script as required.

---

## Output

The workflow generates standardized parent compound structures after each processing stage, allowing inspection of every transformation from the original SMILES (S0) through the final standardized parent structure (S5).

These standardized structures may subsequently be used for compound annotation, generation of International Chemical Identifier (InChI) and InChIKey identifiers, molecular descriptor calculation, and virtual screening workflows.

---

## References

This workflow implements a modified RDKit-based version of the AutoMID pipeline described by:

- Bento AP, et al. (2020)
- Stathias V, et al. (2020)

---

## Citation

If you use this workflow in your research, please cite:

- The ODO-DP Data Descriptor publication.
- The **chemical structure curation pipeline** references shown above.


