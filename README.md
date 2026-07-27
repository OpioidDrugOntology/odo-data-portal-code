# ODO Data Portal Code  

This repository provides the reproducible code and workflows supporting construction of the **Opioid Drug Ontology (ODO) Data Portal (ODO-DP)**.

The repository contains the data curation and processing workflows used to generate the ODO-DP dataset described in the accompanying publication. It includes reproducible Python workflows, demonstration datasets, and environment specifications to support local execution and reuse.

Release **v1.0.0** corresponds to the ODO Data Descriptor publication (DOI will be added automatically once the GitHub Release is archived in Zenodo).

Future releases will correspond to additional ODO-related publications and associated workflows.

All code is designed to run locally using Python 3.x through the provided Conda environment (`environment.yml`).


---

## 🚀 Quickstart  

Clone the repository and create the Conda environment.

```bash
git clone https://github.com/OpioidDrugOntology/odo-data-portal-code.git
cd odo-data-portal-code

conda env create -f environment.yml
conda activate odo-chem
```

Each workflow contains:

- Python scripts
- demonstration input files
- demonstration output files
- workflow-specific documentation

Refer to the README within each workflow directory for execution instructions.


---

## ⚙️ Requirements

- Conda (Miniconda or Anaconda)
- Python 3.10 or later
- All required Python packages are installed automatically using the provided `environment.yml` file.


---

## 📦 Installation

Clone the repository and create the Conda environment.

```bash
git clone https://github.com/OpioidDrugOntology/odo-data-portal-code.git
cd odo-data-portal-code

conda env create -f environment.yml
conda activate odo-chem
```

---

## 📂 Repository Structure

```text
odo-data-portal-code/
│
├── manual_curation_workflow/
│   ├── table_extraction/
│   ├── assay_description_standardization/
│   └── ascii_normalization/
│
├── chembl_curation_workflow/
│   ├── reference_retrieval/
│   └── activity_aggregation/
│
├── compound_processing_workflow/
│   ├── chemical_standardization/
│   ├── compound_annotation/
│   └── descriptor_generation/
│
├── general_data_processing_workflow/
│   ├── qikprop_aggregation/
│   └── final_schema_alignment/
│
├── environment.yml
├── LICENSE
└── README.md

```

---

# 🔬 Workflow Overview

## 1. Manual Curation Workflow

Supports curation of bioactivity data extracted from scientific publications.

Current workflow components include:

- publication table extraction
- assay description standardization
- ASCII normalization

Each processing step is documented within the workflow directory.

---

## 2. ChEMBL Curation Workflow

Processes bioactivity data obtained from ChEMBL.

Current workflow components include:

- document reference retrieval (DOI, PMID, patent identifiers)
- activity aggregation and redundancy processing

Each processing step is documented within the workflow directory.

---

## 3. Compound Processing Workflow

Processes compounds originating from both manually curated publications and ChEMBL-derived datasets.

Current workflow components include:

- chemical standardization
- compound annotation
- molecular descriptor generation

Each processing step is documented within the workflow directory.

---

## 4. General Data Processing Workflow

Contains dataset-level processing utilities applied during construction of the final ODO-DP release.

Current workflow components include:

- QikProp data aggregation
- final schema alignment
- additional data processing utilities

Each processing step is documented within the workflow directory.

---

# 📊 Demonstration Data

Each workflow includes representative demonstration datasets illustrating the expected input and output files for the associated scripts.

The demonstration datasets are intended solely to illustrate workflow execution and are not intended to reproduce the complete ODO-DP dataset.

---


## 📑 Notes

Individual workflow directories contain detailed documentation describing:

- required inputs
- generated outputs
- execution examples
- software dependencies
- workflow-specific implementation details

---

# 📖 Citation

If you use this repository, please cite the accompanying ODO Data Descriptor publication.

(Reference will be added following publication.)

---

# 📄 License

See the LICENSE file for licensing information.
