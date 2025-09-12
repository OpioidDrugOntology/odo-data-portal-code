# ODO Data Portal Code  

This repository provides the reproducible **code and workflows** for the **Opioid Drug Ontology (ODO) Data Portal**.  
It includes scripts, pipelines, and environment specifications to ensure analyses can be run locally by any user.

- **Release v1.0.0** will correspond to the ODO Data Descriptor paper (DOI will be added automatically once the GitHub Release is archived in Zenodo).  
- Future releases will correspond to additional ODO-related publications.  
- All code is designed to run locally with Python 3.x (via conda, see `environment.yml`).  

---

## 🚀 Quickstart  

Run the demo pipeline and descriptor generation locally:  

```bash
# 1. Set up environment (once)
conda env create -f environment.yml
conda activate odo-chem

# 2. Run pipeline (input → standardized output)
python pipelines/AutoMID_pipeline_S1-S5.py \
  --in data/example_input.csv \
  --out data/example_output.csv \
  --smiles-col smiles \
  --id-col odo_id

# 3. Run descriptor generation (output → InChI, InChIKey, MW, MF)
python pipelines/simple_descriptors.py \
  --in data/example_output.csv \
  --out data/example_descriptors.csv

# 4. Run aggregation demo (replicates → single row + median)
python pipelines/aggregate_data.py \
  --infile data/aggregate_demo_input.csv \
  --outfile data/aggregate_demo_output.csv \
  --config config/aggregate.yaml

---

## ⚙️ Requirements

- [Conda (Miniconda/Anaconda)](https://docs.conda.io/en/latest/miniconda.html)  
- Python **3.8+** (tested with Python 3.10/3.11)  
- [RDKit](https://www.rdkit.org/) and [pandas](https://pandas.pydata.org/) (installed via the provided environment file)  
- [PyYAML](https://pyyaml.org/)  

> Inside the conda environment, `python` will already be Python 3, so either `python` or `python3` works.

---

## 📦 Installation

Clone the repo and create the environment:

```bash
git clone https://github.com/OpioidDrugOntology/odo-data-portal-code.git
cd odo-data-portal-code

conda env create -f environment.yml
conda activate odo-chem
``` 

---

## 📂 Repository Structure

```<pre>
odo-data-portal-code/
├── pipelines/
│   ├── AutoMID_pipeline_S1-S5.py     # ⚙️ main standardization pipeline
│   ├── simple_descriptors.py         # ⚛️ descriptor generation script
│   └── aggregate_data.py             # 🧩 redundancy-aware aggregation script
├── data/
│   ├── example_input.csv             # ✏️ demo input (3 compounds)
│   ├── example_output.csv            # 📊 pipeline output (S0→S5 results)
│   ├── example_descriptors.csv       # ⚛️ descriptor output (InChI, InChIKey, MW, MF)
│   ├── aggregate_demo_input.csv      # ✏️ aggregation demo input (replicates/variants)
│   └── aggregate_demo_output.csv     # 📊 aggregation demo output (comma-joined + median)
├── config/
│   ├── config.yaml                   # ⚙️ pipeline defaults
│   └── aggregate.yaml                # ⚙️ aggregation keys & median rule
├── environment.yml                   # 🛠️ conda environment setup
├── LICENSE                           # 📜 license file
└── README.md                         # 📖 project documentation

```


---


## 📊 Demo Data
Workflow:
📄 example_input.csv → ⚙️ Pipeline (S0→S5) → 📊 example_output.csv → ⚛️ Descriptor Generation → 📈 example_descriptors.csv

data/example_input.csv — 3 compounds (neutral, HCl salt, [3H]-labeled)

data/example_output.csv — standardized results after S0→S5 pipeline

data/example_descriptors.csv — computed molecular descriptors (InChI, InChIKey, MW, MF)

data/aggregate_demo_input.csv — demo input with replicates and variants (20 rows)

data/aggregate_demo_output.csv — aggregated output: all values preserved as comma-separated lists plus median of endpoint in endpoint_value_nM_median


🔑 Highlights
---
⚙️ Pipeline (S0→S5)

S2: Removes counter-ions (e.g., .Cl) and retains the largest organic fragment

S5: Clears isotope labels (e.g., [3H]) for clean standardization

⚛️ Descriptor Generation

Converts each standardized SMILES into InChI and InChIKey

Computes molecular properties: Molecular Weight (MW) and Molecular Formula (MF)

🧩 Aggregation 

Identifies redundant experimental rows across sources

Computes median endpoint values (endpoint_value_nM_median)

Concatenates metadata (e.g., PMIDs, notes) into comma-separated lists for transparency


⚙ Re-run Pipeline Locally
---
conda activate odo-chem
python pipelines/AutoMID_pipeline_S1-S5.py \
  --in data/example_input.csv \
  --out data/example_output.csv \
  --smiles-col smiles \
  --id-col odo_id

After running the S1–S5 pipeline, you can compute simple molecular descriptors:


⚛ Run Descriptor Generation

conda activate odo-chem
python pipelines/simple_descriptors.py \
  --in data/example_output.csv \
  --out data/example_descriptors.csv


🧩 Run Aggregation Demo

conda activate odo-chem
python pipelines/aggregate_data.py \
  --in data/aggregate_demo_input.csv \
  --out data/aggregate_demo_output.csv \
  --config config/aggregate.yaml

---


## 📑 Notes

- **Document metadata** (PMID, DOI, patent IDs) were retrieved directly from the [ChEMBL API](https://www.ebi.ac.uk/chembl/) (release 34).  

- **Because these are simple API lookups, no custom code is required here — just cite the ChEMBL API as the source.

- **All non-trivial processing code (pipeline, descriptors) is deposited in this repository for transparency.
