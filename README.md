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


---
## 📂 Repository Structure

</code></pre>
odo-data-portal-code/
├── pipelines/
│ ├── AutoMID_pipeline_S1–S5.py 📗 pipeline (S1–S5 standardization)
│ └── simple_descriptors.py 📑 InChI, InChIKey, MW, MF
├── data/
│ ├── example_input.csv ✏️ demo input (3 compounds)
│ ├── example_output.csv 📊 S0→S5 transformations
│ └── example_descriptors.csv 📈 computed descriptors
├── environment.yml 🛠️ conda environment
├── README.md 📖 quickstart & usage
├── LICENSE 📜 MIT license
└── CITATION.cff 📝 citation metadata
</code></pre>


---


## 📊 Demo Data

- `data/example_input.csv` — 3 compounds (neutral, HCl salt, [3H]-labeled)  
- `data/example_output.csv` — expected S0→S5 results  
- `data/example_descriptors.csv` — simple descriptors (InChI, InChIKey, MW, MF) generated from `example_output.csv`  

**Highlights**  
- **S2:** removes counter-ions (e.g., `.Cl`) and keeps the largest organic fragment  
- **S5:** clears isotope labels (e.g., `[3H]`)  

**Re-run locally (pipeline):**  
```bash
conda activate odo-chem
python pipelines/AutoMID_pipeline_S1-S5.py \
  --in data/example_input.csv \
  --out data/example_output.csv \
  --smiles-col smiles \
  --id-col odo_id


---

## 🧪 Descriptor Generation

After running the S1–S5 pipeline, you can compute simple molecular descriptors:

```bash
conda activate odo-chem
python pipelines/simple_descriptors.py \
  --in data/example_output.csv \
  --out data/example_descriptors.csv

---

## 📑 Notes

- **Document metadata** (PMID, DOI, patent IDs) were retrieved directly from the [ChEMBL API](https://www.ebi.ac.uk/chembl/) (release 34).  

