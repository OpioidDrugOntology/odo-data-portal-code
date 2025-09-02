# ODO Data Portal Code

This repository provides the reproducible **code and workflows** for the **Opioid Drug Ontology (ODO) Data Portal**.  
It includes scripts, pipelines, and environment specifications to ensure analyses can be run locally by any user.

- **Release v1.0.0** will correspond to the ODO Data Descriptor paper (DOI will be added automatically once the GitHub Release is archived in Zenodo).  
- Future releases will correspond to additional ODO-related publications.  
- All code is designed to run locally with Python 3.x (via conda, see `environment.yml`).  

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

```
odo-data-portal-code/
├── pipelines/
│ └── AutoMID_pipeline_S1–S5.py 🧪 generalized CLI script (S1–S5 pipeline)
├── environment.yml 🔧 reproducible conda environment (Python 3 + RDKit + pandas)
├── README.md 📖 quickstart instructions & usage
├── .gitignore 🚫 keep junk out of the repo
├── data/
│ ├── example_input.csv ✏️ tiny demo input file (3 compounds: neutral, salt, isotope)
│ └── example_output.csv 📊 expected output (S0→S5 transformations)
├── LICENSE 📜 MIT license file
└── CITATION.cff 📝 citation metadata (links to Zenodo DOI later)
```
