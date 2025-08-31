# ODO Data Portal Code

This repository contains the reproducible code for the **Opioid Drug Ontology (ODO) Data Portal**.  
It provides scripts, workflows, and environment specifications to reproduce analyses associated with ODO publications.

- **Release v1.0.0** corresponds to the ODO Data Descriptor paper (DOI will be added automatically once the GitHub Release is archived in Zenodo).  
- Future releases will correspond to additional ODO-related publications.  
- All code is designed to run locally with Python 3.x (via conda, see `environment.yml`).  

👉 See the main [README.md](README.md) for installation and usage instructions.

---

## 📂 Repository Structure

```
odo-data-portal-code/
├── pipelines/
│ └── AutoMID_pipeline_S1-S5.py 🧪 generalized CLI script (S1–S5 pipeline)
├── environment.yml 🔧 reproducible conda environment (Python 3 + RDKit + pandas)
├── README.md 📖 quickstart instructions & usage
├── .gitignore 🚫 keep junk out of the repo
├── data/
│ └── example_input.csv 🧪 tiny demo file users can test with
│ └── .gitkeep
├── LICENSE 📜 MIT license file
├── CITATION.cff 📝 citation metadata (links to Zenodo DOI later)
```
