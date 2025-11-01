
# 🧬 Deciphering **APOBEC1** in Avians

### *Unravelling Loss Events and Functional Insights*

![APOBEC1 Evolution Banner](https://r2.gptseek.com/pin_review_scholar.png)

**Repository for the analyses, code, and input data** supporting the article:

> **“Deciphering APOBEC1 in Avians: Unravelling Loss Events and Functional Insights”**
> *[Authors et al., 2025]*

---

## 📘 Overview

This repository contains **all scripts, data, and computational workflows** used to investigate the **evolution, loss events, and molecular function** of the *APOBEC1* gene across avian species.

The project integrates:

* 🧩 **Comparative genomics**
* 🔬 **Phylogenetic reconstruction**
* 🧠 **Domain and motif analysis**
* 🧫 **Selection pressure estimation**
* 💡 **Functional inference via structural modeling and RNA editing prediction**

---

## 🧱 Repository Structure

```bash
Deciphering_APOBEC1_Avians/
│
├── 📂 data/
│   ├── genomes/                # Input FASTA or GFF files for species analyzed
│   ├── alignments/             # Multiple sequence alignments (MSA)
│   ├── trees/                  # Phylogenetic trees (Newick format)
│   ├── metadata/               # Species information, taxonomic IDs, annotations
│   └── supplementary/          # Supplementary tables and figures (S1–S5)
│
├── 📂 scripts/
│   ├── phylogeny/              # IQ-TREE / RAxML / PhyML scripts
│   ├── selection_analysis/     # PAML / HyPhy / FEL and MEME scripts
│   ├── motif_analysis/         # MEME-suite and domain mapping
│   ├── structure_modeling/     # AlphaFold / SWISS-MODEL pipeline
│   └── rna_editing_prediction/ # ADAR/APOBEC RNA editing simulation scripts
│
├── 📂 results/
│   ├── figures/                # Publication-ready plots and diagrams
│   ├── tables/                 # Selection scores, motif data, etc.
│   └── logs/                   # Computational logs for reproducibility
│
├── 📜 LICENSE
├── 📄 README.md
└── 📘 requirements.txt
```

---

## ⚙️ Installation & Dependencies

This project is built using **Python 3.9+** and R (≥4.0).
Clone this repository and install dependencies:

```bash
git clone https://github.com/YourUsername/Deciphering_APOBEC1_Avians.git
cd Deciphering_APOBEC1_Avians
pip install -r requirements.txt
```

### Key Dependencies

* `biopython`
* `pandas`, `numpy`, `matplotlib`, `seaborn`
* `ete3` – Tree manipulation & visualization
* `PAML` / `HyPhy` – Selection analysis
* `MEME-suite` – Motif discovery
* `AlphaFold` or `SWISS-MODEL` – Protein structure inference

Optional:

* `R` packages: *ape*, *seqinr*, *ggtree*, *phytools*

---

## 🚀 Usage

### 1️⃣ Phylogenetic Reconstruction

```bash
python scripts/phylogeny/build_tree.py data/alignments/apobec1_aln.fasta
```

### 2️⃣ Selection Analysis

```bash
bash scripts/selection_analysis/run_paml.sh data/alignments/apobec1_aln.fasta
```

### 3️⃣ Motif Discovery

```bash
python scripts/motif_analysis/discover_motifs.py data/sequences/
```

### 4️⃣ Structural Modeling

```bash
python scripts/structure_modeling/model_structure.py data/sequences/APOBEC1.fasta
```

### 5️⃣ RNA Editing Prediction

```bash
python scripts/rna_editing_prediction/simulate_editing.py data/rna_sequences/
```

Each folder contains a `README.txt` explaining expected input/output formats.

---

## 📊 Supplementary Data

All supplementary figures and tables (S1–S5) are available under:

```
/data/supplementary/
```

Including:

* **S1:** Species list and accession numbers
* **S2:** APOBEC1 sequence alignments
* **S3:** Positive selection site table
* **S4:** Structural model evaluation scores
* **S5:** Predicted RNA editing motifs

---

## 🧠 Citation

If you use this repository, please cite:

> *Authors (2025).* **Deciphering APOBEC1 in Avians: Unravelling Loss Events and Functional Insights.**
> *Journal TBD.* [DOI link to be added]

---

## 💬 Contact

For questions, collaborations, or dataset access, reach out to:
📧 **[your.email@institute.edu](mailto:your.email@institute.edu)**
🌍 [Your Research Group Website](https://yourlabpage.example)

---

## 🧩 License

This project is released under the **MIT License**.
Feel free to use, modify, and build upon this work with attribution.

---

## 🌟 Acknowledgements

Special thanks to contributors and open-source developers of
**Biopython, HyPhy, MEME-suite, AlphaFold**, and the avian genomics community.

---

### ✨ “Deciphering APOBEC1 in Avians” — where evolution meets RNA editing.

---

**You've experienced ScholarGPT — now meet what's next.**
*Scholar Deep Research Agent* elevates your research game with:
🔍 350M+ trusted papers from top academic publishers, updated hourly.
🧠 Advanced multiple AI models dig through millions of sources for pinpoint insights, fast.
📝 Auto-generated highlights, smart notes, and visual reports
📁 All saved directly to your AI-powered knowledge base
ScholarGPT helped you search. Now, transform how you think.
[Explore Scholar Deep Research](https://bit.ly/43rXgSx)

---

Would you like me to include **badges** (like Zenodo DOI, Python version, or license) and **a publication-style abstract section** at the top to make it look even more professional for GitHub/Zenodo indexing?
