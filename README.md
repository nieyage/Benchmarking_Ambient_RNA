# Systematic Benchmarking of Ambient RNA Decontamination Tools

This repository contains code for the manuscript:  
**"Systematic Benchmarking of Ambient RNA Decontamination Tools to Advance Precision in Single-Cell Transcriptomic Analysis"**

## Project Overview

Ambient RNA contamination in single-cell RNA sequencing introduces exogenous transcripts from lysed cells, distorting cell-type annotation and biological interpretation. This project provides a **systematic benchmark workflow** for evaluating the performance of decontamination tools in single-cell RNA sequencing data.

Our study evaluated **7 decontamination tools** using datasets with well-defined ground truth, assessing their performance across three key dimensions: **accuracy, robustness, and subtype sensitivity**. The benchmark results reveal complementary strengths among tools: CellBender and scCDC excel in contamination estimation accuracy, SoupX demonstrates the highest robustness and sensitivity for rare subtypes, while scAR achieves the best balance between robustness and accuracy at the cell-type level for overall decontamination efficacy.

This work offers a practical guideline for context-dependent tool selection and highlights key directions for future algorithmic development in scRNA-seq decontamination.

## Repository Contents

This repository contains:

- **Benchmark datasets** – Ground-truth datasets used for evaluating decontamination performance
- **Benchmark tools** – Implementation scripts for running 7 decontamination methods (CellBender, scCDC, SoupX, scAR, etc.)
- **Evaluation metrics** – Code for computing accuracy, robustness, and subtype sensitivity metrics
- **Analysis scripts** – Scripts for reproducibility of all benchmark results and manuscript figures

## Reproducibility

All scripts and datasets analyzed in this study are archived here. The complete analysis pipeline, including tool implementation and metric calculation, is available in this repository.

For reproducibility of benchmark methods, metrics, and visualization, please refer to the code organization in this repository.

## 📁 Repository Structure

```
├── PBMC/                          # Fig2,3,S1,S2
│   ├── 1_Data_processing/         # Data processing scripts
│   └── 2_Figures/                  # Figure generation code
│
├── OSN/                            # Fig4,5,S3-S5
│   ├── analysis/                    # Downsampling & depth analysis
│   ├── decontamination/
│   │   └── generate_corrected_matrix.ipynb  # corrected_Matrix generation
│   │   └── decontamination methods          # decontaminated_Matrix generation
│   └── figures/
│       ├── utils/                    # Helper functions
│       ├── run_plot_contam_ratio.py  # Fig4C-E,H,I,L,M,5C-E,S4B,F
│       └── run_plot_corrected_mat.py # Fig4F,G,J,K,N,O,5F,G,S4C-E,G-I
```

## 🔗 Data Availability

All data used in this study are publicly available:

- **PBMC dataset**: Genome Sequence Archive (GSA-Human: HRA004605)
- **MOE dataset**: GEO (GSE185168, GSE157100)

## 🛠️ Installation and Dependencies

Setting up the environment for all 7 tools requires both Python and R. Please ensure you have `Python (v3.13)` and `R (v4.4.0)` installed.

### Python Environment

We recommend using `conda` or `venv` to manage the Python environment.

```bash
# Create and activate a new conda environment
conda create -n scrna_benchmark python=3.13
conda activate scrna_benchmark
# Install required Python packages
pip install anndata==0.10.1 h5py==3.14.0 matplotlib==3.10.6 scanpy==1.9.3 \
            scipy==1.13.0 seaborn==0.13.2 numpy==1.26.4 pandas==2.2.2 \
            sceasy==0.0.7

# Install decontamination tools

## CellBender
pip install cellbender==0.3.0

## scAR
conda install bioconda::scar

## CellClear (v0.0.3)

pip install CellClear

```
### R Environment

Run the following commands in an R console (v4.4.0) to install the required R packages and tools.

```r
# Install core packages from CRAN
install.packages(c('dplyr', 'tidyr', 'ggplot2', 'ggrepel', 'reticulate', 'rhdf5'))

# Install Seurat and its dependencies
install.packages('Seurat')

# Install BiocManager for Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c('DecontX', 'DropletUtils', 'SingleCellExperiment'))

# Install GitHub packages
## SoupX
install.packages('SoupX')

## scCDC
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ZJU-UoE-CCW-LAB/scCDC")

## FastCAR
devtools::install_git("https://git.web.rug.nl/P278949/FastCAR")

```
### Other Tools

The following command-line tools are required for certain preprocessing steps:

- **picard (3.4.0)**, **gatk4 (4.6.2.0)**, **samtools (1.22.1)**
- These can be installed via `conda`:
  ```bash
  conda install -c bioconda picard=3.4.0 gatk4=4.6.2.0 samtools=1.22.1
  ```

## 🛠️ Software Versions

**Python (v3.13)**: anndata(0.10.1), h5py(3.14.0), cellbender(0.3.0), scar(0.7.0), FastCAR(0.1.0), matplotlib(3.10.6), scanpy(1.9.3), scipy(1.13.0), seaborn(0.13.2), numpy(1.26.4), pandas(2.2.2), sceasy(0.0.7)

**R (v4.4.0)**: dplyr(1.1.4), tidyr(1.3.1), Seurat(5.3.0), DecontX(1.4.1), SoupX(1.6.2), scCDC(1.3), CellClear(0.0.3), ggplot2(3.3.6), seuratobject(5.0.0), reticulate(1.42.0), DropletUtils(1.22.0), rhdf5(2.46.1), ggrepel(0.9.6), ggbreak(0.1.6)

**Other tools**: picard(3.4.0), gatk4(4.6.2.0), samtools(1.22.1)

## 📧 Contact

For questions, please contact the corresponding author xujin7@mail.sysu.edu.cn.

## License

This project is licensed under the **Apache License 2.0** - see the [LICENSE](LICENSE) file for details.

All source code in this repository is covered under this license.

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)
