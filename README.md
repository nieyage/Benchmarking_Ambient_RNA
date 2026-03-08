# Systematic Benchmarking of Ambient RNA Decontamination Tools

This repository contains code for the manuscript:  
**"Systematic Benchmarking of Ambient RNA Decontamination Tools to Advance Precision in Single-Cell Transcriptomic Analysis"**

## 📁 Repository Structure

```
├── PBMC/                          # Fig2,3,S1,S2
│   ├── 1_Data_processing/         # Data processing scripts
│   └── 2_Figures/                  # Figure generation code
│
├── OSN/                            # Fig4,5,S3-S5
│   ├── analysis/                    # Downsampling & depth analysis
│   ├── decontamination/
│   │   └── generate_corrected_matrix.ipynb  # Matrix generation
│   └── figures/
│       ├── utils/                    # Helper functions
│       ├── run_plot_contam_ratio.py  # Fig4C-E,H,I,L,M,5C-E,S4B,F
│       └── run_plot_corrected_mat.py # Fig4F,G,J,K,N,O,5F,G,S4C-E,G-I
```

## 🔗 Data Availability

All data used in this study are publicly available:

- **PBMC dataset**: Genome Sequence Archive (GSA-Human: HRA004605)
- **MOE dataset**: GEO (GSE185168, GSE157100)

## 🛠️ Software Versions

**Python (v3.13)**: anndata(0.10.1), h5py(3.14.0), cellbender(0.3.0), scar(0.7.0), FastCAR(0.1.0), matplotlib(3.10.6), scanpy(1.9.3), scipy(1.13.0), seaborn(0.13.2), numpy(1.26.4), pandas(2.2.2), sceasy(0.0.7)

**R (v4.4.0)**: dplyr(1.1.4), tidyr(1.3.1), Seurat(5.3.0), DecontX(1.4.1), SoupX(1.6.2), scCDC(1.3), CellClear(0.0.3), ggplot2(3.3.6), seuratobject(5.0.0), reticulate(1.42.0), DropletUtils(1.22.0), rhdf5(2.46.1), ggrepel(0.9.6), ggbreak(0.1.6)

**Other tools**: picard(3.4.0), gatk4(4.6.2.0), samtools(1.22.1)


## 📧 Contact

For questions, please contact the corresponding author.
