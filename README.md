# BladderCancer_LineageDynamics

## Overview

This repository contains the computational analysis code for the research paper "Recurrence in the chemotherapy regimen of bladder carcinoma originates from quiescent epidermoid-like cells" .

## Abstract

Cancer relapse upon chemotherapy remains the primary challenge in cancer treatment. By integrating high-resolution CRISPR/Cas9-based evolving lineage tracing with single-cell RNA sequencing in mice, we tracked the evolution of bladder carcinoma in the presence or absence of chemotherapy. We found that a cell population exhibiting an epidermoid-like cell state was able to survive chemotherapy.

## Key Findings

- Epidermoid-like cells survive chemotherapy and initiate tumor relapse
- Plasticity of epidermoid-like cells shapes relapsed tumor phenotypes  
- IGFBP5 (Insulin-like growth factor-binding protein 5) is a key co-regulator of aggressive progression
- Targeting IGFBP5 can overcome chemoresistance in relapsed bladder carcinoma

## Repository Structure

```
BladderCancer_LineageDynamics/
├── analysis/                    # Main analysis notebooks and scripts
├── bulk/                        # Bulk sequencing data analysis
│   ├── RNA/                     # Bulk RNA-seq analysis
│   └── WES/                     # Whole exome sequencing analysis
├── config/                      # Configuration files and references
├── lineage/                     # Lineage tracing analysis
│   ├── upstream/                # Raw lineage data processing
│   └── downstream/              # Phylogenetic reconstruction and analysis
├── notebooks/                   # Jupyter notebooks for analysis
├── resources/                   # Additional resources
├── results/                     # Analysis results
├── sc_multiomics/              # Single-cell analysis
└── scripts/                    # Utility scripts and functions
```

## Data Types

This study integrates multiple data modalities:

1. **CRISPR/Cas9 Lineage Tracing**: High-resolution evolving lineage tracing
2. **Single-cell RNA-seq (scRNA-seq)**: Transcriptomic profiling of individual cells
3. **Single-cell ATAC-seq (scATAC-seq)**: Chromatin accessibility profiling
4. **Bulk RNA-seq**: Population-level transcriptomic analysis
5. **Whole Exome Sequencing (WES)**: Genomic variant analysis

## Key Analysis Components

### 1. Lineage Tracing Analysis
- Target site processing and phylogenetic reconstruction using Cassiopeia
- Clonal expansion detection and fitness inference using Jungle
- Evolutionary trajectory analysis

### 2. Single-cell Multi-omics
- scRNA-seq processing: Quality control, normalization, and clustering
- scATAC-seq analysis: Peak calling and chromatin accessibility
- Multi-modal integration

### 3. Comparative Analysis
- Treatment vs. control comparisons
- Primary vs. relapsed tumor analysis
- Cell state identification and characterization

## Sample Information

The study includes samples from:
- **MGHU3-LT**: Engineered MGHU3 bladder cancer cell line with lineage tracing
- **RT4-LT**: Engineered RT4 bladder cancer cell line with lineage tracing
- **Treatment conditions**: Vehicle control vs. chemotherapy treatment
- **Time points**: Primary tumors and relapsed tumors after treatment

## Usage

1. **Data preprocessing**: Process raw sequencing data and lineage information
2. **Quality control**: Filter cells and features based on quality metrics
3. **Lineage reconstruction**: Build phylogenetic trees from target site data
4. **Single-cell analysis**: Perform clustering, trajectory, and differential expression analysis
5. **Integration**: Combine lineage and expression information for evolutionary analysis

## Citation

If you use this code, please cite:

```
Zhou, D., Li, Y., Yan, Z. et al. Recurrence in the chemotherapy regimen of bladder carcinoma originates from quiescent epidermoid-like cells. Under review.
```

## Contact

For questions about the analysis code, please contact the corresponding authors.
