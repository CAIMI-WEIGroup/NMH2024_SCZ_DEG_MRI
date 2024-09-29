# NMH2024_SCZ_DEG_MRI
This repository contains scripts and data used in Cui et al., (2024).

## Data
DESeq2_deg_results_all.csv includes differentiated gene expression levels for all genes between SCZ and HC.

data_all_lausanne120.mat includes tabulated imaging metrics in SCZ and HC.

## Scripts
To perform DESeq2 analysis on differential gene expression in SCZ, use R script:
```
DESeq2analysis.R
```

To perform DESeq2 analysis on differential transcript expression in SCZ, use R script:
```
DESeq2isoform.R
```

To perform gene-set enrichment analysis, use MATLAB scripts:
```
scripts_enrichment.m
```

To perform MAGMA gene property analysis, use BASH scripts:
```
magma_gene_property.sh
```

To perform PLS correlation analysis between gene expression and GMV, use MATLAB scripts:
```
scripts_assoc_morph.m
```

To perform PLS correlation analysis between gene expression, GMV, and connectivity, use MATLAB scripts:
```
scripts_assoc_edge_FC_PLS.m
scripts_assoc_edge_SC_PLS.m
```

To perform PLS correlation analysis between gene expression, GMV, and cognitions, use MATLAB scripts:
```
scripts_assoc_clinical.m
```

## if you use any data or script from this repository, please cite:
Cui, LB., Zhao, SW., Zhang, YH. et al. Associated transcriptional, brain and clinical variations in schizophrenia. Nat. Mental Health (2024). https://doi.org/10.1038/s44220-024-00306-1