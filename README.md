# Differential Expression Analysis 

## Description
A personal project of exploring a differential analysis of gene expression workflow by following along a practical guide from a conference workshop. 
The dataset analyzed contains microarray data of genes and different experimental conditions in hTNFTg mouse model of inflammatory polyarthritis [Karagianni et al. 2019](https://doi.org/10.1371/journal.pcbi.1006933).

Part of the project is an implementation of the same workflow using python libraries as an endeavor to potentially create an application that takes in any .tsv file uploaded by the user, runs the analysis automatically and displays the results in a Streamlit dashboard.

## Main Repo Files
### - DEG.qmd:

  Quarto markdown file containing the differential expression analysis workflow in R script code blocks, as well as personal notes on the pipeline's theoretical background.
  
  Contents include:
  - Exploratory Data Analysis
  - Data Distribution
  - Dimensionality Reduction Algorithms (UMAP & PCA)
  - Statistical Analysis (ANOVA & post-hoc tests)
  - Volcano Plot
  - Hierarchical Clustering
  - Functional Enrichment Analysis

  Libraries:

  `preprocessCore` 
  `umap` 
  `ggplot2` 
  `multcomp` 
  `gplots`
  `factoextra`
  `cowplot`
  `RColorBrewer`
  `kableExtra`
  `randomForest`
  `dplyr`
  `plotly`
  `gprofiler2` 
  `caret`
  
---

### - DEG.html:
  
  Rendered DEG.qmd file in HTML format.
  
---

### - DEG.py:
  
  Pending python implementation of the differential expression analysis workflow [_with a Streamlit application potential*_] exploring the following python libraries:
  `streamlit`
  `pandas`
  `seaborn`
  `numpy`
  `umap`
  `scipy`
  `statsmodels`
  `sklearn`
  `bokeh`
  `matplotlib`
  `sys`

  *_**Note**: Needs adjustments to take in any .tsv file!_
