![](https://github.com/GeNeHetX/SCONE/blob/main/readme_logo.png)


**Date:** Decembre 2025  
**SCONE Authors :** [Camille PIGNOLET](https://github.com/CamillePignolet) [(camille.pignolet@inserm.fr)](mailto:camille.pignolet@inserm.fr) 
<br> From GeNeHetX Team [website](https://www.genehetx.fr/)<br>
If you use it in your work, please cite it as: SCONE © 2026 by Camille Pignolet.
<br>


## Overview
SCONE is a tool for analyzing Single Cell RNA-Seq data. It guides the user through quality control (QC), normalization, clustering, UMAP visualization, marker gene identification, and cluster annotation to assign each cluster a cell type label. There is also an automated cell-type annotation (ScType) and a pseudobulk correlation analysis<br><br>

## Prerequisites : 
SCONE requires the R language (at least version 4.0).<br>

If R is installed, you can launch the application directly via a command terminal or work on Rstudio.

- **REQUIRED** install R: [download](https://cran.r-project.org/)
- **REQUIRED** install the Rtools compiler : [download](https://cran.r-project.org/bin/windows/Rtools/)
- Not required, but if you want an IDE (integrated development environment), you can install Rstudio : [download](https://posit.co/download/rstudio-desktop/)
<br><br>

## Installation 

1 - Install all the tools of SCONE, run this command in a R terminal :
```R
    install.packages(c("shiny", "shinyjs", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools", "Seurat"))
```

2- Then, run this command to open the interactive application :
```R
   shiny::runGitHub("SCONE", "GeNeHetX")
```
<br>

## Input File Format

SCONE expects a single `.zip` file containing a standard 10X Genomics–style matrix:

- `matrix.mtx` → sparse gene expression matrix (genes × cells)  
- `genes.tsv` or `features.tsv` → gene annotations  
- `barcodes.tsv` → cell barcodes  

### Metadata

Additional metadata columns (e.g., `sampleID`, `condition`, `patient`, `batch`, etc.)  
must be added as extra columns in `barcodes.tsv`.

- The first column must contain the cell barcode.
- All additional columns are automatically imported as Seurat metadata.

Example structure:

```
barcode        sampleID    condition    patient
CellA-1        S1          Tumor        P01
CellB-1        S1          Tumor        P01
CellC-1        S2          Normal       P02
```

---

### Automated cell-type annotation (ScType)

Ianevski et al., 2022, Nature Communications  
https://doi.org/10.1038/s41467-022-28803-w  

ScType database and functions are obtained from:  
https://github.com/IanevskiAleksandr/sc-type  

**License**: ScType is distributed under the MIT License (see the original repository for details).  
Users are responsible for complying with the ScType license when using the database and scoring functions.
