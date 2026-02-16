![](https://github.com/GeNeHetX/SCONE/blob/main/readme_logo.png)


**Date:** Decembre 2025  
**SCONE Authors :** [Camille PIGNOLET](https://github.com/CamillePignolet) [(camille.pignolet@inserm.fr)](mailto:camille.pignolet@inserm.fr) 
<br> From GeNeHetX Team [website](https://www.genehetx.fr/)<br>
If you use it in your work, please cite it as: SCONE © 2026 by Camille Pignolet.
<br>


## Overview
SCONE is a tool for analyzing Single Cell RNA-Seq data. It guides the user through quality control (QC), normalization, clustering, UMAP visualization, marker gene identification, and cluster annotation to assign each cluster a cell type label.<br><br>

## Prerequisites : 
SCONE requires the R language (at least version 4.0).<br>

If R is installed, you can launch the application directly via a command terminal or work on Rstudio.

- **REQUIRED** install R: [download](https://cran.r-project.org/)
- **REQUIRED** install the Rtools compiler : [download](https://cran.r-project.org/bin/windows/Rtools/)
- Not required, but if you want an IDE (integrated development environment), you can install Rstudio : [download](https://posit.co/download/rstudio-desktop/)
<br>

## Installation 

1 - Install all the tools of SCONE, run this command in a R terminal :
```R
    install.packages(c("shiny", "shinyjs", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools", "Seurat"))
```

2- Then, run this command to open the interactive application :
```R
   shiny::runGitHub("SCONE", "GeNeHetX")
```
