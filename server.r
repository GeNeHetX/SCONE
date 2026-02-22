### ---- GLOBAL OPTIONS ----
options(timeout = 1000)
options(download.file.method = "libcurl")

cran_packages <- c(
  "shiny", "shinydashboard", "shinycssloaders", "magrittr", "shinyjs",
  "Seurat", "Matrix", "utils", "ggplot2", "gridExtra", "jsonlite", "DT", 
  "HGNChelper", "igraph", "ggraph", "scCustomize", "dplyr", 
  "circlize", "vegan"
)
cran_new <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if(length(cran_new)) {
  install.packages(cran_new, dependencies = TRUE)
}
if(!"grr" %in% installed.packages()[, "Package"]) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
    repos = NULL,
    type = "source"
  )
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("ComplexHeatmap")

bioc_new <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if(length(bioc_new)) {
  BiocManager::install(bioc_new, ask = FALSE, update = FALSE)
}

github_packages <- c(
  "GeNeHetX/CancerRNASig",
  "bnprks/BPCells/r",
  "cole-trapnell-lab/monocle3"
)

if(length(github_packages)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  
  for(repo in github_packages) {
    remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
  }
}


## LOAD LIBRARIES
library(shiny)
library(shinydashboard)
library(Seurat)
library(Matrix)
library(utils)
library(ggplot2)
library(shinycssloaders)
library(magrittr)
library(shinyjs)
library(gridExtra)
library(jsonlite)
library(DT)
library(HGNChelper)
library(igraph)
library(ggraph)
library(scCustomize)
library(dplyr)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(CancerRNASig)
library(viridis)
library(colorspace)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

options(shiny.maxRequestSize = 100*1024^3)
DT::datatable(data.frame(), extensions = "Buttons")


server <- function(input, output, session) {
  
  rv <- reactiveValues(seu = NULL, subset_cells = NULL, meta_choices = NULL, log = character(), files = NULL)
  
  append_log <- function(txt) {
    rv$log <- c(rv$log, txt)
  }
  output$tutorial_content <- renderUI({
          HTML("
          <h2>Single-Cell RNA-Seq Analysis with Seurat</h2>

          <p>This tutorial explains step-by-step how to analyze single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. Each step includes explanations so beginners can understand what is happening.</p>

          <hr>

          <h3>1. Create a Seurat Object</h3>
          <p>First, we create a Seurat object. This object stores your expression data and metadata and is the foundation for downstream analysis.</p>
          <pre><code class='r'>
          library(Seurat)
          # 'expressionMatrix' is your raw count matrix (genes x cells)
          seu <- CreateSeuratObject(counts = expressionMatrix)
          </code></pre>

          <hr>

          <h3>2. Quality Control (QC)</h3>
          <p>QC removes poor-quality cells. Cells with too few or too many genes, or high mitochondrial content, are filtered out.</p>
          <pre><code class='r'>
          seu <- subset(seu, subset = nFeature_RNA &gt; 200 &amp;&amp; nFeature_RNA &lt; 5000 &amp;&amp; percent.mt &lt; 10)
          </code></pre>

          <hr>

          <h3>3. Normalize the Data</h3>
          <p>Normalization adjusts for differences in sequencing depth between cells.</p>
          <pre><code class='r'>
          seu <- NormalizeData(seu)
          </code></pre>

          <hr>

          <h3>4. Identify Highly Variable Genes</h3>
          <p>We select genes with high variation across cells; these are most informative for clustering.</p>
          <pre><code class='r'>
          seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
          </code></pre>

          <hr>

          <h3>5. Scale the Data</h3>
          <p>Scaling centers and scales gene expression values; necessary for PCA and clustering.</p>
          <pre><code class='r'>
          seu <- ScaleData(seu)
          </code></pre>

          <hr>

          <h3>6. Perform PCA (Dimensionality Reduction)</h3>
          <p>PCA summarizes variation into principal components (PCs) capturing the main patterns in the data.</p>
          <pre><code class='r'>
          seu <- RunPCA(seu, npcs = 20)
          </code></pre>

          <hr>

          <h3>7. Clustering</h3>
          <p>Group similar cells. Seurat supports multiple clustering methods:</p>
          <ul>
          <li><strong>Leiden clustering</strong> detects communities on a nearest-neighbor graph (recommended).</li>
          <li><strong>Louvain clustering</strong> is an alternative modularity-based method.</li>
          <li><strong>K-means clustering</strong> partitions cells into a pre-defined number of clusters.</li>
          </ul>
          <pre><code class='r'>
          # Compute nearest neighbors graph
          seu <- FindNeighbors(seu, dims = 1:20)

          # Leiden clustering
          seu <- FindClusters(seu, resolution = 0.5, algorithm = 4)  # algorithm = 4 → Leiden

          # Optional Louvain clustering
          # seu <- FindClusters(seu, resolution = 0.5, algorithm = 1)

          # K-means clustering (optional)
          pca_emb <- Embeddings(seu, 'pca')[, 1:20]
          km <- kmeans(pca_emb, centers = 10)
          seu$kmeans_clusters <- as.factor(km$cluster)
          Idents(seu) <- seu$kmeans_clusters
          </code></pre>

          <hr>

          <h3>8. Visualize with UMAP</h3>
          <p>UMAP projects cells into 2D space. Points represent cells; clusters are visually identifiable.</p>
          <pre><code class='r'>
          seu <- RunUMAP(seu, dims = 1:20)
          DimPlot(seu, reduction = 'umap', label = TRUE)
          </code></pre>

          <hr>

          <h3>9. Find marker genes per cluster</h3>
          <p>Identify genes highly expressed in each cluster for annotation or pathway analysis.</p>
          <pre><code class='r'>
          markers <- FindAllMarkers(
              seu,
              only.pos = TRUE,
              min.pct = 0.25,
              logfc.threshold = 0.25
          )
          head(markers)
          </code></pre>

          <hr>

          <h3>10. ScType annotation</h3>
          <p>
          <ScType> is a computational method for automated marker gene selection based only on scRNA-seq data. 
          The open-source portal <a href=\"http://sctype.app\" target=\"_blank\">http://sctype.app</a> provides an interactive web implementation. 
          Reference: <a href=\"https://doi.org/10.1038/s41467-022-28803-w\" target=\"_blank\">Ianevski et al., 2022, Nature Communications</a>.
          </p>
          <pre><code class='r'>
          # Ensure cluster identities are set
          Idents(seu) <- seu$seurat_clusters

          # Load ScType database
          db_url <- 'https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx'
          db_local <- file.path(tempdir(), 'ScTypeDB_full.xlsx')
          if (!file.exists(db_local)) download.file(db_url, db_local, mode = 'wb')

          gs_list <- gene_sets_prepare(db_local, tissue = 'Pancreas')  # change tissue if needed

          # Prepare scaled RNA matrix
          if(is.null(seu[['RNA']]$scale.data)) {
              seu <- ScaleData(seu, features = rownames(seu))
          }
          scRNAseqData_scaled <- as.matrix(seu[['RNA']]$scale.data)

          # Run ScType scoring
          es.max <- sctype_score(
              scRNAseqData = scRNAseqData_scaled,
              scaled = TRUE,
              gs = gs_list$gs_positive,
              gs2 = gs_list$gs_negative
          )

          # Aggregate scores per cluster
          clusters <- unique(as.character(Idents(seu)))
          cL_results <- do.call('rbind', lapply(clusters, function(cl) {
              cells_cl <- colnames(seu)[Idents(seu) == cl]
              scores_cluster <- rowSums(es.max[, cells_cl, drop = FALSE])
              scores_cluster <- sort(scores_cluster, decreasing = TRUE)
              head(data.frame(
                  cluster = cl,
                  type = names(scores_cluster),
                  scores = scores_cluster,
                  ncells = length(cells_cl),
                  stringsAsFactors = FALSE
              ), 10)
          }))

          # Pivot wide and assign label_max
          score_table <- tidyr::pivot_wider(
              cL_results,
              names_from = type,
              values_from = scores,
              values_fill = 0,
              values_fn = list(scores = max)
          )
          score_cols <- setdiff(colnames(score_table), c('cluster', 'ncells'))
          score_matrix <- as.matrix(score_table[, score_cols])
          max_scores <- apply(score_matrix, 1, max)
          max_labels <- score_cols[max.col(score_matrix, ties.method = 'first')]
          score_table$label_max <- ifelse(max_scores < score_table$ncells / 4, 'Unknown', max_labels)

          # Save results
          rv$sctype_annotation <- score_table
          </code></pre>
          <hr>

          <h3>11. Pseudobulk Correlation Analysis</h3>

          <p>
          Pseudobulk analysis aggregates gene expression across cells belonging to the same cluster.
          This reduces single-cell noise and produces bulk-like profiles that are more stable
          for correlation and similarity analysis.
          </p>

          <pre><code class='r'>
          # Ensure cluster identities are set
          Idents(seu) <- seu$seurat_clusters

          # Aggregate expression per cluster (pseudobulk)
          bulk <- AggregateExpression(
              seu,
              group.by = 'seurat_clusters',
              return.seurat = TRUE
          )

          # Extract normalized expression matrix
          mat <- GetAssayData(bulk, assay = 'RNA', slot = 'data')

          # Select top 1000 most variable genes
          top1000 <- names(sort(apply(mat, 1, sd), decreasing = TRUE))[1:1000]

          # Center expression values (improves correlation structure)
          centered_mat <- scale(mat[top1000, ], scale = FALSE)

          # Compute Pearson correlation between clusters
          mat_corr <- cor(centered_mat, method = 'pearson')

          # Visualize correlation heatmap
          library(ComplexHeatmap)
          library(circlize)

          Heatmap(
              mat_corr,
              name = 'Correlation',
              col = colorRamp2(c(0.3, 0.6, 1), c('blue', 'white', 'red')),
              clustering_method_rows = 'ward.D2',
              clustering_method_columns = 'ward.D2'
          )
          </code></pre>

          <p>
          Clusters with high correlation values share similar transcriptional programs,
          while low correlation indicates distinct biological states.
          </p>

          <hr>
          

          <h3>Notes for Beginners</h3>
          <ul>
          <li>Always inspect your data after each step.</li>
          <li>Quality control is crucial: low-quality or doublet cells can bias results.</li>
          <li>Variable features should capture biologically meaningful differences.</li>
          <li>PCA reduces noise and speeds up clustering.</li>
          <li>UMAP or t-SNE plots help visually explore clusters.</li>
          <li>After clustering, marker genes guide interpretation and ScType annotation.</li>
          <li>Aggregating cells (Pseudobulk) based on clustering within each sample allows generation of cell-type–specific pseudobulk profiles, providing more robust and less noisy expression signals for downstream comparison.</li>
          </ul>
          ")

    })

  observeEvent(input$show_tutorial, {
    shinyjs::toggle("tutorial_panel")

  })

# -------------------------
# LOAD FILES AND CREATE SEURAT OBJECT AUTOMATICALLY
# -------------------------
observeEvent(input$input_zip, {
  req(input$input_zip)
  append_log("Starting Seurat object creation...")

  tmpdir <- tempdir()

  tryCatch({
    append_log("Unzipping input...")
    unzip(input$input_zip$datapath, exdir = tmpdir)

    # Find matrix / features / barcodes
    mtx_file <- list.files(tmpdir, pattern = "matrix.*\\.mtx$", full.names = TRUE, ignore.case = TRUE)[1]
    genes_file <- list.files(tmpdir, pattern = "(^genes|features).*\\.tsv$", full.names = TRUE, ignore.case = TRUE)[1]
    barcodes_file <- list.files(tmpdir, pattern = "barcodes.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)[1]

    append_log("Reading 10X matrix...")
    expressionMatrix <- ReadMtx(
      mtx = mtx_file,
      cells = barcodes_file,
      features = genes_file,
      feature.column = 1,
      skip.cell = 1,
      cell.sep = "\t"
    )

    append_log("Creating Seurat object...")
    seu <- CreateSeuratObject(counts = expressionMatrix)

    # -------------------------
    # Read barcode metadata (once)
    # -------------------------
    barcode_df <- read.table(barcodes_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    if (ncol(barcode_df) > 1) {
      colnames(barcode_df)[1] <- "barcode"

      # garder uniquement les barcodes présents dans Seurat
      barcode_df <- barcode_df[barcode_df$barcode %in% colnames(seu), , drop = FALSE]

      # mettre dans le même ordre que Seurat
      barcode_df <- barcode_df[match(colnames(seu), barcode_df$barcode), , drop = FALSE]

      if (ncol(barcode_df) > 1) {
        meta_to_add <- barcode_df[, -1, drop = FALSE]  # tout sauf la 1ère colonne
        rownames(meta_to_add) <- barcode_df$barcode

        # retirer les colonnes déjà présentes pour éviter les doublons
        meta_to_add <- meta_to_add[, !colnames(meta_to_add) %in% colnames(seu@meta.data), drop = FALSE]

        if(ncol(meta_to_add) > 0) {
          seu <- AddMetaData(seu, metadata = meta_to_add)
        }
      }
    }

    append_log("Seurat object created with additional barcode metadata if available.")

    # -------------------------
    # QC metrics
    # -------------------------
    mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes)

    if(length(mito_genes) > 0){
      seu$percent_mito <- Matrix::colSums(seu[mito_genes, ]) / Matrix::colSums(seu) * 100
    } else {
      seu$percent_mito <- 0
    }

    # Stocker Seurat object
    rv$seu <- seu
    rv$files <- list(mtx = mtx_file, genes = genes_file, barcodes = barcodes_file)
    append_log("Seurat object created.")

    # -------------------------
    # QC plots
    # -------------------------
    observe({
      req(rv$seu)
      
      updateSelectInput(
        session,
        "qc_split_meta",
        choices = colnames(rv$seu@meta.data),
        selected = colnames(rv$seu@meta.data)[1]
      )
    })


    output$qc_violin <- renderPlot({
      req(rv$seu)
      df <- rv$seu@meta.data
      split_var <- input$qc_split_meta

      p1 <- ggplot(df,  aes(x = .data[[split_var]],y = nFeature_RNA)) +
        geom_violin(fill = "lightblue") +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
        geom_hline(yintercept = c(600, 5500), linetype = "dashed", color = "red") +
        labs(y = "nFeature_RNA", x = "", title = 'Gene per Cell/Nucleus') +
        theme_classic()

      p2 <- ggplot(df,  aes(x = .data[[split_var]],y = nCount_RNA)) +
        geom_violin(fill = "lightgreen") +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
        geom_hline(yintercept = c(1200, 45000), linetype = "dashed", color = "red") +
        labs(y = "nCount_RNA", x = "", title= 'UMI per Cell/Nucleus') +
        theme_classic()

      # p3 <- QC_Histogram(seurat_object = rv$seu, features = "percent_mito", low_cutoff = 15)
      p3 <- ggplot(df,  aes(x = .data[[split_var]], y = percent.mt)) + 
       geom_violin(fill = "salmon") + 
       geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) + 
       geom_hline(yintercept = 20, linetype = "dashed", color = "red") + 
       labs(y = "percent.mt", x = "", title='Mito Gene % per Cell/Nucleus') + theme_classic()

      df <- df %>% mutate(complexity = nFeature_RNA / nCount_RNA)
      p4 <- ggplot(df,  aes(x = .data[[split_var]],y = complexity)) +
        geom_violin(fill = "lightgrey") +
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
        geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
        labs(y = "Complexity (nFeature/nCount)", x = "", title="Cell complexity") +
        theme_classic()

      grid.arrange(p1, p2, p3, p4, ncol = 4)
    })

    output$qc_scatter1 <- renderPlot({
      QC_Plot_UMIvsGene(seurat_object = rv$seu, low_cutoff_gene = 600, high_cutoff_gene = 5500,
                        low_cutoff_UMI = 500, high_cutoff_UMI = 50000)
    })

    output$qc_scatter2 <- renderPlot({
      QC_Plot_GenevsFeature(seurat_object = rv$seu, feature1 = "percent_mito",
                            low_cutoff_gene = 600, high_cutoff_gene = 5500, high_cutoff_feature = 20)
    })

    output$qc_scatter3 <- renderPlot({
      QC_Plot_UMIvsGene(seurat_object = rv$seu, meta_gradient_name = "percent_mito",
                        low_cutoff_gene = 600, high_cutoff_gene = 5500, high_cutoff_UMI = 45000)
    })

    output$qc_scatter4 <- renderPlot({
      QC_Plot_UMIvsGene(seurat_object = rv$seu, meta_gradient_name = "percent_mito",
                        low_cutoff_gene = 600, high_cutoff_gene = 5500, high_cutoff_UMI = 45000,
                        meta_gradient_low_cutoff = 20)
    })

    # -------------------------
    # SHOW SUMMARY / OVERVIEW
    # -------------------------
    output$seu_summary <- renderPrint({
      req(rv$seu)
      rv$seu
    })

    # Stocker colonnes existantes pour la table DT
    rv$meta_base <- colnames(rv$seu@meta.data)

   output$meta_table <- DT::renderDataTable({
        req(rv$seu)
        
        # toutes les colonnes dynamiques + base
        all_cols <- c(rv$meta_base, grep("sampleID|^published|^(louvain_|kmeans_|leiden_)", 
                                        colnames(rv$seu@meta.data), value = TRUE))
        
        base_names <- sub("\\.\\d+$", "", all_cols)
        keep_idx <- match(unique(base_names), base_names)
        meta_to_show <- rv$seu@meta.data[, all_cols[keep_idx], drop = FALSE]
        DT::datatable(meta_to_show, options = list(pageLength = 10, scrollX = TRUE))
    })


      # -------------------------
      # SUBSET UI (threshold inputs)
      # -------------------------
      output$subset_ui <- renderUI({
        req(rv$seu)
        tagList(
          numericInput("nFeature_min", "Min nFeature_RNA", value = 200),
          numericInput("nFeature_max", "Max nFeature_RNA", value = 5000),
          numericInput("percent_mt_max", "Max percent.mt", value = 10)
        )
      })

      subset_done <- reactiveVal(FALSE)
      rv$seu_original <- seu   # copie intacte
      rv$seu <- seu            # objet de travail
      rv$subset_cells <- NULL


      observeEvent(input$apply_subset, {
          subset_done(TRUE)
          req(rv$seu_original)
          
          # Appliquer le subset
          cells_keep <- WhichCells(rv$seu_original, expression = 
                                    nFeature_RNA >= input$nFeature_min &
                                    nFeature_RNA <= input$nFeature_max &
                                    percent.mt <= input$percent_mt_max)
          
          rv$subset_cells <- subset(rv$seu_original, cells = cells_keep)
          rv$seu <- rv$subset_cells

          # Log
          append_log(paste("Subset applied:", length(cells_keep), "cells kept"))
          
          # Rendu du Seurat object subset avec nb de cellules et features
          output$subset_summary <- renderPrint({
            req(rv$subset_cells)
            cat("Subset Seurat object summary:\n")
            cat("Number of cells:", ncol(rv$subset_cells), "\n")
            cat("Number of features:", nrow(rv$subset_cells), "\n")
          })
        })
        
        seu <- if (!is.null(rv$subset_cells)) rv$subset_cells else rv$seu_original
            
        output$subset_message <- renderUI({
            req(subset_done())

            tags$span(
              "→ Run analysis again with you subsetted Seurat object",
              style = "margin-left:10px; font-style:italic; color:#666;"
            )
          })


    }, error = function(e){
      append_log(paste("ERROR:", e$message))
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
    
  })

  
  # -------------------------
  # RUN ANALYSIS (subset + normalize + PCA + clustering + UMAP)
  # -------------------------
  observeEvent(input$run_analysis, {
    req(rv$seu)
    shinyjs::disable("run_analysis")
    on.exit({ shinyjs::enable("run_analysis") })

    tryCatch({
      seu <- rv$seu

      # Apply subset if thresholds defined
      if (!is.null(input$nFeature_min) && !is.null(input$nFeature_max) && !is.null(input$percent_mt_max)) {
        cells_keep <- WhichCells(seu, expression = 
                                   nFeature_RNA >= input$nFeature_min &
                                   nFeature_RNA <= input$nFeature_max &
                                   percent.mt <= input$percent_mt_max)
        seu <- subset(seu, cells = cells_keep)
        append_log(paste("Subset applied:", length(cells_keep), "cells kept"))
      } else {
        append_log("No subset thresholds applied; running on all cells")
      }

      # Normalize, HVG, Scale, PCA
      append_log("Normalizing data...")
      seu <- NormalizeData(seu, verbose = FALSE)
      append_log(paste("Finding variable features (n=", input$n_hvg, ")", sep=""))
      seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=input$n_hvg, verbose=FALSE)
      append_log("Scaling data...")
      seu <- ScaleData(seu, verbose = FALSE)
      append_log(paste("Running PCA (npcs=", input$n_pcs, ")", sep=""))
      seu <- RunPCA(seu, npcs = input$n_pcs, verbose = FALSE)

      # Clustering
      if (input$clust_method == "louvain") {
        cluster_name <- paste0("louvain_n", input$knn, "_res", input$resolution)
        append_log(paste("Running Louvain clustering:", cluster_name))
        seu <- FindNeighbors(seu, dims = 1:input$n_pcs, k.param = input$knn)
        seu <- FindClusters(seu, resolution = input$resolution)
        seu@meta.data[[cluster_name]] <- Idents(seu)
        Idents(seu) <- seu@meta.data[[cluster_name]]

      } else if (input$clust_method == "leiden") {
        cluster_name <- paste0("leiden_n", input$leiden_k, "_res", input$resolution)
        append_log(paste("Running Leiden clustering:", cluster_name))
        seu <- FindNeighbors(seu, dims = 1:input$n_pcs, k.param = input$leiden_k)
        seu <- FindClusters(seu, resolution = input$resolution, algorithm = 4) # 4 = Leiden
        seu@meta.data[[cluster_name]] <- Idents(seu)
        Idents(seu) <- seu@meta.data[[cluster_name]]

      } else { # kmeans
        cluster_name <- paste0("kmeans_n", input$kmeans_centers, "_res", input$resolution)
        append_log(paste("Running K-means clustering:", cluster_name))
        pca_emb <- Embeddings(seu, "pca")[, 1:input$n_pcs, drop = FALSE]
        km <- kmeans(pca_emb, centers = input$kmeans_centers)
        seu@meta.data[[cluster_name]] <- as.factor(km$cluster)
        Idents(seu) <- seu@meta.data[[cluster_name]]
      }

      # UMAP
      if (!"umap" %in% names(seu@reductions)) {
        append_log("Running UMAP...")
        seu <- RunUMAP(seu, dims = 1:input$n_pcs, verbose = FALSE)
      }

      rv$seu <- seu
      rv$meta_choices <- grep("sampleID|^published|^(louvain_|kmeans_|leiden_)", names(seu@meta.data), value = TRUE)
      append_log(paste("Clustering done:", cluster_name))
      append_log("Analysis completed.")

      updateSelectInput(session, "meta_color_by", choices = rv$meta_choices, selected = cluster_name)
      output$metadata_selector <- renderUI({
        selectInput("meta_choice", "Select metadata for plotting", choices = rv$meta_choices)
      })

    }, error = function(e){
      append_log(paste("ERROR:", e$message))
      showNotification(paste("Error:", e$message), type="error", duration = 8)
    })
  })

  # -------------------------
  # UMAP plots
  # -------------------------
  output$umapPlot <- renderPlot({
    req(rv$seu)
    tryCatch({
      DimPlot(rv$seu, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP - clusters")
    }, error = function(e) ggplot() + ggtitle("Unable to plot UMAP"))
  })

  output$umapMetaPlot <- renderPlot({
    req(rv$seu)
    req(input$meta_color_by)
    if (!(input$meta_color_by %in% colnames(rv$seu@meta.data))) {
      ggplot() + ggtitle("Selected metadata not found")
    } else {
      DimPlot(rv$seu, reduction = "umap", label = TRUE, group.by = input$meta_color_by) +
        ggtitle(paste("UMAP colored by", input$meta_color_by))
    }
  })

  # -------------------------
  # Log output
  # -------------------------
  output$status <- renderText({
    if (length(rv$log) == 0) return("Idle")
    paste(rev(rv$log), collapse = "\n")
  })
 
  # -------------------------
  # UMAP plots
  # -------------------------
  # output$umapPlot <- renderPlot({
  #   req(rv$seu)
  #   tryCatch({
  #     DimPlot(rv$seu, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP - clusters")
  #   }, error = function(e) ggplot() + ggtitle("Unable to plot UMAP"))
  # })
    # UMAP pour Violin
  output$umapViolinPlot <- renderPlot({
    req(rv$seu)
    req(input$split_meta)
    DimPlot(rv$seu, reduction = "umap", group.by = input$split_meta, label = TRUE)
  })

  # UMAP pour DotPlot
  output$umapDotPlot <- renderPlot({
    req(rv$seu)
    req(input$dot_split_meta)
    DimPlot(rv$seu, reduction = "umap", group.by = input$dot_split_meta, label = TRUE)
  })

  
  output$umapMetaPlot <- renderPlot({
    req(rv$seu)
    req(input$meta_color_by)
    if (!(input$meta_color_by %in% colnames(rv$seu@meta.data))) {
      ggplot() + ggtitle("Selected metadata not found")
    } else {
      DimPlot(rv$seu, reduction = "umap", group.by = input$meta_color_by, label = TRUE) 
        # +ggtitle(paste("UMAP colored by", input$meta_color_by))
    }
  })
  
  # -------------------------
  # Violin
  # -------------------------
  gene_families <- reactive({
      path <- "markers_celltypes.json"
      req(file.exists(path))
      jsonlite::fromJSON(path)
    })

  output$violin_gene_ui <- renderUI({
      req(rv$seu)
      genes <- rownames(rv$seu)
      selectizeInput(
        inputId = "violin_gene",
        label = "Gene(s)",
        choices = genes,
        selected = NULL,
        multiple = TRUE,
        options = list(
          placeholder = 'Start typing gene...',
          maxOptions = 1000,
          create = FALSE
        )
      )
    })

    output$violin_gene_family_ui <- renderUI({
      req(gene_families())   # s’assure que le JSON est chargé
      selectInput(
        "violin_gene_family", 
        "Gene family (optional)", 
        choices = c("None" = "", names(gene_families())),
        selected = ""
      )
    })

    observe({
      req(rv$seu)
      updateSelectizeInput(session, "violin_gene", choices = rownames(rv$seu), server = TRUE)
    })

    # Split meta UI
    output$split_meta_ui <- renderUI({
      req(rv$meta_choices)
      selectInput("split_meta", "Split by (metadata)", choices = rv$meta_choices, selected = rv$meta_choices[1])
    })

    observeEvent(input$plot_violin, {
      req(rv$seu)

      genes_manual <- input$violin_gene
      genes_family <- if(input$violin_gene_family != "") gene_families()[[input$violin_gene_family]] else character(0)
      genes_to_plot <- unique(c(genes_manual, genes_family))
      genes_exist <- genes_to_plot[genes_to_plot %in% rownames(rv$seu)]
      
      if (length(genes_exist) == 0) {
        showNotification("None of the selected genes are found in the Seurat object", type="error")
        return()
      }

      output$violinPlot <- renderPlot({
        df <- rv$seu@meta.data
        expr_data <- GetAssayData(rv$seu, layer="data")  # Seurat v5

        plots <- lapply(genes_exist, function(gene) {
          df[[gene]] <- expr_data[gene, ]
          if (input$split_by && !is.null(input$split_meta) && input$split_meta %in% colnames(df)) {
            ggplot(df, aes_string(x=input$split_meta, y=gene)) +
              geom_violin(fill="lightblue") +
              geom_jitter(width=0.2, size=0.5) +
              theme_classic() +
              ggtitle(gene)
          } else {
            ggplot(df, aes_string(x="1", y=gene)) +
              geom_violin(fill="lightblue") +
              geom_jitter(width=0.2, size=0.5) +
              theme_classic() +
              ggtitle(gene)
          }
        })
        gridExtra::grid.arrange(grobs=plots, ncol=3)
      })
    })

  # -------------------------
  # DotPlot
  # -------------------------
  observe({
    req(rv$seu)
    updateSelectizeInput(session, "dot_genes", choices = rownames(rv$seu), server = TRUE )
  })

  output$dot_gene_family_ui <- renderUI({
      req(gene_families())  # JSON déjà chargé
      selectInput(
        "dot_gene_family", 
        "Gene family (optional)", 
        choices = c("None" = "", names(gene_families())),
        selected = ""
      )
    })
    
    output$dot_split_meta_ui <- renderUI({
        req(rv$meta_choices)

        selectInput(
          "dot_split_meta",
          "Group DotPlot by:",
          choices = rv$meta_choices,
          selected = "seurat_clusters"  # valeur par défaut
        )
      })


    observeEvent(input$plot_dot, {
        req(rv$seu)

        genes_manual <- input$dot_genes
        genes_family <- if(input$dot_gene_family != "")
          gene_families()[[input$dot_gene_family]]
        else character(0)

        genes_to_plot <- unique(c(genes_manual, genes_family))
        genes_exist <- genes_to_plot[genes_to_plot %in% rownames(rv$seu)]

        if (length(genes_exist) == 0) {
          showNotification(
            "None of the selected genes are found in the Seurat object",
            type="error"
          )
          return()
        }

        output$dotPlot <- renderPlot({
          group_var <- if (
            !is.null(input$dot_split_meta) &&
            input$dot_split_meta %in% colnames(rv$seu@meta.data)
          ) {
            input$dot_split_meta
          } else {
            "seurat_clusters"
          }

          DotPlot(rv$seu,
            features = genes_exist, 
            group.by = group_var,
            cols = c("blue", "red")
          ) + coord_flip()

        })

      })



    # -------------------------
    # FIND MARKERS
    # -------------------------

    # UI dynamique pour choisir metadata
    output$marker_split_meta_ui <- renderUI({
      req(rv$meta_choices)

      selectInput(
        "marker_split_meta",
        "Group cells by:",
        choices = rv$meta_choices,
        selected = "seurat_clusters"
      )
    })

    # UI dynamique pour choisir clusters selon le mode
    output$cluster_select_ui <- renderUI({
      req(rv$seu)
      req(input$marker_split_meta)
      clusters <- sort(unique(rv$seu@meta.data[[input$marker_split_meta]]))

      if (input$marker_mode == "vsrest") {

        selectInput(
          "cluster_a",
          "Cluster of interest:",
          choices = clusters,
          selected = clusters[1]
        )

      } else { # vscluster

        tagList(
          selectInput(
            "cluster_a",
            "Cluster A:",
            choices = clusters,
            selected = clusters[1]
          ),
          selectInput(
            "cluster_b",
            "Cluster B:",
            choices = clusters,
            selected = if (length(clusters) > 1) clusters[2] else clusters[1]
          )
        )
      }
    })

    output$umapMarkerPlot <- renderPlot({
      req(rv$seu)
      req(input$marker_split_meta)
      DimPlot(rv$seu, reduction = "umap", group.by = input$marker_split_meta, label = TRUE) 
    })

    # calcul des markers
    observeEvent(input$run_markers, {

      req(rv$seu)
      req(input$marker_split_meta)
      seu <- rv$seu

      # mettre la metadata choisie comme ident
      Idents(seu) <- seu[[input$marker_split_meta]][,1]

      mode <- input$marker_mode

      # Progress bar
      withProgress(message = "Finding markers...", {

        # 1 vs 1 : vérif clusters différents
        if (mode == "vscluster" && input$cluster_a == input$cluster_b) {
          showNotification("Please select two different clusters", type = "error")
          return()
        }

        # 1 vs 1
        if (mode == "vscluster") {

          req(input$cluster_a, input$cluster_b)

          res <- FindMarkers(
            seu,
            ident.1 = input$cluster_a,
            ident.2 = input$cluster_b,
            logfc.threshold = 0.25
          )

          if (nrow(res) > 0) {
            res$cluster <- paste(input$cluster_a, "vs", input$cluster_b)
          }

        } else { # vs rest

          req(input$cluster_a)

          res <- FindMarkers(
            seu,
            ident.1 = input$cluster_a,
            logfc.threshold = 0.25
          )

          if (nrow(res) > 0) {
            res$cluster <- input$cluster_a
          }

        } # fin mode

      }) # fin withProgress

      # ajouter colonne gene
      res$gene <- rownames(res)

      # protection si aucun marker
      if (nrow(res) == 0) {

        showNotification("No markers found with current thresholds", type = "warning")

        output$marker_table <- DT::renderDataTable({
          DT::datatable(data.frame())
        })

        output$top_marker_table <- DT::renderDataTable({
          DT::datatable(data.frame())
        })

        return()
      }

      # -------------------------
      # Full table
      # -------------------------
      output$marker_table <- DT::renderDataTable({
        DT::datatable(
          res,
          extensions = "Buttons",
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
              list(extend = "copy",  text = "Copy"),
              list(extend = "csv",   text = "CSV",   filename = "markers_table"),
              list(extend = "excel", text = "Excel", filename = "markers_table")
            )
          )
        )
      })


      # -------------------------
      # Top markers : juste la liste des gènes par cluster
      # -------------------------
      # -------------------------
      # Top markers : liste des gènes par cluster
      # -------------------------
      if (input$marker_mode == "vscluster") {
        
        # Cluster A : logFC positif → enrichi dans cluster A
        top_cluster_a <- res |>
          dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) |>
          dplyr::slice_max(avg_log2FC, n = input$top_n) |>
          dplyr::mutate(top_for = input$cluster_a) |>
          dplyr::select(top_for, gene, avg_log2FC)
        
        # Cluster B : logFC négatif → enrichi dans cluster B
        top_cluster_b <- res |>
          dplyr::filter(p_val_adj < 0.05 & avg_log2FC < 0) |>
          dplyr::slice_min(avg_log2FC, n = input$top_n) |>
          dplyr::mutate(top_for = input$cluster_b, avg_log2FC = abs(avg_log2FC)) |>
          dplyr::select(top_for, gene, avg_log2FC)
        
        top_markers <- dplyr::bind_rows(top_cluster_a, top_cluster_b)
        
      } else {
        
        # vs rest ou FindAllMarkers
        top_markers <- res |>
          dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) |>
          dplyr::group_by(cluster) |>
          dplyr::slice_max(avg_log2FC, n = input$top_n) |>
          dplyr::ungroup() |>
          dplyr::select(cluster, gene, avg_log2FC)
      }

      output$top_marker_table <- DT::renderDataTable({
        DT::datatable(
          top_markers,
          extensions = "Buttons",
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
              list(extend = "copy",  text = "Copy"),
              list(extend = "csv",   text = "CSV",   filename = "top_markers"),
              list(extend = "excel", text = "Excel", filename = "top_markers")
            )
          )
        )
      })




    }) 

    output$marker_table_title <- renderUI({
    HTML("
    <b>Columns explanation:</b><br>
    <ul>
      <li><b>avg_log2FC</b>: average log2 fold change of expression between the two groups</li>
      <li><b>pct.1</b>: fraction of cells expressing the gene in the cluster of interest (ident.1)</li>
      <li><b>pct.2</b>: fraction of cells expressing the gene in the other cluster (ident.2 or rest)</li>
      <li><b>p_val_adj</b>: adjusted p-value for differential expression</li>
      <li><b>cluster</b>: the cluster in which the gene is differential (cluster of interest)</li>
    </ul>
    ")
  })

  ##--------------------------
  # ScType Annotation Server
  ##--------------------------

  # UI dynamique pour choisir metadata
  output$type_split_meta_ui <- renderUI({
    req(rv$meta_choices)
    selectInput(
      "sctype_split_meta",
      "Group cells by (metadata)",
      choices = rv$meta_choices,
      selected = "seurat_clusters"
    )
  })

  # UMAP initial coloré par metadata choisie
  output$sctype_umap_plot <- renderPlot({

    req(rv$seu)
    req(input$sctype_split_meta)

    seu <- rv$seu
    cluster_col <- input$sctype_split_meta

    if (!(cluster_col %in% colnames(seu@meta.data))) {
      showNotification(
        paste("Column", cluster_col, "not found in Seurat object"),
        type = "error"
      )
      return()
    }

    Idents(seu) <- seu@meta.data[[cluster_col]]

    if (!is.null(rv$sctype_annotation)) {

      label_map <- setNames(rv$sctype_annotation$label_max, rv$sctype_annotation$cluster)
      seu$ScType_label <- as.character(seu@meta.data[[cluster_col]])
      for (cl in names(label_map)) {
        seu$ScType_label[seu@meta.data[[cluster_col]] == cl] <- label_map[cl]
      }
      DimPlot(seu, reduction = "umap", group.by = "ScType_label", label = TRUE
      ) +       theme(
              legend.position = "bottom",       # place la légende en bas
              legend.title = element_text(size = 12),
              legend.text  = element_text(size = 10),
              plot.title = element_text(hjust = 0.5)
    )
    } else {

      DimPlot(
        seu,
        reduction = "umap",
        group.by = cluster_col,
        label = TRUE
      )
    }
  })


  # ==============================
  # Lancer ScType (subset + full)
  # ==============================
  observeEvent(input$run_sctype_annotation, {

   req(rv$seu)
    req(input$sctype_split_meta)

    seu <- rv$seu
    cluster_col <- input$sctype_split_meta

    if (!(cluster_col %in% colnames(seu@meta.data))) {
      showNotification(
        paste("Column", cluster_col, "not found in Seurat object"),
        type = "error"
      )
      return()
    }

    Idents(seu) <- seu@meta.data[[cluster_col]]

    withProgress(message = "Running ScType annotation...", value = 0, {

    DefaultAssay(seu) <- "RNA"
    seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seu[["RNA"]])))

    append_log(paste("Seurat object", ifelse(seurat_package_v5, "v5", "v4"), "detected"))
    incProgress(0.05)

    if (is.null(seu[["RNA"]]$scale.data) || ncol(seu[["RNA"]]$scale.data) == 0) {
      append_log("scale.data absent ou vide — recalcul avec ScaleData()")
      seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
    }

    scRNAseqData_scaled <- as.matrix(seu[["RNA"]]$scale.data)

    if (nrow(scRNAseqData_scaled) == 0 || ncol(scRNAseqData_scaled) == 0) {
      stop("Erreur : scale.data vide même après ScaleData()")
    }
    incProgress(0.1)

    db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    db_local <- file.path(tempdir(), "ScTypeDB_full.xlsx")
    if (!file.exists(db_local)) download.file(db_url, db_local, mode = "wb")
    gs_list <- gene_sets_prepare(db_local, input$sctype_tissue)
    incProgress(0.2)

    es.max <- sctype_score(
      scRNAseqData = scRNAseqData_scaled,
      scaled = TRUE,
      gs = gs_list$gs_positive,
      gs2 = gs_list$gs_negative
    )
    # Ne garder que les cellules du subset
    es.max <- es.max[, colnames(seu), drop = FALSE]
    incProgress(0.5)

    clusters <- unique(as.character(Idents(seu)))

    cL_results <- do.call("rbind", lapply(clusters, function(cl) {
      cells_cl <- colnames(seu)[Idents(seu) == cl]
      if (length(cells_cl) == 0) return(NULL)
      scores_cluster <- rowSums(es.max[, cells_cl, drop = FALSE])
      scores_cluster <- sort(scores_cluster, decreasing = TRUE)
      head(data.frame(
        cluster = cl,
        type = names(scores_cluster),
        scores = scores_cluster,
        ncells = length(cells_cl),
        stringsAsFactors = FALSE
      ), 10)
    }))

    if (is.null(cL_results) || nrow(cL_results) == 0) {
      stop("No ScType scores generated.")
    }
    incProgress(0.75)

    score_table <- tidyr::pivot_wider(
      cL_results,
      names_from = type,
      values_from = scores,
      values_fill = 0,
      values_fn = list(scores = max)
    )

    score_cols <- setdiff(colnames(score_table), c("cluster", "ncells"))
    if (length(score_cols) == 0) stop("No score columns detected after pivot.")

    # label max
    score_matrix <- as.matrix(score_table[, score_cols])
    max_scores <- apply(score_matrix, 1, max)
    max_labels <- score_cols[max.col(score_matrix, ties.method = "first")]
    score_table$label_max <- ifelse(
      max_scores < score_table$ncells / 4,
      "Unknown",
      max_labels
    )

    rv$sctype_annotation <- score_table
    incProgress(1)
    append_log("ScType annotation completed.")
  })
})

    # ==============================
    # Table ScType
    # ==============================
    output$sctype_table <- DT::renderDataTable({
      req(rv$sctype_annotation)
      df <- rv$sctype_annotation
      score_cols <- setdiff(colnames(df), c("cluster", "ncells", "label_max"))

      # Colorer la colonne max par ligne
      max_per_row <- apply(df[, score_cols], 1, function(x) names(x)[which.max(x)])

      datatable(
        df,
        extensions = "Buttons",
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = list(
            list(extend = "copy", text = "Copy"),
            list(extend = "csv", text = "CSV", filename = "sctype_scores"),
            list(extend = "excel", text = "Excel", filename = "sctype_scores")
          )
        )
      ) %>%
        formatStyle(
          columns = score_cols,
          backgroundColor = styleEqual(max_per_row, rep('rgba(255,0,0,0.3)', length(max_per_row)))
        )
    })


    ## -----------------------------
    ## Circlepack pour ScType scores
    ## -----------------------------
    output$sctype_circlepack_plot <- renderPlot({
      req(rv$sctype_annotation)  # nécessite que l'annotation ScType ait été calculée

      # edges : cluster → cell type
      edges <- rv$sctype_annotation
      edges <- edges[order(edges$cluster), ]
      edges_long <- tidyr::pivot_longer(edges, cols = -c(cluster, ncells, label_max),
                                        names_to = "type", values_to = "scores")
      edges_long <- edges_long[edges_long$scores > 0, ]
      edges_long$from <- paste0("cluster ", edges_long$cluster)
      edges_long$to <- paste0(edges_long$type, "_", edges_long$cluster)
      edges_final <- edges_long[, c("from", "to", "scores")]

      # nodes : cluster + cell types
      nodes_lvl1 <- unique(edges_final$from)
      nodes_lvl1_df <- data.frame(
        cluster = nodes_lvl1,
        ncells = sapply(nodes_lvl1, function(x) sum(edges_final$scores[edges_final$from == x])),
        Colour = "#f1f1ef",
        ord = 1,
        realname = nodes_lvl1,
        stringsAsFactors = FALSE
      )

      nodes_lvl2_df <- data.frame(
        cluster = edges_final$to,
        ncells = edges_final$scores,
        Colour = rep(RColorBrewer::brewer.pal(12, "Set3"), length.out = nrow(edges_final)),
        ord = 2,
        realname = edges_final$to,
        stringsAsFactors = FALSE
      )

      nodes <- rbind(nodes_lvl1_df, nodes_lvl2_df)
      mygraph <- graph_from_data_frame(d = edges_final[, c("from", "to")], vertices = nodes)

      # plot
      ggraph(mygraph, layout = 'circlepack', weight = I(ncells)) + 
        geom_node_circle(aes(filter = ord == 1, fill = I("#F5F5F5"), colour = I("#D3D3D3")), alpha = 0.9) +
        geom_node_circle(aes(filter = ord == 2, fill = I(Colour), colour = I("#D3D3D3")), alpha = 0.9) +
        geom_node_text(aes(filter = ord == 2, label = realname), colour = "black", size = 5) +
        geom_node_label(aes(filter = ord == 1, label = realname), colour = "black", size = 6, fill = "white") +
        theme_void()
    })


    # ----------------------------------
    # proj signature with AddModuleScore
    # ----------------------------------
      observeEvent(input$run_signature_projection, {

        req(rv$seu)
        req(input$signature_select)
        seu_full <- rv$seu
        cells_subset <- if(!is.null(rv$subset_cells)) Cells(rv$subset_cells) else Cells(seu_full)

        seu <- subset(
          x = seu_full,
          cells = cells_subset
        )
        for (red in Reductions(seu_full)) {
          if (red %in% Reductions(seu_full)) {
            seu[[red]] <- subset(seu_full[[red]], cells = Cells(seu))
          }
        }
        for (g in names(seu_full@graphs)) {
          seu@graphs[[g]] <- seu_full@graphs[[g]][Cells(seu), Cells(seu)]
        }
        selected_sigs <- CancerRNASig::signatures$geneset[input$signature_select]
        genes_present <- rownames(seu)
        selected_sigs <- lapply(selected_sigs, function(sig) intersect(sig, genes_present))
        selected_sigs <- selected_sigs[sapply(selected_sigs, length) > 0]

        if(length(selected_sigs) == 0){
          showNotification("No genes from selected signature found in Seurat object.", type = "error")
          return()
        }
        append_log(paste("Running AddModuleScore for:", paste(names(selected_sigs), collapse = ", ")))
        seu <- AddModuleScore(
          object = seu,
          features = selected_sigs,
          name = "SigScore_"
        )
        rv$signature_seu <- seu
      })


   output$signature_umap_plots <- renderUI({

      req(rv$signature_seu)
      req(input$signature_select)

      n <- length(input$signature_select)

      # max 2 colonnes pour garder de gros plots
      ncol <- min(2, n)
      col_width <- 12 / ncol

      plot_list <- lapply(seq_len(n), function(i) {
        plotname <- paste0("sig_plot_", i)

        column(
          width = col_width,
          plotOutput(plotname, height = "600px")
        )
      })

      rows <- split(plot_list, ceiling(seq_along(plot_list) / ncol))

      do.call(tagList, lapply(rows, fluidRow))
    })


    observe({

    req(rv$signature_seu)
    req(input$signature_select)

    for (i in seq_along(input$signature_select)) {

      local({
        ii <- i
        plotname <- paste0("sig_plot_", ii)

        output[[plotname]] <- renderPlot({

        feature_name <- paste0("SigScore_", ii)

        p <- FeaturePlot(
          rv$signature_seu,
          features = feature_name,
          reduction = "umap",
          pt.size=1.3, 
        )
        p +
          scale_color_distiller(palette = "RdBu") +
          ggtitle(input$signature_select[ii])
      })

      })
    }
  })

    # -----------------------
    # Annotate manual clusters
    # -----------------------

    output$annot_split_meta_ui <- renderUI({
        req(rv$meta_choices)
        selectInput(
          "annot_split_meta",
          "Group cells by:",
          choices = rv$meta_choices,
          selected = rv$meta_choices[1]
        )
      })

      output$cluster_labels_ui <- renderUI({
        req(rv$seu)
        req(input$annot_split_meta)
        clusters <- sort(unique(rv$seu@meta.data[[input$annot_split_meta]]))
        tagList(
          lapply(clusters, function(cl){
            fluidRow(
              column(5, textInput(paste0("label_", cl), 
                                  label = paste("Cluster", cl, "label"), 
                                  value = cl)),
              column(5, textInput(paste0("color_", cl), 
                                  label = "Color (hex)", 
                                  value = rainbow(length(clusters))[which(clusters==cl)]))
            )
          })
        )
      })

    annot_colors <- reactiveVal(NULL)

    observeEvent(input$apply_labels, {
      req(rv$seu)
      req(input$annot_split_meta)
      clusters <- sort(unique(rv$seu@meta.data[[input$annot_split_meta]]))
      labels <- sapply(clusters, function(cl) input[[paste0("label_", cl)]])
      colors <- sapply(clusters, function(cl) input[[paste0("color_", cl)]])
      unique_labels <- unique(labels)
      final_colors <- setNames(colors[match(unique_labels, labels)], unique_labels)
      annot_colors(list(clusters = clusters, labels = labels, colors = final_colors))
    })

    output$annot_umap <- renderPlot({
      req(rv$seu)
      req(input$annot_split_meta)
      seu <- rv$seu
      umap_df <- as.data.frame(Embeddings(seu, "umap"))
      colnames(umap_df) <- c("UMAP1","UMAP2")
      cluster_vec <- as.character(seu@meta.data[[input$annot_split_meta]])

      if (!is.null(annot_colors())) {
        clusters <- annot_colors()$clusters
        labels <- annot_colors()$labels
        colors_input <- annot_colors()$colors
        cluster_to_label <- setNames(labels, clusters)
        cell_labels <- cluster_to_label[cluster_vec]
        umap_df$label <- cell_labels
        label_colors <- colors_input

        ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label)) +
        geom_point(size = 1) +
        scale_color_manual(values = label_colors) +
        theme_minimal() +
        theme(
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        ) +
        ggtitle("UMAP annotated by user labels")

      } else {
        DimPlot(seu, reduction = "umap", group.by = input$annot_split_meta, label = TRUE) +
          theme_minimal()
      }
    })

    output$download_umap_pdf <- downloadHandler(
      filename = function() {
        paste0("UMAP_annotated_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        req(rv$seu)
        req(input$annot_split_meta)
        seu <- rv$seu
        umap_df <- as.data.frame(Embeddings(seu, "umap"))
        colnames(umap_df) <- c("UMAP1","UMAP2")
        cluster_vec <- as.character(seu@meta.data[[input$annot_split_meta]])

        if (!is.null(annot_colors())) {
          clusters <- annot_colors()$clusters
          labels <- annot_colors()$labels
          colors_input <- annot_colors()$colors
          cluster_to_label <- setNames(labels, clusters)
          cell_labels <- cluster_to_label[cluster_vec]
          umap_df$label <- cell_labels
          label_colors <- colors_input
        } else {
          umap_df$label <- cluster_vec
          label_colors <- rainbow(length(unique(cluster_vec)))
        }

        pdf(file, width = 10, height = 8)
        print(
          ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label)) +
            geom_point(size = 1) +
            scale_color_manual(values = label_colors) +
            theme_minimal() +
            theme(
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14)
            ) +
            ggtitle("UMAP annotated by user labels")
        )
        dev.off()
      }
    )
    output$download_metadata_csv <- downloadHandler(
      filename = function() { paste0("metadata_with_labels_", Sys.Date(), ".csv") },
      content = function(file) {
        req(rv$seu)
        req(input$annot_split_meta)
        seu <- rv$seu

        meta_df <- seu@meta.data
        cluster_vec <- as.character(meta_df[[input$annot_split_meta]])

        if (!is.null(annot_colors())) {
          clusters <- annot_colors()$clusters
          labels <- annot_colors()$labels
          cluster_to_label <- setNames(labels, clusters)
          meta_df$custom_label <- cluster_to_label[cluster_vec]
        } else {
          meta_df$custom_label <- cluster_vec
        }

        write.csv(meta_df, file, row.names = TRUE)
      }
    ) 

    ##-----------------------------
    ## Pseudobulk
    ## ----------------------------

    observe({
      req(rv$seu)
      # split_meta : colonnes contenant "sample"
      updateSelectInput(
        session, "split_meta",
        choices = colnames(rv$seu@meta.data)[grepl("sample", colnames(rv$seu@meta.data), ignore.case = TRUE)],
        selected = "sample_id"
      )
      
      # annot_meta : colonnes contenant "sctype" ou "published"
      updateSelectInput(
        session, "annot_meta",
        choices = colnames(rv$seu@meta.data)[grepl("sctype|published", colnames(rv$seu@meta.data), ignore.case = TRUE)],
        selected = "sctype_label"
      )
    })

    rv$pseudobulk_data <- reactiveVal(NULL)

observeEvent(input$run_pseudobulk, {

  req(rv$seu)
  req(input$split_meta)
  req(input$annot_meta)
  req(input$k_monocle)

  seu <- rv$seu
  split_meta <- input$split_meta
  annot_meta <- input$annot_meta
  Kmonocle <- as.numeric(input$k_monocle)

  withProgress(message = "Running Pseudobulk Analysis...", value = 0, {

    # Split Seurat object par sample
    obj.list <- SplitObject(seu, split.by = split_meta)
    moclust_list <- list()
    n_samples <- length(obj.list)
    inc <- 0.9 / n_samples

    for (i in seq_along(obj.list)) {
      incProgress(0.05, detail = paste("Processing sample:", names(obj.list)[i]))

      # Expression matrix compatible Seurat v5
      seurat_v5 <- !("counts" %in% names(slotNames(obj.list[[i]]@assays$RNA)))
      if(seurat_v5){
        expr_mat <- GetAssayData(obj.list[[i]], assay = "RNA", layer = "counts")
      } else {
        expr_mat <- obj.list[[i]]@assays$RNA@counts
      }

      cell_annot <- obj.list[[i]]@meta.data
      gene_annot <- data.frame(gene_short_name = rownames(expr_mat))
      rownames(gene_annot) <- rownames(expr_mat)

      # Créer CellDataSet Monocle3
      cds <- new_cell_data_set(expr_mat, cell_metadata = cell_annot, gene_metadata = gene_annot)
      cds <- preprocess_cds(cds)
      cds <- reduce_dimension(cds)
      cds <- cluster_cells(cds, k = Kmonocle)

      # Ajouter clusters Monocle à Seurat metadata
      cluster_name <- paste0("Monocle_k", Kmonocle)
      obj.list[[i]] <- AddMetaData(obj.list[[i]], metadata = clusters(cds), col.name = cluster_name)

      sample_name <- names(obj.list)[i]
      cluster_vec <- clusters(cds)  
      barcodes <- colnames(obj.list[[i]])  # ce sont les vrais barcodes comme Monc_S01_1
      cluster_vec_prefixed <- paste0(sample_name, "_", cluster_vec)
      moclust_list[[sample_name]] <- setNames(cluster_vec_prefixed, barcodes)
      incProgress(inc)
    }

    # Fusionner tous les clusters
    combined_clusters <- unlist(moclust_list)
    df_clusters <- data.frame(
      cluster = combined_clusters,
      barcode = names(combined_clusters),
      stringsAsFactors = FALSE
    )
    # df_clusters$barcode contient Monc_S01.Monc_S01_1 etc.
    df_clusters$barcode <- sub(".*\\.", "", df_clusters$barcode)
    df_clusters$split_meta <- sub("\\..*$", "", rownames(df_clusters))

    seu@meta.data$MonocleCluster <- df_clusters$cluster[match(rownames(seu@meta.data), df_clusters$barcode)]

    if(all(is.na(seu$MonocleCluster))){
      showNotification("No cells assigned to clusters. Check your split_meta and K.", type = "error")
      return()
    }

    # Pseudobulk
    bulk <- AggregateExpression(seu, group.by = "MonocleCluster", return.seurat = TRUE)
    mat <- GetAssayData(bulk, assay = "RNA", layer = "data")

    # Top 1000 genes les plus variables
    top1000 <- names(sort(apply(mat, 1, sd, na.rm = TRUE), decreasing = TRUE))[1:1000]
    centered_mat <- scale(mat[top1000,], scale = FALSE)
    mat_corr <- cor(centered_mat, use = "pairwise.complete.obs")

    df_clusters$barcode <- gsub("_", "-", df_clusters$barcode)

    # Top annotation (metadata choisi)
    top_annot_df <- data.frame(
      annot = df_clusters$split_meta[match(colnames(mat_corr), df_clusters$barcode)],
      stringsAsFactors = FALSE
    )
    top_vals <- unique(top_annot_df$annot)
    top_vals <- top_vals[!is.na(top_vals)]
    top_col <- setNames(viridis(length(top_vals), option = "C"), top_vals)
    top_annot <- HeatmapAnnotation(
      top_meta = top_annot_df$annot,
      col = list(top_meta = top_col),
      simple_anno_size = unit(1.2, "cm")
    )

     # Bottom annotation : histogramme des proportions
    prop_tab <- prop.table(table(seu@meta.data[[annot_meta]], seu@meta.data$MonocleCluster), 2)
    cell_types <- rownames(prop_tab)
    n <- length(cell_types)
    cols <- qualitative_hcl(n, palette = "Dark3") 
    bottom_col <- setNames(cols, cell_types)
    bottom_annot <- HeatmapAnnotation(
      bar = anno_barplot(t(prop_tab), gp = gpar(fill = bottom_col), height = unit(6, "cm")),
      simple_anno_size = unit(1, "cm")
    )

    # Stocker les données réactives
    rv$pseudobulk_data(list(
      mat_corr = mat_corr,
      top_annot = top_annot,
      bottom_annot = bottom_annot,
      bottom_col = bottom_col
    ))

    incProgress(1)
  })
})


    # Render Plot réactif
    output$pseudobulk_heatmap <- renderPlot({
      data <- rv$pseudobulk_data()
      req(data)

    # Heatmap + annotation
    ht <- Heatmap(
      data$mat_corr,
      col = colorRamp2(c(0.3, 0.6, 1), c("blue", "white", "red")),
      top_annotation = data$top_annot,
      bottom_annotation = data$bottom_annot,
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      show_row_names = TRUE,
      show_column_names = TRUE,
      name = "Corr"
    )

    # Légende du barplot
    barplot_legend <- Legend(
      labels = names(data$bottom_col),
      title = "Cell types",
      legend_gp = gpar(fill = data$bottom_col)
    )

    draw(ht, annotation_legend_list = list(barplot_legend))
    })

    output$download_pseudobulk_pdf <- downloadHandler(
      filename = function() {
        paste0("pseudobulk_heatmap_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        data <- rv$pseudobulk_data()
        req(data)

        # Ouvrir le PDF
        pdf(file, width = 12, height = 15)  # ajuste la taille si besoin

        # Créer la heatmap
        ht <- Heatmap(
          data$mat_corr,
          col = colorRamp2(c(0.3, 0.6, 1), c("blue", "white", "red")),
          top_annotation = data$top_annot,
          bottom_annotation = data$bottom_annot,
          clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson",
          clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2",
          show_row_names = TRUE,
          show_column_names = TRUE,
          name = "Corr"
        )

        # Légende barplot si nécessaire
        if (!is.null(data$bottom_col)) {
          barplot_legend <- Legend(
            labels = names(data$bottom_col),
            title = "Cell types",
            legend_gp = gpar(fill = data$bottom_col)
          )
          draw(ht, annotation_legend_list = list(barplot_legend))
        } else {
          draw(ht)
        }

        dev.off()
      }
    )
}