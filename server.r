packages <- c(
  "shiny", "shinydashboard", "shinycssloaders", "magrittr", "shinyjs",
  "Seurat", "Matrix", "utils", "ggplot2", "gridExtra", "jsonlite", "DT", 
  "HGNChelper", "igraph", "ggraph"
)
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) {
  # Bioconductor manager pour certains packages si nécessaire
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for(pkg in new_packages) {
    # Certains packages sont sur CRAN, d'autres sur Bioconductor
    if(pkg %in% c("Seurat", "Matrix")) {
      install.packages(pkg, dependencies = TRUE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

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

<h3>Notes for Beginners</h3>
<ul>
<li>Always inspect your data after each step.</li>
<li>Quality control is crucial: low-quality or doublet cells can bias results.</li>
<li>Variable features should capture biologically meaningful differences.</li>
<li>PCA reduces noise and speeds up clustering.</li>
<li>UMAP or t-SNE plots help visually explore clusters.</li>
<li>After clustering, marker genes guide interpretation and ScType annotation.</li>
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

      # QC mitochondrial
      mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes)

      # Store Seurat object
      rv$seu <- seu
      rv$files <- list(mtx = mtx_file, genes = genes_file, barcodes = barcodes_file)
      append_log("Seurat object created.")

      # -------------------------
      # QC PLOTS
      # -------------------------
      output$qc_violin <- renderPlot({
        req(rv$seu)
        df <- rv$seu@meta.data
        p1 <- ggplot(df, aes(x = "", y = nFeature_RNA)) + geom_violin(fill = "lightblue") + geom_jitter(width = 0.2, size = 0.5) + theme_classic()
        p2 <- ggplot(df, aes(x = "", y = nCount_RNA)) + geom_violin(fill = "lightgreen") + geom_jitter(width = 0.2, size = 0.5) + theme_classic()
        p3 <- ggplot(df, aes(x = "", y = percent.mt)) + geom_violin(fill = "salmon") + geom_jitter(width = 0.2, size = 0.5) + theme_classic()
        grid.arrange(p1, p2, p3, ncol = 3)
      })

      output$qc_scatter1 <- renderPlot({ FeatureScatter(rv$seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") })
      output$qc_scatter2 <- renderPlot({ FeatureScatter(rv$seu, feature1 = "nCount_RNA", feature2 = "percent.mt") })
 

      # -------------------------
      # SHOW SUMMARY / OVERVIEW
      # -------------------------
      output$seu_summary <- renderPrint({ req(rv$seu); rv$seu })

      rv$meta_base <- colnames(rv$seu@meta.data)  # Stocke les colonnes existantes

      output$meta_table <- DT::renderDataTable({
        req(rv$seu)
        dynamic_cols <- grep("^(louvain_|kmeans_|leiden_)", colnames(rv$seu@meta.data), value = TRUE)
        meta_to_show <- rv$seu@meta.data[, c(rv$meta_base, dynamic_cols), drop = FALSE]
        
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
      rv$meta_choices <- grep("^(louvain_|kmeans_|leiden_)", names(seu@meta.data), value = TRUE)
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

    # ==============================
    # Assay RNA
    # ==============================
    DefaultAssay(seu) <- "RNA"
    seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seu[["RNA"]])))

    append_log(paste("Seurat object", ifelse(seurat_package_v5, "v5", "v4"), "detected"))
    incProgress(0.05)

    # ==============================
    # Extraire scale.data (ou recalculer si vide)
    # ==============================
    if (is.null(seu[["RNA"]]$scale.data) || ncol(seu[["RNA"]]$scale.data) == 0) {
      append_log("scale.data absent ou vide — recalcul avec ScaleData()")
      seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
    }

    scRNAseqData_scaled <- as.matrix(seu[["RNA"]]$scale.data)

    if (nrow(scRNAseqData_scaled) == 0 || ncol(scRNAseqData_scaled) == 0) {
      stop("Erreur : scale.data vide même après ScaleData()")
    }
    incProgress(0.1)

    # ==============================
    # Charger base ScType
    # ==============================
    db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    db_local <- file.path(tempdir(), "ScTypeDB_full.xlsx")
    if (!file.exists(db_local)) download.file(db_url, db_local, mode = "wb")
    gs_list <- gene_sets_prepare(db_local, input$sctype_tissue)
    incProgress(0.2)

    # ==============================
    # ScType scoring
    # ==============================
    es.max <- sctype_score(
      scRNAseqData = scRNAseqData_scaled,
      scaled = TRUE,
      gs = gs_list$gs_positive,
      gs2 = gs_list$gs_negative
    )
    # Ne garder que les cellules du subset
    es.max <- es.max[, colnames(seu), drop = FALSE]
    incProgress(0.5)

    # ==============================
    # Agrégation par cluster
    # ==============================
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

    # ==============================
    # Pivot wide
    # ==============================
    score_table <- tidyr::pivot_wider(
      cL_results,
      names_from = type,
      values_from = scores,
      values_fill = 0,
      values_fn = list(scores = max)
    )

    score_cols <- setdiff(colnames(score_table), c("cluster", "ncells"))
    if (length(score_cols) == 0) stop("No score columns detected after pivot.")

    # ==============================
    # Label max
    # ==============================
    score_matrix <- as.matrix(score_table[, score_cols])
    max_scores <- apply(score_matrix, 1, max)
    max_labels <- score_cols[max.col(score_matrix, ties.method = "first")]

    score_table$label_max <- ifelse(
      max_scores < score_table$ncells / 4,
      "Unknown",
      max_labels
    )

    # ==============================
    # Sauvegarde
    # ==============================
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

      seu <- if (!is.null(rv$subset_cells)) rv$subset_cells else rv$seu

      selected_sigs <- CancerRNASig::signatures$geneset[input$signature_select]

      # filtre gènes présents
      genes_present <- rownames(seu)

      selected_sigs <- lapply(selected_sigs, function(sig) {
        intersect(sig, genes_present)
      })

      append_log(paste(
        "Running AddModuleScore for:",
        paste(input$signature_select, collapse = ", ")
      ))

      seu <- AddModuleScore(
        seu,
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
          reduction = "umap"
        )

        p +
          scale_color_distiller(palette = "RdYlBu") +
          coord_fixed() +
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


  
}