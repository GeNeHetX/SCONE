packages <- c(
  "shiny", "shinydashboard", "shinycssloaders", "magrittr", "shinyjs",
  "Seurat", "Matrix", "utils", "ggplot2", "gridExtra", "jsonlite", "DT"
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
        <p>First, we need to create a Seurat object. This object stores all of your expression data and associated metadata. It is the foundation for all downstream analysis.</p>
        <pre><code class='r'>
        library(Seurat)
        # 'expressionMatrix' is your raw count matrix (genes x cells)
        seu <- CreateSeuratObject(counts = expressionMatrix)
        </code></pre>

        <hr>

        <h3>2. Quality Control (QC)</h3>
        <p>QC is critical to remove poor-quality cells. Cells with very few genes detected might be dead or broken. Cells with too many genes could be doublets (two cells counted as one). We also often remove cells with high mitochondrial gene content because they indicate stressed or dying cells.</p>
        <pre><code class='r'>
        # Filter out low-quality cells
        seu <- subset(seu, subset = nFeature_RNA &gt; 200 &amp;&amp; nFeature_RNA &lt; 5000 &amp;&amp; percent.mt &lt; 10)
        </code></pre>

        <hr>

        <h3>3. Normalize the Data</h3>
        <p>Normalization adjusts for differences in sequencing depth between cells. It ensures that the gene expression values are comparable across all cells.</p>
        <pre><code class='r'>
        seu <- NormalizeData(seu)
        </code></pre>

        <hr>

        <h3>4. Identify Highly Variable Genes</h3>
        <p>We select genes that show high variation across cells. These are the most informative genes for understanding cell types and states. Using only these genes reduces noise and speeds up computations.</p>
        <pre><code class='r'>
        seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
        </code></pre>

        <hr>

        <h3>5. Scale the Data</h3>
        <p>Scaling centers and scales each gene's expression. This step is important because PCA and clustering assume data are scaled. Without scaling, highly expressed genes would dominate the analysis.</p>
        <pre><code class='r'>
        seu <- ScaleData(seu)
        </code></pre>

        <hr>

        <h3>6. Perform PCA (Dimensionality Reduction)</h3>
        <p>PCA (Principal Component Analysis) reduces the dataset's complexity by summarizing the variation into principal components (PCs). These PCs capture the most important patterns in the data.</p>
        <pre><code class='r'>
        seu <- RunPCA(seu, npcs = 20)
        </code></pre>

        <hr>

        <h3>7. Clustering</h3>
        <p>We group similar cells together. Seurat allows two methods:</p>
        <ul>
        <li><strong>Louvain clustering</strong> finds clusters based on nearest-neighbor graphs of cells.</li>
        <li><strong>K-means clustering</strong> divides cells into a pre-defined number of clusters based on PCA embeddings.</li>
        </ul>
        <pre><code class='r'>
        # Louvain clustering
        seu <- FindNeighbors(seu, dims = 1:20)
        seu <- FindClusters(seu, resolution = 0.5)

        # OR K-means clustering
        pca_emb <- Embeddings(seu, 'pca')[, 1:20]
        km <- kmeans(pca_emb, centers = 10)
        seu$kmeans_clusters <- as.factor(km$cluster)
        Idents(seu) <- seu$kmeans_clusters
        </code></pre>

        <hr>

        <h3>8. Visualize with UMAP</h3>
        <p>UMAP projects cells into 2D space so we can visually inspect clusters. Each point represents a single cell, and cells in the same cluster often group together in this plot.</p>
        <pre><code class='r'>
        seu <- RunUMAP(seu, dims = 1:20)
        DimPlot(seu, reduction = 'umap', label = TRUE)
        </code></pre>

        <hr>

        <h3>Notes for Beginners</h3>
        <ul>
        <li>Always inspect your data after each step.</li>
        <li>QC is crucial: low-quality or doublet cells can bias your results.</li>
        <li>Variable features should capture the most biologically meaningful differences between cells.</li>
        <li>PCA reduces noise and speeds up clustering.</li>
        <li>UMAP or t-SNE plots are useful for visually exploring your data and understanding clusters.</li>
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

      DimPlot(rv$seu, reduction = "umap", group.by = input$marker_split_meta, label = TRUE) +
        theme_minimal()
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
