library(shiny)
library(shinydashboard)
library(Seurat)
library(Matrix)
library(utils)
library(ggplot2)
library(shinycssloaders)
library(magrittr)
library(shinyjs)

options(shiny.maxRequestSize = 500*1024^2)

server <- function(input, output, session) {
  rv <- reactiveValues(seu = NULL, meta_choices = NULL, log = character(), files = NULL)
  
  append_log <- function(txt) {
    rv$log <- c(rv$log, txt)
  }
  
  # Show uploaded file list
  output$file_list <- renderText({
    if (is.null(input$input_zip)) return("No file uploaded")
    paste("Uploaded:", input$input_zip$name)
  })
  
  output$loaded_info <- renderUI({
    if (is.null(rv$seu)) {
      tagList(tags$b("No Seurat object loaded"))
    } else {
      ncell <- ncol(rv$seu)
      ngene <- nrow(rv$seu)
      tagList(
        tags$b("Seurat object ready:"),
        tags$ul(
          tags$li(paste("Cells:", ncell)),
          tags$li(paste("Genes:", ngene)),
          tags$li(paste("Assay:", DefaultAssay(rv$seu)))
        )
      )
    }
  })
  
  output$status <- renderText({
    if (length(rv$log) == 0) return("Idle")
    paste(rev(rv$log), collapse = "\n")
  })
  
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

  # -------------------------
  # Main Seurat pipeline
  # -------------------------
  observeEvent(input$run, {
    req(input$input_zip)
    append_log("Starting analysis...")
    shinyjs::disable("run")
    on.exit({ shinyjs::enable("run") })
    
    tmpdir <- tempdir()
    tryCatch({
      append_log("Unzipping input...")
      unzip(input$input_zip$datapath, exdir = tmpdir)
      
      # Find matrix / features / barcodes
      mtx_file <- list.files(tmpdir, pattern = "matrix.*\\.mtx$", full.names = TRUE, ignore.case = TRUE)
      genes_file <- list.files(tmpdir, pattern = "(^genes|features).*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
      barcodes_file <- list.files(tmpdir, pattern = "barcodes.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
      
      if (length(mtx_file) == 0 || length(genes_file) == 0 || length(barcodes_file) == 0) {
        all_files <- list.files(tmpdir, recursive = TRUE, full.names = TRUE)
        if (length(mtx_file) == 0) mtx_file <- all_files[grepl("matrix.*\\.mtx$", all_files, ignore.case = TRUE)]
        if (length(genes_file) == 0) genes_file <- all_files[grepl("(^|/)features.*\\.tsv$|(^|/)genes.*\\.tsv$", all_files, ignore.case = TRUE)]
        if (length(barcodes_file) == 0) barcodes_file <- all_files[grepl("(^|/)barcodes.*\\.tsv$", all_files, ignore.case = TRUE)]
      }
      
      if (length(mtx_file) < 1) stop("matrix.mtx not found inside the zip.")
      if (length(genes_file) < 1) stop("genes.tsv / features.tsv not found inside the zip.")
      if (length(barcodes_file) < 1) stop("barcodes.tsv not found inside the zip.")
      
      mtx_file <- mtx_file[1]
      genes_file <- genes_file[1]
      barcodes_file <- barcodes_file[1]
      
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
        # QC METRICS OUTPUT
        # -------------------------
        observe({
        req(rv$seu)
        seu <- rv$seu

        # Calculate percent mitochondrial if missing
        if (!"percent.mt" %in% colnames(seu@meta.data)) {
            mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
            seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes)
        }

        output$qc_violin <- renderPlot({
            VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    pt.size = 0.1) + ggtitle("QC Metrics")
        })

        output$qc_scatter1 <- renderPlot({
            FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        })

        output$qc_scatter2 <- renderPlot({
            FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
        })
        })


      append_log("Normalizing data...")
      seu <- NormalizeData(seu, verbose = FALSE)
      append_log(paste("Finding variable features (n=", input$n_hvg, ")", sep=""))
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = input$n_hvg, verbose = FALSE)
      append_log("Scaling data...")
      seu <- ScaleData(seu, verbose = FALSE)
      append_log(paste("Running PCA (npcs=", input$n_pcs, ")", sep=""))
      seu <- RunPCA(seu, npcs = input$n_pcs, verbose = FALSE)
      
      # Clustering
      if (input$clust_method == "louvain") {
        append_log("Finding neighbors...")
        seu <- FindNeighbors(seu, dims = 1:input$n_pcs, k.param = input$knn, verbose = FALSE)
        append_log(paste("Finding clusters (resolution=", input$resolution, ")", sep=""))
        seu <- FindClusters(seu, resolution = input$resolution, verbose = FALSE)
      } else {
        append_log("Running Kmeans on PCA embeddings...")
        pca_emb <- Embeddings(seu, "pca")[, 1:input$n_pcs, drop = FALSE]
        km <- kmeans(pca_emb, centers = input$kmeans_centers)
        seu$kmeans_clusters <- as.factor(km$cluster)
        Idents(seu) <- seu$kmeans_clusters
      }
      
      append_log("Running UMAP...")
      seu <- RunUMAP(seu, dims = 1:input$n_pcs, verbose = FALSE)
      
      # -------------------------
      # Incorporate external metadata
      # -------------------------
      if (!is.null(input$meta_file)) {
        append_log("Loading external metadata...")
        meta_ext <- read.table(input$meta_file$datapath, sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)
        # Ensure barcodes match
        common_cells <- intersect(colnames(seu), rownames(meta_ext))
        if (length(common_cells) == 0) {
          append_log("WARNING: no matching barcodes in external metadata.")
        } else {
          seu <- AddMetaData(seu, metadata = meta_ext[common_cells, , drop=FALSE])
          append_log(paste("External metadata merged:", ncol(meta_ext), "columns"))
        }
      }
      
      rv$seu <- seu
      rv$files <- list(mtx = mtx_file, genes = genes_file, barcodes = barcodes_file)
      append_log("Seurat pipeline completed.")
      
      # Update metadata choices
      meta_df <- seu@meta.data
      meta_choices <- names(meta_df)
      rv$meta_choices <- meta_choices
      updateSelectInput(session, "meta_color_by", choices = meta_choices, selected = "seurat_clusters")
      output$metadata_selector <- renderUI({
        selectInput("meta_choice", "Select metadata for splitting / viewing", choices = meta_choices)
      })
      
    }, error = function(e) {
      append_log(paste("ERROR:", e$message))
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
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
      DimPlot(rv$seu, reduction = "umap", group.by = input$meta_color_by, label = FALSE) +
        ggtitle(paste("UMAP colored by", input$meta_color_by))
    }
  })
  
  # -------------------------
  # Violin
  # -------------------------
  output$split_meta_ui <- renderUI({
    req(rv$meta_choices)
    selectInput("split_meta", "Split by (metadata)", choices = rv$meta_choices, selected = rv$meta_choices[1])
  })
  
  observeEvent(input$plot_violin, {
    req(rv$seu)
    gene <- input$violin_gene
    if (gene == "") {
      showNotification("Please provide a gene symbol", type="warning")
      return()
    }
    output$violinPlot <- renderPlot({
      seu <- rv$seu
      if (!(gene %in% rownames(seu))) {
        gene2 <- toupper(gene)
        if (gene2 %in% rownames(seu)) gene <- gene2
      }
      if (!(gene %in% rownames(seu))) {
        plot.new(); title(paste("Gene", gene, "not found"))
      } else {
        if (input$split_by && !is.null(input$split_meta) && input$split_meta %in% colnames(seu@meta.data)) {
          VlnPlot(seu, features = gene, group.by = "seurat_clusters", split.by = input$split_meta, pt.size = 0.1) +
            ggtitle(paste("Violin:", gene))
        } else {
          VlnPlot(seu, features = gene, group.by = "seurat_clusters", pt.size = 0.1) + ggtitle(paste("Violin:", gene))
        }
      }
    })
  })
  
  # -------------------------
  # DotPlot
  # -------------------------
  observeEvent(input$plot_dot, {
    req(rv$seu)
    genes <- unlist(strsplit(input$dot_genes, ","))
    genes <- trimws(genes[genes != ""])
    if (length(genes) == 0) {
      showNotification("Provide at least one gene", type="warning")
      return()
    }
    output$dotPlot <- renderPlot({
      genes_in <- genes[genes %in% rownames(rv$seu)]
      if (length(genes_in) == 0) {
        plot.new(); title("No provided genes found in data")
      } else {
        DotPlot(rv$seu, features = genes_in, group.by = "seurat_clusters") + RotatedAxis() +
          ggtitle("DotPlot by cluster")
      }
    })
  })

    # -------------------------
    # FIND MARKERS
    # -------------------------
    observeEvent(input$run_markers, {
    req(rv$seu)
    seu <- rv$seu

    cluster1 <- input$cluster_a
    cluster2 <- input$cluster_b
    mode <- input$marker_mode

    # Identify identities
    Idents(seu) <- seu$seurat_clusters

    if (mode == "all") {
        res <- FindAllMarkers(seu, only.pos = TRUE)
    } else if (mode == "vscluster") {
        res <- FindMarkers(seu, ident.1 = cluster1, ident.2 = cluster2)
    } else {
        # vs rest
        res <- FindMarkers(seu, ident.1 = cluster1)
    }

    output$marker_table <- renderDataTable({
        res
    })
    })
    

}
