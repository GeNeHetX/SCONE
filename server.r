### ---- GLOBAL OPTIONS ----
options(timeout = 1000)
options(download.file.method = "libcurl")

cran_packages <- c(
  "shiny", "shinydashboard", "shinycssloaders", "magrittr", "shinyjs",
  "Seurat", "Matrix", "utils", "ggplot2", "gridExtra", "jsonlite", "DT",
  "HGNChelper", "igraph", "ggraph", "scCustomize", "dplyr",
  "circlize", "vegan", "colourpicker", "openxlsx", "tidyr", "tibble"
)
cran_new <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(cran_new)) install.packages(cran_new, dependencies = TRUE)

if (!"grr" %in% installed.packages()[, "Package"]) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
    repos = NULL, type = "source"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("ComplexHeatmap")
bioc_new <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(bioc_new)) BiocManager::install(bioc_new, ask = FALSE, update = FALSE)

github_packages <- c(
  "GeNeHetX/CancerRNASig",
  "bnprks/BPCells/r",
  "cole-trapnell-lab/monocle3",
  "immunogenomics/presto",
  "must-bioinfo/fastCNV"
)
if (length(github_packages)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  
  for (repo in github_packages) {
    pkg_name <- basename(repo)  # ex: "CancerRNASig" depuis "GeNeHetX/CancerRNASig"
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      tryCatch(
        remotes::install_github(repo, dependencies = TRUE, upgrade = "never"),
        error = function(e) message("Could not install ", repo, ": ", e$message)
      )
    } else {
      message(pkg_name, " already installed — skipping.")
    }
  }
}

library(shiny); library(shinydashboard); library(Seurat); library(Matrix)
library(utils); library(ggplot2); library(shinycssloaders); library(magrittr)
library(shinyjs); library(gridExtra); library(jsonlite); library(DT)
library(presto); library(HGNChelper); library(igraph); library(ggraph)
library(scCustomize); library(dplyr); library(monocle3); library(ComplexHeatmap)
library(circlize); library(vegan); library(CancerRNASig); library(viridis)
library(colorspace); library(colourpicker); library(openxlsx); library(tidyr)
library(tibble); library(fastCNV)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

options(shiny.maxRequestSize = 100 * 1024^3)
DT::datatable(data.frame(), extensions = "Buttons")

subset_done          <- reactiveVal(FALSE)
celltype_subset_done <- reactiveVal(FALSE)
cnv_correction_done  <- reactiveVal(FALSE)

# ── Helper : choisir la bonne reduction UMAP ──────────────────────────────────
get_umap_reduction <- function(seu) {
  if ("umap_cnv_corrected" %in% names(seu@reductions)) "umap_cnv_corrected" else "umap"
}

server <- function(input, output, session) {

  rv <- reactiveValues(
    seu               = NULL,
    seu_original      = NULL,
    subset_cells      = NULL,
    meta_choices      = NULL,
    meta_base         = NULL,
    log               = character(),
    files             = NULL,
    sctype_annotation = NULL,
    selected_sig_names = NULL,
    signature_seu     = NULL,
    marker_res        = NULL,
    pseudobulk_matrix = NULL,
    pseudobulk_data   = NULL
  )

  append_log <- function(txt) rv$log <- c(rv$log, txt)

  # ==============================
  # TUTORIAL
  # ==============================
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
          <p>QC identifies and removes low-quality cells before analysis. Four metrics are visualized as violin plots, each with a red dashed line indicating the filtering threshold. Thresholds should always be set by visually inspecting the violin plots for your dataset.</p>

          <h4>nFeature_RNA - Number of detected genes per cell</h4>
          <p>Counts the number of unique genes detected in each cell, meaning genes with at least 1 UMI count. This is the most direct measure of cell quality.</p>
          <ul>
            <li><b>Too low</b>: empty droplet, dead cell, or poor RNA capture.</li>
            <li><b>Too high</b>: likely a doublet, two cells captured in the same droplet.</li>
          </ul>
          <pre><code class='r'>seu &lt;- subset(seu, subset = nFeature_RNA &gt; 600 &amp;&amp; nFeature_RNA &lt; 5500)</code></pre>

          <h4>nCount_RNA - Total UMI count per cell</h4>
          <p>A UMI (Unique Molecular Identifier) is a short random barcode attached to each RNA molecule before PCR amplification. UMI counts represent the number of original molecules captured, correcting for PCR duplicates. nCount_RNA is the total number of UMIs detected in a cell and reflects overall transcriptional activity.</p>
          <ul>
            <li><b>Too low</b>: poor RNA capture, empty droplet, or a dying cell with degraded RNA.</li>
            <li><b>Too high</b>: likely a doublet, or a very large and highly active cell.</li>
          </ul>
          <pre><code class='r'>seu &lt;- subset(seu, subset = nCount_RNA &gt; 1200 &amp;&amp; nCount_RNA &lt; 45000)</code></pre>

          <h4>percent.mt - Mitochondrial gene percentage</h4>
          <p>Proportion of UMIs mapping to mitochondrial genes (prefixed MT-). Mitochondria have their own genome and produce their own RNA independently of the nucleus.</p>
          <p>High mitochondrial content can have two very different causes:</p>
          <ul>
            <li><b>Cell damage or lysis</b>: when a cell membrane is disrupted, cytoplasmic RNA leaks out of the droplet. Mitochondria, being membrane-bound organelles, are retained longer and their RNA becomes over-represented. This is a technical artifact.</li>
            <li><b>Organelle burst</b>: the nucleus is intact but the cytoplasm has been lost, leaving only mitochondrial RNA.</li>
            <li><b>Tumor or metabolically active cells</b>: some cancer cells, particularly those relying on oxidative phosphorylation, genuinely express high levels of mitochondrial transcripts. In this case, high percent.mt is a biological signal and not an artifact. A strict threshold may inadvertently remove real cancer cells.</li>
          </ul>
          <p>The threshold must be adapted to your tissue and biological context. A cutoff of 10% is standard for healthy tissue but may need to be relaxed to 20-30% in tumor samples.</p>
          <pre><code class='r'>seu &lt;- subset(seu, subset = percent.mt &lt; 20)</code></pre>

          <h4>Complexity - nFeature_RNA divided by nCount_RNA</h4>
          <p>This ratio measures transcriptional complexity: how many distinct genes are detected per UMI captured. A complex cell expresses many different genes at moderate levels, while a low-complexity cell expresses very few genes at very high levels.</p>
          <ul>
            <li><b>High complexity (close to 1)</b>: each UMI maps to a different gene, indicating a rich and diverse transcriptome.</li>
            <li><b>Low complexity (close to 0)</b>: many UMIs mapping to very few genes. Typical of red blood cells, empty droplets, or highly specialized cells.</li>
          </ul>
          <pre><code class='r'>seu &lt;- subset(seu, subset = nFeature_RNA &gt; 600 &amp;&amp; nFeature_RNA &lt; 5500 &amp;&amp; nCount_RNA &gt; 1200 &amp;&amp; nCount_RNA &lt; 45000 &amp;&amp; percent.mt &lt; 20)</code></pre>

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

  observeEvent(input$show_tutorial, shinyjs::toggle("tutorial_panel"))

  # ==============================
  # FONCTION COMMUNE init_seu
  # ==============================
  init_seu <- function(seu) {
    rv$seu          <- seu
    rv$seu_original <- seu
    rv$subset_cells <- NULL
    rv$meta_choices <- sort(grep(
      "sampleID|^published|sctype_|manual_|^(louvain_|kmeans_|leiden_)",
      colnames(seu@meta.data), value = TRUE
    ))
    rv$meta_base    <- colnames(seu@meta.data)

    existing <- sort(grep(
      "sampleID|sctype_|published_|manual_|^(louvain_|kmeans_|leiden_)",
      colnames(seu@meta.data), value = TRUE
    ))
    if (length(existing) > 0)
      updateSelectInput(session, "meta_color_by", choices = existing, selected = existing[1])

    updateSelectInput(session, "qc_split_meta",
      choices = colnames(seu@meta.data), selected = colnames(seu@meta.data)[1])
  }

  # ==============================
  # LOAD ZIP
  # ==============================
  observeEvent(input$input_zip, {
    req(input$input_zip)
    append_log("Starting Seurat object creation...")
    tmpdir <- file.path(tempdir(), paste0("upload_", as.integer(Sys.time())))
    dir.create(tmpdir, showWarnings = FALSE)

    tryCatch({
      unzip(input$input_zip$datapath, exdir = tmpdir)
      mtx_file      <- list.files(tmpdir, pattern = "matrix.*\\.mtx$",          full.names = TRUE, ignore.case = TRUE)[1]
      genes_file    <- list.files(tmpdir, pattern = "(^genes|features).*\\.tsv$", full.names = TRUE, ignore.case = TRUE)[1]
      barcodes_file <- list.files(tmpdir, pattern = "barcodes.*\\.tsv$",         full.names = TRUE, ignore.case = TRUE)[1]

      expressionMatrix <- ReadMtx(mtx = mtx_file, cells = barcodes_file, features = genes_file,
        feature.column = 1, skip.cell = 1, skip.feature = 1, cell.sep = "\t")
      seu <- CreateSeuratObject(counts = expressionMatrix)

      barcode_df <- read.table(barcodes_file, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
      if (ncol(barcode_df) > 1) {
        colnames(barcode_df)[1] <- "barcode"
        barcode_df <- barcode_df[barcode_df$barcode %in% colnames(seu), , drop = FALSE]
        barcode_df <- barcode_df[match(colnames(seu), barcode_df$barcode), , drop = FALSE]
        if (ncol(barcode_df) > 1) {
          meta_to_add <- barcode_df[, -1, drop = FALSE]
          rownames(meta_to_add) <- barcode_df$barcode
          meta_to_add <- meta_to_add[, !colnames(meta_to_add) %in% colnames(seu@meta.data), drop = FALSE]
          if (ncol(meta_to_add) > 0) seu <- AddMetaData(seu, metadata = meta_to_add)
        }
      }

      mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
      if (length(mito_genes) > 0) {
        mat_seu <- seu@assays$RNA$counts
        seu$percent_mito <- Matrix::colSums(mat_seu[mito_genes, ]) / Matrix::colSums(mat_seu) * 100
      } else seu$percent_mito <- 0

      rv$files <- list(mtx = mtx_file, genes = genes_file, barcodes = barcodes_file)
      init_seu(seu)
      append_log("Seurat object created.")

    }, error = function(e) {
      append_log(paste("ERROR:", e$message))
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  })

  # ==============================
  # LOAD RDS
  # ==============================
  observeEvent(input$input_rds, {
    req(input$input_rds)
    append_log("Loading Seurat RDS object...")
    withProgress(message = "Loading Seurat object...", {
      tryCatch({
        seu <- readRDS(input$input_rds$datapath)
        if (!inherits(seu, "Seurat")) {
          showNotification("Not a valid Seurat object.", type = "error"); return()
        }
        if (!"percent_mito" %in% colnames(seu@meta.data)) {
          mito_genes <- grep("^MT-", rownames(seu), value = TRUE)
          if (length(mito_genes) > 0) {
            mat_seu <- seu@assays$RNA$counts
            seu$percent_mito <- Matrix::colSums(mat_seu[mito_genes, ]) / Matrix::colSums(mat_seu) * 100
          } else seu$percent_mito <- 0
        }
        init_seu(seu)
        existing <- grep("^(louvain_|kmeans_|leiden_)", colnames(seu@meta.data), value = TRUE)
        append_log(paste("Loaded:", ncol(seu), "cells,", nrow(seu), "features.",
          if ("umap_cnv_corrected" %in% names(seu@reductions)) "CNV-corrected UMAP detected."
          else if ("umap" %in% names(seu@reductions)) "UMAP detected." else "No UMAP.",
          if (length(existing) > 0) paste("Clusterings:", paste(existing, collapse = ", ")) else "No clustering."))
        showNotification(paste("Loaded:", ncol(seu), "cells,", nrow(seu), "genes."), type = "message")
      }, error = function(e) {
        append_log(paste("ERROR loading RDS:", e$message))
        showNotification(paste("Error:", e$message), type = "error", duration = 8)
      })
    })
  })

  # ==============================
  # SAVE RDS
  # ==============================
  output$download_seurat_rds <- downloadHandler(
    filename = function() paste0("SeuratObject_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content  = function(file) {
      req(rv$seu)
      withProgress(message = "Saving...", saveRDS(rv$seu, file))
      append_log("Seurat object saved.")
    }
  )

  # ==============================
  # OUTPUTS GLOBAUX
  # ==============================
  output$seu_summary <- renderPrint({ req(rv$seu); rv$seu })

  output$meta_table <- DT::renderDataTable({
    req(rv$seu)
    DT::datatable(rv$seu@meta.data, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$subset_ui <- renderUI({
    req(rv$seu)
    tagList(
      numericInput("nFeature_min", "Min nFeature_RNA", value = 200),
      numericInput("nFeature_max", "Max nFeature_RNA", value = 5000),
      numericInput("high_cutoff_UMI", "Max UMI",       value = 45000)
    )
  })

  output$subset_message <- renderUI({
    req(subset_done())
    tags$span("→ Run analysis again with your subsetted Seurat object.",
      style = "margin-left:10px; font-style:italic; color:#666;")
  })

  # ==============================
  # QC PLOTS
  # ==============================
  observe({
    req(rv$seu)
    updateSelectInput(session, "qc_split_meta",
      choices = colnames(rv$seu@meta.data), selected = colnames(rv$seu@meta.data)[1])
  })

  output$qc_violin <- renderPlot({
    req(rv$seu, input$qc_split_meta, input$qc_split_meta %in% colnames(rv$seu@meta.data))
    df        <- rv$seu@meta.data
    split_var <- input$qc_split_meta

    p1 <- ggplot(df, aes(x = .data[[split_var]], y = nFeature_RNA)) +
      geom_violin(fill = "lightblue") + geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
      geom_hline(yintercept = c(600, 5500), linetype = "dashed", color = "red") +
      labs(y = "nFeature_RNA", x = "", title = "Gene per Cell/Nucleus") + theme_classic()

    p2 <- ggplot(df, aes(x = .data[[split_var]], y = nCount_RNA)) +
      geom_violin(fill = "lightgreen") + geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
      geom_hline(yintercept = c(1200, 45000), linetype = "dashed", color = "red") +
      labs(y = "nCount_RNA", x = "", title = "UMI per Cell/Nucleus") + theme_classic()

    p3 <- if ("percent_mito" %in% colnames(df)) {
      ggplot(df, aes(x = .data[[split_var]], y = percent_mito)) +
        geom_violin(fill = "salmon") + geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
        geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
        labs(y = "percent_mito", x = "", title = "Mito Gene %") + theme_classic()
    } else ggplot() + ggtitle("percent_mito not found") + theme_classic()

    df <- df %>% mutate(complexity = nFeature_RNA / nCount_RNA)
    p4 <- ggplot(df, aes(x = .data[[split_var]], y = complexity)) +
      geom_violin(fill = "lightgrey") + geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      labs(y = "Complexity", x = "", title = "Cell complexity") + theme_classic()

    grid.arrange(p1, p2, p3, p4, ncol = 4)
  })

  output$qc_scatter1 <- renderPlot({
    req(rv$seu)
    QC_Plot_UMIvsGene(seurat_object = rv$seu, low_cutoff_gene = input$nFeature_min,
      high_cutoff_gene = input$nFeature_max, low_cutoff_UMI = 500, high_cutoff_UMI = input$high_cutoff_UMI)
  })
  output$qc_scatter2 <- renderPlot({
    req(rv$seu)
    QC_Plot_GenevsFeature(seurat_object = rv$seu, feature1 = "percent_mito",
      low_cutoff_gene = input$nFeature_min, high_cutoff_gene = input$nFeature_max,
      high_cutoff_feature = input$percent_mt_max)
  })
  output$qc_scatter3 <- renderPlot({
    req(rv$seu)
    QC_Plot_UMIvsGene(seurat_object = rv$seu, meta_gradient_name = "percent_mito",
      low_cutoff_gene = input$nFeature_min, high_cutoff_gene = input$nFeature_max,
      high_cutoff_UMI = input$high_cutoff_UMI, meta_gradient_low_cutoff = input$percent_mt_max)
  })
  output$qc_scatter4 <- renderPlot({
    req(rv$seu)
    QC_Plot_UMIvsGene(seurat_object = rv$seu, meta_gradient_name = "percent_mito",
      low_cutoff_gene = input$nFeature_min, high_cutoff_gene = input$nFeature_max,
      high_cutoff_UMI = input$high_cutoff_UMI)
  })

  # ==============================
  # SUBSET QC
  # ==============================
  observeEvent(input$apply_subset, {
    subset_done(TRUE)
    req(rv$seu_original)
    cells_keep <- WhichCells(rv$seu_original, expression =
      nFeature_RNA >= input$nFeature_min & nFeature_RNA <= input$nFeature_max)
    rv$subset_cells <- subset(rv$seu_original, cells = cells_keep)
    rv$seu          <- rv$subset_cells
    append_log(paste("Subset applied:", length(cells_keep), "cells kept"))
    output$subset_summary <- renderPrint({
      req(rv$subset_cells)
      cat("Subset summary:\n")
      cat("Cells kept:", ncol(rv$subset_cells), "\n")
      cat("Features:", nrow(rv$subset_cells), "\n")
    })
  })

  # ==============================
  # SUBSET BY CELL TYPE
  # ==============================
  annot_col_available <- reactive({
    req(rv$seu)
    found <- grep("sctype|manual_annot|published", colnames(rv$seu@meta.data),
                  value = TRUE, ignore.case = TRUE)
    if (length(found) == 0) return(NULL)
    return(sort(found))
  })

  output$celltype_subset_ui <- renderUI({
    req(input$subset_by_celltype == "yes")
    if (is.null(annot_col_available())) {
      return(tags$p(style = "color:#E53935; font-style:italic;",
        "No annotation found. Please run ScType, Manual Annotation, or add a published annotation."))
    }
    tagList(
      selectInput("celltype_annot_col", "Select annotation to use:",
        choices = annot_col_available(), selected = annot_col_available()[1]),
      uiOutput("celltype_table_ui"),
      actionButton("apply_celltype_subset", "Apply Cell Type Subset",
        icon = icon("filter"), style = "background:#4CAF50; color:white; margin-top:10px;")
    )
  })

  output$celltype_table_ui <- renderUI({
    req(input$celltype_annot_col, rv$seu)
    cell_types <- sort(unique(rv$seu@meta.data[[input$celltype_annot_col]]))
    cell_types <- cell_types[!is.na(cell_types)]
    checkboxGroupInput("celltype_keep", "Select cell types to keep:",
      choices = cell_types, selected = cell_types)
  })

  observeEvent(input$apply_celltype_subset, {
    req(rv$seu_original, input$celltype_annot_col, input$celltype_keep)
    annot_col <- input$celltype_annot_col
    if (!(annot_col %in% colnames(rv$seu_original@meta.data))) {
      rv$seu_original@meta.data[[annot_col]] <- rv$seu@meta.data[[annot_col]][
        match(rownames(rv$seu_original@meta.data), rownames(rv$seu@meta.data))]
    }
    cells_keep <- rownames(rv$seu_original@meta.data)[
      rv$seu_original@meta.data[[annot_col]] %in% input$celltype_keep]
    if (length(cells_keep) == 0) {
      showNotification("No cells found.", type = "error"); return()
    }
    seu_sub <- subset(rv$seu_original, cells = cells_keep)
    new_col <- paste0(annot_col, "_subset")
    seu_sub@meta.data[[new_col]] <- seu_sub@meta.data[[annot_col]]
    rv$seu               <- seu_sub
    rv$subset_cells      <- seu_sub
    rv$meta_choices      <- colnames(seu_sub@meta.data)
    rv$sctype_annotation <- NULL
    celltype_subset_done(TRUE)
    append_log(paste("Cell type subset:", length(cells_keep), "cells —",
      paste(input$celltype_keep, collapse = ", ")))
    output$celltype_subset_summary <- renderPrint({
      cat("Cell type subset summary:\n")
      cat("Cells kept:", ncol(seu_sub), "\n")
      cat("Cell types:", paste(input$celltype_keep, collapse = ", "), "\n")
      cat("Column added:", new_col, "\n")
    })
  })

  output$celltype_subset_message <- renderUI({
    req(celltype_subset_done())
    tags$span("→ Cell type subset applied. Run analysis again.",
      style = "margin-left:10px; font-style:italic; color:#666;")
  })

  # ==============================
  # CNV CORRECTION
  # ==============================
  output$cnv_annot_col_ui <- renderUI({
    req(input$cnv_correction == "yes")
    if (is.null(annot_col_available())) {
      return(tags$p(style = "color:#E53935; font-style:italic;",
        "No annotation found. Please run ScType or Manual Annotation first."))
    }
    selectInput("cnv_annot_col", "Select annotation to use:",
      choices = annot_col_available(), selected = annot_col_available()[1])
  })

  output$cnv_ref_cells_ui <- renderUI({
    req(input$cnv_correction == "yes", input$cnv_annot_col, rv$seu)
    cell_types <- sort(unique(rv$seu@meta.data[[input$cnv_annot_col]]))
    cell_types <- cell_types[!is.na(cell_types)]
    tagList(
      checkboxGroupInput("cnv_ref_celltypes", "Reference cell types:", choices = cell_types, selected = cell_types)
    )
  })

  observeEvent(input$run_cnv_correction, {
    req(rv$seu, input$cnv_annot_col, input$cnv_ref_celltypes)
    if (is.null(annot_col_available())) {
      showNotification("No annotation found.", type = "error"); return()
    }
    if (length(input$cnv_ref_celltypes) == 0) {
      showNotification("Please select at least one reference cell type.", type = "error"); return()
    }
    annot_col <- input$cnv_annot_col
    seu       <- rv$seu
    if (!(annot_col %in% colnames(seu@meta.data))) {
      showNotification(paste("Column", annot_col, "not found."), type = "error"); return()
    }

    withProgress(message = "Running CNV correction...", value = 0, {

      # ── 1. fastCNV per sample ────────────────────────────────────────────
      append_log("Step 1: Running fastCNV per sample...")

      save_path <- file.path(tempdir(), paste0("CNV_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
      append_log(paste("CNV outputs will be saved to:", save_path))
      rv$cnv_save_path <- save_path  # stocker pour le download

      sample_col <- colnames(seu@meta.data)[grepl("sample", colnames(seu@meta.data), ignore.case = TRUE)][1]

      if (!is.null(sample_col) && !is.na(sample_col)) {
        sample_list  <- SplitObject(seu, split.by = sample_col)
        sample_names <- names(sample_list)
        sample_list  <- fastCNV(seuratObj = sample_list, sampleName = sample_names,
          referenceVar = annot_col, referenceLabel = input$cnv_ref_celltypes,
          printPlot = FALSE, doPlot = TRUE, savePath = save_path, outputType = "png")
        seu <- merge(sample_list[[1]], sample_list[-1])
        DefaultAssay(seu) <- "RNA"
        seu <- JoinLayers(seu)
      } else {
        append_log("No sample column — running fastCNV on full object.")
        seu <- fastCNV(seuratObj = seu, referenceVar = annot_col,
          referenceLabel = input$cnv_ref_celltypes,
          printPlot = FALSE, doPlot = TRUE, savePath = save_path, outputType = "png")
      }
      incProgress(0.3)

      # ── 2. Extraire matrice CNV + sauvegarder + garder seulement cnv_fraction ─────────
      append_log("Step 2: Extracting CNV arm matrix...")
      chr_arms <- grep("_CNV$", colnames(seu@meta.data), value = TRUE)
      chr_arms <- chr_arms[!grepl("classification|fraction|cluster", chr_arms)]
      if (length(chr_arms) == 0) {
        showNotification("No CNV arm columns found.", type = "error"); return()
      }

      cnv_arm_matrix <- as.matrix(seu@meta.data[, chr_arms])

      # Sauvegarder la matrice CNV bras chromosomiques
      cnv_df <- data.frame(
        cell_barcode = rownames(seu@meta.data),
        cnv_arm_matrix,
        stringsAsFactors = FALSE
      )

      # Ajouter aussi cnv_fraction et cnv_clusters si disponibles
      extra_cols <- grep("cnv_fraction|cnv_clusters|_CNV_classification", colnames(seu@meta.data), value = TRUE)
      if (length(extra_cols) > 0) {
        cnv_df <- cbind(cnv_df, seu@meta.data[, extra_cols, drop = FALSE])
      }

      write.table(cnv_df,
        file      = file.path(save_path, "cnv_arm_matrix.tsv"),
        sep       = "\t",
        row.names = FALSE,
        quote     = FALSE
      )
      append_log(paste("CNV arm matrix saved to:", file.path(save_path, "cnv_arm_matrix.tsv")))

      # Retirer les colonnes bras — garder uniquement cnv_fraction
      cols_to_remove <- c(chr_arms,
        grep("_CNV_cluster|_CNV_classification", colnames(seu@meta.data), value = TRUE))
      seu@meta.data <- seu@meta.data[, !colnames(seu@meta.data) %in% cols_to_remove, drop = FALSE]
      incProgress(0.4)

      # ── 3. PCA sur matrice CNV ────────────────────────────────────────────
      append_log("Step 3: PCA on CNV matrix...")
      cnv_pca    <- prcomp(cnv_arm_matrix, center = TRUE, scale. = FALSE, rank. = 10)
      cnv_scores <- cnv_pca$x
      cnv_vars   <- paste0("cnv_pc", 1:ncol(cnv_scores))
      for (i in 1:ncol(cnv_scores)) seu[[cnv_vars[i]]] <- cnv_scores[colnames(seu), i]
      incProgress(0.5)

      # ── 4. PCA transcriptomique standard ──────────────────────────────────
      append_log("Step 4: Standard transcriptomic PCA...")
      seu <- seu |>
        NormalizeData(verbose = FALSE) |>
        FindVariableFeatures(nfeatures = input$n_hvg, verbose = FALSE) |>
        ScaleData(verbose = FALSE) |>
        RunPCA(npcs = input$n_pcs, verbose = FALSE)
      expr_pcs <- Embeddings(seu, "pca")
      incProgress(0.6)

      # ── 5. Regression CNV PCs ─────────────────────────────────────────────
      append_log("Step 5: Regressing CNV signal...")
      cnv_mat <- as.matrix(seu@meta.data[colnames(seu), cnv_vars])
      expr_pcs_corrected <- apply(expr_pcs, 2, function(pc_col) residuals(lm(pc_col ~ cnv_mat)))
      rownames(expr_pcs_corrected) <- colnames(seu)
      colnames(expr_pcs_corrected) <- paste0("PC_", 1:input$n_pcs)

      # Retirer les colonnes cnv_pc du metadata
      seu@meta.data <- seu@meta.data[, !colnames(seu@meta.data) %in% cnv_vars, drop = FALSE]
      incProgress(0.7)

      # ── 6. Stocker PCA corrigee + UMAP ───────────────────────────────────
      append_log("Step 6: Storing corrected PCA and running UMAP...")
      seu[["pca_cnv_corrected"]] <- CreateDimReducObject(
        embeddings = expr_pcs_corrected, key = "PCACORR_", assay = "RNA")

      seu <- seu |>
        FindNeighbors(reduction = "pca_cnv_corrected", dims = 1:input$n_pcs, verbose = FALSE) |>
        FindClusters(resolution = input$resolution, verbose = FALSE) |>
        RunUMAP(reduction = "pca_cnv_corrected", dims = 1:input$n_pcs,
                reduction.name = "umap_cnv_corrected", verbose = FALSE)
      incProgress(0.9)

      rv$seu          <- seu
      rv$seu          <- seu
      rv$meta_choices <- sort(grep(
        "sampleID|^published|sctype_|manual_|cnv_fraction|^(louvain_|kmeans_|leiden_)",
        names(seu@meta.data), value = TRUE
      ))
      cnv_correction_done(TRUE)
      incProgress(1)
      append_log("CNV correction completed — umap_cnv_corrected available.")

      output$cnv_correction_summary <- renderPrint({
        cat("CNV correction summary:\n")
        cat("Annotation used:", annot_col, "\n")
        cat("Reference cell types:", paste(input$cnv_ref_celltypes, collapse = ", "), "\n")
        cat("CNV arms detected:", length(chr_arms), "\n")
        cat("CNV PCs computed:", ncol(cnv_scores), "\n")
        cat("Transcriptomic PCs corrected:", input$n_pcs, "\n")
        cat("New reduction: pca_cnv_corrected\n")
        cat("New UMAP: umap_cnv_corrected\n")
        cat("Metadata kept: cnv_fraction only\n")
      })
    })
  })

  output$cnv_correction_message <- renderUI({
    req(cnv_correction_done())
    tags$span("→ CNV correction applied. UMAP computed on corrected PCA (umap_cnv_corrected).",
      style = "margin-left:10px; font-style:italic; color:#666;")
  })

  output$cnv_feature_plot_ui <- renderUI({
  req(cnv_correction_done())
  req("cnv_fraction" %in% colnames(rv$seu@meta.data))
  tagList(
    hr(),
    h4("CNV Fraction per cell"),
    plotOutput("cnv_fraction_plot", height = "500px") %>% withSpinner()
  )
})

output$cnv_fraction_plot <- renderPlot({
  req(rv$seu)
  req(cnv_correction_done())
  req("cnv_fraction" %in% colnames(rv$seu@meta.data))
  red <- get_umap_reduction(rv$seu)
  FeaturePlot(rv$seu, features = "cnv_fraction", reduction = red, pt.size = 1) +
    scale_color_viridis_c(option = "C") +
    ggtitle("CNV Fraction")
})

output$download_cnv_results <- downloadHandler(
  filename = function() paste0("CNV_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
  content  = function(file) {
    req(rv$cnv_save_path)
    files_to_zip <- list.files(rv$cnv_save_path, full.names = TRUE)
    zip(file, files = files_to_zip, flags = "-j")
  }
)

output$cnv_download_ui <- renderUI({
  req(cnv_correction_done())
  downloadButton(
    "download_cnv_results",
    "Download CNV results (.zip)",
    icon  = icon("download"),
    style = "background:#7B1FA2; color:white; margin-top:10px;"
  )
})

  # ==============================
  # RUN ANALYSIS
  # ==============================
  observeEvent(input$run_analysis, {
    req(rv$seu)
    shinyjs::disable("run_analysis")
    on.exit(shinyjs::enable("run_analysis"))

    tryCatch({
      seu <- rv$seu

      if (!is.null(input$nFeature_min) && !is.null(input$nFeature_max)) {
        cells_keep <- WhichCells(seu, expression =
          nFeature_RNA >= input$nFeature_min & nFeature_RNA <= input$nFeature_max)
        seu <- subset(seu, cells = cells_keep)
        append_log(paste("Subset applied:", length(cells_keep), "cells kept"))
      }

      append_log("Normalizing...")
      seu <- NormalizeData(seu, verbose = FALSE)
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = input$n_hvg, verbose = FALSE)
      seu <- ScaleData(seu, verbose = FALSE)
      seu <- RunPCA(seu, npcs = input$n_pcs, verbose = FALSE)
      append_log(paste("PCA done (npcs =", input$n_pcs, ")"))

      if (input$clust_method == "louvain") {
        cluster_name <- paste0("louvain_n", input$knn, "_res", input$resolution)
        seu <- FindNeighbors(seu, dims = 1:input$n_pcs, k.param = input$knn)
        seu <- FindClusters(seu, resolution = input$resolution)
        seu@meta.data[[cluster_name]] <- Idents(seu)
      } else if (input$clust_method == "leiden") {
        cluster_name <- paste0("leiden_n", input$leiden_k, "_res", input$resolution)
        seu <- FindNeighbors(seu, dims = 1:input$n_pcs, k.param = input$leiden_k)
        seu <- FindClusters(seu, resolution = input$resolution, algorithm = 4)
        seu@meta.data[[cluster_name]] <- Idents(seu)
      } else {
        cluster_name <- paste0("kmeans_n", input$kmeans_centers, "_res", input$resolution)
        pca_emb <- Embeddings(seu, "pca")[, 1:input$n_pcs, drop = FALSE]
        km <- kmeans(pca_emb, centers = input$kmeans_centers)
        seu@meta.data[[cluster_name]] <- as.factor(km$cluster)
      }
      Idents(seu) <- seu@meta.data[[cluster_name]]

      if (!"umap" %in% names(seu@reductions)) {
        append_log("Running UMAP...")
        seu <- RunUMAP(seu, dims = 1:input$n_pcs, verbose = FALSE)
      }

      rv$seu          <- seu
      rv$meta_choices <- sort(grep(
          "sampleID|^published|sctype_|manual_|cnv_fraction|^(louvain_|kmeans_|leiden_)",
          names(seu@meta.data), value = TRUE
        ))
      append_log(paste("Analysis completed:", cluster_name))
      updateSelectInput(session, "meta_color_by", choices = rv$meta_choices, selected = cluster_name)

    }, error = function(e) {
      append_log(paste("ERROR:", e$message))
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  })

  # ==============================
  # UMAP PLOTS — tous utilisent get_umap_reduction
  # ==============================
  output$umapPlot <- renderPlot({
    req(rv$seu)
    red <- get_umap_reduction(rv$seu)
    tryCatch(
      DimPlot(rv$seu, reduction = red, label = TRUE, repel = TRUE) +
        ggtitle(if (red == "umap_cnv_corrected") "UMAP - CNV corrected" else "UMAP - clusters"),
      error = function(e) ggplot() + ggtitle("Unable to plot UMAP")
    )
  })

  output$umapMetaPlot <- renderPlot({
    req(rv$seu, input$meta_color_by)
    red <- get_umap_reduction(rv$seu)
    if (!(input$meta_color_by %in% colnames(rv$seu@meta.data)))
      ggplot() + ggtitle("Selected metadata not found")
    else
      DimPlot(rv$seu, reduction = red, label = TRUE, group.by = input$meta_color_by)
  })

  output$umapViolinPlot <- renderPlot({
    req(rv$seu, input$split_meta)
    DimPlot(rv$seu, reduction = get_umap_reduction(rv$seu), group.by = input$split_meta, label = TRUE)
  })

  output$umapDotPlot <- renderPlot({
    req(rv$seu, input$dot_split_meta)
    DimPlot(rv$seu, reduction = get_umap_reduction(rv$seu), group.by = input$dot_split_meta, label = TRUE)
  })

  output$umapMarkerPlot <- renderPlot({
    req(rv$seu, input$marker_split_meta)
    DimPlot(rv$seu, reduction = get_umap_reduction(rv$seu), group.by = input$marker_split_meta, label = TRUE)
  })

  # ==============================
  # LOG
  # ==============================
  output$status <- renderText({
    if (length(rv$log) == 0) return("Idle")
    paste(rev(rv$log), collapse = "\n")
  })

  # ==============================
  # VIOLIN
  # ==============================
  gene_families <- reactive({
    path <- "markers_celltypes.json"
    req(file.exists(path))
    jsonlite::fromJSON(path)
  })

  output$violin_gene_ui <- renderUI({
    req(rv$seu)
    selectizeInput("violin_gene", "Gene(s)", choices = rownames(rv$seu), selected = NULL,
      multiple = TRUE, options = list(placeholder = "Start typing gene...", maxOptions = 1000))
  })

  output$violin_gene_family_ui <- renderUI({
    req(gene_families())
    selectInput("violin_gene_family", "Gene family (optional)",
      choices = names(gene_families()), multiple = TRUE, selected = "")
  })

  observe({
    req(rv$seu)
    updateSelectizeInput(session, "violin_gene", choices = rownames(rv$seu), server = TRUE)
  })

  output$split_meta_ui <- renderUI({
    req(rv$meta_choices)
    sorted_choices <- sort(rv$meta_choices)
    selectInput("split_meta", "Split by (metadata)", choices = sorted_choices, selected = sorted_choices[1])
  })

  observeEvent(input$plot_violin, {
    req(rv$seu)
    genes_manual <- input$violin_gene
    genes_family <- character(0)
    if (!is.null(input$violin_gene_family) && length(input$violin_gene_family) > 0 &&
        !("" %in% input$violin_gene_family))
      genes_family <- unlist(gene_families()[input$violin_gene_family])

    genes_exist <- unique(c(genes_manual, genes_family))
    genes_exist <- genes_exist[genes_exist %in% rownames(rv$seu)]
    if (length(genes_exist) == 0) {
      showNotification("No genes found in Seurat object.", type = "error"); return()
    }

    output$violinPlot <- renderPlot({
      df        <- rv$seu@meta.data
      expr_data <- GetAssayData(rv$seu, layer = "data")
      plots <- lapply(genes_exist, function(gene) {
        df[[gene]] <- expr_data[gene, ]
        if (input$split_by && !is.null(input$split_meta) && input$split_meta %in% colnames(df)) {
          ggplot(df, aes_string(x = input$split_meta, y = gene, fill = input$split_meta)) +
            geom_violin() + geom_jitter(width = 0.1, size = 0.1) +
            theme_classic() + theme(legend.position = "none") + ggtitle(gene)
        } else {
          ggplot(df, aes_string(x = "1", y = gene)) +
            geom_violin() + geom_jitter(width = 0.1, size = 0.5) +
            theme_classic() + theme(legend.position = "none") + ggtitle(gene)
        }
      })
      gridExtra::grid.arrange(grobs = plots, ncol = 3)
    })
  })

  # ==============================
  # DOTPLOT
  # ==============================
  observe({
    req(rv$seu)
    updateSelectizeInput(session, "dot_genes", choices = rownames(rv$seu), server = TRUE)
  })

  output$dot_gene_family_ui <- renderUI({
    req(gene_families())
    selectInput("dot_gene_family", "Gene family (optional)",
      choices = c("None" = "", names(gene_families())), multiple = TRUE, selected = "")
  })

  output$dot_split_meta_ui <- renderUI({
    req(rv$meta_choices)
    sorted_choices <- sort(rv$meta_choices)
    selectInput("dot_split_meta", "Group DotPlot by:", choices = sorted_choices, selected = sorted_choices[1])
  })

  observeEvent(input$plot_dot, {
    req(rv$seu)
    genes_manual     <- input$dot_genes
    selected_families <- input$dot_gene_family[input$dot_gene_family != ""]
    genes_family <- if (length(selected_families) > 0)
      unlist(lapply(selected_families, function(f) gene_families()[[f]])) else character(0)

    genes_exist <- unique(c(genes_manual, genes_family))
    genes_exist <- genes_exist[genes_exist %in% rownames(rv$seu)]
    if (length(genes_exist) == 0) {
      showNotification("No genes found.", type = "error"); return()
    }

    output$dotPlot <- renderPlot({
      group_var <- if (!is.null(input$dot_split_meta) &&
                       input$dot_split_meta %in% colnames(rv$seu@meta.data))
        input$dot_split_meta else colnames(rv$seu@meta.data)[1]

      dp   <- DotPlot(rv$seu, features = genes_exist, group.by = group_var)
      data <- dp$data

      exp_mat <- data %>%
        dplyr::select(features.plot, id, avg.exp) %>%
        tidyr::pivot_wider(names_from = id, values_from = avg.exp) %>%
        tibble::column_to_rownames("features.plot") %>% as.matrix()

      percent_mat <- data %>%
        dplyr::select(features.plot, id, pct.exp) %>%
        tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>%
        tibble::column_to_rownames("features.plot") %>% as.matrix()
      percent_mat <- percent_mat / 100

      col_dist   <- as.dist(1 - cor(exp_mat, use = "pairwise.complete.obs", method = "pearson"))
      col_hclust <- hclust(col_dist, method = "ward.D2")
      col_dend   <- as.dendrogram(col_hclust)

      grp_cols  <- setNames(colorspace::qualitative_hcl(ncol(exp_mat), palette = "Dark3"), colnames(exp_mat))
      column_ha <- HeatmapAnnotation(
        Group = anno_simple(colnames(exp_mat), col = grp_cols, gp = gpar(col = NA)),
        annotation_label = group_var, show_annotation_name = TRUE, simple_anno_size = unit(4, "cm"))

      col_fun <- circlize::colorRamp2(
        c(0, quantile(exp_mat, 0.90, na.rm = TRUE) / 2, quantile(exp_mat, 0.90, na.rm = TRUE)),
        viridis::viridis(100)[c(1, 50, 100)])

      unit_size <- 4
      layer_fun <- function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = NA, fill = NA))
        grid.circle(x = x, y = y, r = pindex(percent_mat, i, j) * unit(unit_size, "mm"),
          gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))
      }

      lgd_list <- list(
        Legend(title = "Avg\nexpression", col_fun = col_fun, direction = "vertical",
          title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14),
          legend_height = unit(6, "cm"), grid_width = unit(8, "mm")),
        Legend(labels = c("0%", "25%", "50%", "75%", "100%"), title = "% expressing",
          title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 14),
          grid_height = unit(10, "mm"), grid_width = unit(10, "mm"),
          graphics = list(
            function(x, y, w, h) grid.circle(x=x, y=y, r=unit(0.01,             "mm"), gp=gpar(fill="black")),
            function(x, y, w, h) grid.circle(x=x, y=y, r=unit(0.25*unit_size,   "mm"), gp=gpar(fill="black")),
            function(x, y, w, h) grid.circle(x=x, y=y, r=unit(0.5 *unit_size,   "mm"), gp=gpar(fill="black")),
            function(x, y, w, h) grid.circle(x=x, y=y, r=unit(0.75*unit_size,   "mm"), gp=gpar(fill="black")),
            function(x, y, w, h) grid.circle(x=x, y=y, r=unit(1   *unit_size,   "mm"), gp=gpar(fill="black"))))
      )

      hp <- Heatmap(exp_mat, show_heatmap_legend = FALSE, rect_gp = gpar(type = "none"),
        layer_fun = layer_fun, top_annotation = column_ha, cluster_columns = col_dend,
        cluster_column_slices = FALSE, cluster_rows = TRUE,
        clustering_distance_rows = "pearson", clustering_method_rows = "ward.D2",
        row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 14),
        column_names_rot = 45, border = "black", name = "Expression")

      draw(hp, annotation_legend_list = lgd_list,
        padding = unit(c(5, 5, 5, 5), "mm"), legend_gap = unit(5, "mm"))
    })
  })

  # ==============================
  # FIND MARKERS
  # ==============================
  output$marker_split_meta_ui <- renderUI({
    req(rv$meta_choices)
    sorted_choices <- sort(rv$meta_choices)
    selectInput("marker_split_meta", "Group cells by:", choices = sorted_choices, selected = sorted_choices[1])
  })

  output$cluster_select_ui <- renderUI({
    req(rv$seu, input$marker_split_meta)
    clusters <- sort(unique(rv$seu@meta.data[[input$marker_split_meta]]))
    if (input$marker_mode == "vsrest") {
      selectInput("cluster_a", "Cluster of interest:", choices = clusters, selected = clusters[1])
    } else {
      tagList(
        selectInput("cluster_a", "Cluster A:", choices = clusters, selected = clusters[1]),
        selectInput("cluster_b", "Cluster B:", choices = clusters,
          selected = if (length(clusters) > 1) clusters[2] else clusters[1])
      )
    }
  })

  observeEvent(input$run_markers, {
    req(rv$seu, input$marker_split_meta)
    seu  <- rv$seu
    Idents(seu) <- seu[[input$marker_split_meta]][, 1]
    mode <- input$marker_mode

    withProgress(message = "Finding markers...", {
      if (mode == "vscluster" && input$cluster_a == input$cluster_b) {
        showNotification("Please select two different clusters.", type = "error"); return()
      }
      res <- if (mode == "vscluster") {
        req(input$cluster_a, input$cluster_b)
        r <- FindMarkers(seu, ident.1 = input$cluster_a, ident.2 = input$cluster_b, logfc.threshold = 0.25)
        if (nrow(r) > 0) r$cluster <- paste(input$cluster_a, "vs", input$cluster_b)
        r
      } else {
        req(input$cluster_a)
        r <- FindMarkers(seu, ident.1 = input$cluster_a, logfc.threshold = 0.25)
        if (nrow(r) > 0) r$cluster <- input$cluster_a
        r
      }
    })

    res$gene      <- rownames(res)
    rv$marker_res <- res

    if (nrow(res) == 0) {
      showNotification("No markers found.", type = "warning")
      output$marker_table     <- DT::renderDataTable(DT::datatable(data.frame()))
      output$top_marker_table <- DT::renderDataTable(DT::datatable(data.frame()))
      return()
    }

    output$marker_table <- DT::renderDataTable(DT::datatable(res, extensions = "Buttons",
      options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
        buttons = list(list(extend="copy",text="Copy"),
          list(extend="csv",text="CSV",filename="markers"),
          list(extend="excel",text="Excel",filename="markers")))))

    top_markers <- if (mode == "vscluster") {
      top_a <- res |> dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) |>
        dplyr::slice_max(avg_log2FC, n = input$top_n) |>
        dplyr::mutate(top_for = input$cluster_a) |> dplyr::select(top_for, gene, avg_log2FC)
      top_b <- res |> dplyr::filter(p_val_adj < 0.05 & avg_log2FC < 0) |>
        dplyr::slice_min(avg_log2FC, n = input$top_n) |>
        dplyr::mutate(top_for = input$cluster_b, avg_log2FC = abs(avg_log2FC)) |>
        dplyr::select(top_for, gene, avg_log2FC)
      dplyr::bind_rows(top_a, top_b)
    } else {
      res |> dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) |>
        dplyr::group_by(cluster) |> dplyr::slice_max(avg_log2FC, n = input$top_n) |>
        dplyr::ungroup() |> dplyr::select(cluster, gene, avg_log2FC)
    }

    output$top_marker_table <- DT::renderDataTable(DT::datatable(top_markers, extensions = "Buttons",
      options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip",
        buttons = list(list(extend="copy",text="Copy"),
          list(extend="csv",text="CSV",filename="top_markers"),
          list(extend="excel",text="Excel",filename="top_markers")))))
  })

  output$marker_table_title <- renderUI({
    HTML("<b>Columns:</b><br><ul>
      <li><b>avg_log2FC</b>: log2 fold change</li>
      <li><b>pct.1</b>: fraction expressing in cluster of interest</li>
      <li><b>pct.2</b>: fraction expressing in other cluster</li>
      <li><b>p_val_adj</b>: adjusted p-value</li>
      <li><b>cluster</b>: cluster of interest</li></ul>")
  })

  output$download_top_markers_tsv <- downloadHandler(
    filename = function() paste0("top_markers_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"),
    content  = function(file) {
      req(rv$seu)
      seu <- rv$seu
      Idents(seu) <- seu[[input$marker_split_meta]][, 1]
      all_clusters <- sort(unique(Idents(seu)))
      withProgress(message = "Computing markers for all clusters...", {
        df_export <- lapply(all_clusters, function(cl) {
          res_cl      <- FindMarkers(seu, ident.1 = cl, logfc.threshold = 0.25, only.pos = TRUE)
          res_cl$gene <- rownames(res_cl)
          genes <- res_cl |> dplyr::filter(p_val_adj < 0.05) |>
            dplyr::slice_max(avg_log2FC, n = input$top_n) |> dplyr::pull(gene)
          c(genes, rep(NA, input$top_n - length(genes)))
        })
      })
      df_export           <- as.data.frame(df_export)
      colnames(df_export) <- as.character(all_clusters)
      write.table(df_export, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    }
  )

  # ==============================
  # SCTYPE
  # ==============================
  output$type_split_meta_ui <- renderUI({
    req(rv$meta_choices)
    sorted_choices <- sort(rv$meta_choices)
    selectInput("sctype_split_meta", "Group cells by (metadata)",
      choices = sorted_choices, selected = sorted_choices[1])
  })

  output$sctype_umap_plot <- renderPlot({
    req(rv$seu, input$sctype_split_meta)
    seu         <- rv$seu
    cluster_col <- input$sctype_split_meta
    red         <- get_umap_reduction(seu)
    if (!(cluster_col %in% colnames(seu@meta.data))) return()
    Idents(seu) <- seu@meta.data[[cluster_col]]
    if (!is.null(rv$sctype_annotation)) {
      label_map        <- setNames(rv$sctype_annotation$label_max, rv$sctype_annotation$cluster)
      seu$ScType_label <- as.character(seu@meta.data[[cluster_col]])
      for (cl in names(label_map)) seu$ScType_label[seu@meta.data[[cluster_col]] == cl] <- label_map[cl]
      DimPlot(seu, reduction = red, group.by = "ScType_label", label = TRUE) +
        theme(legend.position = "bottom")
    } else {
      DimPlot(seu, reduction = red, group.by = cluster_col, label = TRUE)
    }
  })

  observeEvent(input$run_sctype_annotation, {
    req(rv$seu, input$sctype_split_meta)
    seu         <- rv$seu
    cluster_col <- input$sctype_split_meta
    if (!(cluster_col %in% colnames(seu@meta.data))) return()
    Idents(seu) <- seu@meta.data[[cluster_col]]

    withProgress(message = "Running ScType annotation...", value = 0, {
      DefaultAssay(seu) <- "RNA"
      if (is.null(seu[["RNA"]]$scale.data) || ncol(seu[["RNA"]]$scale.data) == 0)
        seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
      scRNAseqData_scaled <- as.matrix(seu[["RNA"]]$scale.data)
      incProgress(0.1)

      db_local <- file.path(tempdir(), "ScTypeDB_full.xlsx")
      if (!file.exists(db_local))
        download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
          db_local, mode = "wb")
      gs_list <- gene_sets_prepare(db_local, input$sctype_tissue)
      incProgress(0.2)

      es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE,
        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
      es.max <- es.max[, colnames(seu), drop = FALSE]
      incProgress(0.5)

      clusters   <- unique(as.character(Idents(seu)))
      cL_results <- do.call("rbind", lapply(clusters, function(cl) {
        cells_cl <- colnames(seu)[Idents(seu) == cl]
        if (length(cells_cl) == 0) return(NULL)
        scores_cluster <- sort(rowSums(es.max[, cells_cl, drop = FALSE]), decreasing = TRUE)
        head(data.frame(cluster = cl, type = names(scores_cluster),
          scores = scores_cluster, ncells = length(cells_cl), stringsAsFactors = FALSE), 10)
      }))
      if (is.null(cL_results) || nrow(cL_results) == 0) stop("No ScType scores generated.")
      incProgress(0.75)

      score_table  <- tidyr::pivot_wider(cL_results, names_from = type, values_from = scores,
        values_fill = 0, values_fn = list(scores = max))
      score_cols   <- setdiff(colnames(score_table), c("cluster", "ncells"))
      score_matrix <- as.matrix(score_table[, score_cols])
      max_scores   <- apply(score_matrix, 1, max)
      max_labels   <- score_cols[max.col(score_matrix, ties.method = "first")]
      score_table$label_max <- ifelse(max_scores < score_table$ncells / 4, "Unknown", max_labels)

      cluster_to_label           <- setNames(score_table$label_max, score_table$cluster)
      seu@meta.data$sctype_annot <- cluster_to_label[as.character(seu@meta.data[[cluster_col]])]

      rv$sctype_annotation <- score_table
      rv$seu               <- seu
      rv$meta_choices      <- colnames(seu@meta.data)
      incProgress(1)
      append_log("ScType annotation completed.")
    })
  })

  output$sctype_table <- DT::renderDataTable({
    req(rv$sctype_annotation)
    df         <- rv$sctype_annotation
    score_cols <- setdiff(colnames(df), c("cluster", "ncells", "label_max"))
    max_per_row <- apply(df[, score_cols], 1, function(x) names(x)[which.max(x)])
    datatable(df, extensions = "Buttons", rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip",
        buttons = list(list(extend="copy",text="Copy"),
          list(extend="csv",text="CSV",filename="sctype"),
          list(extend="excel",text="Excel",filename="sctype")))) %>%
      formatStyle(columns = score_cols,
        backgroundColor = styleEqual(max_per_row, rep("rgba(255,0,0,0.3)", length(max_per_row))))
  })

  output$sctype_circlepack_plot <- renderPlot({
    req(rv$sctype_annotation)
    edges      <- rv$sctype_annotation[order(rv$sctype_annotation$cluster), ]
    edges_long <- tidyr::pivot_longer(edges, cols = -c(cluster, ncells, label_max),
      names_to = "type", values_to = "scores")
    edges_long <- edges_long[edges_long$scores > 0, ]
    edges_long$from <- paste0("cluster ", edges_long$cluster)
    edges_long$to   <- paste0(edges_long$type, "_", edges_long$cluster)
    edges_final     <- edges_long[, c("from", "to", "scores")]

    nodes_lvl1 <- unique(edges_final$from)
    nodes <- rbind(
      data.frame(cluster = nodes_lvl1,
        ncells = sapply(nodes_lvl1, function(x) sum(edges_final$scores[edges_final$from == x])),
        Colour = "#f1f1ef", ord = 1, realname = nodes_lvl1, stringsAsFactors = FALSE),
      data.frame(cluster = edges_final$to, ncells = edges_final$scores,
        Colour = rep(RColorBrewer::brewer.pal(12, "Set3"), length.out = nrow(edges_final)),
        ord = 2, realname = edges_final$to, stringsAsFactors = FALSE)
    )
    mygraph <- graph_from_data_frame(d = edges_final[, c("from", "to")], vertices = nodes)
    ggraph(mygraph, layout = "circlepack", weight = I(ncells)) +
      geom_node_circle(aes(filter = ord == 1, fill = I("#F5F5F5"), colour = I("#D3D3D3")), alpha = 0.9) +
      geom_node_circle(aes(filter = ord == 2, fill = I(Colour),   colour = I("#D3D3D3")), alpha = 0.9) +
      geom_node_text(aes(filter  = ord == 2, label = realname), colour = "black", size = 5) +
      geom_node_label(aes(filter = ord == 1, label = realname), colour = "black", size = 6, fill = "white") +
      theme_void()
  })

  # ==============================
  # SIGNATURE PROJECTION
  # ==============================
  custom_signatures <- reactiveVal(list())

  observeEvent(input$custom_sig_file, {
    req(input$custom_sig_file)
    df       <- openxlsx::read.xlsx(input$custom_sig_file$datapath)
    sig_list <- lapply(df, function(col) col[!is.na(col) & col != ""])
    custom_signatures(sig_list)
    updateSelectizeInput(session, "custom_sig_select", choices = names(sig_list))
    append_log(paste("Custom signatures loaded:", paste(names(sig_list), collapse = ", ")))
  })

  observeEvent(input$run_signature_projection, {
    req(rv$seu)
    all_selected <- c(input$signature_select, input$custom_sig_select)
    if (length(all_selected) == 0) {
      showNotification("Please select at least one signature.", type = "error"); return()
    }
    seu_full     <- rv$seu
    cells_subset <- if (!is.null(rv$subset_cells)) Cells(rv$subset_cells) else Cells(seu_full)
    seu          <- subset(x = seu_full, cells = cells_subset)
    for (red in Reductions(seu_full)) seu[[red]] <- subset(seu_full[[red]], cells = Cells(seu))
    for (g in names(seu_full@graphs)) seu@graphs[[g]] <- seu_full@graphs[[g]][Cells(seu), Cells(seu)]

    all_sigs      <- c(CancerRNASig::signatures$geneset, custom_signatures())
    selected_sigs <- all_sigs[all_selected]
    genes_present <- rownames(seu)
    selected_sigs <- lapply(selected_sigs, function(sig) intersect(sig, genes_present))
    selected_sigs <- selected_sigs[sapply(selected_sigs, length) > 0]
    if (length(selected_sigs) == 0) {
      showNotification("No genes found in Seurat object.", type = "error"); return()
    }
    append_log(paste("AddModuleScore for:", paste(names(selected_sigs), collapse = ", ")))
    seu <- AddModuleScore(object = seu, features = selected_sigs, name = "SigScore_")
    rv$selected_sig_names <- names(selected_sigs)
    rv$signature_seu      <- seu
  })

  output$signature_umap_plots <- renderUI({
    req(rv$signature_seu, rv$selected_sig_names)
    n         <- length(rv$selected_sig_names)
    ncol      <- min(2, n)
    col_width <- 12 / ncol
    plot_list <- lapply(seq_len(n), function(i)
      column(width = col_width, plotOutput(paste0("sig_plot_", i), height = "600px")))
    rows <- split(plot_list, ceiling(seq_along(plot_list) / ncol))
    do.call(tagList, lapply(rows, fluidRow))
  })

  observe({
    req(rv$signature_seu, rv$selected_sig_names)
    red <- get_umap_reduction(rv$signature_seu)
    for (i in seq_along(rv$selected_sig_names)) {
      local({
        ii       <- i
        sig_name <- rv$selected_sig_names[ii]
        output[[paste0("sig_plot_", ii)]] <- renderPlot({
          p <- FeaturePlot(rv$signature_seu, features = paste0("SigScore_", ii),
            reduction = red, pt.size = 1.3)
          p + scale_color_distiller(palette = "RdBu") + ggtitle(sig_name)
        })
      })
    }
  })

  output$dot_split_meta_sig_ui <- renderUI({
    req(rv$signature_seu)
    sorted_choices <- sort(grep("^(louvain_|kmeans_|leiden_|sampleID)",
      colnames(rv$signature_seu@meta.data), value = TRUE))
    selectInput("dot_split_meta_sig", "Group DotPlot by:",
      choices = sorted_choices, selected = sorted_choices[1])
  })

  output$signature_dotplot <- renderPlot({
    req(rv$signature_seu, rv$selected_sig_names, input$dot_split_meta_sig)
    seu        <- rv$signature_seu
    n_sigs     <- length(rv$selected_sig_names)
    score_cols <- paste0("SigScore_", seq_len(n_sigs))
    score_cols <- score_cols[score_cols %in% colnames(seu@meta.data)]
    if (length(score_cols) == 0) return()

    meta         <- seu@meta.data
    meta$cluster <- meta[[input$dot_split_meta_sig]]
    df_long <- do.call(rbind, lapply(seq_along(score_cols), function(i)
      data.frame(cluster = meta$cluster, signature = rv$selected_sig_names[i],
        score = meta[[score_cols[i]]], stringsAsFactors = FALSE)))

    df_summary <- df_long %>% dplyr::group_by(cluster, signature) %>%
      dplyr::summarise(avg_score = mean(score, na.rm = TRUE),
        pct_exp = mean(score > 0, na.rm = TRUE) * 100, .groups = "drop")

    ggplot(df_summary, aes(x = cluster, y = signature)) +
      geom_point(aes(size = pct_exp, color = avg_score)) +
      scale_color_distiller(palette = "RdBu", name = "Avg score") +
      scale_size_continuous(range = c(1, 10), name = "% cells > 0") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 14), axis.title = element_blank()) +
      labs(title = "Signature scores per cluster")
  })

  output$download_sig_scores <- downloadHandler(
    filename = function() paste0("signature_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"),
    content  = function(file) {
      req(rv$signature_seu, rv$selected_sig_names)
      n_sigs     <- length(rv$selected_sig_names)
      score_cols <- paste0("SigScore_", seq_len(n_sigs))
      score_cols <- score_cols[score_cols %in% colnames(rv$signature_seu@meta.data)]
      df_export  <- rv$signature_seu@meta.data[, score_cols, drop = FALSE]
      colnames(df_export) <- rv$selected_sig_names
      df_export$cell_barcode <- rownames(df_export)
      df_export <- df_export[, c("cell_barcode", rv$selected_sig_names)]
      write.table(df_export, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # ==============================
  # MANUAL ANNOTATION
  # ==============================
  output$annot_split_meta_ui <- renderUI({
    req(rv$meta_choices)
    sorted_choices <- sort(grep("^(louvain_|kmeans_|leiden_)", names(rv$seu@meta.data), value = TRUE))
    selectInput("annot_split_meta", "Group cells by:", choices = sorted_choices, selected = sorted_choices[1])
  })

  output$cluster_labels_ui <- renderUI({
    req(rv$seu, input$annot_split_meta)
    clusters <- sort(unique(rv$seu@meta.data[[input$annot_split_meta]]))
    tagList(lapply(clusters, function(cl) {
      fluidRow(
        column(5, textInput(paste0("label_", cl), label = paste("Cluster", cl),
          value = isolate(input[[paste0("label_", cl)]]) %||% cl)),
        column(5, colourInput(paste0("color_", cl), label = "Color",
          value = isolate(input[[paste0("color_", cl)]]) %||% rainbow(length(clusters))[which(clusters == cl)],
          showColour = "both"))
      )
    }))
  })

  annot_colors <- reactiveVal(NULL)

  observeEvent(input$apply_labels, {
    req(rv$seu, input$annot_split_meta)
    clusters      <- sort(unique(rv$seu@meta.data[[input$annot_split_meta]]))
    labels        <- sapply(clusters, function(cl) input[[paste0("label_", cl)]])
    colors        <- sapply(clusters, function(cl) input[[paste0("color_", cl)]])
    unique_labels <- unique(labels)
    final_colors  <- setNames(colors[match(unique_labels, labels)], unique_labels)
    annot_colors(list(clusters = clusters, labels = labels, colors = final_colors))
    cluster_to_label              <- setNames(labels, clusters)
    rv$seu@meta.data$manual_annot <- cluster_to_label[as.character(rv$seu@meta.data[[input$annot_split_meta]])]
    rv$meta_choices               <- colnames(rv$seu@meta.data)
    append_log("Manual annotation added as manual_annot.")
  })

  output$annot_umap <- renderPlot({
    req(rv$seu, input$annot_split_meta)
    seu         <- rv$seu
    red         <- get_umap_reduction(seu)
    umap_df     <- as.data.frame(Embeddings(seu, red))
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    cluster_vec <- as.character(seu@meta.data[[input$annot_split_meta]])
    if (!is.null(annot_colors())) {
      cluster_to_label <- setNames(annot_colors()$labels, annot_colors()$clusters)
      umap_df$label    <- cluster_to_label[cluster_vec]
      label_colors     <- annot_colors()$colors
    } else {
      umap_df$label <- cluster_vec
      label_colors  <- rainbow(length(unique(cluster_vec)))
    }
    ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label)) +
      geom_point(size = 1) + scale_color_manual(values = label_colors) +
      theme_minimal() + theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      ggtitle("UMAP annotated by user labels")
  })

  output$download_umap_pdf <- downloadHandler(
    filename = function() paste0("UMAP_annotated_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content  = function(file) {
      req(rv$seu, input$annot_split_meta)
      seu         <- rv$seu
      red         <- get_umap_reduction(seu)
      umap_df     <- as.data.frame(Embeddings(seu, red))
      colnames(umap_df) <- c("UMAP1", "UMAP2")
      cluster_vec <- as.character(seu@meta.data[[input$annot_split_meta]])
      if (!is.null(annot_colors())) {
        cluster_to_label <- setNames(annot_colors()$labels, annot_colors()$clusters)
        umap_df$label    <- cluster_to_label[cluster_vec]
        label_colors     <- annot_colors()$colors
      } else {
        umap_df$label <- cluster_vec
        label_colors  <- rainbow(length(unique(cluster_vec)))
      }
      pdf(file, width = 10, height = 8)
      print(ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label)) +
        geom_point(size = 1) + scale_color_manual(values = label_colors) +
        theme_minimal() + ggtitle("UMAP annotated by user labels"))
      dev.off()
    }
  )

  output$download_metadata_csv <- downloadHandler(
    filename = function() paste0("metadata_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      req(rv$seu, input$annot_split_meta)
      meta_df     <- rv$seu@meta.data
      cluster_vec <- as.character(meta_df[[input$annot_split_meta]])
      if (!is.null(annot_colors())) {
        cluster_to_label     <- setNames(annot_colors()$labels, annot_colors()$clusters)
        meta_df$custom_label <- cluster_to_label[cluster_vec]
      } else meta_df$custom_label <- cluster_vec
      write.csv(meta_df, file, row.names = TRUE)
    }
  )

  # ==============================
  # PSEUDOBULK
  # ==============================
  observe({
    req(rv$seu)
    updateSelectInput(session, "split_meta",
      choices  = colnames(rv$seu@meta.data)[grepl("sample", colnames(rv$seu@meta.data), ignore.case = TRUE)],
      selected = "sample_id")
    updateSelectInput(session, "annot_meta",
      choices  = colnames(rv$seu@meta.data)[grepl(
        "^(louvain_|kmeans_|leiden_|sctype_|manual_|published)",
        colnames(rv$seu@meta.data), ignore.case = TRUE)],
      selected = "sctype_label")
  })

  rv$pseudobulk_matrix <- reactiveVal(NULL)
  rv$pseudobulk_data   <- reactiveVal(NULL)

  observeEvent(input$run_pseudobulk, {
    req(rv$seu, input$split_meta, input$annot_meta, input$k_monocle)
    seu      <- rv$seu
    Kmonocle <- as.numeric(input$k_monocle)

    withProgress(message = "Running Pseudobulk Analysis...", value = 0, {
      obj.list     <- SplitObject(seu, split.by = input$split_meta)
      moclust_list <- list()
      inc          <- 0.9 / length(obj.list)

      for (i in seq_along(obj.list)) {
        incProgress(0.05, detail = paste("Processing:", names(obj.list)[i]))
        seurat_v5 <- !("counts" %in% names(slotNames(obj.list[[i]]@assays$RNA)))
        expr_mat  <- if (seurat_v5)
          GetAssayData(obj.list[[i]], assay = "RNA", layer = "counts")
        else obj.list[[i]]@assays$RNA@counts

        cell_annot <- obj.list[[i]]@meta.data
        gene_annot <- data.frame(gene_short_name = rownames(expr_mat))
        rownames(gene_annot) <- rownames(expr_mat)

        cds <- new_cell_data_set(expr_mat, cell_metadata = cell_annot, gene_metadata = gene_annot)
        cds <- preprocess_cds(cds)
        cds <- reduce_dimension(cds)
        cds <- cluster_cells(cds, k = Kmonocle)

        obj.list[[i]] <- AddMetaData(obj.list[[i]], metadata = clusters(cds),
          col.name = paste0("Monocle_k", Kmonocle))

        sample_name <- names(obj.list)[i]
        cluster_vec <- clusters(cds)
        moclust_list[[sample_name]] <- setNames(
          paste0(sample_name, "_", cluster_vec), colnames(obj.list[[i]]))
        incProgress(inc)
      }

      combined_clusters <- unlist(moclust_list)
      df_clusters <- data.frame(cluster = combined_clusters, barcode = names(combined_clusters),
        stringsAsFactors = FALSE)
      df_clusters$barcode    <- sub(".*\\.", "", df_clusters$barcode)
      df_clusters$split_meta <- sub("\\..*$", "", rownames(df_clusters))

      seu@meta.data$MonocleCluster <- df_clusters$cluster[
        match(rownames(seu@meta.data), df_clusters$barcode)]

      if (all(is.na(seu$MonocleCluster))) {
        showNotification("No cells assigned to clusters.", type = "error"); return()
      }

      bulk    <- AggregateExpression(seu, group.by = "MonocleCluster", return.seurat = TRUE)
      mat     <- GetAssayData(bulk, assay = "RNA", layer = "data")
      top1000 <- names(sort(apply(mat, 1, sd, na.rm = TRUE), decreasing = TRUE))[1:1000]
      centered_mat <- scale(mat[top1000, ], scale = FALSE)
      mat_corr     <- cor(centered_mat, use = "pairwise.complete.obs")
      df_clusters$barcode <- gsub("_", "-", df_clusters$barcode)

      rv$pseudobulk_matrix(list(mat_corr = mat_corr, seu = seu, df_clusters = df_clusters))
      incProgress(1)
      append_log("Pseudobulk matrix computed.")
    })
  })

  observeEvent(list(rv$pseudobulk_matrix(), input$annot_meta), {
    req(rv$pseudobulk_matrix())
    data        <- rv$pseudobulk_matrix()
    mat_corr    <- data$mat_corr
    seu         <- data$seu
    df_clusters <- data$df_clusters

    cluster_to_sample <- tapply(df_clusters$split_meta, df_clusters$cluster,
      function(x) names(sort(table(x), decreasing = TRUE))[1])
    names(cluster_to_sample) <- gsub("_", "-", names(cluster_to_sample))

    top_annot_df <- data.frame(annot = cluster_to_sample[colnames(mat_corr)], stringsAsFactors = FALSE)
    top_vals     <- unique(top_annot_df$annot[!is.na(top_annot_df$annot)])
    top_col      <- setNames(viridis(length(top_vals), option = "C"), top_vals)
    top_annot    <- HeatmapAnnotation(top_meta = top_annot_df$annot,
      col = list(top_meta = top_col), simple_anno_size = unit(1.2, "cm"))

    req(input$annot_meta %in% colnames(seu@meta.data))
    prop_tab     <- prop.table(table(seu@meta.data[[input$annot_meta]], seu@meta.data$MonocleCluster), 2)
    cell_types   <- rownames(prop_tab)
    bottom_col   <- setNames(qualitative_hcl(length(cell_types), palette = "Dark3"), cell_types)
    bottom_annot <- HeatmapAnnotation(
      bar = anno_barplot(t(prop_tab), gp = gpar(fill = bottom_col), height = unit(6, "cm")),
      simple_anno_size = unit(1, "cm"))

    rv$pseudobulk_data(list(mat_corr = mat_corr, top_annot = top_annot,
      bottom_annot = bottom_annot, bottom_col = bottom_col))
  })

  output$pseudobulk_heatmap <- renderPlot({
    req(rv$pseudobulk_data())
    data <- rv$pseudobulk_data()
    ht <- Heatmap(data$mat_corr,
      col = colorRamp2(c(0.3, 0.6, 1), c("blue", "white", "red")),
      top_annotation = data$top_annot, bottom_annotation = data$bottom_annot,
      clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
      clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
      show_row_names = TRUE, show_column_names = TRUE, name = "Corr")
    barplot_legend <- Legend(labels = names(data$bottom_col), title = "Cell types",
      legend_gp = gpar(fill = data$bottom_col))
    draw(ht, annotation_legend_list = list(barplot_legend))
  })

  output$download_pseudobulk_pdf <- downloadHandler(
    filename = function() paste0("pseudobulk_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content  = function(file) {
      req(rv$pseudobulk_data())
      data <- rv$pseudobulk_data()
      pdf(file, width = 12, height = 15)
      ht <- Heatmap(data$mat_corr,
        col = colorRamp2(c(0.3, 0.6, 1), c("blue", "white", "red")),
        top_annotation = data$top_annot, bottom_annotation = data$bottom_annot,
        clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
        clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
        show_row_names = TRUE, show_column_names = TRUE, name = "Corr")
      barplot_legend <- Legend(labels = names(data$bottom_col), title = "Cell types",
        legend_gp = gpar(fill = data$bottom_col))
      draw(ht, annotation_legend_list = list(barplot_legend))
      dev.off()
    }
  )

} # fin server