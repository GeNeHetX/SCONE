library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(magrittr)
library(shinyjs)
library(shinyBS)

# Palette Sakura
header_bg <- "#FCE4EC"
header_text <- "#EC407A"
box_accent <- "#F48FB1"
box_dark <- "#EC407A"
btn_bg <- "#F48FB1"
btn_text <- "#fff"

ui <- dashboardPage(
  skin = "black",

  ## ---------------- Header ----------------
  dashboardHeader(
    titleWidth = 300,
    title = tags$span(
      tags$span("SCONE", style = "font-family:'Arial'; font-weight:bold;"), tags$img( src = "logo.jpg", height = "60px", style ="margin-left:100px; margin-top:-1px;" ), 
      style = "display:flex; align-items:center;" 
    ),
    tags$li(
      class = "dropdown header-menu-left",
      actionLink("show_tutorial", "Tutorial", icon = icon("book-open"))
    )
  ),
  

  ## ---------------- Sidebar ----------------
  dashboardSidebar(
    width = 300,
    useShinyjs(),
    br(),
    br(),
    br(),
    br(),
    fileInput("input_zip", "Upload scRNA zipfolder :\n 1 matrix.mtx + 1 barcodes.tsv & 1 features.tsv", accept = ".zip"),
    tags$hr(),
    h4("Seurat Settings"),
    numericInput("n_hvg", "High Variable Genes (HVG)", value = 1000),
    numericInput("n_pcs", "PCA components", value = 20),

    radioButtons("clust_method",label = div("Clustering method",span(icon("circle-info"), id="clust_info")), choices = list("Louvain" = "louvain", "Leiden"  = "leiden", "Kmeans"  = "kmeans"),inline=TRUE, selected = "leiden"),
    bsTooltip( "clust_info", title = HTML(
        "Louvain: community detection based on modularity<br/>
        Leiden: improved Louvain, better partitions<br/>
        K-means: partition cells in k clusters based on PCA")),
        
    # Paramètres conditionnels
    conditionalPanel("input.clust_method == 'louvain'", numericInput("knn", "Neighbors (k)", value = 20)),
    conditionalPanel("input.clust_method == 'leiden'", numericInput("leiden_k", "Neighbors (k)", value = 20)),
    conditionalPanel("input.clust_method == 'kmeans'", numericInput("kmeans_centers", "Centers (k)", value = 10)),
    sliderInput( "resolution", "Resolution of clustering", min = 0.1, max = 2, value = 0.5, step = 0.1),
    actionButton("run_analysis", "RUN Seurat Analysis", icon = icon("play"), style = paste0("background:", btn_bg, "; color:", btn_text, ";")),

    # tags$hr(),
    # uiOutput("loaded_info"),

    tags$hr(),
    h4("Logs"),
    verbatimTextOutput("status", placeholder = TRUE),

    tags$hr(),
    p(style="font-size:12px;color:#555;", "SCONE © — Single-Cell Omics Navigation & Exploration \n please cite authors : Camille Pignolet (GeNeHetX Team, INSERM, 2026)")
  ),

  ## ---------------- Body (CSS moved here) ----------------
  dashboardBody(
  tags$head(
    tags$style(HTML(sprintf("

      /*******************************
       * HEADER (title + navbar)
       *******************************/

      /* Couleur + style du titre */
      .main-header .logo { 
        color:%s !important; 
        font-weight:bold; 
        font-size:20px;
      }

      /* Couleur du bandeau du header */
      .main-header .navbar { 
        background:%s !important; 
      }

      /* --- AUGMENTER LA HAUTEUR DU HEADER --- */
      .main-header .logo {
        height: 90px !important;          /* hauteur du logo/title */
        line-height: 90px !important;
      }
      .main-header .navbar {
        height: 90px !important;          /* hauteur de la barre */
        min-height: 90px !important;
      }  

        /* Ajuster la hauteur du header pour ne laisser que le logo */
        .main-header .logo {
        height: auto !important;
        line-height: normal !important;
        padding-top: 10px;
        padding-bottom: 10px;
        }

        /* Ajuster le padding du contenu pour que le body commence juste après le header */
        .content-wrapper,
        .right-side,
        .main-sidebar {
        padding-top: 10px !important;  /* même hauteur que le header */
        }

    
      /* Liens du header en rose */
      .main-header .navbar .nav > li.header-menu-left > a {
        color: %s !important;
        font-weight: bold;
        font-size: 15px;
      }

      /* Survol des liens */
      .main-header .navbar .nav > li.header-menu-left > a:hover {
        color: black !important;
        background: #F8BBD0 !important;
      }

      /*******************************
       * BOXES (style rose)
       *******************************/
      .box.box-solid.box-info > .box-header {
        background:%s !important;
        color:white !important;
      }

      .box.box-solid.box-info {
        border-color:%s !important;
      }

      .box.box-solid.box-success > .box-header {
        background:#F9EAF1 !important;
        color:#333 !important;
      }

      /*******************************
       * Titre + logo: marges
       *******************************/
      .main-header .logo {
        padding-left: 15px;
        padding-right: 30px;
      }

       /**********************************
        * FILE INPUT (Upload) — rose
        **********************************/
        .btn-file {
        background: #F48FB1 !important;   /* Rose */
        color: white !important;
        border-color: #F48FB1 !important;
        }

        .btn-file:hover {
        background: #EC407A !important;   /* Rose plus foncé */
        color: white !important;
        }


        /**********************************
        * RADIO BUTTONS — rose
        **********************************/
        .radio input[type='radio'] + span:before,
        .radio-inline input[type='radio'] + span:before {
        border: 2px solid #F48FB1 !important;
        }

        .radio input[type='radio']:checked + span:after,
        .radio-inline input[type='radio']:checked + span:after {
        background-color: #EC407A !important;
        border-color: #EC407A !important;
        }


        /**********************************
        * SLIDER INPUT — curseur rose
        **********************************/
        .irs-bar,
        .irs-bar-edge {
        background: #F48FB1 !important;
        border-color: #F48FB1 !important;
        }

        .irs-slider {
        background: #EC407A !important;
        border-color: #EC407A !important;
        }

        .irs-line {
        border-color: #F48FB1 !important;
        }

        .irs-grid-pol {
        background: #F48FB1 !important;
        }


    ", header_text, header_bg, header_text, box_dark, box_accent)))
  ),
      div(
      id = "tutorial_panel",
      style = "display:none;",   # caché par défaut
      fluidRow(
        column(12,
          box(
            width = 12,
            status = "info",
            solidHeader = TRUE,
            title = "Tutorial",
            htmlOutput("tutorial_content")
          )
        )
      )
    ),

### Overview (séparé)
fluidRow(
  column(12,
    tabBox(
      width = 12, id = "overview_tabBox",
      tabPanel("Overview",
        fluidRow(
          box(width = 6, status = "info", solidHeader = TRUE,
              title = "Seurat Object",
              verbatimTextOutput("seu_summary") %>% withSpinner()
          ),
          box(width = 6, status = "info", solidHeader = TRUE,
              title = "Metadata Table",
              DT::dataTableOutput("meta_table") %>% withSpinner()
          )
        )
      )
    )
  )
),

### Main Analysis Tabs (QC, UMAP, Violin, DotPlot…)
fluidRow(
  column(12,
    tabBox(
      width = 12, id = "analysis_tabBox",

      # QC Metrics + Subset
      tabPanel("QC Metrics",
        fluidRow(
          box(width = 12, status = "info", solidHeader = TRUE,
              title = "QC Metrics",
              plotOutput("qc_violin", height = "600px") %>% withSpinner()
          ),
          box(width = 6, status = "info", solidHeader = TRUE,
              title = "nCount vs nFeature",
              plotOutput("qc_scatter1", height = "500px") %>% withSpinner()
          ),
          box(width = 6, status = "info", solidHeader = TRUE,
              title = "nCount vs percent.mt",
              plotOutput("qc_scatter2", height = "500px") %>% withSpinner()
          )
        ),
        fluidRow(
          box(width = 12, status = "info", solidHeader = TRUE,
              title = "Subset Cells Based on QC",
              fluidRow(
                column(4, numericInput("nFeature_min", "Min nFeature_RNA", value = 200)),
                column(4, numericInput("nFeature_max", "Max nFeature_RNA", value = 5000)),
                column(4, numericInput("percent_mt_max", "Max percent.mt", value = 10))
              ),
              fluidRow(
                column(12,
                       actionButton("apply_subset", "Apply Subset", icon = icon("filter"),
                                    style = paste0("background:", btn_bg, "; color:", btn_text, ";")),
                      uiOutput("subset_message"),
                      verbatimTextOutput("subset_summary") %>% withSpinner()
                )
              )
          )
        )
      ),

      # # UMAP
      # tabPanel("UMAP",
      #   fluidRow(
      #     box(width = 12, status = "info", solidHeader = TRUE,
      #         title = "UMAP plot",
      #         plotOutput("umapPlot", height = "1200px") %>% withSpinner()
      #     )
      #   )
      # ),

      # UMAP Metadata
      tabPanel("UMAP",
        fluidRow(
          box(width = 4, status = "info", solidHeader = TRUE, title = "Metadata selection",
              # uiOutput("metadata_selector"),
              selectInput("meta_color_by", "Color by metadata", choices = NULL)
          ),
          box(width = 8, status = "info", solidHeader = TRUE, title = "UMAP by metadata",
              plotOutput("umapMetaPlot", height = "900px") %>% withSpinner()
          )
        )
      ),

      # Violin
      tabPanel("Violin Plot",
        fluidRow(
          box(width = 4, status = "info", solidHeader = TRUE, title = "Violin settings",
            selectizeInput("violin_gene", "Gene",
                          choices = NULL,    # rempli côté serveur
                          multiple = TRUE,
                          options = list(placeholder='Type gene name...', maxOptions=1000)
            ),
            uiOutput("violin_gene_family_ui"),   # <- ici le select des familles
            checkboxInput("split_by", "Split by metadata", FALSE),
            uiOutput("split_meta_ui"),
            actionButton("plot_violin", "Plot Violin", icon = icon("chart-area")),
            plotOutput("umapViolinPlot", height = "600px") %>% withSpinner()
          ),
          box(width = 8, status = "info", solidHeader = TRUE, title = "Violin plot",
              plotOutput("violinPlot", height = "900px") %>% withSpinner()
          )
        )
      ),

      # DotPlot
      tabPanel("DotPlot",
        fluidRow(
          box(width = 4, status = "info", solidHeader = TRUE, title = "DotPlot settings",
              # selectizeInput pour les gènes
              selectizeInput(
                "dot_genes",
                "Gene(s)",
                choices = NULL,  # rempli côté serveur
                multiple = TRUE,
                options = list(
                  placeholder='Type gene name...',
                  maxOptions=1000
                )
              ),
              # Famille de gènes
              uiOutput("dot_gene_family_ui"),
              uiOutput("dot_split_meta_ui"),
              actionButton("plot_dot", "Plot DotPlot", icon = icon("dot-circle")),
              plotOutput("umapDotPlot", height = "600px") %>% withSpinner()
          ),
          box(width = 8, status = "info", solidHeader = TRUE, title = "DotPlot",
              plotOutput("dotPlot", height = "900px") %>% withSpinner()
          )
        )
      ),
      tabPanel("Find Markers",
        fluidRow(
          box(width = 4, status = "info", solidHeader = TRUE,
              title = "Marker Settings",
              uiOutput("marker_split_meta_ui"),
              selectInput("marker_mode","Comparison mode",
                choices = list(
                  "Cluster vs rest" = "vsrest",
                  "Cluster vs cluster" = "vscluster"
                )
              ),
              uiOutput("cluster_select_ui"),
              numericInput("top_n", "Top markers per cluster",
                          value = 10, min = 5, max = 100),
              actionButton("run_markers", "Run Marker Analysis", icon = icon("dna"),
                style = "background:#EC407A;color:white;"
              ),
              plotOutput("umapMarkerPlot", height = "600px") %>% withSpinner()
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Marker Results",
              h4("Top markers"),
              dataTableOutput("top_marker_table"),
              hr(),
              uiOutput("marker_table_title"),
              h4("Full results"),
              dataTableOutput("marker_table")
          )
        )
      ),
        tabPanel("Cluster Annotation",
              fluidRow(
                # Settings
                box(width = 4, status = "info", solidHeader = TRUE,
                    title = "Annotation Settings",
                    
                    # Choix metadata
                    uiOutput("annot_split_meta_ui"),

                    # Cluster labels + couleur dynamique
                    uiOutput("cluster_labels_ui"),

                    actionButton("apply_labels", "Label UMAP manually",
                                icon = icon("palette"),
                                style = "background:#4CAF50; color:white;"),
                    downloadButton("download_metadata_csv", "Download metadata + customize label")
                ),
                
                # UMAP
                box(width = 8, status = "info", solidHeader = TRUE,
                    title = "UMAP with annotated clusters",
                    downloadButton("download_umap_pdf", "Download customize UMAP"),
                    plotOutput("annot_umap", height = "600px") %>% withSpinner()
                )
              )
            

                    )
                )
            )
        )
    )
)