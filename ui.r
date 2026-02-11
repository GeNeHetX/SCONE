library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(magrittr)
library(shinyjs)

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

    radioButtons("clust_method", "Clustering method",
                 choices = list("Louvain" = "louvain", "Kmeans" = "kmeans")),

    conditionalPanel("input.clust_method == 'louvain'", numericInput("knn", "Neighbors (k)", value = 20)),
    conditionalPanel("input.clust_method == 'kmeans'", numericInput("kmeans_centers", "Centers (k)", value = 10)),

    sliderInput("resolution", "Resolution of clustering", 0.1, 2, 0.5, step = 0.1),

    actionButton("run", "RUN - Seurat Analysis", icon = icon("play"),
                 style = paste0("background:", btn_bg, "; color:", btn_text, ";")),

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


  ### Analysis

   fluidRow(
      column(12,
             tabBox(
               width = 12, id = "tabs",
              #  tabPanel("Tutorial",
              #           fluidRow(
              #               column(width = 12,
              #               htmlOutput("tutorial_content")
              #               )
              #           )
              #   ),
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
                        )
               ),
               tabPanel("UMAP",
                        fluidRow(
                          box(width = 12, status = "info", solidHeader = TRUE,
                              title = "UMAP plot",
                              plotOutput("umapPlot", height = "1200px") %>% withSpinner()
                          )
                        )
               ),
               tabPanel("UMAP Metadata",
                        fluidRow(
                          box(width = 4, status = "info", solidHeader = TRUE, title = "Metadata selection",
                              uiOutput("metadata_selector"),
                              selectInput("meta_color_by", "Color by metadata", choices = NULL)
                          ),
                          box(width = 8, status = "info", solidHeader = TRUE, title = "UMAP by metadata",
                              plotOutput("umapMetaPlot", height = "900px") %>% withSpinner()
                          )
                        )
               ),
               tabPanel("Violin Plot",
                        fluidRow(
                          box(width = 4, status = "info", solidHeader = TRUE, title = "Violin settings",
                              textInput("violin_gene", "Gene"),
                              checkboxInput("split_by", "Split by metadata", FALSE),
                              uiOutput("split_meta_ui"),
                              actionButton("plot_violin", "Plot Violin", icon = icon("chart-area"))
                          ),
                          box(width = 8, status = "info", solidHeader = TRUE, title = "Violin plot",
                              plotOutput("violinPlot", height = "900px") %>% withSpinner()
                          )
                        )
               ),
               tabPanel("DotPlot",
                        fluidRow(
                          box(width = 4, status = "info", solidHeader = TRUE, title = "DotPlot settings",
                              tags$label("Genes (comma-separated)"),
                              tags$textarea(id="dot_genes", rows=6, placeholder="CD3D,MS4A1,GNLY"),
                              br(),
                              actionButton("plot_dot", "Plot DotPlot", icon = icon("dot-circle"))
                          ),
                          box(width = 8, status = "info", solidHeader = TRUE, title = "DotPlot",
                              plotOutput("dotPlot", height = "900px") %>% withSpinner()
                          )
                        )
               ),
               tabPanel("Find Markers",
                        fluidRow(
                            box(width = 4, status = "info", solidHeader = TRUE, title = "Marker Settings",
                                selectInput("marker_mode", "Mode",
                                            choices = list(
                                            "Find all markers (FindAllMarkers)" = "all",
                                            "Cluster vs cluster" = "vscluster",
                                            "Cluster vs rest" = "vsrest"
                                            )),
                                
                                uiOutput("cluster_select_ui")  # Defined below
                                ,
                                actionButton("run_markers", "Run Marker Analysis",
                                            icon = icon("dna"),
                                            style = "background:#EC407A;color:white;")
                            ),

                            box(width = 8, status = "info", solidHeader = TRUE,
                                title = "Marker Results",
                                dataTableOutput("marker_table") %>% withSpinner()
                            )
                        )
               ),
               tabPanel("Annotate Clusters",
                        fluidRow(
                         
                        )
               )
             )
      )
    )
  )
)