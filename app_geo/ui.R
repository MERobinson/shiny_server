library(shiny)
library(ggplot2)
library(ggvis)
library(dplyr)
library(DT)
library(heatmap3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotly)

shinyUI(fluidPage(
  
  theme = "bootstrap.flatly.css",

  titlePanel(tags$h3(tags$a(
    imageOutput("icon", height = "100px", width = "100px", inline = TRUE),
    href="http://10.197.211.94:3838"), 
    "GEOdude")),

  sidebarLayout(
    sidebarPanel(
      div(style="display:inline-block", textInput("gse", "GSE Acession Number")),
      div(style="display:inline-block", actionButton("go_button", "GO")),
      div(tags$small("Enter an accession number from gene expression omnibus, this should be a series ID starting GSE followed by a unique number (e.g. GSE55976). Press GO to update.")),
      hr(),
      uiOutput("geo_link"),
      hr(),
      uiOutput("choose_gpl"),
      tags$small("Some data series will contain multiple different platform types, e.g. different microarrays used, select the platform you would like to analyse (see GEO page for details)"),
      hr(),
      uiOutput("upload_counts"),
      uiOutput("upload_text"),
      hr(),
      uiOutput("group_select"),
      tags$small("This will determine groupings for plots - you can select a column from the phenotype data provided by the authors, or create your own groupings by selecting CUSTOM"),
      hr(),
      uiOutput("custom_groupings")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Experiment Info",
                 tags$h3("Abstract"),
                 htmlOutput("abstract_text")),
        tabPanel("Phenotype Info",
                 DT::dataTableOutput("pdata_table")),
        tabPanel("Features Info",
                 DT::dataTableOutput("fdata_table")),
        tabPanel("Heatmap",
                 uiOutput("max_genes"),
                 verbatimTextOutput("temp_value"),
                 plotOutput("heatmap"),
                 plotOutput("hm_leg")),
        tabPanel("PCA",
                 uiOutput("choose_pc1"),
                 uiOutput("choose_pc2"),
                 plotlyOutput("pca"),
                 verbatimTextOutput("pca_hover")),
        tabPanel("GSEA",
                 uiOutput("choose_numer"),
                 uiOutput("choose_denom"),
                 selectInput("gs_type", "Gene set input type:",
                             c("MSigDB: H - hallmark" = "H",
                               "MSigDB: C1 - positional" = "c1",
                               "MSigDB: C2 - curated" = "c2",
                               "MSigDB: C3 - motif" = "c3",
                               "MSigDB: C4 - computational" = "c4",
                               "MsigDB: C5 - GO" = "c5",
                               "MSigDB: C6 - oncogenic" = "c6",
                               "MSigDB: C7 - immunological" = "c7",
                               "CUSTOM: GMT file" = "GMT")),
                 uiOutput("gmt_upload"),
                 uiOutput("choose_id_type"),
                 selectInput("org", "Organism:",
                             c("Homo sapiens" = "Hs",
                               "Mus musculus" = "Mm")),
                 plotOutput("gsea_plot"))
      )
    )
  )
))
