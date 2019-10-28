library(shiny)
library(tidyverse)
library(MultiAssayExperiment)
library(plotly)
library(DT)
library(org.Hs.eg.db)
library(shinyWidgets)
library(officer)
library(rvg)
library(colourpicker)
library(randomcoloR)
options(stringsAsFactors = F)

# load CCLE data
ccledat <- readRDS("/poolio/public_data/CCLE/CCLE_MAE_processed.2019-10-27.rds")

# just make a summary table for the assay info
# assay_info <- data.frame("Assay Type" = names(ccledat$mae),
#                          "Sample #" = sapply(names(ccledat$mae), function(x) {
#                            ncol(assay(ccledat$mae[,,x])) }),
#                          "Feature #" = sapply(names(ccledat$mae), function(x) {
#                            nrow(assay(ccledat$mae[,,x])) }),
#                          check.names = F)
# saveRDS(assay_info, "assay_info.rds")
assay_info <- readRDS("assay_info.rds")

# organise tsne dim data
tsne_dims <- lapply(names(ccledat$dimred), function(x)  {
  y <- cbind(ccledat$dimred[[x]], colData(ccledat$mae)[ccledat$cc,]) %>%
    as.data.frame()
  y$assay_type <- x
  return(y)
}) %>% bind_rows()

# function to generate pptx
gen_pptx <- function(plot, file) {
  read_pptx() %>% 
    add_slide(layout = "Title and Content", master = "Office Theme") %>% 
    ph_with_vg(doc, ggobj = plot, type = 'body') %>% 
    print(target = file)
}

############### ui ############### 
ui <- fluidPage(
    
  theme = "bootstrap.flatly.css",

  titlePanel(tags$h3(tags$a(
    imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
                href="http://10.197.211.94:3838"), 
    "Cancer Cell Line Encyclopedia")),
  
  tabsetPanel(
    tabPanel("Sample Info",
             plotOutput("pheno_tsne_plot"),
             br(),
             DT::DTOutput('sample_info_dt')),
    tabPanel("Assay Info",
             DT::DTOutput('assay_info_dt'),
             br(), br(),
             DT::DTOutput('feature_info_dt'),
             br(), br(),
             plotOutput("violin_plot")),
    # tabPanel("Tissue-level",
    #          sidebarLayout(
    #            sidebarPanel(
    #              
    #            ),
    #            mainPanel(
    #              plotOutput("tissue_level_plot"),
    #              div(style="display:inline-block", 
    #                  downloadButton("dl_tissue_level_ppt", label = "Download PPT")),
    #              div(style="display:inline-block", 
    #                  downloadButton("dl_tissue_level_png", label = "Download PNG")),
    #              br(),
    #              DT::DTOutput('tisse_level_dt')
    #            )
    #          )
    #         ),
    tabPanel("tSNE",
             # sidebarLayout(
             #   # sidebarPanel(
             #   #   selectInput("assay_type", "Choose Assay:",
             #   #               c("RPPA","RNAseq","Metabolites","Chromatin"),
             #   #               selected = "Protein", multiple = FALSE,
             #   #               selectize = TRUE, width = NULL, size = NULL)#,
             #   #   # div(style="display:inline-block", textInput("target_id", "Target Name",
             #   #   #                                             value = "beta-Catenin")),
             #   #   # div(style="display:inline-block", actionButton("go_button", "GO")),
             #   #   # htmlOutput("gene_check"),
             #   #   # br(),
             #   #   # div(style="display:inline-block",
             #   #   #     colourInput("col_expr_low","Low expr", value = "#3498dbff",
             #   #   #                 allowTransparent = TRUE)),
             #   #   # div(style="display:inline-block",
             #   #   #     colourInput("col_expr_high","High expr", value = "#C23C3C",
             #   #   #                 allowTransparent = TRUE)) ),
             #   # ),
               mainPanel(width = "100%",
                 plotOutput("expr_tsne_plot"),
                 # div(style="display:inline-block",
                 #     downloadButton("dl_sample_level_ppt", label = "Download PPT")),
                 br(),
                 div(style="display:inline-block",
                     downloadButton("dl_expr_tsne_png", label = "Download PNG")),
                 # br(),
                 # DT::DTOutput('per_sample_expr_dt')
             )
          # )
    )
  )
)

############### server ############### 
server <- function(input, output) {
  
  # icon
  output$icon <- renderImage(list(src = "www/pigeon_icon.jpg",
                                  height = "50px", width = "50px"), 
                             deleteFile = F)
  
  # assay summary table
  output$assay_info_dt <- DT::renderDT(
    DT::datatable(assay_info, editable = "cell", escape = F, rownames = F,
                  selection = "single")
  )
  
  # get assay type
  get_assay_type <- reactive({
    return(names(ccledat$mae)[input$assay_info_dt_rows_selected])
  })
  
  # features table
  output$feature_info_dt <- DT::renderDT({
    assay_type <- get_assay_type()
    if (is_empty(assay_type)) {
      return()
    } else {
      DT::datatable(ccledat$fdat[[assay_type]], rownames = F,
                    selection = "single")
    }
  })
  
  # get feature
  get_feature <- reactive({
    assay_type <- get_assay_type()
    if (is_empty(assay_type)) {
      return()
    } else {
      return(ccledat$fdat[[assay_type]][input$feature_info_dt_rows_selected,1])
    }
  })
  
  # samples table
  output$sample_info_dt <- DT::renderDT(
    DT::datatable(as.data.frame(colData(ccledat$mae)),
                  editable = "cell", escape = F, rownames = F)
  )
  
  # plot simple violin plot
  output$violin_plot <- renderPlot({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) {
      return()
    }
    plotdat <- tsne_dims
    plotdat$expr <- scale(as.numeric(assay(ccledat$mae[feature,ccledat$cc,assay_type])))[,1]
    plotdat_summary <- plotdat %>% group_by(Group) %>%
      summarise(mean = mean(expr)) %>%
      arrange(-mean)
    print(head(plotdat_summary))
    plotdat$Group <- factor(plotdat$Group, levels = plotdat_summary$Group)
    print(head(plotdat))
    ggplot(plotdat, aes(x = Group, y = expr, fill = Group)) +
      geom_violin(draw_quantiles = T, trim = F, alpha = .5) +
      theme_bw(base_size = 18) +
      xlab("Sample Type") + ylab(paste0(feature," Level")) +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # plot tSNE of pheno dat
  output$pheno_tsne_plot <- renderPlot({
    palette <- distinctColorPalette(length(unique(tsne_dims$Group)))
    ggplot(tsne_dims, aes(x = V1, y = V2)) +
      geom_point(aes(fill = Group), col = "black", shape = 21, stroke = .2) +
      facet_grid(cols = vars(assay_type)) +
      scale_fill_manual(values = palette, name = "Sample Type") +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank())
  })
  
  # plot tSNE of expr dat
  get_expr_tsne_plot <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) {
      return()
    }
    plotdat <- tsne_dims
    plotdat$expr <- scale(as.numeric(assay(ccledat$mae[feature,ccledat$cc,assay_type])))[,1]
    if (any(plotdat$expr > 3)) {
      plotdat[plotdat$expr > 3,]$expr <- 3
    }
    if (any(plotdat$expr < -3)) {
      plotdat[plotdat$expr < -3,]$expr <- -3
    }
    colnames(plotdat)[c(1,2)] <- c("dim1","dim2")
    ggplot(plotdat, aes(x = dim1, y = dim2, label = Name)) +
      geom_point(aes(fill = expr), col = "black", shape = 21, stroke = .2) +
      facet_grid(cols = vars(assay_type)) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           name = paste0(feature, " expr\n[z-score]"), limits = c(-3,3)) +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank())
  })
  
  # render tSNE of expr dat
  output$expr_tsne_plot <- renderPlot({
    get_expr_tsne_plot()
  })

  # download expression tSNE as png
  output$dl_expr_tsne_png <- downloadHandler(
    filename = function() {
      paste0("CCLE_tSNE_", get_feature(), "_levels.", Sys.Date(), ".png")
    },
    content = function(file) {
      plot <- get_expr_tsne_plot()
      ggsave(plot = plot, filename = file, height = 4, width = 18)
    }
  )
}

############### run ############### 
shinyApp(ui = ui, server = server)
