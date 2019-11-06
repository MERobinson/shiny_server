library(shiny)
library(tidyverse)
library(DT)
options(stringsAsFactors = F)

# load data
dat <- readRDS("/poolio/internal_data/surf_proteomics/GD_surfprot_BCR_ABL1_4h.2019-11-05.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   title = "SurfProt", 
   theme = "bootstrap.flatly.css",
   
   titlePanel(tags$h3(tags$a(
      imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
      href="http://10.197.211.94:80"), "Surface Proteomics")),
   
   # Show a plot of the generated distribution
   mainPanel(width = '100%',
      plotOutput('MAplot', 
                 click = clickOpts(id = "plot_click"),
                 hover = hoverOpts(id = "plot_hover"),
                 brush = brushOpts(id = "plot_brush")),
      br(),
      DT::DTOutput('protein_info_dt'),
      br(),
      textInput('id_input', 'Enter a list of genes or proteins of interest', 
                value = "", width = NULL, placeholder = NULL)
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   # plot MA
   output$MAplot <- renderPlot({
      ggplot(dat, aes(x = av_LFQ, y = l2fc, label = protein)) + 
         geom_point(aes(col = p_value < 0.05), alpha = .5, size = .8) +
         scale_color_manual(values = c("grey20","firebrick"),
                            labels = c("non-significant","p < 0.05", " ")) +
         scale_y_continuous(limits = c(-10,20)) +
         theme_bw(base_size = 16) +
         theme(panel.grid = element_blank(),
               legend.background = element_blank(),
               legend.title = element_blank(), 
               legend.key = element_blank()) 
   })
   
   # get selected proteins from graph
   sel_prot_graph <- reactive({
      brushed <- row.names(brushedPoints(dat, input$plot_brush))
      clicked <- row.names(nearPoints(dat, input$plot_click, addDist = TRUE))
      return(unique(brushed, clicked))
   })
   
   # datatable of info on selected proteins
   output$protein_info_dt <- DT::renderDT({
      sel_prot <- sel_prot_graph()
      if (is_empty(sel_prot)) {
         return()
      } else {
         seldat <- dat[sel_prot, c(1:8,10:12,14:17)]
         print(head(seldat))
         colnames(seldat) <- c("Protein ID", "UniProt ID", "Gene Symbols", "Total UniqPep",
                               "av LFQ", "av iBAQ", "av HF", "L2FC", "P-value", "Adjusted P-value",
                               "Specificity", "TM Domain", "# TM Domain", "GO Membrane Anno", "Protein Name")
         DT::datatable(seldat, rownames = F,
                       extensions = 'Buttons',
                       options = list(dom = 'Bfrtip',
                                      buttons = list(list(extend = 'copy',
                                                          title = NULL),
                                                     list(extend = 'print',
                                                          title = "Surface Proteomics - Selected Proteins"),
                                                     list(extend = 'csv',
                                                          filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                          title = NULL),
                                                     list(extend = 'excel',
                                                          filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                          title = NULL),
                                                     list(extend = 'pdf',
                                                          filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                          title = "Surface Proteomics - Selected Protein"))))
      }
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

