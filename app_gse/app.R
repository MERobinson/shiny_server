library(shiny)
library(tidyverse)
library(clusterProfiler)
options(stringsAsFactors = F)

# load msigdb gene sets
rdata_files <- list.files(path = "/poolio/resources/genesets", pattern = "*.rdata",
                          recursive = T, full.names = T)
lapply(rdata_files, load, .GlobalEnv)

# function to generate pptx
gen_pptx <- function(plot, file) {
   read_pptx() %>% 
      add_slide(layout = "Title and Content", master = "Office Theme") %>% 
      ph_with_vg(doc, ggobj = plot, type = 'body') %>% 
      print(target = file)
}

ui <- fluidPage(
   
   title = "ssPEA", 
   theme = "bootstrap.flatly.css",
   
   titlePanel(tags$h3(tags$a(
      imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
      href="http://10.197.211.94:3838"), 
      "Pathway Enrichment Anlysis")),
   
   tabsetPanel(
      tabPanel("Data Upload",
               selectInput('input_type', "Select Dataset or Upload:", 
                           choices = list("")),
               DT::DTOutput('sample_info_dt'))
   )
)

server <- function(input, output) {
   output$MAplot <- renderPlotly({
      ggplot(dat, aes(x = AveExpr, y = logFC, label = protein)) + 
         geom_point(aes(col = P.Value < 0.05), alpha = .5, size = .8) +
         scale_color_manual(values = c("grey20","firebrick"),
                            labels = c("non-significant","p < 0.05", " ")) +
         scale_y_continuous(limits = c(-10,20)) +
         theme_bw(base_size = 16) +
         theme(panel.grid = element_blank(),
               legend.background = element_blank(),
               legend.title = element_blank(), 
               legend.key = element_blank(), legend) 
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

