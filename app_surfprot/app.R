library(shiny)
library(tidyverse)
library(plotly)
options(stringsAsFactors = F)

# load data
dat <- readRDS("/poolio/internal_data/surf_proteomics/GD_surfprot_BCR_ABL1_4h.2019-10-30.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Surface Proteomics"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotlyOutput("MAplot")
      )
   )
)

# Define server logic required to draw a histogram
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

