library(shiny)
library(tidyverse)
library(org.Hs.eg.db)
options(stringsAsFactors = F)

genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")


############### ui ############### 
ui <- fluidPage(
    
    theme = "bootstrap.flatly.css",
    
    titlePanel(div(img(src='geodude.jpg',
                       style=c("width:110px","height:110px")),
                   "SurviveR"),

    sidebarLayout(
        sidebarPanel(
            div(style="display:inline-block", seletizeInput("gene", "Gene", genes, selected = "TP53")),
            div(style="display:inline-block", actionButton("go_button", "GO")),
            div(tags$small("Enter gene symbol and press go")),
            hr(),
            uiOutput("choose_dataset"),
            div(tags$small("Select a pre-processed dataset from the list")),
            hr(),
            uiOutput("choose_datatype"),
            tags$small("Some datasets have both survival and regression data, choose which to display")
        ),

        mainPanel(
           plotOutput("survplot")
        )
    )
)

############### server ############### 
server <- function(input, output) {

    # parse gene input
    
    
    # update GPL input choices
    output$choose_dataset <- renderUI({
        gse_data <- get_gse_data()
        gpl_list <- unname(sapply(gse_data, annotation))
        gpl_list <- unname(sapply(gpl_list, function(x) {
            gpl_data <- getGEO(x)
            paste0(x, ": ", gpl_data@header$title)
        }))
        selectInput("gpl", "Platform", gpl_list)
    })
    
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

############### run ############### 
shinyApp(ui = ui, server = server)
