library(shiny)
library(DT)
library(dplyr)
data("starwars")

ui <- fluidPage(title = "Ratios",
                sidebarLayout(
                    sidebarPanel(width = 2,
                                 actionButton("button", "Go"), # Emulates data loading
                                 sliderInput("seuil_j", "Threshold J",
                                             min = 0,  max = 80, value = 35, step = 0.5)),
                    mainPanel( 
                        fluidRow(column(width = 12,
                                        DT::dataTableOutput("ratios"))))
                )
)

server <- function(input, output, session) {
    donnees_ratios <- reactive({
        req(input$button)
        set.seed(14)
        starwars %>% 
            select(1:10) %>% # DataTables is not happy with list columns
            mutate(signe = sample(c(1, -1), replace = TRUE, size = nrow(.)),
                   ratio_j = signe * mass / height) %>% 
            select(name, mass, height, signe, ratio_j, everything())
    })
    
    output$ratios <- DT::renderDataTable({
        donnees_ratios() %>% 
            creer_DT() %>% 
            formatter_DT(input)
    })
}

creer_DT <- function(donnees) {
    datatable(donnees, 
              rownames = FALSE, 
              class = 'cell-border stripe compact hover',
              extensions = c("Buttons"),
              options = list(
                  dom = 'Blfrtip',
                  buttons = list(
                      list(extend = "pdf", 
                           exportOptions = list(stripHtml = FALSE,
                                                columns = ':visible'),
                           orientation = 'landscape'),
                      list(extend = "print", 
                           exportOptions = list(stripHtml = FALSE,
                                                columns = ':visible')),
                      "excel", "csv", "colvis"),
                  language = list(
                      decimal = ",",
                      thousands = "&#8239;"  # small unbreakable space
                  )
              )
    )
}

formatter_DT <- function(table, input) {
    table %>% 
        formatPercentage(columns = c("ratio_j"),
                         digits = 1L, dec.mark = ",", mark = "&#8239;") %>%
        formatRound(columns = c("height", "mass"),
                    digits = 1L, dec.mark = ",", mark = "&#8239;") %>%
        format_seuil("ratio_j", input$seuil_j)
}

format_seuil <- function(table, column, seuil) {
    # Threshold for the aboslute value, and different coloring if higher or lower
    formatStyle(table, column, 
                fontWeight = styleInterval(
                    c(-seuil / 100, seuil / 100), c("bold", "normal", "bold")),
                color = styleInterval(
                    c(-seuil / 100, seuil / 100), c("red", "black", "orange")
                ))
}

shinyApp(ui, server)