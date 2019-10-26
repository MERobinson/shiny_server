library(shiny)
library(tidyverse)
library(plotly)
library(org.Hs.eg.db)
library(shinyWidgets)
library(survival)
library(survminer)
options(stringsAsFactors = F)

genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
aliases <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = genes,
                                 keytype = "SYMBOL",
                                 column = "ALIAS")

# tmp until i sort out db
# datasets <- list.files("/poolio/public_data/GEO_expr/rds", pattern = '.rds', full.names = T)
# names(datasets) <- sub("^.+/(.+)\\.[0-9\\-]+\\.rds", "\\1", datasets)
# datasets <- lapply(datasets, function(x) readRDS(x))
# saveRDS(datasets, "/poolio/public_data/GEO_expr/rds/combined_datsets.2019-10-24.rds")
datasets <- readRDS("/poolio/public_data/GEO_expr/rds/combined_datsets.2019-10-24.rds")

############### ui ############### 
ui <- fluidPage(
    
    theme = "bootstrap.flatly.css",

    titlePanel(img(src='pigeon_icon.jpg', style=c("width:90px","height:90px")),
               "PIGEONS! - PatIent Gene ExpressiON & Survival"),

    sidebarLayout(
        sidebarPanel(
            div(style="display:inline-block", textInput("gene", "Gene", value = "TP53")),
            div(style="display:inline-block", actionButton("go_button", "GO")),
            htmlOutput("gene_check"),
            uiOutput("choose_dataset"),
            br(),
            uiOutput("choose_probes"),
            tags$small("If multiple probes selected, average will be used"),
            br(), br(),
            checkboxGroupButtons("percentile", "Percentile:", c("25%","50%"), justified = T, status = "success", selected = "50%")
        ),

        mainPanel(
          plotOutput("survplot_median"),
          plotOutput("survplot_quart")
        )
    )
)

############### server ############### 
server <- function(input, output) {

    # get gene
    get_gene <- reactive({
      input$go_button
      return(isolate(input$gene))
    })

    # check gene
    output$gene_check <- renderUI({
      gene <- get_gene()
      if (is.null(gene) | gene == "") {
        return(tags$small("Gene field is empty - enter a gene name"))
      } else {
        if (gene %in% genes) {
          return(NULL)
        } else {
          palias <- aliases[which(aliases$ALIAS == gene),]$SYMBOL
          pgenes <- genes[grep(gene, genes, ignore.case = T)]
          pdist <- adist(gene, pgenes)
          pgenes <- pgenes[which(pdist == min(pdist))]
          pmatch <- unique(c(palias, pgenes))
          pmatch <- ifelse(is_empty(pmatch), "none found", paste0(pmatch, collapse = ", "))
          return(tags$small(paste0("Gene symbol is not valid, potential matches: ", pmatch)))
        }
      }
    })

    # dataset choice ui
    output$choose_dataset <- renderUI({
      selectInput("dataset", "Dataset", names(datasets))
    })

    # get dataset
    get_dat <- reactive({
      datname <- input$dataset
      return(datasets[[datname]])
    })

    # probe choice ui
    output$choose_probes <- renderUI({
      dat <- get_dat()
      gene <- get_gene()
      if (is_empty(dat) | is_empty(gene)) {
        return()
      } else {
        probes <- fData(dat)[fData(dat)$symbol == gene,]
      }
      validate(need(nrow(probes) > 0, "Gene does not exist in selected dataset"))
      checkboxGroupButtons("probes", "Probe(s):", rownames(probes), justified = T,
                           status = "info", selected = rownames(probes)[1],
                           checkIcon = list(yes = icon("ok", lib = "glyphicon"),
                                            no = icon("remove", lib = "glyphicon")))
    })

    # get probe ids
    get_expr <- reactive({
      dat <- get_dat()
      gene <- get_gene()
      probes <- input$probes
      if (length(probes) > 1) {
        expr <- colMeans(assay(dat)[rownames(dat) %in% probes,])
      } else {
        expr <- assay(dat)[rownames(dat) %in% probes,]
      }
      return(expr)
    })

    # get exprcats
    get_exprcats <- reactive({
      percs <- input$percentile
      expr <- get_expr()
      cats <- list()
      if ("50%" %in% percs) {
        th <- quantile(expr, .5)
        cats$median <- ifelse(expr > th, "Upper", "Lower")
      }
      if ("25%" %in% percs) {
        lwth <- quantile(expr, .25)
        upth <- quantile(expr, .75)
        cats$quart <- ifelse(expr > upth, "Upper", NA)
        cats$quart[expr < lwth] <- "Lower"
      }
      return(cats)
    })

    # organise survival data
    get_osdat <- reactive({
      dat <- get_dat()
      cats <- get_exprcats()
      lapply(cats, function(cat) {
        osdat <- data.frame(time = dat$years_to_death,
                            events = dat$censor_death,
                            expr = cat)
        osdat$expr <- factor(osdat$expr, levels = c("Lower", "Upper"))
        return(osdat)
      })
    })

    # median survival plot
    output$survplot_median <- renderPlot({
      osdat <- get_osdat()
      if ("median" %in% names(osdat)) {
        fit <- survfit(Surv(time, events) ~ expr, data = osdat[["median"]])
        ggsurvplot(fit,
                   data = osdat[["median"]],
                   pval = T,
                   ylab = "OS probability",
                   xlab = "Time [years]",
                   censor = F,
                   legend.title = "",
                   legend.labs = c("Lower (50%)","Upper (50%)"),
                   legend = "bottom",
                   risk.table = F,
                   palette = c("forestgreen","firebrick"),
                   axes.offset = F)
      } else {
        return()
      }
    })

    # quartile survival plot
    output$survplot_quart <- renderPlot({
      osdat <- get_osdat()
      if ("quart" %in% names(osdat)) {
        fit <- survfit(Surv(time, events) ~ expr, data = osdat[["quart"]])
        ggsurvplot(fit,
                   data = osdat[["quart"]],
                   pval = T,
                   ylab = "OS probability",
                   xlab = "Time [years]",
                   censor = F,
                   legend.title = "",
                   legend.labs = c("Lower (25%)","Upper (25%)"),
                   legend = "bottom",
                   risk.table = F,
                   palette = c("forestgreen","firebrick"),
                   axes.offset = F)
      } else {
        return()
      }
    })
}

############### run ############### 
shinyApp(ui = ui, server = server)
