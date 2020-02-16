library(shiny)
library(tidyverse)
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(officer)
library(rvg)
library(pzfx)
library(ggpubr)
library(colourpicker)
library(ggnewscale)
library(readxl)
options(stringsAsFactors = F)

################################# Setup #######################################

# function to generate pptx
gen_pptx <- function(plot, file, height = 5, width = 5, left = 1, top = 1) {
    read_pptx() %>% 
        add_slide(layout = "Title and Content", master = "Office Theme") %>% 
        ph_with_vg_at(doc, ggobj = plot, type = 'body', 
                      height = height, width = width,
                      left = left, top = top) %>% 
        print(target = file)
}

# load internal datasets
internal_dat <- list.files("data", pattern = ".xlsx", full.names = T)
names(internal_dat) <- sub("^data/(.+)\\.xlsx", "\\1", internal_dat)

#################################### UI #######################################
ui <- function(request) {
    fluidPage(
        
        # title
        title = "scWB", 
        theme = shinytheme("flatly"),
        
        titlePanel(tags$h3(tags$a(
            imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
            href="http://10.197.211.94:3838"), 
            "Prism")),
        
        sidebarLayout(
            sidebarPanel(width = 3,
                         selectInput(inputId, label, choices, selected = NULL, multiple = FALSE,
                                     selectize = TRUE, width = NULL, size = NULL)
                         selectInput(inputId = "existing_dat",
                                     label = "Choose a dataset:",
                                     choice = names(internal_df)),
                         fileInput(inputId = "input_file",
                                   label = "Upload an excel file:",
                                   multiple = F,
                                   accept = c(".xlsx")),
                         div(style="display:inline-block;width:45%;", 
                             colourInput("col_selected","Low", value = "#202080",
                                         allowTransparent = TRUE)),
                         div(style="display:inline-block;width:45%;", 
                             colourInput("col_other","High", value = "#000000",
                                         allowTransparent = TRUE)),
                         textInput("ylab_text", "Y-axis label:"),
                         textInput("xlab_text", "X-axis label:")
            ),
            mainPanel(width = 9,
                      plotOutput("main_plot") %>% withSpinner(color="#0dc5c1"),
                      div(style="display:inline-block;",
                          downloadButton("dl_main_plot_ppt", label = "PPT",
                                         style = "font-size:12px;height:30px;padding:5px;"),
                          downloadButton("dl_main_plot_png", label = "PNG",
                                         style = "font-size:12px;height:30px;padding:5px;"),
                          bookmarkButton(label = "Save Session", 
                                         style = "font-size:13px;height:30px;padding:5px;")),
                      br(), br(),
                      tags$h4("Summary Stats:"),
                      DT::DTOutput('stats_dt'),
                      br(), br(),
                      tags$h4("Input Data:"),
                      DT::DTOutput('data_dt')
            )
        )
        
    )
}


################################# Server ######################################
server <- function(input, output) {
    
    # icon
    output$icon <- renderImage(list(src = "../hexagons/prism.png",
                                    height = "84px", width = "72px"), 
                               deleteFile = F)
    
    # load prism data
    load_data <- reactive({
        if (is_empty(input$input_file)) return(NULL)
        ext <- sub("^(.+)\\.(xlsx)","\\2",input$input_file$name)
        if(ext == "xlsx") {
            df <- as.data.frame(read_excel(path = input$input_file$datapath))
        }
        return(df)
    })
    
    # get name from file
    get_name <- reactive({
        sub("^(.+)\\.(pzfx|csv|xlsx)","\\1",input$input_file$name)
    })
    
    # show data
    output$data_dt <- DT::renderDT(
        DT::datatable(as.data.frame(load_data()), 
                      escape = F, 
                      rownames = F,
                      options = list(dom = 't',
                                     pagelength = 5,
                                     lengthMenu = list(c(5, 10, -1),
                                                       c('5', '10', 'All'))))
    )
    
    # run aov
    get_stats <- reactive({
        df <- load_data()
        if (is_empty(df)) return(NULL)
        statdat <- na.omit(gather(df, "condition", "value"))
        statdat$condition <- factor(statdat$condition, levels = colnames(df))
        
        return(stat_res)
    })
    
    # show data
    output$stats_dt <- DT::renderDT(
        DT::datatable(get_stats(), 
                      escape = F, 
                      rownames = F,
                      extensions = 'Buttons',
                      options = list(dom = 'frtipB',
                                     buttons = list('pageLength',
                                                    list(extend = 'copy',
                                                         title = NULL),
                                                    list(extend = 'csv',
                                                         filename = paste0(get_name(), "_stat_res.csv"),
                                                         title = NULL),
                                                    list(extend = 'excel',
                                                         filename = paste0(get_name(), "_stat_res.xlsx"),
                                                         title = NULL))),
                      caption = "Select rows to display significance bar.")
    )
    
    # generate dotplot
    gen_main_plot <- reactive({
        df <- load_data()
        if (is_empty(df)) return(NULL)
        plotdat <- na.omit(gather(df, "category", "value"))
        plotdat$category <- factor(plotdat$category, 
                                   levels = colnames(df))
        select_col <- colnames(df)[(input$data_dt_columns_selected+1)]
        plotdat$selected <- plotdat$category %in% select_col
        
        # scaling
        if (input$scale_opt == T & input$control_type == "single") {
            plotdat$value <- scale(plotdat$value)
        } 
        if (input$scale_opt == T & input$control_type == "paired") {
            plotdat <- lapply(seq(1, ncol(df), by = 2), function(x) {
                conds <- colnames(df)[c(x,(x+1))]
                tmpdat <- plotdat[plotdat$category %in% conds,]
                tmpdat$value <- tmpdat$value / mean(tmpdat$value)
                return(tmpdat)
            }) %>% bind_rows()
        } 
        
        # basic plot outline
        ymax <- max(plotdat$value, na.rm = T)
        ymin <- min(plotdat$value, na.rm = T)
        yrange <- ymax-ymin
        p <- ggplot(plotdat, aes(x = category, y = value)) +
            coord_cartesian(clip = "off", expand = c(1,0), 
                            ylim = c(floor(ymin - (yrange/20)),
                                     ceiling(ymax + (yrange/20)))) +
            scale_fill_manual(values = c(input$col_other, input$col_selected),
                              name = input$selected_cond_label) +
            ylab(input$ylab_text) +
            theme(panel.grid = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_line(size = .8, color = "black"),
                  axis.ticks.length = unit(.25, "cm"),
                  axis.line.y = element_line(size = .8, color = "black"),
                  axis.line.x = element_line(size = .8, color = "black"),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.title.x = element_blank(),
                  axis.text = element_text(color = "black", family = "Arial", size = 18),
                  text = element_text(color = "black", family = "Arial", size = 18),
                  plot.margin = unit(c(5,15,2,2), "lines"),
                  legend.position = c(1.2,.5))
        
        # different plot type elements
        if (input$plot_type == "dotplot") {
            p <- p + 
                geom_dotplot(aes(fill = selected),
                             bins = 50,
                             binaxis='y', stackdir='center', 
                             dotsize = (1.2 + 1/yrange),
                             lwd = .6, alpha = .8) +
                stat_summary(fun.y = mean, geom = "tile",
                             width = .3, height = (yrange/75), 
                             fill = "black") +
                stat_summary(fun.y = mean,
                             fun.ymin = function(x) mean(x) - sd(x),
                             fun.ymax = function(x) mean(x) + sd(x),
                             geom = "errorbar", width = .2,  height = .15) 
        } else if (input$plot_type == "violin") {
            p <- p + 
                geom_violin(aes(fill = selected))
        } else if (input$plot_type == "boxplot") {
            p <- p + 
                geom_boxplot(aes(fill = selected), width = .75)
        }
        
        # different significance label elements
        if (input$sig_style == "bars") {
            stat_res <- get_stats()[input$stats_dt_rows_selected,]
            comps <- as.list(split(stat_res[,1:2], seq(nrow(stat_res))))
            comps <- unname(lapply(comps, function(x) unname(as.character(x))))
            p <- p + ggpubr::stat_compare_means(comparisons = comps, 
                                                label = "p.format",
                                                size = 6,
                                                bracket.size = 1.1)
        } else if (input$sig_style == "heatmap") {
            sigdat <- get_stats()
            sigdat$p_value <- as.numeric(sigdat$p_value)
            sigdat$nl10pval <- -log10(sigdat$p_value) * sign(sigdat$Difference)
            if (input$control_type == "single") {
                denom_lvl <- colnames(df)[1]
                sigdat <- sigdat[sigdat$Denominator == denom_lvl,]
                sigdat <- rbind(data.frame(Numerator = denom_lvl,
                                           Denominator = denom_lvl,
                                           Difference = 0,
                                           CI_lower = 0,
                                           CI_upper = 0,
                                           p_value = 1,
                                           nl10pval = 0),
                                sigdat)
                sigdat$y <- ceiling(ymax + (yrange/20)) + (yrange/20)
                p <- p + 
                    new_scale_fill() +
                    geom_tile(data = sigdat, color = "black", lwd = (yrange/20), height = (yrange/12),
                              aes(x = Numerator, y = y, fill = nl10pval)) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                         midpoint = 0, name = "Signed\np-value\n[-log10]")
            } else if (input$control_type == "paired") {
                denom_lvl <- colnames(df)[seq(1, ncol(df), 2)]
                numer_lvl <- colnames(df)[seq(2, ncol(df), 2)]
                matched_idx <- sapply(seq_along(numer_lvl), function(x) {
                    which(sigdat$Denominator == denom_lvl[x] &
                              sigdat$Numerator == numer_lvl[x]) })
                sigdat <- sigdat[matched_idx,]
                sigdat$x <- sapply(sigdat$Numerator, function(x) which(colnames(df) == x) - .5)
                sigdat$y <- ceiling(ymax + (yrange/20)) + (yrange/20)
                print(sigdat)
                p <- p + 
                    new_scale_fill() +
                    geom_tile(data = sigdat, color = "black", lwd = (yrange/20), height = (yrange/12),
                              aes(x = x, y = y, fill = nl10pval), width = 2) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                         midpoint = 0, name = "Signed\np-value\n[-log10]")
            }
        } else if (input$sig_style == "effect_size") {
            if (input$control_type == "single") {
                denom_lvl <- colnames(df)[1]
                effsizes <- lapply(colnames(df)[-1], function(x) {
                    tmp <- effsize::cohen.d(plotdat[plotdat$category == x,]$value,
                                            plotdat[plotdat$category == denom_lvl,]$value)
                    data.frame(condition = x, d = tmp$estimate) }) %>% bind_rows()
                effsizes <- rbind(data.frame(condition = denom_lvl, d = 0), effsizes)
                effsizes$y <- ceiling(ymax + (yrange/20)) + (yrange/20)
                p <- p + 
                    new_scale_fill() +
                    geom_tile(data = effsizes, color = "black", lwd = (yrange/20), height = (yrange/12),
                              aes(x = condition, y = y, fill = d)) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                         midpoint = 0, name = "Cohen's d")
            } else if (input$control_type == "paired") {
                denom_lvl <- colnames(df)[seq(1, ncol(df), 2)]
                numer_lvl <- colnames(df)[seq(2, ncol(df), 2)]
                effsizes <- lapply(seq_along(numer_lvl), function(x) {
                    tmp <- effsize::cohen.d(plotdat[which(plotdat$category == numer_lvl[x]),]$value,
                                            plotdat[which(plotdat$category == denom_lvl[x]),]$value)
                    data.frame(x = (which(colnames(df) == numer_lvl[x]) - 0.5),
                               d = tmp$estimate) }) %>% bind_rows()
                effsizes$y <- ceiling(ymax + (yrange/20)) + (yrange/20)
                p <- p + 
                    new_scale_fill() +
                    geom_tile(data = effsizes, color = "black", lwd = (yrange/20), height = (yrange/12),
                              aes(x = x, y = y, fill = d), width = 2) +
                    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                         midpoint = 0, name = "Cohen's d")
            }
        }
        return(p)
    })
    
    # output main plot
    output$main_plot <- renderPlot({
        plot <- gen_main_plot()
        validate(need(!is_empty(plot), "Upload a file."))
        return(plot)
    })
    
    # download main plot as png
    output$dl_main_plot_png <- downloadHandler(
        filename = function() {
            name <- get_name()
            paste0(name, "_levels_",input$plot_type, "_", input$sig_style,  ".png")
        },
        content = function(file) {
            plot <- gen_main_plot()
            df <- load_data()
            height <- ifelse(input$sig_style == "bars", 6.5, 6)
            # width <- ifelse(input$sig_style == "bars", (ncol(df) * 0.8)+5,
            #                 (ncol(df) * 0.8)+5)
            ggsave(plot = plot, filename = file, 
                   height = height, width = (ncol(df) * 0.8)+5)
        }
    )
    
    # download main plot as ppt
    output$dl_main_plot_ppt <- downloadHandler(
        filename = function() {
            name <- get_name()
            paste0(name, "_levels_",input$plot_type, "_", input$sig_style, ".pptx")
        },
        content = function(file) {
            plot <- gen_main_plot()
            df <- load_data()
            file_pptx <- tempfile(fileext = ".pptx")
            height <- ifelse(input$sig_style == "bars", 6.5, 6)
            # width <- ifelse(input$sig_style == "bars", (ncol(df) * 0.8)+5,
            #                 (ncol(df) * 0.8)+5)
            gen_pptx(plot, file_pptx, 
                     height = height, width = (ncol(df) * 0.8)+5)
            file.rename(from = file_pptx, to = file)
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking = "server")
