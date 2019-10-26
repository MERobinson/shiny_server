library(shiny)
library(GEOquery)
library(ggplot2)
library(ggvis)
library(dplyr)
library(DT)
library(heatmap3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotly)

# load msigdb gene sets
rdata_files <- list.files(path = "~/Resources/Gene_sets", pattern = "*.rdata",
                          recursive = T, full.names = T)
lapply(rdata_files, load, .GlobalEnv)

# function to generate palette
select_colours <- function(n, palette = "default") {
  palettes <- list(
    "default" = c('chartreuse3', 'cornflowerblue',
                  'darkgoldenrod1', 'peachpuff3',
                  'mediumorchid2', 'turquoise3',
                  'wheat4', 'slategray2')
  )
  if (n > length(palettes[[palette]])) {
    print(paste0("Warning: insufficient colours in selected ",
          "palette, returning randomly chosen colours"))
    return(sample(colors(), n))
  }
  return(palettes[[palette]][1:n])
}

# function to read .gmt gene set files
read_gmt <- function(gmt_file) {
  gene_set <- scan(gmt_file, what="", sep="\n")
  gene_set <- strsplit(gene_set, "\\t+")
  names(gene_set) <- sapply(gene_set, `[[`, 1)
  gene_set <- lapply(gene_set, `[`, -c(1,2))
  gene_set <- lapply(gene_set, function(x) { unique(x) })
  return(gene_set)
}

# function to convert a set list into a T2G for clusterProfiler input
term_to_gene <- function(gene_set_list) {
  term2genes <- data.frame(term = c(rep(names(gene_set_list), lengths(gene_set_list))),
                           gene = unlist(gene_set_list))
}

# function to create seperate legend
leg_func <- function (legend, lwd = 3, cex = 1.1, col) {
  plot(0, xaxt = "n", bty = "n", yaxt = "n",
       type = "n", xlab = "", ylab = "")
  legend("center", legend = legend, lwd = lwd, col = col,
          bty = "n", cex = cex)
}

# function for calc ratio of geometric means
gm_ratios <- function(exprs, idx_A, idx_B) {
  exprs <- exprs + abs(min(exprs))
  apply(exprs, 1, function(x) {
    exp(mean(log(x[idx_A]), na.rm = T)) / exp(mean(log(x[idx_B]), na.rm = T))
  })
}


############## server ############### 
shinyServer(function(input, output) {
  
  # get GSE data
  get_gse_data <- reactive({
    input$go_button
    gse <- isolate(input$gse)
    if (is.null(gse) | gse == "") {
      return()
    } else {
      withProgress(message = 'Fetching data from GEO', value = .3, { 
        getGEO(gse) 
      })
    }
  })
  
  # link to GEO
  output$geo_link <- renderUI({
    if (input$gse == "") {
      return()
    } else {
      tags$a(href=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$gse),
             "Link to GEO record")
    }
  })
  
  # update GPL input choices
  output$choose_gpl <- renderUI({
    gse_data <- get_gse_data()
    gpl_list <- unname(sapply(gse_data, annotation))
    gpl_list <- unname(sapply(gpl_list, function(x) {
      gpl_data <- getGEO(x)
      paste0(x, ": ", gpl_data@header$title)
    }))
    selectInput("gpl", "Platform", gpl_list)
  })
  
  # get GPL data
  get_gpl_data <- reactive({
    gse_data <- get_gse_data()
    gpl_list <- unname(sapply(gse_data, annotation))
    if (grepl("^GPL", input$gpl)) {
      gpl_input <- sub("^([GPL0-9]+):.+$", "\\1", input$gpl)
      gpl_data <- gse_data[[match(gpl_input, gpl_list)]]
    } else {
      return()
    }
    # if custom data given 
    if (!is.null(input$counts_file)) {
      count_data <- read.delim(input$counts_file$datapath, row.names = 1)
      pd <- new("AnnotatedDataFrame", data=pData(gpl_data))
      fd <- bitr(rownames(count_data), 
                 fromType = "ENTREZID",
                 toType = c("SYMBOL","ENSEMBL","REFSEQ",
                            "ACCNUM","UNIPROT"),
                 OrgDb = org_db)
      fd <- new("AnnotatedDataFrame", data=fd)
    }
    return(gpl_data)
  })
  
  # add ui for counts upload if no expr data
  output$upload_counts <- renderUI({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      expr_data <- exprs(gpl_data)
      if (typeof(expr_data) == "logical") {
        fileInput("counts_file", "Upload count matrix:",
                  accept=c('text/plain')) 
      } else {
        return()
      }
    }
  })
  
  # text for counts upload
  output$upload_text <- renderUI({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      expr_data <- exprs(gpl_data)
      if (typeof(expr_data) == "logical") {
        tags$small("If this panel is displayed it indicates that the selected dataset does not provide expression data as part of the GEO database. Count data may be manually uploaded if available (e.g. from supplementary files on GEO page), but make sure that the rownames (i.e. sample names) of your count matrix match the rownames of the phenotype data from the GEO database (as shown in left-most column of table under the Phenotype Info tab).")
      } else {
        return()
      }
    }
  })
  
  ### Change these to inputs for ID type and org for custom count table upload
  # # gmt id types
  # output$choose_id_type <- renderUI({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     selectInput("gmt_id_type", "Choose ID type:",
  #                 keytypes(org.Hs.eg.db))
  #   } else {
  #     return()
  #   }
  # })
  # 
  # # match anno columns
  # output$match_anno1 <- renderUI({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     selectInput("gmt_id_type", "Choose ID type:",
  #                 ketypes(org.Hs.eg.db))
  #   } else {
  #     return()
  #   }
  # })
  
  # select pdata col to group data by (or custom)
  output$group_select <- renderUI({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      selectInput("group_cname", "Groupings:", "") 
    } else {
      gpl_data <- pData(gpl_data)
      selectInput("group_cname", "Groupings:",
                  c("CUSTOM",colnames(gpl_data))) 
    }
  })
  
  # dynamic number of text inputs for custom groupings
  output$custom_groupings <- renderUI({
    if (input$group_cname == "CUSTOM") {
      gpl_data <- get_gpl_data()
      gpl_data <- pData(gpl_data)
      lapply(seq_along(rownames(gpl_data)), function(i) {
        textInput(paste0("group",i), 
                  rownames(gpl_data)[i],
                  "")
      })
    } else {
      return()
    }
  })
  
  # abstract text
  output$abstract_text <- renderUI({
    gpl_dat <- get_gpl_data()
    if (is.null(gpl_dat)) {
      return(tags$h5("Enter GEO accession number in left panel"))
    } else {
      return(tags$h5(experimentData(gpl_dat)@abstract))
    }
  })
  
  # pheno table
  output$pdata_table <- DT::renderDataTable({
    gpl_dat <- get_gpl_data()
    if (is.null(gpl_dat)) {
      return()
    } else {
      gpl_pdat <- pData(gpl_dat)
      # remove columns of no interest
      useless_cols <- apply(gpl_pdat, 2, function(x) length(unique(x)) == 1)
      gpl_pdat <- gpl_pdat[,!useless_cols]
      useless_cols <- colnames(gpl_pdat) %in% c("supplementary_file","data_row_count")
      gpl_pdat <- gpl_pdat[,!useless_cols]
      DT::datatable(
        data = gpl_pdat,
        options = list(columnDefs = list(list(
          targets = 1:ncol(gpl_pdat),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 30 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
            "}")
        ))),
        callback = JS('table.page(3).draw(false);'))
    }
  })
  
  # summary table
  output$fdata_table <- DT::renderDataTable({
    gpl_dat <- get_gpl_data()
    if (is.null(gpl_dat)) {
      return()
    } else {
      gpl_fdat <- fData(gpl_dat)
      DT::datatable(
        data = gpl_fdat,
        options = list(columnDefs = list(list(
          targets = 1:ncol(gpl_fdat),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 20 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
            "}")
        ))),
        callback = JS('table.page(3).draw(false);'))
    }
  })
  
  # order genes by var
  var_order_data <- reactive({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      plot_data <- exprs(gpl_data)
      vars <- apply(plot_data, 1, var)
      vars <- vars[order(-vars)]
      plot_data <- plot_data[match(names(vars), rownames(plot_data)),] 
    }
  })
  
  # input for n of genes to plot
  output$max_genes <- renderUI({
    plot_data <- var_order_data()
    if (is.null(plot_data)) {
      tags$h5("Enter GSE ID, GPL ID and groupings info before attempting clustering.")
    } else {
      max <- nrow(na.omit(plot_data))
      output$temp_value <- renderText(paste0("Maximum: ", max, " (number of complete cases)"))
      numericInput('n_genes', 'Number of genes:', 100,
                   min = 1, max = max)
    }
  })
  
  # grouping levels
  get_cond_levels <- reactive({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else if (input$group_cname == "CUSTOM") {
      gpl_data <- pData(gpl_data)
      cond_levels <- sapply(1:nrow(gpl_data), function(i) {
        if(input[[paste0("group", i)]] == "") {
          return("NA")
        } else {
          return(input[[paste0("group", i)]])
        }
      })
    } else {
      cond_levels <- pData(gpl_data)[,input$group_cname]
    }
    cond_levels <- as.factor(cond_levels)
    cols_levels <- cond_levels
    levels(cols_levels) <- select_colours(nlevels(cols_levels))
    list(cond_levels = cond_levels, 
         cols_levels = cols_levels)
  })
  
  # output heatmap
  output$heatmap <- renderPlot({
    plot_data <- var_order_data()
    if (is.null(plot_data)) {
      return()
    } else {
      cond_levels <- get_cond_levels()[[1]]
      cols_levels <- get_cond_levels()[[2]]
      cols <- as.character(cols_levels)
      heatmap3(plot_data[1:input$n_genes,],
               cexRow = 0.5, scale = "row", labRow = NA, labCol = NA,
               ColSideColors = cols, ColSideLabs = NA, 
               na.rm = T, col = topo.colors(100))
    }
  })
  
  # legend for heatmap
  output$hm_leg <- renderPlot({
    cond_levels <- get_cond_levels()[[1]]
    if (is.null(cond_levels)) {
      return()
    } else {
      cols_levels <- get_cond_levels()[[2]]
      leg_func(legend = levels(cond_levels),
               col = levels(cols_levels))
    }
  })
  
  # perform pca
  get_pca_data <- reactive({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      pca_data <- t(na.omit(exprs(gpl_data)))
      pca_data <- as.data.frame(prcomp(pca_data, scale = T)$x)
      cond_levels <- get_cond_levels()[[1]] 
      pca_data$group <- as.character(cond_levels)
      return(pca_data)
    }
  })
  
  # input ui for choosing PCs
  output$choose_pc1 <- renderUI({
    pca_data <- get_pca_data()
    if (is.null(pca_data)) {
      tags$h5("Enter GSE ID, GPL ID and groupings info before attempting PCA.")
    } else {
      selectInput('pc1', 'First PC:', colnames(pca_data), "PC1")
    }
  })
  output$choose_pc2 <- renderUI({
    pca_data <- get_pca_data()
    if (is.null(pca_data)) {
      return()
    } else {
      selectInput('pc2', 'Second PC:', colnames(pca_data), "PC2")
    }
  })
  
  # plot PCA
  output$pca <- renderPlotly({
    pca_data <- get_pca_data()
    if (is.null(pca_data)) {
      return()
    } else {
      col_name <- "group"
      ggplot(data = pca_data, 
             aes_string(x = input$pc1, y = input$pc2, col = col_name)) +
        geom_point() +
        theme_bw(base_size = 18) +
        theme(legend.position = "bottom", 
              legend.title=element_blank(),
              panel.grid = element_blank()) +
        xlab(input$pc1) + ylab(input$pc2)
    }
  })
  
  # # hover info PCA
  # output$pca_hover <- renderPrint({
  #   d <- event_data("plotly_hover")
  #   if (is.null(d)) "Hover on a point!" else d
  # })
  
  # input ui for choosing contrast levels for GSEA
  output$choose_numer <- renderUI({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      tags$h5("Enter GSE ID, GPL ID and groupings info before attempting GSEA.")
    } else {
      cond_levels <- get_cond_levels()
      selectInput('numer', 'Numerator:', levels(cond_levels[[1]]))
    }
  })
  output$choose_denom <- renderUI({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      cond_levels <- get_cond_levels()[[1]]
      lvl2 <- ifelse(length(levels(cond_levels)) > 1,
                     levels(cond_levels)[2],
                     levels(cond_levels)[1])
      selectInput('denom', 'Denominator:', levels(cond_levels), lvl2)
    }
  })
  
  # file upload for custom
  output$gmt_upload <- renderUI({
    gs_type <- input$gs_type
    if (gs_type == "GMT") {
      fileInput("gmt_file", "Choose GMT file:",
                accept=c('text/plain'))
    } else {
      return()
    }
  })
  
  # gmt id types
  output$choose_id_type <- renderUI({
    gs_type <- input$gs_type
    if (gs_type == "GMT") {
      selectInput("gmt_id_type", "Choose ID type:",
                  keytypes(org.Hs.eg.db))
    } else {
      return()
    }
  })
  
  # match anno columns
  output$match_anno1 <- renderUI({
    gs_type <- input$gs_type
    if (gs_type == "GMT") {
      selectInput("gmt_id_type", "Choose ID type:",
                  ketypes(org.Hs.eg.db))
    } else {
      return()
    }
  })
  
  # annotation data
  anno <- reactive({
    gpl_data <- get_gpl_data()
    if (input$org == "Hs") {
      org_db <- "org.Hs.eg.db"
    } else if (input$org == "Mm") {
      org_db <- "org.Mm.eg.db"
    }
    bitr(fData(gpl_data)$'GENE SYMBOL', 
         fromType = "SYMBOL",
         toType = c("ENTREZID","ENSEMBL","REFSEQ",
                    "ACCNUM","UNIPROT"),
         OrgDb = org_db)
  })
    
  # getting gene set
  get_gene_sets <- reactive({
    gs_type <- input$gs_type
    if (gs_type == "GMT") {
      gmt_file <- input$gmt_file
      if (is.null(gmt_file)) return()
      gene_set <- read_gmt(gmt_file$datapath)
      return(gene_set)
    } else {
      gene_set <- eval(parse(text = paste0(input$org, ".", gs_type)))
      return(gene_set)
    }
  })
  
  # temp
  output$temp <- renderTable({
    gene_set <- get_gene_sets()
    if (is.null(gene_set)) return()
    summaryTable <- data.frame("Set name" = names(gene_set),
                               "N genes" = lengths(gene_set))
    summaryTable
  })
  
  # perform gsea test
  gsea_test <- reactive({
    gpl_data <- get_gpl_data()
    if (is.null(gpl_data)) {
      return()
    } else {
      cond_levels <- get_cond_levels()
      gene_sets <- get_gene_sets()
      gene_sets <- term_to_gene(gene_sets)
      numer_idx <- which(cond_levels$cond_levels == input$numer)
      denom_idx <- which(cond_levels$cond_levels == input$denom)
      exprs_mat <- exprs(gpl_data)
      rownames(exprs_mat) <- fData(gpl_data)$"Entrez_Gene_ID"
      gm_values <- gm_ratios(exprs_mat, numer_idx, denom_idx)
      gm_values <- gm_values[order(-gm_values)]
      gsea_res <- GSEA(gm_values, TERM2GENE = gene_sets,
                       minGSSize = 10, pvalueCutoff = 0.1)
      return(gsea_res)
    }
  })
  
  # plot PCA
  output$gsea_plot <- renderPlot({
    print("mooooooo")
    gsea_res <- gsea_test()
    print("booooooo")
    if (is.null(gsea_res)) {
      return()
    } else {
      if (is.null(chosen_gene_set)) {
        chosen_gene_set <- gsea_res$ID[1]
      }
      print("gooooooo")
      gseaplot(gsea_res, geneSetID = chosen_gene_set)
    }
  })
  
})
