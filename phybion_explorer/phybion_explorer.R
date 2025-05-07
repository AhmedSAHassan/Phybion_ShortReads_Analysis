# setup and loading packages ----------------------------------------------

library("shiny")
library("ggplot2")
library("DESeq2")
library("mosdef")
library("dplyr")
library("GeneTonic")
library("shinythemes")
library("shinycssloaders")
library("plotly")
library("patchwork")

dds_phybion <- readRDS("PhyBioN_SR_Analysis/dds.RDS")
dds_phybion <- estimateSizeFactors(dds_phybion)


dds_list <- list(
  dds_t1 = dds_phybion[, dds_phybion$Time_point == "2hrs"],
  dds_t2 = dds_phybion[, dds_phybion$Time_point == "6hrs"]
)

anno_df <- mosdef::get_annotation_orgdb(dds_phybion, "org.Hs.eg.db", id_type = "ENSEMBL")

# defining helper functions -----------------------------------------------

plot_gene_modified <- function(dds_list, genes, intgroup = "condition", assay = "counts",
                               annotation_obj = NULL, normalized = TRUE, transform = TRUE,
                               color_by = "Donors", plot_titles = c("2hrs", "6hrs")) {

  if (!is.null(plot_titles) && length(plot_titles) != length(dds_list)) {
    stop("The length of `plot_titles` must match the number of elements in `dds_list`.")
  }

  df_list <- lapply(seq_along(dds_list), function(i) {
    dds <- dds_list[[i]]

    data_list <- lapply(genes, function(gene) {
      df <- mosdef::get_expr_values(de_container = dds, gene = gene, intgroup = intgroup,
                                    assay = assay, normalized = normalized)
      df$sample_id <- rownames(df)
      df$gene_id <- gene
      df$source <- plot_titles[i]

      if (color_by %in% colnames(colData(dds))) {
        df[[color_by]] <- colData(dds)[, color_by]
      } else {
        df[[color_by]] <- NA
      }

      return(df)
    })

    bind_rows(data_list)
  })

  df <- bind_rows(df_list)
  df$plotby <- interaction(df[[intgroup]])

  if (!is.null(annotation_obj)) {
    df$gene_name <- annotation_obj$gene_name[match(df$gene_id, annotation_obj$gene_id)]
  } else {
    df$gene_name <- df$gene_id
  }

  plots <- lapply(unique(df$gene_name), function(gene_name) {
    gene_df <- df[df$gene_name == gene_name, ]

    ggplot(gene_df, aes(x = .data$plotby, y = .data$exp_value, color = .data[[color_by]],
                        group = .data$source,
                        text = paste0("Sample: ", sample_id,
                                      "<br>Value: ", round(exp_value, 2),
                                      "<br>Time: ", source,
                                      "<br>Dose: ", .data$plotby))) +
      geom_point(position = position_jitter(width = 0.2), shape = 16, size = 3, alpha = 0.8) +
      stat_summary(fun = mean, geom = "line", aes(group = source), color = "grey70", linewidth = 1) +
      facet_wrap(~source) +
      theme_minimal(base_size = 14) +
      labs(title = gene_name, x = "", y = ifelse(transform, "Normalized counts (log10)", "Normalized counts"),
           color = color_by) +
      scale_color_brewer(palette = "Set2") +
      {
        if (transform) scale_y_log10() else scale_y_continuous()
      }
  })

  return(plots)
}


# here comes the app definition ------------------------------------------------

## ui definition ---------------------------------------------------------------

phybion_ui <- fluidPage(
  theme = shinytheme("flatly"),
  title = "phybion_explorer: Exploring Gene Dynamics in Response to X-ray",
  # Application title
  titlePanel(div(
    h2("📈 Exploring Gene Dynamics in Response to X-ray"),
    style = "margin-bottom: 20px"
  )),

  # Sidebar with input element(s)
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        selectizeInput("gene_name", label = "Select gene(s)", choices = NULL,
                       selected = c("FDXR", "MDM2", "TICAM1", "NFKBIZ", "GPN1", "PMAIP1"), multiple = TRUE),
        selectInput("color_by", label = "Color by", choices = c("Sex", "Donors", "Dose")),
        downloadButton("download_plot", "Download Plot"),
        uiOutput("genecard_links")
      )
    ),

    mainPanel(
      h4("Here you can visualize normalized expression across doses and timepoints."),
      withSpinner(uiOutput("dynamic_plot_ui"), type = 4)
    )
  )
)

## server definition -----------------------------------------------------------
phybion_server <- function(input, output, session) {

  # to speed up the part to populate the entries of the selectizeInput element
  observe({
    updateSelectizeInput(
      session = session,
      inputId = "gene_name",
      choices = c(Choose = "", anno_df$gene_name),
      selected = c("FDXR", "MDM2", "TICAM1", "NFKBIZ", "GPN1", "PMAIP1"),
      server = TRUE
    )
  })

  # Store ggplot objects
  gene_plot_data <- reactive({
    req(input$gene_name)

    gene_ids <- anno_df$gene_id[match(input$gene_name, anno_df$gene_name)]

    plot_gene_modified(
      dds_list = dds_list,
      genes = gene_ids,
      intgroup = "Dose",
      annotation_obj = anno_df,
      color_by = input$color_by
    )
  })


  gene_plots <- reactive({
    lapply(gene_plot_data(), ggplotly, tooltip = "text")
  })

  output$dynamic_plot_ui <- renderUI({
    req(input$gene_name)
    plot_ids <- seq_along(input$gene_name)

    plot_outputs <- lapply(plot_ids, function(i) {
      plotlyOutput(outputId = paste0("gene_plot_", i), height = "500px")
    })

    do.call(tagList, plot_outputs)
  })

  observe({
    req(input$gene_name)
    plots <- gene_plots()

    lapply(seq_along(plots), function(i) {
      output[[paste0("gene_plot_", i)]] <- renderPlotly({
        plots[[i]]
      })
    })
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(paste(input$gene_name, collapse = "_"), "_expression_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = wrap_plots(gene_plot_data(), ncol = 1),
             device = "png", width = 12, height = 4 * length(input$gene_name))
    }
  )

  output$genecard_links <- renderUI({
    req(input$gene_name)
    links <- lapply(input$gene_name, function(gene) {
      gene_link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene)
      tags$a(href = gene_link, target = "_blank", paste("GeneCards:", gene))
    })
    tagList(tags$h5("GeneCards Links:"), tags$ul(lapply(links, tags$li)))
  })
}

# run the application ----------------------------------------------------------
shinyApp(ui = phybion_ui, server = phybion_server)

