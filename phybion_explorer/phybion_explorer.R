# setup and loading packages ----------------------------------------------

library("shiny")
library("ggplot2")
library("DESeq2")
library("mosdef")
library("dplyr")
library("GeneTonic")

dds_phybion <- readRDS("PhyBioN_SR_Analysis/dds.RDS")
dds_phybion <- estimateSizeFactors(dds_phybion)


dds_list <- list(
  dds_t1 = dds_phybion[, dds_phybion$Time_point == "2hrs"],
  dds_t2 = dds_phybion[, dds_phybion$Time_point == "6hrs"]
)

anno_df <- mosdef::get_annotation_orgdb(dds_phybion, "org.Hs.eg.db", id_type = "ENSEMBL")

# defining helper functions -----------------------------------------------

plot_gene_modified <- function (dds_list, gene, intgroup = "condition", assay = "counts",
                                annotation_obj = NULL, normalized = TRUE, transform = TRUE,
                                labels_display = FALSE, labels_repel = TRUE, plot_type = "auto",
                                return_data = FALSE, color_by = "Donors", plot_titles = NULL) {
  plot_type <- match.arg(plot_type, c("auto", "jitteronly", "boxplot", "violin", "sina"))

  # Ensure plot_titles matches the length of dds_list
  if (!is.null(plot_titles) && length(plot_titles) != length(dds_list)) {
    stop("The length of `plot_titles` must match the number of elements in `dds_list`.")
  }

  # Collect data from both dds objects
  df_list <- lapply(seq_along(dds_list), function(i) {
    dds <- dds_list[[i]]
    if (!(all(intgroup %in% colnames(colData(dds))))) {
      stop("`intgroup` not found in the colData slot of the dds object for dds ", i)
    }

    df <- get_expression_values(dds = dds, gene = gene, intgroup = intgroup,
                                assay = assay, normalized = normalized)
    df$sample_id <- rownames(df)

    # Use custom titles if provided
    df$source <- if (!is.null(plot_titles)) {
      plot_titles[i]
    } else {
      paste0("dds_", i)  # Default identifier for the dds object
    }

    # Add Donors column if it exists
    if ("Donors" %in% colnames(colData(dds))) {
      df$Donors <- colData(dds)[, "Donors"]
    } else {
      stop("The `Donors` column is missing in the colData of the dds object for dds ", i)
    }

    return(df)
  })
  df <- bind_rows(df_list)

  # Add gene annotation if provided
  if (!is.null(annotation_obj)) {
    genesymbol <- annotation_obj$gene_name[match(gene, annotation_obj$gene_id)]
  } else {
    genesymbol <- gene
  }

  onlyfactors <- df[, match(intgroup, colnames(df))]
  df$plotby <- interaction(onlyfactors)
  min_by_groups <- min(table(df$plotby))

  if (return_data) {
    return(df)
  }

  # Create ggplot
  p <- ggplot(df, aes(x = .data$plotby, y = .data$exp_value, col = .data[[color_by]])) +
    scale_x_discrete(name = "") +
    scale_color_discrete(name = "Experimental\ngroup") +
    theme_bw() +
    facet_wrap(~source, scales = "fixed")  # Use custom titles for facets

  jit_pos <- position_jitter(width = 0.2, height = 0, seed = 42)

  # Plot individual points for samples/donors
  p <- p + geom_point(position = jit_pos, shape = 16, size = 3)

  # Display labels if required
  if (labels_display) {
    if (labels_repel) {
      p <- p + ggrepel::geom_text_repel(aes(label = .data$sample_id),
                                        min.segment.length = 0, position = jit_pos)
    } else {
      p <- p + geom_text(aes(label = .data$sample_id), hjust = -0.1, vjust = 0.1, position = jit_pos)
    }
  }

  # Add a connected mean line across all groups
  p <- p + stat_summary(fun = mean, geom = "line", aes(group = 1),
                        color = "grey80", linewidth = 0.8)

  y_label <- if (assay == "counts" & normalized) {
    "Normalized counts"
  } else if (assay == "counts" & !normalized) {
    "Counts"
  } else if (assay == "abundance") {
    "TPM - Transcripts Per Million"
  } else {
    assay
  }

  if (transform) {
    p <- p + scale_y_log10(name = paste0(y_label, " (log10 scale)"))
  } else {
    p <- p + scale_y_continuous(name = y_label)
  }

  p <- p + labs(title = paste0(genesymbol, " - ", gene))

  return(p)
}


# here comes the app definition ------------------------------------------------

## ui definition ---------------------------------------------------------------

phybion_ui <- fluidPage(

  # Application title
  titlePanel("Exploring the Phybion dataset"),

  # Sidebar with input element(s)
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene_name", label = "select a gene", choices = NULL, selected = "FDXR")
    ),

    # Show a plot of the gene selected
    mainPanel(
      plotOutput("gene_plot")
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
      selected = "FDXR",
      server = TRUE
    )
  })

  output$gene_plot <- renderPlot({
    # check that something is selected, otherwise "throw a gentle error message"
    validate(
      need({input$gene_name != ""}, message = "Please select a gene")
    )

    # fetch the gene id "again"
    gene_id <- anno_df$gene_id[match(input$gene_name, anno_df$gene_name)]

    p <- plot_gene_modified(dds_list,
                            intgroup = "Dose",
                            gene = gene_id,
                            annotation_obj = anno_df,
                            plot_titles = c("2h", "6h"))
    p
  })
}

# run the application ----------------------------------------------------------
shinyApp(ui = phybion_ui, server = phybion_server)
