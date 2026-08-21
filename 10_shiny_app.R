library(shiny)
library(plotly)
library(dplyr)
library(DT)

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

scatter_df <- readRDS(
  here::here(
    "results",
    "07-additional-vis",
    "scatter_df.rds"
  )
)



# ---------------------------------------------------------
# UI
# ---------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("DESeq2 Exposure × Time:Exposure Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Significance"),
      
      checkboxGroupInput(
        inputId = "sig_groups",
        label = "Show:",
        choices = c(
          "Neither" = "Neither",
          "Exposure" = "Exposure",
          "Time:Exposure" = "Time:Exposure",
          "Both" = "Both"
        ),
        selected = c(
          "Neither",
          "Exposure",
          "Time:Exposure",
          "Both"
        )
      ),
      
      hr(),
      
      h4("Significance threshold"),
      
      numericInput(
        inputId = "padj_cutoff",
        label = "Adjusted p-value cutoff:",
        value = 0.05,
        min = 0.0001,
        max = 1,
        step = 0.01
      ),
      
      hr(),
      
      h4("Effect size"),
      
      numericInput(
        inputId = "lfc_cutoff",
        label = "Minimum absolute log2FC:",
        value = 0,
        min = 0,
        step = 0.1
      ),
      
      hr(),
      
      h4("Table columns"),
      
      strong("Functional annotations"),
      
      selectizeInput(
        inputId = "annotation_cols",
        label = NULL,
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Select annotation columns..."
        )
      ),
      
      fluidRow(
        
        column(
          6,
          actionButton(
            inputId = "select_all_annotations",
            label = "Select all",
            width = "100%"
          )
        ),
        
        column(
          6,
          actionButton(
            inputId = "clear_annotations",
            label = "Clear all",
            width = "100%"
          )
        )
        
      ),
      
      br(),
      
      strong("Genome information"),
      
      selectizeInput(
        inputId = "genome_cols",
        label = NULL,
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Select genome columns..."
        )
      ),
      
      fluidRow(
        
        column(
          6,
          actionButton(
            inputId = "select_all_genome",
            label = "Select all",
            width = "100%"
          )
        ),
        
        column(
          6,
          actionButton(
            inputId = "clear_genome",
            label = "Clear all",
            width = "100%"
          )
        )
      ),
      
      
      hr(),
      
      actionButton(
        inputId = "reset_selection",
        label = "Reset brush selection"
      ),
      
      width = 3
    ),
    
    mainPanel(
      
      h4("Exposure vs Time:Exposure"),
      
      plotlyOutput(
        outputId = "scatter",
        height = "700px"
      ),
      
      hr(),
      
      h4("Selected genes"),
      
      DTOutput(
        outputId = "gene_table"
      ),
      
      hr(),
      
      h4("Gene details"),
      
      uiOutput(
        outputId = "gene_details"
      )
      
    )
    
  )
)


# ---------------------------------------------------------
# Server
# ---------------------------------------------------------

server <- function(input, output, session) {
  
    # Columns that aren't annotations
    
    core_cols <- c(
      "gene_id",
      "exposure_log2FC",
      "exposure_padj",
      "interaction_log2FC",
      "interaction_padj",
      "significance"
    )
    
    genome_cols <- c(
      "locus_tag",
      "operon",
      "Start",
      "Stop",
      "Strand",
      "Type"
    )
    
    annotation_cols <- c(
      "Gene",
      "Product",
      "COG",
      "COG_def",
      "COG_cat",
      "COG_cat_def",
      "KO",
      "KO_def",
      "KEGG_Pathway",
      "KEGG_path_def",
      "KEGG_Module",
      "KEGG_module_def",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "PFAMs",
      "GO_bak",
      "EC_bak",
      "Other_bak",
      "EC_egg",
      "RefSeq_bak",
      "UniRef100_bak",
      "UniRef90_bak",
      "UniRef50_bak"
    )
    
    # Only offer columns that actually exist
    annotation_cols <- intersect(
      annotation_cols,
      names(scatter_df)
    )
    
    genome_cols <- intersect(
      genome_cols,
      names(scatter_df)
    )
    
    observe({
      
      updateSelectizeInput(
        session,
        "annotation_cols",
        choices = annotation_cols,
        selected = intersect(
          c(
            "Gene",
            "Product",
            "COG",
            "COG_def",
            "COG_cat"
          ),
          annotation_cols
        ),
        server = TRUE
      )
      
      updateSelectizeInput(
        session,
        "genome_cols",
        choices = genome_cols,
        selected = intersect(
          c(
            "locus_tag",
            "operon",
            "Start",
            "Stop",
            "Strand",
            "Type"
          ),
          genome_cols
        ),
        server = TRUE
      )
      
    })
    
   

  
  
  # -------------------------------------------------------
  # Apply significance and effect-size filters
  # -------------------------------------------------------
  
  filtered_data <- reactive({
    
    df <- scatter_df
    
    # Recalculate significance using current cutoff
    df <- df %>%
      mutate(
        significance = case_when(
          
          exposure_padj < input$padj_cutoff &
            interaction_padj < input$padj_cutoff ~ "Both",
          
          exposure_padj < input$padj_cutoff ~ "Exposure",
          
          interaction_padj < input$padj_cutoff ~ "Time:Exposure",
          
          TRUE ~ "Neither"
        )
      )
    
    # Significance group filter
    df <- df %>%
      filter(significance %in% input$sig_groups)
    
    # Effect-size filter
    df <- df %>%
      filter(
        abs(exposure_log2FC) >= input$lfc_cutoff |
          abs(interaction_log2FC) >= input$lfc_cutoff
      )
    
    df
  })
  
  
  # -------------------------------------------------------
  # Hover text
  # -------------------------------------------------------
  
  plot_data <- reactive({
    
    filtered_data() %>%
      mutate(
        hover_text = paste0(
          "<b>", gene_id, "</b>",
          
          "<br>Gene: ",
          ifelse(is.na(Gene), "NA", Gene),
          
          "<br>Product: ",
          ifelse(is.na(Product), "NA", Product),
          
          "<br><br><b>Exposure</b>",
          
          "<br>log2FC: ",
          round(exposure_log2FC, 3),
          
          "<br>padj: ",
          signif(exposure_padj, 3),
          
          "<br><br><b>Time:Exposure</b>",
          
          "<br>log2FC: ",
          round(interaction_log2FC, 3),
          
          "<br>padj: ",
          signif(interaction_padj, 3),
          
          "<br><br>COG: ",
          ifelse(is.na(COG), "NA", COG),
          
          "<br>COG description: ",
          ifelse(is.na(COG_def), "NA", COG_def),
          
          "<br>KO: ",
          ifelse(is.na(KO), "NA", KO)
        )
      )
  })
  
  
  # -------------------------------------------------------
  # Plot
  # -------------------------------------------------------
  
  output$scatter <- renderPlotly({
    
    df <- plot_data()
    
    # Background genes
    p <- ggplot() +
      
      geom_point(
        data = df %>%
          filter(significance == "Neither"),
        
        aes(
          x = exposure_log2FC,
          y = interaction_log2FC,
          text = hover_text,
          key = gene_id
        ),
        
        colour = "grey85",
        size = 1.2,
        alpha = 0.3
      ) +
      
      # Significant genes
      geom_point(
        data = df %>%
          filter(significance != "Neither"),
        
        aes(
          x = exposure_log2FC,
          y = interaction_log2FC,
          colour = significance,
          text = hover_text,
          key = gene_id
        ),
        
        size = 1.8,
        alpha = 0.8
      ) +
      
      scale_colour_manual(
        values = c(
          "Exposure" = "#edae49",
          "Time:Exposure" = "#d1495b",
          "Both" = "#00798c"
        )
      ) +
      
      geom_hline(
        yintercept = 0,
        colour = "grey60",
        linewidth = 0.4
      ) +
      
      geom_vline(
        xintercept = 0,
        colour = "grey60",
        linewidth = 0.4
      ) +
      
      theme_bw(base_size = 11) +
      
      labs(
        x = "Exposure log2 Fold Change",
        y = "Time:Exposure log2 Fold Change",
        colour = "Significance"
      )
    
    ggplotly(
      p,
      tooltip = "text",
      source = "scatter"
    ) %>%
      highlight(
        on = "plotly_click",
        off = "plotly_doubleclick",
        persistent = TRUE,
        dynamic = TRUE,
        color = "black",
        opacityDim = 0.25,
        size = 8
      )%>%
      layout(
        dragmode = "select"
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToAdd = c(
          "select2d",
          "lasso2d"
        )
      )
  })
  
  
  # -------------------------------------------------------
  # Selected genes from brush
  # -------------------------------------------------------
  
  selected_data <- reactive({
    
    event_data(
      event = "plotly_selected",
      source = "scatter"
    )
  })
  
  
  # -------------------------------------------------------
  # Gene table
  # -------------------------------------------------------
  
  output$gene_table <- renderDT({
    
    df <- filtered_data()
    
    selection <- selected_data()
    
    # Restrict to brushed genes
    if (!is.null(selection) && nrow(selection) > 0) {
      
      df <- df %>%
        filter(gene_id %in% selection$key)
      
    }
    
    # Core DESeq2 columns
    core_cols <- c(
      "gene_id",
      "exposure_log2FC",
      "exposure_padj",
      "interaction_log2FC",
      "interaction_padj",
      "significance"
    )
    
    # Add selected annotation columns
    table_cols <- c(
      core_cols,
      input$annotation_cols
    )
    
    # Only retain columns that actually exist
    table_cols <- intersect(
      table_cols,
      names(df)
    )
    
    df <- df %>%
      select(all_of(table_cols))
    
    datatable(
      df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "500px",
        autoWidth = TRUE
      ),
      selection = "single",
      rownames = FALSE
    )
    
  }, server = TRUE)
  
  
  
  #
  selected_gene <- reactive({
    
    row <- input$gene_table_rows_selected
    
    if (length(row) == 0) {
      return(NULL)
    }
    
    df <- filtered_data()
    
    selection <- selected_data()
    
    if (!is.null(selection) && nrow(selection) > 0) {
      
      df <- df %>%
        filter(gene_id %in% selection$key)
      
    }
    
    df[row, , drop = FALSE]
    
  })
  ####
  
  output$gene_details <- renderUI({
    
    gene <- selected_gene()
    
    if (is.null(gene)) {
      
      return(
        p(
          "Brush or lasso genes in the plot, then click a gene in the table.",
          style = "color: grey;"
        )
      )
      
    }
    
    tagList(
      
      h4(gene$gene_id),
      
      h5("DESeq2 results"),
      
      fluidRow(
        
        column(
          3,
          strong("Exposure log2FC"),
          br(),
          round(gene$exposure_log2FC, 3)
        ),
        
        column(
          3,
          strong("Exposure padj"),
          br(),
          signif(gene$exposure_padj, 4)
        ),
        
        column(
          3,
          strong("Time:Exposure log2FC"),
          br(),
          round(gene$interaction_log2FC, 3)
        ),
        
        column(
          3,
          strong("Time:Exposure padj"),
          br(),
          signif(gene$interaction_padj, 4)
        )
        
      ),
      
      hr(),
      
      h5("Genome annotation"),
      
      fluidRow(
        
        column(
          4,
          strong("Gene"),
          br(),
          gene$Gene
        ),
        
        column(
          8,
          strong("Product"),
          br(),
          gene$Product
        )
        
      ),
      
      br(),
      
      fluidRow(
        
        column(
          4,
          strong("COG"),
          br(),
          gene$COG
        ),
        
        column(
          8,
          strong("COG description"),
          br(),
          gene$COG_def
        )
        
      ),
      
      br(),
      
      fluidRow(
        
        column(
          4,
          strong("KO"),
          br(),
          gene$KO
        ),
        
        column(
          8,
          strong("KO description"),
          br(),
          gene$KO_def
        )
        
      ),
      
      br(),
      
      fluidRow(
        
        column(
          6,
          strong("KEGG pathway"),
          br(),
          gene$KEGG_Pathway
        ),
        
        column(
          6,
          strong("KEGG module"),
          br(),
          gene$KEGG_Module
        )
        
      ),
      
      hr(),
      
      h5("Genomic location"),
      
      fluidRow(
        
        column(
          3,
          strong("Start"),
          br(),
          gene$Start
        ),
        
        column(
          3,
          strong("Stop"),
          br(),
          gene$Stop
        ),
        
        column(
          3,
          strong("Strand"),
          br(),
          gene$Strand
        ),
        
        column(
          3,
          strong("Type"),
          br(),
          gene$Type
        )
        
      )
      
    )
    
  })
  
  
  
  # -------------------------------------------------------
  # Reset selection
  # -------------------------------------------------------
  
  observeEvent(input$reset_selection, {
    
    plotlyProxy(
      "scatter",
      session
    ) %>%
      plotlyProxyInvoke(
        "relayout",
        list(
          selections = list()
        )
      )
    
  })
  
}


# ---------------------------------------------------------
# Run app
# ---------------------------------------------------------

shinyApp(
  ui = ui,
  server = server
)
