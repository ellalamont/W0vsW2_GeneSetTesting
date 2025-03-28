library(shiny)
library(gmodels)
library(grid)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?


source("Import_data.R") # for list_dfs_2

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20))



# Define UI ----
ui <- fluidPage(
  titlePanel("Sputum Pheatmap"),
  
  fluidRow(
    
    column(width = 2.5,

           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           # Dropdown for selecting the gene set within the chosen rda file.
           selectInput("my_GeneSet",
                       label = "Gene Set",
                       choices = NULL),
           # Add checkbox to toggle gene set points
           checkboxInput("show_gene_types", label = "Show gene types", value = FALSE),
    ),
    
    column(width = 5,
           uiOutput("dynamic_pheatmap")
           # plotOutput("pheatmap", width = "200%", height = "600")
           # plotOutput("pheatmap", width = "100%", height = "600px")
    )
    
  )
  
)

# Define server logic ----
server <- function(input, output, session) {
  
  # When a new gene set source is selected, update the gene set dropdown
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(session, "my_GeneSet",
                      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
                      selected = NULL)
  })
  
  # Dynamic UI for heatmap height
  output$dynamic_pheatmap <- renderUI({
    req(input$my_GeneSetSource, input$my_GeneSet)
    
    # Count the number of genes in the selected set
    num_genes <- sum(rownames(my_tpm) %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Dynamically set plot height (base height + extra space per gene)
    plot_height <- max(400, min(2500, num_genes * 40))  # Adjust as needed
    
    plotOutput("pheatmap", height = paste0(plot_height, "px"))
  })
  
  
  # Render the pheatmap
  output$pheatmap <- renderPlot({
    # Ensure both inputs are available
    req(input$my_GeneSetSource, input$my_GeneSet)
    
    # Subset the dataframe using genes in the selected gene set
    # my_tpm has genes as rownames and samples as columns
    my_data <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Determine annotation_row based on checkbox input
    annotation_row_data <- if(input$show_gene_types) gene_annot["Product"] else NULL
    
  
    p <- pheatmap(my_data, 
                  annotation_col = my_pipeSummary["Week"], 
                  annotation_row = annotation_row_data,  # Conditional annotation
                  annotation_colors = my_annotation_colors,
                  scale = "row", 
                  fontsize = 18)
    p
    
    # pheatmap returns a complex grid object; use grid.draw() to render it in Shiny.
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))