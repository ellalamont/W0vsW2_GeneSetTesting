library(shiny)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?


source("Import_data.R") # list_dfs and allGeneSetList

# Define UI ----
ui <- fluidPage(
  titlePanel("Forest Plots"),
  
  fluidRow(
    column(width = 3,
           
           # Dropdown for selecting with DEG file to use
           selectInput("my_DEG_file",
                       label = "DEG file",
                       choices = df_names),
           
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           
           # Adjust text size
           numericInput("my_text_size", 
                        label = "Text size", 
                        value = 1.1, min = 0.5, max = 2, step = 0.1),
           
           ),
    
    column(width = 9, # Max is 12...
           plotOutput("forest_plot", height = "650px"))
           # uiOutput("dynamic_forest_plot")) # To make the Forest plot change length based on number of gene sets
    )
  )


# Define server logic ----
server <- function(input, output, session){
  
  # To make the Forest plot change length based on number of gene sets
  # output$dynamic_forest_plot <- renderUI({
  #   n_GeneSets <- length(allGeneSetList[[input$my_GeneSetSource]])
  #   base_height <- 200 # minimum height (pixels)
  #   height_per_gene <- 1.5 # how many pixels per gene set
  #   plot_height <- max(base_height, n_GeneSets * height_per_gene)
  #   plotOutput("forest_plot", height = paste0(plot_height, "px"))
  # })
  
  # Make the Forest Plot
  output$forest_plot <- renderPlot({
    
    split_name <- strsplit(input$my_DEG_file, "\\.")[[1]] # Separate the df_name by the .
    Right_DataSet <- split_name[1]
    Left_DataSet  <- split_name[3]
    
    plotGeneSetForest(file = list_dfs[[input$my_DEG_file]],
                      geneSets = allGeneSetList[[input$my_GeneSetSource]],
                      main = input$my_DEG_file,
                      left.label = paste0(Left_DataSet, " (n=3)"), 
                      right.label = paste0(Right_DataSet, " (n=3)"),
                      xRange = 4, # Changes how far out the log2fold change axis goes
                      text.cex = input$my_text_size, pt.cex = 1.25, lwd = 3.5) 
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)


