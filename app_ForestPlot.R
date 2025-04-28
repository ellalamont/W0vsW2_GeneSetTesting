library(shiny)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?


source("Import_data.R") # list_dfs and allGeneSetList

# Define UI ----
ui <- fluidPage(
  titlePanel("Forest Plots"),
  
  fluidRow(
    column(width = 4,
           
           # Dropdown for selecting with DEG file to use
           selectInput("my_DEG_file",
                       label = "DEG file",
                       choices = df_names),
           
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList))
           ),
    
    column(width = 8, # Max is 12...
           plotOutput("forest_plot"))
    )
  )


# Define server logic ----
server <- function(input, output, session){
  
  # Make the Forest Plot
  output$forest_plot <- renderPlot({
    
    split_name <- strsplit(input$my_DEG_file, "\\.")[[1]] # Separate the df_name by the .
    Right_DataSet <- split_name[1]
    Left_DataSet  <- split_name[3]
    
    plotGeneSetForest(file = list_dfs[[input$my_DEG_file]],
                      geneSets =allGeneSetList[[input$my_GeneSetSource]],
                      main = input$my_DEG_file,
                      left.label = paste0(Left_DataSet, " (n=3)"), 
                      right.label = paste0(Right_DataSet, " (n=3)"),
                      xRange = 4, # Changes how far out the log2fold change axis goes
                      text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)


