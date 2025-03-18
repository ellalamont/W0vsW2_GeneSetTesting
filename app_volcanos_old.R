library(shiny)
library(gmodels)

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
  titlePanel("MICROM 431 Volcanos"),
  
  fluidRow(
    
    column(width = 2.5,
           # h3("Choose the condition"),
           selectInput("my_comparison",
                       label = "Volcano plots",
                       choices = df_names),
           textInput("my_GeneID", 
                     label = "Gene ID",
                     value = "Rv..."), 
    ),
    
    column(width = 5,
           plotlyOutput("volcano_plot",
                        width = 1000, height = 600),
    ),
  )
  
)

# Define server logic ----
server <- function(input, output) {
  
  # Volcano Plot
  output$volcano_plot <- renderPlotly({
    
    single_gene <- list_dfs_2[[input$my_comparison]] %>% 
      filter(GENE_ID == input$my_GeneID)
    
    my_volcano <- list_dfs_2[[input$my_comparison]] %>%
      ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_ID)) + 
      geom_point() + 
      
      # Add a differently colored point
      geom_point(data = single_gene, color = "yellow", aes(col = DE, label = DE_labels, text = GENE_ID)) + 
      
      labs(title = input$my_comparison) + 
      geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
      geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
      scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
    # geom_label_repel(max.overlaps = 10) # Can do geom_text_repel or geom_label_rebel
    
    # Determine the max and min axes values for labeling 
    plot_build <- ggplot_build(my_volcano)
    y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
    x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
    x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
    
    # Add the gene number annotations
    text_up <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant up") %>% nrow()
    text_down <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant down") %>% nrow()
    my_volcano_annotated <- my_volcano +
      annotate("text", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
      annotate("text", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
    
    final_plot <- my_volcano_annotated + my_plot_themes 
    final_plot
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)