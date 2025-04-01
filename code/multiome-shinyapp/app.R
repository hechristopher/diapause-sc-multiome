library(shiny)
library(bslib)

source("rna_plots.R")
library(Seurat)
library(SeuratData)
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "diapause multiome",
  sidebar = sidebar(
    selectInput(
      "assay",
      label = "choose an assay",
      choices = 
        list(
          "RNA",
          "ATAC",
          "gene activity",
          "motif activity"
        )
    ),
    textInput("feature", label = "Feature"),
    actionButton("plot", label = "Plot")
  ),
  textOutput("selected_feature"),
  plotOutput("feature_umap", width = "100%")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$selected_feature <- renderText({
    paste0("Plotting: [", input$feature, "] in assay ", input$assay)
  })
  output$feature_umap <- renderPlot({
    FeaturePlot(pbmc3k.final, input$feature)
  }, height = 300, width = 300)
}

# Run the application 
shinyApp(ui = ui, server = server)
