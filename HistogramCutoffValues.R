library(shiny)
library(bslib)
library(mclust)
library(tidyverse)
set.seed(123)
expr <- readRDS("expression_values.rds")


# Function to calculate cutoff value for a gene and cancer type
calculate_cutoff <- function(data, GENE, CANCER) {
  filtered_data <- data %>% filter(name == GENE & `tcga code` == CANCER) %>% pull(value)
  
  fit <- Mclust(filtered_data, G = 2)
  means <- fit$parameters$mean
  v <- mean(means)
  
  hist_GMM <- ggplot(data %>% filter(name == GENE & `tcga code` == CANCER), aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), fill = "blue", color = "black", alpha = 0.5, bins = 30) +
    geom_density(color = "red", size = 1) +
    geom_vline(xintercept = v, color = "red", linetype = "dashed", linewidth = 1.5) +
    annotate("text", x = Inf, y = Inf, label = paste("Cutoff:", round(v, 2)), vjust = 1,hjust = 1, color = "black", size = 8) + 
    labs(y = "Frequency", x = "Expression (log2+1)") +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))  
  
  return(list(plot = hist_GMM, cutoff = v))
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Histogram: Cutoff Value"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("GENE", "Select Gene:", choices = sort(unique(expr$name))),
      selectInput("CANCER", "Select Cancer Type:", choices = sort(unique(expr$`tcga code`))),
      actionButton("update", "Update"),
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      plotOutput("histPlot"),
      textOutput("CutoffValue")
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  
  observeEvent(input$update, {
    result <- calculate_cutoff(expr, input$GENE, input$CANCER)
    
    output$histPlot <- renderPlot({
      result$plot
    })
    
    output$CutoffValue <- renderText({
      paste("Cutoff Value:", round(as.numeric(result$cutoff, 2)))
    })
    output$downloadPlot <- downloadHandler(
      filename = function(){
        paste("histogramcutoff_", input$CANCER, "_", input$GENE, ".pdf", sep = "")
      },content = function(file){
        ggsave(file, plot = result$plot, device = "pdf")
      }
    )
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

