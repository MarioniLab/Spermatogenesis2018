library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Stylize front-end
  tags$style(type='text/css', ".select-input { font-size: 12px; line-height: 12px;} .select-dropdown { font-size: 12px; line-height: 12px; }"),
  
  # Application title
  titlePanel("Tsne representation of spermatogenesis"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", label = "Select gene", choices = NULL, options =
                       list(placeholder = 'Select a gene name', maxItems = 1,
                            maxOptions = 10)),
      helpText("Please select the gene symbol you would like to visualize.",
               "Visualization of the expression of whole chromosomes is achieved when",
               "entering the chromosome name (e.g. ChrX)."),
      selectizeInput("dataset", label = "Select dataset",
                     choices = c("B6", "P5", "P10", "P15", "P20", "P25", "P30", "P35"),
                     options = list(placeholder = 'Select one or multiple datasets',
                                    maxItems = 10)),
      helpText("Please select one or multiple datasets for visualization.",
               "The cells were mapped to adult B6 samples as reference."),
      selectInput("group", label = "Select groups",
                  choices = NULL, multiple = TRUE),
      helpText("Please select the somatic or germ cell types for visualization.\n",
               "eP: early Pachytene spermatocytes, mP: mid Pachytene spermatocytes,",
               "lP: late Pachytene spermatocytes, D: Diplotene spermatocytes, M: Meiosis,",
               "S: Spermatid, PTM: Peritubular Myoid Cells, tMg: testicular Macrophages"),
      actionButton("goButton", "Go!")
    ),
    
    mainPanel(
      tabPanel("Plot",
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"), plotOutput("tsne"), plotOutput("sample"))
               ),
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"), plotOutput("cluster"), plotOutput("boxplot"))
               ))
      
    )
  )
))