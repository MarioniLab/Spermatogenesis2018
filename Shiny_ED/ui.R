library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Stylize front-end
  tags$style(type='text/css', ".select-input { font-size: 12px; line-height: 12px;} .select-dropdown { font-size: 12px; line-height: 12px; }"),
  
  # Application title
  titlePanel("Tsne representation of spermatognesis"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", label = "Select gene", choices = NULL, options = 
                       list(placeholder = 'Select a gene name', maxItems = 1, 
                            maxOptions = 10)),
      selectizeInput("dataset", label = "Select dataset", 
                  choices = c("B6", "P10", "P15", "P20", "P30", "P35", "Tc0", "Tc1"),
                  options = list(placeholder = 'Select one or multiple datasets',
                                 maxItems = 8)),
      selectInput("group", label = "Select groups", 
                     choices = NULL, multiple = TRUE),
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