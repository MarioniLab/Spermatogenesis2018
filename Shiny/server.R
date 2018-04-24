library(shiny)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scater)
load(paste("/Users/", 
                     system("whoami", intern=TRUE),
                     "/Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce.RData",
                     sep = ""))

shinyServer(function(input, output, session) {
  
  updateSelectizeInput(session, 'gene', 
                       choices = c(rowData(sce)$Symbol, unique(mouse.genes$Chromosome.scaffold.name)), 
                       server = TRUE)
  
  #### Visualize gene expression on tSNE
  createPlot.tsne <- eventReactive(input$goButton, {
    
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)]
    gene <- input$gene  # read in gene name
    
    if(grepl("Chr", gene)){
      cur_data <- logcounts(sce)[rowData(sce)$ID %in% 
                                       mouse.genes[mouse.genes$Chromosome.scaffold.name == gene,1],]
      Gene <- Matrix::colSums(cur_data)/Matrix::colSums(logcounts(sce))
    }
    else{
      Gene <- logcounts(sce)[rowData(sce)$Symbol == gene,]
    }
    
    ggplot(data.frame(tSNE1 = reducedDims(sce)$TSNE[,1],
                      tSNE2 = reducedDims(sce)$TSNE[,2],
                      Gene = Gene,
                      shown = ifelse(grepl(paste(input$dataset,collapse="|"), 
                                          colData(sce)$Sample),
                                     TRUE, FALSE))) +
      geom_point(aes(tSNE1, tSNE2, colour = Gene, alpha = shown)) + theme_minimal() + 
      scale_color_viridis() + scale_alpha_manual(values = c(0,1)) + 
      guides(alpha=FALSE)
    })
  
  #### Visualize sample info on tSNE
  createPlot.sample <- eventReactive(input$goButton, {
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)]
    
    ggplot(data.frame(tSNE1 = reducedDims(sce)$TSNE[,1],
                      tSNE2 = reducedDims(sce)$TSNE[,2],
                      sample = colData(sce)$Sample,
                      shown = ifelse(grepl(paste(input$dataset,collapse="|"), 
                                           colData(sce)$Sample),
                                     TRUE, FALSE))) +
      geom_point(aes(tSNE1, tSNE2, colour = sample, alpha = shown)) + theme_minimal() + 
      scale_color_brewer(palette = "Set1") + scale_alpha_manual(values = c(0,1)) + 
      guides(alpha=FALSE)
  })
  
  #### Visualize cluster info on tSNE
  createPlot.cluster <- eventReactive(input$goButton, {
    
    # Read in current dataset
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)]
    
    ggplot(data = data.frame(tSNE1 = reducedDims(sce)$TSNE[,1],
                             tSNE2 = reducedDims(sce)$TSNE[,2],
                             group = colData(sce)$Cluster,
                             shown = ifelse(grepl(paste(input$dataset,collapse="|"), 
                                                  colData(sce)$Sample),
                                            TRUE, FALSE))) +
      geom_point(aes(tSNE1, tSNE2, colour = group, alpha = shown)) +
      scale_color_manual(values = col_vector) + theme_minimal() + 
      scale_alpha_manual(values = c(0,1)) + 
      guides(alpha=FALSE)
    })
  
  #### Visualize gene expression across boxplots
  createBoxplot <- eventReactive(input$goButton, {
    
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)] 
    gene <- input$gene  # read in gene name
    
    if(grepl("Chr", gene)){
      cur_data <- logcounts(sce)[rowData(sce)$ID %in% 
                                       mouse.genes[mouse.genes$Chromosome.scaffold.name == gene,1],]
      Gene <- Matrix::colSums(cur_data)/Matrix::colSums(logcounts(sce))
      limits = c(0,max(Gene))
    }
    else{
      Gene <- logcounts(sce)[rowData(sce)$Symbol == gene,]
      limits = c(0,15)
    }
    
    ggplot(data.frame(value = Gene[grepl(paste(input$dataset,collapse="|"), 
                                         colData(sce)$Sample)],
                      cluster = factor(sce$Cluster[grepl(paste(input$dataset,collapse="|"), 
                                                         colData(sce)$Sample)],
                                     levels = names(col_vector)),
                      sample = factor(colData(sce)$Sample[grepl(paste(input$dataset,collapse="|"), 
                                                         colData(sce)$Sample)], 
                                      levels = input$dataset))) +
      geom_boxplot(aes(x = cluster, y = value, fill = cluster, 
                       group = interaction(cluster, sample))) + 
      scale_fill_manual(values = col_vector) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.background = element_blank()) +
      ylim(limits)
  })
  
  
  # Create the output of the event 
  output$tsne <- renderPlot({
    createPlot.tsne()
  })
  output$sample <- renderPlot({
    createPlot.sample()
  })
  output$cluster <- renderPlot({
    createPlot.cluster()
  })
  output$boxplot <- renderPlot({
    createBoxplot()
  })
  
})