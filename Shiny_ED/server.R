library(shiny)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scater)
load(paste("/Users/", 
                     system("whoami", intern=TRUE),
                     "/Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce_ED.RData",
                     sep = ""))
to.show <- c(0,1)
names(to.show) <- c("Excluded", "Included")

shinyServer(function(input, output, session) {
  
  updateSelectizeInput(session, 'gene', 
                       choices = c(rowData(sce)$Symbol, unique(mouse.genes$Chromosome.scaffold.name)), 
                       server = TRUE)
  updateSelectInput(session, 'group', 
                       choices = levels(colData(sce)$AnnotatedClusters), 
                       selected = levels(colData(sce)$AnnotatedClusters))
  
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
                      shown = factor(ifelse(grepl(paste(input$dataset, collapse="|"), 
                                          colData(sce)$Sample) & 
                                            colData(sce)$AnnotatedClusters %in% input$group,
                                          "Included", "Excluded"), levels = c("Excluded", "Included")))) +
      geom_point(aes(tSNE1, tSNE2, colour = Gene, alpha = shown)) + theme_minimal() + 
      scale_color_viridis() + scale_alpha_manual(values = to.show) + 
      guides(alpha=FALSE)
    })
  
  #### Visualize sample info on tSNE
  createPlot.sample <- eventReactive(input$goButton, {
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)]
    
    cur_color_vector <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Set3"))
    
    ggplot(data.frame(tSNE1 = reducedDims(sce)$TSNE[,1],
                      tSNE2 = reducedDims(sce)$TSNE[,2],
                      sample = colData(sce)$Sample,
                      shown = factor(ifelse(grepl(paste(input$dataset, collapse="|"), 
                                                  colData(sce)$Sample) & 
                                              colData(sce)$AnnotatedClusters %in% input$group,
                                            "Included", "Excluded"), levels = c("Excluded", "Included")))) +
      geom_point(aes(tSNE1, tSNE2, colour = sample, alpha = shown)) + theme_minimal() + 
      scale_color_manual(values = cur_color_vector) + scale_alpha_manual(values = to.show) + 
      guides(alpha=FALSE)
  })
  
  #### Visualize cluster info on tSNE
  createPlot.cluster <- eventReactive(input$goButton, {
    
    # Read in current dataset
    #cur_sce <- sce[,grepl(input$select, colData(sce)$Sample)]
    
    ggplot(data = data.frame(tSNE1 = reducedDims(sce)$TSNE[,1],
                             tSNE2 = reducedDims(sce)$TSNE[,2],
                             group = colData(sce)$AnnotatedClusters,
                             shown = factor(ifelse(grepl(paste(input$dataset, collapse="|"), 
                                                         colData(sce)$Sample) & 
                                                     colData(sce)$AnnotatedClusters %in% input$group,
                                                   "Included", "Excluded"), levels = c("Excluded", "Included")))) +
      geom_point(aes(tSNE1, tSNE2, colour = group, alpha = shown)) +
      scale_color_manual(values = color_vector) + theme_minimal() + 
      scale_alpha_manual(values = to.show) + 
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
    }
    else{
      Gene <- logcounts(sce)[rowData(sce)$Symbol == gene,]
    }
    
    ggplot(data.frame(value = Gene[grepl(paste(input$dataset,collapse="|"), 
                                         colData(sce)$Sample) & 
                                     colData(sce)$AnnotatedClusters %in% input$group],
                      cluster = factor(sce$AnnotatedClusters[grepl(paste(input$dataset,collapse="|"), 
                                                         colData(sce)$Sample) & 
                                            colData(sce)$AnnotatedClusters %in% input$group],
                                     levels = names(color_vector)),
                      sample = factor(colData(sce)$Sample[grepl(paste(input$dataset,collapse="|"), 
                                                         colData(sce)$Sample) & 
                                          colData(sce)$AnnotatedClusters %in% input$group]))) +
      geom_boxplot(aes(x = cluster, y = value, fill = cluster, 
                       group = interaction(cluster, sample))) + 
      scale_fill_manual(values = color_vector) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.background = element_blank()) 
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