
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinysky)
library(ggplot2)
library(grid)
library(gridExtra)
library(Gviz)
library(GenomicRanges)
library(biomaRt)

## source plot functions
source("plots/demoPlot2.R")
source("plots/demoPlot.R")
source("plots/plotting.R")
source("page1.r")
source("page2.r")

#readProjectOverview <- function(){
#  con <- file("www/oncoscape.txt")
#  content <- readLines(con)
#  close(con)
#  content
#}

# Load the TCGA results
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")
#tcgaResults = results
# Load the CCLE results
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step2.rdata")
#ccleResults = results
# Save some memory
#rm(results)

# Get the chromosomal location for all genes 
#sortedGeneLoc = sortGenesByLocation(tcgaResults, ccleResults)

# Load the TCGA data
# Copy number
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna.rdata")
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna_ccle.rdata")
#rm(ccle.cna)
# Gene expression
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs.rdata")
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs_ccle.rdata")
#rm(ccle.exprs)
# DNA methylation
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth.rdata")
# Methylation annotation data
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation.rdata")
# Project Achilles data
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/PROJECT_ACHILLES/20130620/RDATA/achilles.rdata")


# Reformat results for plotting page 1
# TCGA oncogene (OG) scores
#tcgaResultsHeatmapOG = heatmapDataframe(tcgaResults, 
#		 	     	        scores=list(combined="og.score", Meth="og.methylation",
#			     	        CNA="og.cna", Mut="og.mutations",
#			                shRNA="og.achilles", Expr="og.exprs"))
# TCGA tumor suppressor (TS) scores
#tcgaResultsHeatmapTS = heatmapDataframe(tcgaResults, 
#			     	        scores=list(combined="ts.score", Meth="ts.methylation",
#			     	        CNA="ts.cna", Mut="ts.mutations",
#			                shRNA="ts.achilles", Expr="ts.exprs"))
# TCGA combined score
#tcgaResultsHeatmapCombined = heatmapDataframe(tcgaResults) 
# CCLE OG scores
#ccleResultsHeatmapOG = heatmapDataframe(ccleResults, 
#		 	     	        scores=list(combined="og.score", Meth="og.methylation",
#			     	        CNA="og.cna", Mut="og.mutations",
#			                shRNA="og.achilles", Expr="og.exprs"))
# CCLE OG scores
#ccleResultsHeatmapTS = heatmapDataframe(ccleResults, 
#			     	        scores=list(combined="ts.score", Meth="ts.methylation",
#			     	        CNA="ts.cna", Mut="ts.mutations",
#			                shRNA="ts.achilles", Expr="ts.exprs"))
# CCLE combined score
#ccleResultsHeatmapCombined = heatmapDataframe(ccleResults)

shinyServer(function(input, output, session) {
  
  ###################################################################################
  ## load starter rdata object for widgets and app
  ################################################################################### 
  starterWidgets <- load("www/starterWidgets.RData")
  genes <- apply(geness, 1, function(r) r)
  chrms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  ## this will return the genes from the selected pathway in comp4
  ## currently set to false 4 values
  demoGenesByPathway <- c(genes[1],genes[100],genes[1000],genes[2000])
  starterData <- load("www/starter.RData")  

  
#   showshinyalert(session,"welcomeToApp",styleclass="primary",
#     HTML("
#         <h1 style='color:#1C1C1C' >OncoScape</h1>
#         <span style='color:red'><b>beta version</b></span>
#         
#         <p style='color:#1C1C1C' > 
#         OncoScape is a package for gene prioritization in the R statistical programming environment. 
#         The analysis is run in a contrast fashion, i.e. always two groups of samples are compared with each other. 
#         Examples include:
#         1. tumors vs. normals
#         2. cell lines vs. normals
#         3. treatment responders vs resistant
#         4. samples with mutations in gene X vs wild type
#         Currently, analyses of five data types are implemented in OncoScape:
#         1. gene expression
#         2. DNA copy number
#         3. DNA methylation
#         4. mutation
#         5. shRNA knock-down data
#         Aberrations in each gene are called for each data type separately and scored as 0 (no aberration found) or 1 (aberration found). 
#         These scores are summed across data types to give the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) 
#         aberrations and calculates independent scores for both directions. It is possible to run the analysis on any combination of these data types.
#         </p>
#         
#         <button class='btn btn-primary'>
#           <font size='3'><b>Get started</b></font>
#         </button>
#         
#         </br>
#         
#         <hr>
#         <h3 >Publicatoin</h3>
#         <p>some text about about this header</p>
#         <h3 >License, terms of use, privacy </h3>
#         <p>some text about about this header</p>
#         <h3 >Github repository (?)</h3>
#         <p>some text about about this header</p>
#         <h3 >Contact</h3>
#         <p>some text about about this header</p>
#         "
#     )
#   )  
  ###################################################################################  
  


  ###################################################################################
  ## comp0: Project overview/about
  ###################################################################################
  ## html
  #output$projectDescriptionContent <- renderUI({
  #  HTML(paste(projectDescription, collapse=""))
  #})


  ###################################################################################
  ## comp1: Genes over cancer type
  ###################################################################################
  ## widgets
  ## cancer selector comp1
  output$cancerSelectorC1 <- renderUI({
    selectInput("cancerSelectorChoiceC1", "Select a Cancer type", apply(cancers, 1, function(r) r))
  })
  ## sample set selector comp1
  output$sampleSelectorC1 <- renderUI({
    radioButtons("sampleSelectorC1", label = "Select Sample Set",
                 choices = list("Tumors" = "tumors", "Cell lines" = "cell-lines", "Combined" = "combined"),
                 selected = "tumors")
  })
  ## score type selector comp1
  output$scoreSelectInputC1 <- renderUI({
    selectInput("selectScoreTypeC1", label = "Score type", 
                choices = list("Oncogene score" = "og.score", "Tumor suppressor score" = "ts.score",
                               "Combined score" = "combined.score"), selected = "og.score")
  })
  
  ## score cut off selector comp1
  output$scoreCutoffSelectorC1 <- renderUI({
    selectInput("scoreCutoff", label = "Score cut-off", 
              choices = list("2" = 2, "3" = 3, "4" = 4),selected = 2)
  })
  
  
  ## tables
  output$genesResTable <- renderDataTable({
    input$refreshPlot
    if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
      geneDataFrameResultSet(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
    }else{
      geneDataFrameResultSet(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
    }
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5))
  
  ## plots
  ## view 1
  output$distPlot <- renderPlot({
    input$refreshPlot
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Generating plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 30
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    
    if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
      comp1view1Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
    }else{
      comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
    }
  })
  ## view 2
  output$distPlot2 <- renderPlot({
    input$refreshPlot
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Generating plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 30
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    
    
    if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
      comp1view2Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
    }else{
      comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
    }
  })
  
  ## downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.tsv', sep="\t")
    },
    content <- function(filename){
      write.csv(geneDataFrameResultSet(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1), filename)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('plot-', Sys.Date(), '.jpeg', sep="")
    },
    content <- function(file){
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = 800,
                       res = 100, units = "px")
      }
      g <- arrangeGrob(comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1), 
                  comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1))
      ggsave(file, plot = g, device = device)
    }
  )
  
  
  ## *******************************************************************************
  ## verbatim debug text
#   output$cancer <- renderText({
#     input$refreshPlot
#     paste("cancer : ", isolate(input$cancerSelectorChoiceC1))
#   })
#   output$score <- renderText({
#     input$refreshPlot
#     paste("score : ", isolate(input$selectScoreTypeC1))
#   })
#   output$number <- renderText({
#     input$refreshPlot    
#     paste("number : ", isolate(input$numberOfGenes))
#   })
#   output$sample <- renderText({
#     input$refreshPlot    
#     paste("sample : ", isolate(input$sampleSelectorC1))
#   })
  
  ###################################################################################
  
  
  
  ###################################################################################
  ## comp2: single gene plots
  ###################################################################################
  ## widgets
  ## cancer selector comp2
  output$cancerSelectorC2 <- renderUI({
    selectInput("cancerSelectorChoiceC2", "Select a Cancer type", apply(cancers, 1, function(r) r))
  })
  ## gene selector comp2
  updateSelectizeInput(session, 'geneSelectorChoiceC2', choices = genes, server = TRUE, selected=genes[100])
  ## sample set selector comp2
  output$sampleSelectorC2 <- renderUI({
    radioButtons("sampleSelectorC2", label = "Select Sample Set",
                 choices = list("Tumors" = 1, "Cell-lines" = 2, "Tumors vs Cell lines" = 3),
                 selected = 1)
  })
  
  ## gene expression
  output$geneExpressionPlot <- renderPlot({
    input$refreshPlotC2
    if (length(isolate(input$cancerSelectorChoiceC2))>0 && length(isolate(input$geneSelectorChoiceC2))>0 && length(isolate(input$sampleSelectorC2))>0){
      getPage2Plots(isolate(input$cancerSelectorChoiceC2), isolate(input$geneSelectorChoiceC2), isolate(input$sampleSelectorC2))[["gene.expression"]]
    }else{
      getPage2Plots(input$cancerSelectorChoiceC2, input$geneSelectorChoiceC2, input$sampleSelectorC2)
    }
    #demoPlot()
  })
  ## copy number
  output$cnvPlot <- renderPlot({
    input$refreshPlotC2
    if (length(isolate(input$cancerSelectorChoiceC2))>0 && length(isolate(input$geneSelectorChoiceC2))>0 && length(isolate(input$sampleSelectorC2))>0){
      getPage2Plots(isolate(input$cancerSelectorChoiceC2), isolate(input$geneSelectorChoiceC2), isolate(input$sampleSelectorC2))[["acgh"]]
    }else{
      getPage2Plots(input$cancerSelectorChoiceC2, input$geneSelectorChoiceC2, input$sampleSelectorC2)
    }
    #demoPlot()
  })
  ## dna methylation
  output$dnaMethPlot <- renderPlot({
    input$refreshPlotC2
    if (length(isolate(input$cancerSelectorChoiceC2))>0 && length(isolate(input$geneSelectorChoiceC2))>0 && length(isolate(input$sampleSelectorC2))>0){
      getPage2Plots(isolate(input$cancerSelectorChoiceC2), isolate(input$geneSelectorChoiceC2), isolate(input$sampleSelectorC2))[["methylation"]]
    }else{
      getPage2Plots(input$cancerSelectorChoiceC2, input$geneSelectorChoiceC2, input$sampleSelectorC2)
    }
    #demoPlot()
  })
  ## achilles
  output$achillesPlot <- renderPlot({
    input$refreshPlotC2
    if (length(isolate(input$cancerSelectorChoiceC2))>0 && length(isolate(input$geneSelectorChoiceC2))>0 && length(isolate(input$sampleSelectorC2))>0){
      getPage2Plots(isolate(input$cancerSelectorChoiceC2), isolate(input$geneSelectorChoiceC2), isolate(input$sampleSelectorC2))[["achilles"]]
    }else{
      getPage2Plots(input$cancerSelectorChoiceC2, input$geneSelectorChoiceC2, input$sampleSelectorC2)
    }
    #demoPlot()
  })
  ###################################################################################
  
  
  
  ###################################################################################
  ## comp3: Genes over chromosomes
  ###################################################################################
  ## widgets
  ## cancer selector comp3
  output$cancerSelectorC3 <- renderUI({
    selectInput("cancerSelectorChoiceC3", "Select a Cancer type", apply(cancers, 1, function(r) r))
  })    
  ## score type selector comp3
  output$scoreSelectInputC3 <- renderUI({
    selectInput("selectScoreTypeC3", label = "Score type", 
                choices = list("Choice 1" = 1, "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)
  })
  ## chr selector comp3
  output$chrSelector <- renderUI({
    selectInput("selectChrType", "Select Chromosome", chrms)
  })
  ## sample type selector comp3
  output$sampleSelectorC3 <- renderUI({
    radioButtons("sampleSelectorC3", label = "Select Sample Set",
                 choices = list("Tumors" = 1, "Cell lines" = 2),
                 selected = 1)
  })
  
  ## view 1
  output$selectedScorePlot <- renderText({
    "Section Under Development"
  })
  ## view 2
  output$perctAffectedSamplesPlot <- renderText({
    "Section Under Development"
  })
  ###################################################################################
  

  
  ###################################################################################
  ## comp4: Pathways 
  ###################################################################################
  ## widgets
  ## pathway selector comp4
  output$pathwaySelector <- renderUI({
    selectInput("selectPathwayType", label = "Select Pathway", 
                choices = list("Choice 1" = 1, "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)
  })
  ## sample set selector comp4
  output$sampleSelectorC4 <- renderUI({
    radioButtons("sampleSelectorC4", label = "Select Sample Set",
                 choices = list("Tumors" = 1, "Cell-lines" = 2),
                 selected = 1)
  })
  ## genes selector comp4
  updateSelectizeInput(session, 'geneSelectorChoiceC4', choices = genes, selected = demoGenesByPathway, server = TRUE)
  
  
  ## view 1
  output$pathwayPlot <- renderText({
    "Section Under Development"
  })
  ###################################################################################
  
  
})
