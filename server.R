
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
library(plyr)
library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(pathview)

## source plot functions
source("plotting.R")
source("page1.r")
source("page2.r")
source("page3.r")
source("page4.r")

## functions
showProgress <- function() {
  # Create 0-row data frame which will be used to store data
  dat <- data.frame(x = numeric(0), y = numeric(0))
  withProgress(message = 'Working ', value = 0, {
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
}


# Get the chromosomal location for all genes
#sortedGeneLoc = sortGenesByLocation(tcgaResults, ccleResults)

# Load the TCGA results
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")
#tcgaResults = results
# Load the CCLE results
#load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step2.rdata")
#ccleResults = results
# Save some memory
#rm(results)


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
#		 	     	        scores=list(Combined="og.score", Meth="og.methylation",
#			     	        CNA="og.cna", Mut="og.mutations",
#			                shRNA="og.achilles", Expr="og.exprs"))
# TCGA tumor suppressor (TS) scores
#tcgaResultsHeatmapTS = heatmapDataframe(tcgaResults, 
#			     	        scores=list(Combined="ts.score", Meth="ts.methylation",
#			     	        CNA="ts.cna", Mut="ts.mutations",
#			                shRNA="ts.achilles", Expr="ts.exprs"))
# TCGA combined score
#tcgaResultsHeatmapCombined = heatmapDataframe(tcgaResults) 
# CCLE OG scores
#ccleResultsHeatmapOG = heatmapDataframe(ccleResults, 
#		 	     	        scores=list(Combined="og.score", Meth="og.methylation",
#			     	        CNA="og.cna", Mut="og.mutations",
#			                shRNA="og.achilles", Expr="og.exprs"))
# CCLE OG scores
#ccleResultsHeatmapTS = heatmapDataframe(ccleResults, 
#			     	        scores=list(Combined="ts.score", Meth="ts.methylation",
#			     	        CNA="ts.cna", Mut="ts.mutations",
#			                shRNA="ts.achilles", Expr="ts.exprs"))
# CCLE combined score
#ccleResultsHeatmapCombined = heatmapDataframe(ccleResults)

shinyServer(function(input, output, session) {
  
  ###################################################################################
  ## load starter rdata object for widgets and app
  ################################################################################### 
  load("starterWidgets.RData")
  #load("starter.RData")
  #load("gloc.RData")
  genes <- apply(geness, 1, function(r) r)
  chrms <- c("1","2","3","4","5","6","7","8","9","10","11","12",
             "13","14","15","16","17","18","19","20","21","22","23","24")
  ## this will return the genes from the selected pathway in comp4
  ## currently set to false 4 values
  demoGenesByPathway <- c(genes[1],genes[100],genes[1000],genes[2000])
  ## read pathways file
  pathtemp <- read.delim("kegg_pathways.tsv",sep="\t",header=FALSE)
  pathways <- paste(pathtemp$V2,':',pathtemp$V3, sep= "")


  ###################################################################################
  ## comp0: Project overview/about
  ###################################################################################
  
  ###################################################################################
  
  

  ###################################################################################
  ## comp1: Genes over cancer type
  ###################################################################################
  ## widgets
  ## cancer selector comp1
  output$cancerSelectorC1 <- renderUI({
    selectInput("cancerSelectorChoiceC1", label = NULL, choices = apply(cancers, 1, function(r) r))
  })
  ## sample set selector comp1
  output$sampleSelectorC1 <- renderUI({
    radioButtons("sampleSelectorC1", label = NULL,
                 choices = list("Tumors" = "tumors", "Cell lines" = "cell-lines"),
                 selected = "tumors")
  })
  ## gene selection criteria
  output$geneSelectionMethodC1 <- renderUI({
    selectInput("geneSelectionMethodC1Value", label = "Choose a method below for selecting genes", 
                choices = c("Use Score Type and Cutoff" = "type1","Upload Your Gene list" = "type2",
                                                                        "Copy Paste Your Gene list" = "type3"))
  })
  
  ## gene selection procedure
  output$geneSelectionPanelC1 <- renderUI({
    if(is.null(input$geneSelectionMethodC1Value))
    {
      return()
    }
    switch(input$geneSelectionMethodC1Value,
           "type1" = c(selectInput("selectScoreTypeC1", label = "Select Score type", 
                                 choices = list("Oncogene score" = "og.score", "Tumor suppressor score" = "ts.score",
                                                "Combined score" = "combined.score"), selected = "og.score"),
                              selectInput("scoreCutoff", label = "Score cut-off", 
                                                             choices = list(">2" = 2, ">3" = 3, ">4" = 4, "<-2" = -2, "<-3" = -3, "<-4" = -4),
                                                             selected = 3)
                              ),
           "type2" = fileInput('geneListUploadC1', 'Upload Gene List File',accept = c(".tsv")),
           "type3" = shiny::tags$textarea(id="geneListValuesC1", rows=10, cols=10, "Copy Paste your genes here separated by comma")
           )
  })
  
#   ## score cut off selector comp1
#   if(input$geneSelectionMethodC1Value == "type1")
#   {
#     output$scoreCutoffSelectorC1 <- renderUI({
#       input$selectScoreTypeC1
#       if (input$selectScoreTypeC1 == 'combined.score'){  
#         selectInput("scoreCutoff", label = "Score cut-off", 
#                     choices = list(">2" = 2, ">3" = 3, ">4" = 4, "<-2" = -2, "<-3" = -3, "<-4" = -4),selected = 3)
#       }else{
#         selectInput("scoreCutoff", label = "Score cut-off", 
#                     choices = list(">2" = 2, ">3" = 3, ">4" = 4),selected = 3)      
#       }
#       
#     })    
#   }
  
  
  ## tables
  output$genesResTable <- renderDataTable({

    input$refreshPlot
    input$input$geneSelectionMethodC1Value
    
    if (input$geneSelectionMethodC1Value == "type1")
    {
      if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
        geneDataFrameResultSet(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
      }else{  
        geneDataFrameResultSet(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
      }        
    }else if (input$geneSelectionMethodC1Value == "type2"){
      userfile <- input$geneListUploadC1
      if (!(is.null(userfile)))
      {
        userdata <- read.delim(userfile$datapath,sep="\t")
        
        if (length(isolate(input$cancerSelectorChoiceC1))>0  && length(isolate(input$sampleSelectorC1))>0){
          geneFileDataFrameResultSet(isolate(input$cancerSelectorChoiceC1),userdata,isolate(input$sampleSelectorC1))
        }else{
          geneFileDataFrameResultSet(input$cancerSelectorChoiceC1,userdata,input$sampleSelectorC1)
        }                
      }else{
        return()
      }
    }else{
      genelist <- as.data.frame(strsplit(input$geneListValuesC1,',')[[1]])
      colnames(genelist) <-c ("uploadedGenes")
      if (nrow(genelist)>0 && genelist[1,1] != "Copy Paste your genes here separated by comma"){
        if (length(isolate(input$cancerSelectorChoiceC1))>0  && length(isolate(input$sampleSelectorC1))>0){
          geneFileDataFrameResultSet(isolate(input$cancerSelectorChoiceC1),genelist,isolate(input$sampleSelectorC1))
        }else{
          geneFileDataFrameResultSet(input$cancerSelectorChoiceC1,genelist,input$sampleSelectorC1)
        }                        
      }else{
        return()
      }
    }
    
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  
  ## plots
  ## view 1
  output$distPlot2 <- renderPlot({
    input$refreshPlot
    input$input$geneSelectionMethodC1Value
    
    if (input$geneSelectionMethodC1Value == "type1")
    {
      if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
        showProgress()
        comp1view2Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
      }else{
        showProgress()
        comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
      }
    }else if(input$geneSelectionMethodC1Value == "type2")
    {
      userfile <- input$geneListUploadC1
      if (!(is.null(userfile)))
      {
        userdata <- read.delim(userfile$datapath,sep="\t")  
        if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
          showProgress()
          comp1view2Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
        }else{
          showProgress()
          comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
        }      
      }else{
        return()
      }
    }else{
      if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
        showProgress()
        comp1view2Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
      }else{
        showProgress()
        comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
      }
    }
    
  })  
  ## view 2
  output$distPlot <- renderPlot({
    input$refreshPlot
    input$input$geneSelectionMethodC1Value
    
    if (input$geneSelectionMethodC1Value == "type1")
    {
      if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
        showProgress()
        comp1view1Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
      }else{
        showProgress()
        comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
      }  
    }else if(input$geneSelectionMethodC1Value == "type2")
    {
      userfile <- input$geneListUploadC1
      if (!(is.null(userfile)))
      {
        userdata <- read.delim(userfile$datapath,sep="\t")  
        if (length(isolate(input$cancerSelectorChoiceC1))>0  && length(isolate(input$sampleSelectorC1))>0){
          comp1view1FilePlot(isolate(input$cancerSelectorChoiceC1),userdata,isolate(input$sampleSelectorC1))
        }else{
          comp1view1FilePlot(input$cancerSelectorChoiceC1,userdata,input$sampleSelectorC1)
        }  
      }else{
        return()
      }
      
    }else{
      if (length(isolate(input$scoreCutoff))>0 && length(isolate(input$cancerSelectorChoiceC1))>0 && length(isolate(input$selectScoreTypeC1))>0 && length(isolate(input$sampleSelectorC1))>0){
        showProgress()
        comp1view1Plot(isolate(input$scoreCutoff),isolate(input$cancerSelectorChoiceC1),isolate(input$selectScoreTypeC1),isolate(input$sampleSelectorC1))
      }else{
        showProgress()
        comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1)
      }  
    }
    
    
  })
  
  ## downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.tsv', sep="\t")
    },
    content <- function(filename){
      
      if (input$geneSelectionMethodC1Value == "type1")
      {
        write.csv(geneDataFrameResultSet(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1), filename) 
        
      }else if (input$geneSelectionMethodC1Value == "type2")
      {
        userfile <- input$geneListUploadC1
        if (!(is.null(userfile)))
        {
          userdata <- read.delim(userfile$datapath,sep="\t")
          write.csv(geneFileDataFrameResultSet(input$cancerSelectorChoiceC1,userdata,input$sampleSelectorC1),filename)          
        }
      }else{
        genelist <- as.data.frame(strsplit(input$geneListValuesC1,',')[[1]])
        colnames(genelist) <-c ("uploadedGenes")
        if (nrow(genelist)>0 && genelist[1,1] != "Copy Paste your genes here separated by comma"){
          write.csv(geneFileDataFrameResultSet(input$cancerSelectorChoiceC1,genelist,input$sampleSelectorC1), filename)                        
        }
      }
      
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
      
      if (input$geneSelectionMethodC1Value == "type1")
      {
        g <- arrangeGrob(comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1), 
                         comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1))
        ggsave(file, plot = g, device = device)
      }else if (input$geneSelectionMethodC1Value == "type2")
      {
        userfile <- input$geneListUploadC1
        if (!(is.null(userfile)))
        {
          userdata <- read.delim(userfile$datapath,sep="\t")
          g <- arrangeGrob(comp1view2FilePlot(input$cancerSelectorChoiceC1,userdata,input$sampleSelectorC1),
                           comp1view1FilePlot(input$cancerSelectorChoiceC1,userdata,input$sampleSelectorC1))
          ggsave(file,plot = g, device = device)
        }
      }else{
        g <- arrangeGrob(comp1view2Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1), 
                         comp1view1Plot(input$scoreCutoff,input$cancerSelectorChoiceC1,input$selectScoreTypeC1,input$sampleSelectorC1))
        ggsave(file, plot = g, device = device)
      }  
      
    }
  )
  
  
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
    radioButtons("sampleSelectorC2", label = "Select Sample type",
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
                choices = list("Oncogene Score" = "OG", "Tumor Suppressor Score" = "TS",
                                 "Combined score" = "CO"), selected = "TS")
  })
  ## chr selector comp3
  output$chrSelector <- renderUI({
    selectInput("selectChrType", "Select Chromosome", chrms)
  })
  
  ## sample type selector comp3
  #output$sampleSelectorC3 <- renderUI({
  #  radioButtons("sampleSelectorC3", label = "Select Sample Set",
  #               choices = list("Tumors" = 1, "Cell lines" = 2),
  #               selected = 1)
  #})
  
  ## plots
  ## view 1
  output$selectedScorePlot <- renderPlot({
    input$refreshPlotC3
    if (length(isolate(input$cancerSelectorChoiceC3))>0 && length(isolate(input$selectScoreTypeC3))>0 && length(isolate(input$selectChrType))>0 ){
      showProgress()
      comp3view1Plot(isolate(input$cancerSelectorChoiceC3),isolate(input$selectScoreTypeC3),isolate(input$selectChrType))
    }else{
      showProgress()
      comp3view1Plot(input$cancerSelectorChoiceC3,input$selectScoreTypeC3,input$selectChrType)
    }
  })
  ## view 2
  output$perctAffectedSamplesPlot <- renderPlot({
    input$refreshPlotC3
    if (length(isolate(input$cancerSelectorChoiceC3))>0 && length(isolate(input$selectScoreTypeC3))>0 && length(isolate(input$selectChrType))>0 ){
      showProgress()
      comp3view2Plot(isolate(input$cancerSelectorChoiceC3),isolate(input$selectScoreTypeC3),isolate(input$selectChrType))
    }else{
      showProgress()
      comp3view2Plot(input$cancerSelectorChoiceC3,input$selectScoreTypeC3,input$selectChrType)
    }
  })
  
  ## downloads
  ## view 1
  output$downloadDataView1C3 <- downloadHandler(
    filename = function() {
      paste('plot-', Sys.Date(), '.jpeg', sep="")
    },
    content <- function(file){
      jpeg(filename=file,width=1000,height=500,units="px")
      comp3view1Plot(input$cancerSelectorChoiceC3,input$selectScoreTypeC3,input$selectChrType)
      dev.off()
    }
  )
  ## view 2
  output$downloadDataView2C3 <- downloadHandler(
    filename = function() {
      paste('plot-', Sys.Date(), '.jpeg', sep="")
    },
    content <- function(file){
      jpeg(filename=file,width=1000,height=500,units="px")
      comp3view2Plot(input$cancerSelectorChoiceC3,input$selectScoreTypeC3,input$selectChrType)
      dev.off()
    }
  )
  
  
  ###################################################################################
  

  
  ###################################################################################
  ## comp4: Pathways 
  ###################################################################################
  ## widgets
  ## cancer selector comp4
  output$cancerSelectorC4 <- renderUI({
    selectInput("cancerSelectorChoiceC4", "Select a Cancer type" , choices = c('All',apply(cancers, 1, function(r) r)), selected ="All")
  })
  ## score type selector comp4
  output$scoreSelectInputC4 <- renderUI({
    selectInput("selectScoreTypeC4", label = "Select Score type", 
                choices = list("Oncogene score" = "og.score", "Tumor suppressor score" = "ts.score",
                               "Combined score" = "combined.score"), selected = "og.score")
  })  
  ## pathway selector comp4
  updateSelectizeInput(session, 'pathwaySelectorChoiceC4', choices = pathways, selected = pathways[1], server = TRUE)
  ## sample set selector comp4
  output$sampleSelectorC4 <- renderUI({
    radioButtons("sampleSelectorC4", "Select Sample type",
                 choices = list("Tumors" = "tcga", "Cell lines" = "ccle", "Tumors vs Cell lines" = "both"),
                 selected = "tcga")
  })
  
  ## plots
  output$pathwayPlot <- renderImage({
    input$refreshPlotC4
    if (length(isolate(input$pathwaySelectorChoiceC4))>0 && length(input$cancerSelectorChoiceC4)>0 && length(input$selectScoreTypeC4)>0 && length(input$sampleSelectorC4)>0)
    {
      showProgress()
      list(src = generatePathview2(isolate(input$pathwaySelectorChoiceC4), isolate(input$cancerSelectorChoiceC4),isolate(input$sampleSelectorC4),
                                   isolate(input$selectScoreTypeC4)),
           contentType = 'image/png',
           width = 1000,
           height = 1000,
           alt = "This is alternate text")
    }else{
      showProgress()
      list(src = generatePathview2(input$pathwaySelectorChoiceC4, input$cancerSelectorChoiceC4,input$sampleSelectorC4,
                                   input$selectScoreTypeC4),
           contentType = 'image/png',
           width = 1000,
           height = 1000,
           alt = "This is alternate text")
      }
    }, deleteFile = TRUE)
  
  
  ## downloads
  output$downloadPlotC4 <- downloadHandler(
    filename = function() { 
      paste('plot-', Sys.Date(), '.png', sep="")
    },
    content <- function(file){
      file.copy(generatePathview2(input$pathwaySelectorChoiceC4, input$cancerSelectorChoiceC4,input$sampleSelectorC4,
                                  input$selectScoreTypeC4),file)
    }
  )
  

  ###################################################################################
  
  ###################################################################################
  ## comp5: Summary Statistics
  ###################################################################################
  
  ## tables
  output$sampleOverviewC5 <- renderDataTable({
    sampleOverview
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  output$genesCutoff1C5 <- renderDataTable({
    genesCutoffOne
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  output$genesCutoff2C5 <- renderDataTable({
    genesCutoffTwo
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  output$genesCutoff3C5 <- renderDataTable({
    genesCutoffThree
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  output$genesCutoff4C5 <- renderDataTable({
    genesCutoffFour
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5)
  )
  
  
  ###################################################################################
  
  
})
