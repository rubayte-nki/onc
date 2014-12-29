
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinysky)

## source plot functions
source("plots/demoPlot2.R")
source("plots/demoPlot.R")



readProjectOverview <- function(){
  con <- file("www/oncoscape.txt")
  content <- readLines(con)
  close(con)
  content
}



shinyServer(function(input, output, session) {
  
  ###################################################################################
  ## load starter rdata object for widgets
  ###################################################################################  
  starter <- load("www/starter.RData")
  genes <- apply(geness, 1, function(r) r)
  chrms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
             "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  ## this will return the genes from the selected pathway in comp4
  ## currently set to false 4 values
  demoGenesByPathway <- c(genes[1],genes[100],genes[1000],genes[2000])
  projectDescription <- readProjectOverview()

  ###################################################################################
  ## comp0: Project overview/about
  ###################################################################################  
  showshinyalert(session,"welcomeToApp",styleclass="primary",
    HTML("
        <h1 style='color:#1C1C1C' >OncoScape</h1>
        <span style='color:red'><b>beta version</b></span>
        
        <p style='color:#1C1C1C' > 
        OncoScape is a package for gene prioritization in the R statistical programming environment. 
        The analysis is run in a contrast fashion, i.e. always two groups of samples are compared with each other. 
        Examples include:
        1. tumors vs. normals
        2. cell lines vs. normals
        3. treatment responders vs resistant
        4. samples with mutations in gene X vs wild type
        Currently, analyses of five data types are implemented in OncoScape:
        1. gene expression
        2. DNA copy number
        3. DNA methylation
        4. mutation
        5. shRNA knock-down data
        Aberrations in each gene are called for each data type separately and scored as 0 (no aberration found) or 1 (aberration found). 
        These scores are summed across data types to give the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) 
        aberrations and calculates independent scores for both directions. It is possible to run the analysis on any combination of these data types.
        </p>
        
        <button class='btn btn-primary'>
          <font size='3'><b>Get started</b></font>
        </button>
        
        </br>
        
        <hr>
        <h3 >Publicatoin</h3>
        <p>some text about about this header</p>
        <h3 >License, terms of use, privacy </h3>
        <p>some text about about this header</p>
        <h3 >Github repository (?)</h3>
        <p>some text about about this header</p>
        <h3 >Contact</h3>
        <p>some text about about this header</p>
        "
    )
  )  
  ###################################################################################  
  
  

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
                 choices = list("Tumors" = "tumors", "Cell lines" = "cell-lines"),
                 selected = "tumors")
  })
  ## score type selector comp1
  output$scoreSelectInputC1 <- renderUI({
    selectInput("selectScoreTypeC1", label = "Score type", 
                choices = list("Choice 1" = "1", "Choice 2" = "2",
                               "Choice 3" = "3"), selected = "1")
  })
  ## number of genes text box
  output$numberOfGenesSelectInput <- renderUI({
    textInput("numberOfGenes", label = "# of Top Genes", 
              value = "50")
  })
  
  ## tables
  output$genesResTable <- renderDataTable({
    input$refreshPlot
    if (length(isolate(input$numberOfGenes))>0){
      data.frame(genes[1:isolate(input$numberOfGenes)])
    }else{
      data.frame(genes[1:input$numberOfGenes])
    }
  }, options = list(lengthMenu = list(c(5, 15, 25, 50, -1), list('5', '15', '25', '50', 'All')), pageLength = 5))
  
  ## plots
  ## view 1
  output$distPlot <- renderPlot({
    input$refreshPlot
    if (length(isolate(input$numberOfGenes))>0){
      demoPlot2(isolate(input$numberOfGenes))
    }else{
      demoPlot2(input$numberOfGenes)
    }
  })
  ## view 2
  output$distPlot2 <- renderPlot({
    input$refreshPlot
    if (length(isolate(input$numberOfGenes))>0){
      demoPlot2(isolate(input$numberOfGenes))
    }else{
      demoPlot2(input$numberOfGenes)
    }
  })
  ## *******************************************************************************
  ## verbatim debug text
  output$cancer <- renderText({
    input$refreshPlot
    paste("cancer : ", isolate(input$cancerSelectorChoiceC1))
  })
  output$score <- renderText({
    input$refreshPlot
    paste("score : ", isolate(input$selectScoreTypeC1))
  })
  output$number <- renderText({
    input$refreshPlot    
    paste("number : ", isolate(input$numberOfGenes))
  })
  output$sample <- renderText({
    input$refreshPlot    
    paste("sample : ", isolate(input$sampleSelectorC1))
  })
  
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
  updateSelectizeInput(session, 'geneSelectorChoiceC2', choices = genes, server = TRUE)
  ## sample set selector comp2
  output$sampleSelectorC2 <- renderUI({
    radioButtons("sampleSelectorC2", label = "Select Sample Set",
                 choices = list("Tumors" = 1, "Cell-lines" = 2, "tumors vs cell-lines" = 3),
                 selected = 1)
  })
  
  ## gene expression
  output$geneExpressionPlot <- renderPlot({
    demoPlot()
  })
  ## copy number
  output$cnvPlot <- renderPlot({
    demoPlot()
  })
  ## dna methylation
  output$dnaMethPlot <- renderPlot({
    demoPlot()
  })
  ## achilles
  output$achillesPlot <- renderPlot({
    demoPlot()
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
  output$selectedScorePlot <- renderPlot({
    demoPlot()
  })
  ## view 2
  output$perctAffectedSamplesPlot <- renderPlot({
    demoPlot()
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
  output$pathwayPlot <- renderPlot({
    demoPlot()
  })
  ###################################################################################
  

  
  ###################################################################################
  ## comp5: About 
  ###################################################################################
  ## html
  output$projectDescriptionContent <- renderUI({
    #paste(projectDescription, collapse="")
    HTML(paste(projectDescription, collapse=""))
  })
  
})
