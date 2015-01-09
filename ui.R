
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinysky)

shinyUI(fluidPage(

  navbarPage("OncoScape",
             windowTitle ="OncoScape",
             icon="www/NKIlogo.png",

             ## comp0
             tabPanel("What is OncoScape?",collapsable = TRUE,
                      mainPanel(
                        ## app about text
                        fluidRow(
                          column(5,
                                 wellPanel(
                                 tags$h4("Publication"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h4("Github repository"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h4("Contacts"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h4("App version"),
                                 tags$p("beta"),
                                 tags$br(),
                                 tags$h4("License, terms of use, privacy"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$hr(),
                                 tags$img(src="NKIlogo.png")
                                 )
                              ),
                          column(7,
                                 tags$h3("What is OncoScape?"),
                                 tags$hr(),
                                 tags$p("OncoScape is a package for gene prioritization in the R statistical programming environment. The analysis is run in a contrast fashion, i.e. always two groups of samples are compared with each other. Examples include:"),
                                 tags$ul(tags$li("Tumors vs. Normals"),
                                         tags$li("Cell lines vs. Normals"),
                                         tags$li("Treatment responders vs Resistant"),
                                         tags$li("Samples with mutations in gene X vs Wild type")),
                                 tags$p("Currently, analyses of five data types are implemented in OncoScape:"),
                                 tags$ul(tags$li("Gene Expression"),
                                         tags$li("DNA Copy Number"),
                                         tags$li("DNA Methylation"),
                                         tags$li("mutation"),
                                         tags$li("shRNA knock-down data")),
                                 tags$p("Aberrations in each gene are called for each data type separately and scored as 0 (no aberration found) or 1 (aberration found). 
                                        These scores are summed across data types to give the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) 
                                        aberrations and calculates independent scores for both directions. It is possible to run the analysis on any combination of these data types.")
                          )
                        )  
                      )
             ),
             
             
             ## comp1
             tabPanel("Top Genes in Cancer types",                      
                      sidebarLayout(
                        column(2,
                          wellPanel(
                            ## cancer type select
                            uiOutput("cancerSelectorC1"),
                          
                            ## score type
                            uiOutput("scoreSelectInputC1"),                          
                          
                            ## number of cutoff genes
                            uiOutput("scoreCutoffSelectorC1"),
                            
                            ## tumors or cell-lines
                            uiOutput("sampleSelectorC1"),
                            tags$hr(),
                          
                            ## action button
                            actionButton("refreshPlot",label="Refresh",class='btn btn-primary')
                            )
                          ),
                        mainPanel(
                          fluidRow(
                            ## result genes
                            column(10,
                                   HTML("<h3> Result Gene(s)"),
                                   downloadButton('downloadData', 'Download Data', class='btn btn-primary'),
                                   HTML("</h3>"),
                                   dataTableOutput('genesResTable')
                            )
                          ),
                          fluidRow(
                            ## plot window 
                            column(10,                            
                                   HTML("<hr><h3> Result Gene(s) plot"),
                                   downloadButton('downloadPlot', 'Download Plot', class='btn btn-primary'),
                                   HTML("</h3>"),
                                   tabsetPanel(
                                     ## heatmap view type 1
                                     tabPanel("view1",plotOutput("distPlot")),
                                     ## heatmap view type 1
                                     tabPanel("view2",plotOutput("distPlot2"))
                                   )
                            )
                          )
                      )
                    )    
             ),
             
             ## comp2
             tabPanel("Single Gene in Cancer types",
                      sidebarLayout(
                        column(2,
                               wellPanel(
                                 ## cancer type select
                                 uiOutput("cancerSelectorC2"),
                                 
                                 ## gene select
                                 selectizeInput('geneSelectorChoiceC2', label = "Select Gene", choices = NULL, options = list(maxItems = 1,placeholder="Select Gene")),
                                 
                                 ## select sample set
                                 uiOutput("sampleSelectorC2"),
                                 tags$hr(),
                               
                                 ## action button
                                 actionButton("refreshPlotC2",label="Refresh",class='btn btn-primary')
                               )   
                        ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (5, HTML("<h3>Gene Expression</h3>"),
                                    plotOutput("geneExpressionPlot")),
                            column (5, HTML("<h3>Copy Numbers</h3>"),
                                    plotOutput("cnvPlot"))
                            ),
                          fluidRow(
                            column (5, HTML("<h3>DNA Methylation</h3>"),
                                    plotOutput("dnaMethPlot")),
                            column (5, HTML("<h3>Achilles</h3>"),
                                    plotOutput("achillesPlot"))
                          )
                        )
                      )
             ),
             
             ## comp3
             tabPanel("Genes across Choromosomes",
                      sidebarLayout(
                          column(2,
                                 wellPanel(
                                   ## cancer type select
                                   uiOutput("cancerSelectorC3"),
                                 
                                   ## score type select
                                   uiOutput("selectScoreTypeC3"),
                                 
                                   ## chr select
                                   uiOutput("chrSelector"),
                                   
                                   ## tumors or cell-lines
                                   uiOutput("sampleSelectorC3"),
                                   tags$hr(),
                                   
                                   ## action button
                                   actionButton("refreshPlotC3",label="Refresh",class='btn btn-primary')                         
                                 )     
                        ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Result Plot</h3>"),
                                    textOutput("selectedScorePlot"))
                          )
                        )
                      )
             ),
             
             ## comp4
             tabPanel("Pathways",
                      sidebarLayout(
                        column(2,
                               wellPanel(
                                 ## pathway type select
                                 uiOutput("pathwaySelector"),
                                 
                                 ## multilpe gene select
                                 selectizeInput('geneSelectorChoiceC4', label = "Selected Gene",  choices = NULL, options = list(maxItems = 500,placeholder="Selected Genes")),
                                 
                                 ## tumors or cell lines
                                 uiOutput("sampleSelectorC4"),
                                 tags$hr(),
                                 
                                 ## action button
                                 actionButton("refreshPlotC4",label="Refresh",class='btn btn-primary')                         
                               )
                      ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Pathway Plot</h3>"),
                                    textOutput("pathwayPlot"))
                          )
                        )
                      )
             )
             
    )
  
))