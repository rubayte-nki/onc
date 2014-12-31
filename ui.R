
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinysky)

shinyUI(fluidPage(

  tags$p(""),
    
  shinyalert("welcomeToApp", click.hide = TRUE, auto.close.after = NULL),
  
  
  navbarPage("OncoScape beta",
             
             
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
                            
                            ## number of cutoff genes
                            uiOutput("numberOfGenesSelectInput"),
                            
                            ## tumors or cell-lines
                            uiOutput("sampleSelectorC1"),
                            tags$hr(),
                          
                            ## action button
                            actionButton("refreshPlot",label="Refresh",class='btn btn-primary')
                          )
                        ),
                        mainPanel(
                          fluidRow(
                            ## plot window 
                            column(10,                            
                              HTML("<h3> Result Gene(s) plot"),
                              downloadButton('downloadPlot', 'Download Plot', class='btn btn-primary'),
                              HTML("</h3>"),
                              tabsetPanel(
                                ## heatmap view type 1
                                tabPanel("view1",plotOutput("distPlot")),
                                ## heatmap view type 1
                                tabPanel("view2",plotOutput("distPlot2"))
                                )
                            )
                          ),
                          fluidRow(
                            ## result genes
                            column(10,
                                   HTML("<hr><h3> Result Gene(s)"),
                                   downloadButton('downloadData', 'Download Data', class='btn btn-primary'),
                                   HTML("</h3>"),
                                   dataTableOutput('genesResTable')
                            )
                          ),
                          fluidRow(
                            column(10,
                                   HTML("<hr><span style='color:red'>debug</span><P></p>"),
                                   verbatimTextOutput("cancer"),
                                   verbatimTextOutput("score"),
                                   verbatimTextOutput("number"),
                                   verbatimTextOutput("sample")
                                   )
                            )
                        )
                      )
             ),
             
             ## comp2
             tabPanel("Single Gene in Cancer types",
                      sidebarLayout(
                        column(2,
                               ## cancer type select
                               uiOutput("cancerSelectorC2"),
                               
                               ## gene select
                               selectizeInput('geneSelectorChoiceC2', label = "Select Gene", choices = NULL, options = list(maxItems = 1,placeholder="Select Gene")),
                               
                               ## select sample set
                                uiOutput("sampleSelectorC2")
                               
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
                                 ## cancer type select
                                 uiOutput("cancerSelectorC3"),
                                 
                                 ## score type select
                                 uiOutput("selectScoreTypeC3"),
                                 
                                 ## chr select
                                 uiOutput("chrSelector"),
                                 
                                 ## tumors or cell-lines
                                 uiOutput("sampleSelectorC3")  
                        ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Selected Scores</h3>"),
                                    plotOutput("selectedScorePlot"))
                          ),
                          fluidRow(
                            column (10, HTML("<h3>Percentage Affected Samples</h3>"),
                                    plotOutput("perctAffectedSamplesPlot"))
                          )
                        )
                      )
             ),
             
             ## comp4
             tabPanel("Pathways",
                      sidebarLayout(
                        column(2,
                               ## pathway type select
                               uiOutput("pathwaySelector"),
                               
                               ## multilpe gene select
                               selectizeInput('geneSelectorChoiceC4', label = "Selected Gene",  choices = NULL, options = list(maxItems = 500,placeholder="Selected Genes")),
                               
                               ## tumors or cell lines
                               uiOutput("sampleSelectorC4")
                        ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Pathway Plot</h3>"),
                                    plotOutput("pathwayPlot"))
                          )
                        )
                      )
             ),
             
             
             ## comp5
             tabPanel("About",
                        mainPanel(
                          ## app about text
                          column(12,
                                 tags$h3("OncoScape"),
                                 uiOutput("projectDescriptionContent"),
                                 tags$br(),
                                 tags$h3("Publication"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h3("License, terms of use, privacy"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h3("Github repository (?)"),
                                 tags$p("some text about about this header."),
                                 tags$br(),
                                 tags$h3("Contact"),
                                 tags$p("some text about about this header."),
                                 tags$br()
                                 )  
                        )
             )
             
             
  ) 
  
))