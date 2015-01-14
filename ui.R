
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
                                 #tags$h4("Publication"),
                                 #tags$p("some text about about this header."),
                                 #tags$br(),
                                 #tags$h4("Github repository"),
                                 #tags$p("some text about about this header."),
                                 #tags$br(),
                                 tags$h4("App version"),
                                 tags$p("beta"),
                                 tags$br(),
                                 tags$h4("License, terms of use, privacy"),
                                 tags$p("needs to be defined"),
                                 tags$br(),
                                 tags$h4("Contacts"),
                                 tags$p("If you have any questions or suggestions regarding OncoScape or this app, please contact Lodewyk Wessels (l.wessels@nki.nl)."),
                                 tags$br(),
                                 tags$hr(),
                                 tags$img(src="NKIlogo.png")
                                 )
                              ),
                          column(7,
                                 tags$h3("What is OncoScape?"),
                                 tags$hr(),
                                 tags$p("OncoScape is a package for cancer gene prioritization for the R statistical programming environment. It compares molecular profiling data of two groups of samples in order to identify genes that show significant differences between these groups. Currently, OncoScape performs an analysis of the following five data types:"),
                                 tags$ul(tags$li("Gene Expression"),
                                         tags$li("DNA Copy Number"),
                                         tags$li("DNA Methylation"),
                                         tags$li("mutation"),
                                         tags$li("shRNA knock-down data")),
                                 tags$p("Aberrations in each gene are called for each data type separately and scored as 0 (no significant difference found) or 1 (significant difference found). These scores are summed across all data types giving the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) aberrations and calculates independent scores for both directions, the oncogene score and tumor suppressor score, respectively. Furthermore, a combined score is calculated as oncogene score minus tumor suppressor score."),
                                 tags$p("OncoScape can be applied to the comparison of any arbitrary groups of samples, such as:"),
                                 tags$ul(tags$li("tumors vs. normals"),
                                         tags$li("cell lines vs. normals"),
                                         tags$li("samples sensitive to treatment vs resistant ones"),
                                         tags$li("samples with mutations in gene X vs wild type ones"),
                                 	 tags$li("different cancer subtypes")),
                                 tags$p("This web page provides access to the results of our comprehensive analysis of tumor samples from 11 cancer types and cell lines from 10 cancer types. Tumor data was obtained from The Cancer Genome Atlas (TCGA) and cell line data from the Cancer Cell Line Encyclopedia (CCLE).")
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
                            tags$h4("Select Cancer type"),
                            uiOutput("cancerSelectorC1"),
                            tags$hr(),
                          
                            ## gene selection criteria
                            tags$h4("Gene Selection"),
                            ## score type
                            uiOutput("scoreSelectInputC1"),                                                    
                            ## number of cutoff genes
                            uiOutput("scoreCutoffSelectorC1"),
                            tags$h6("OR"),
                            fileInput('geneListUploadC1', 'Upload Gene List File',
                                      accept=c('text/csv', 
                                               'text/comma-separated-values,text/plain', 
                                               '.csv')),
                            tags$hr(),
                            
                            ## tumors or cell-lines
                            tags$h4("Select Sample type"),
                            uiOutput("sampleSelectorC1"),
                            tags$hr(),
                          
                            ## action button
                            actionButton("refreshPlot",label="Refresh Results",class='btn btn-primary')
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
                                     ## view type 1
                                     tabPanel("Detailed Aberration Profiles",plotOutput("distPlot2")),
                                     ## view type 2
                                     tabPanel("Summary Heat-Map",plotOutput("distPlot"))
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
                                   uiOutput("scoreSelectInputC3"),
                                 
                                   ## chr select
                                   uiOutput("chrSelector"),
                                   
                                   ## tumors or cell-lines
                                   #uiOutput("sampleSelectorC3"),
                                   tags$hr(),
                                   
                                   ## action button
                                   actionButton("refreshPlotC3",label="Refresh Results",class='btn btn-primary')                         
                                 )     
                        ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Score Plot"),
                                    downloadButton('downloadDataView1C3', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    plotOutput("selectedScorePlot"))
                          ),
                          tags$hr(),
                          fluidRow(
                            column (10, HTML("<h3>Affected Samples Plot"),
                                    downloadButton('downloadDataView2C3', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    plotOutput("perctAffectedSamplesPlot"))
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