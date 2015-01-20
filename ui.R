
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
             tabPanel("What is OncoScape?",
                      sidebarLayout(
                        column(2,
                               wellPanel(
                                 # shiny::tags$h4("Publication"),
                                 # shiny::tags$p("some text about about this header."),
                                 # shiny::tags$br(),
                                 # shiny::tags$h4("Github repository"),
                                 # shiny::tags$p("some text about about this header."),
                                 # shiny::tags$br(),
                                 HTML("<h4>App version</h4>"),
                                 HTML("<p>beta</p>"),
                                 HTML("<br/>"),
                                 HTML("<h4>License, terms of use, privacy</h4>"),
                                 HTML("<p>needs to be defined</p>"),
                                 HTML("<br/>"),
                                 HTML("<h4>Contacts</h4>"),
                                 HTML("<p>If you have any questions or suggestions regarding OncoScape or this app, please contact Lodewyk Wessels (l.wessels@nki.nl).</p>"),
                                 HTML("<br/>"),
                                 HTML("<hr/>"),
                                 HTML("<img src='NKIlogo.png'/>")
                               )
                        ),
                        mainPanel(
                          ## app about text
                          fluidRow(
                            column(9,
                                   HTML("<h3>What is OncoScape?</h3>"),
                                   HTML("<p>OncoScape is a package for cancer gene prioritization for the R statistical programming environment. It compares molecular profiling data of two groups of samples in order to identify genes that show significant differences between these groups. Currently, OncoScape performs an analysis of the following five data types:</p>"),
                                   shiny::tags$ul(shiny::tags$li("Gene Expression"),
                                                  shiny::tags$li("DNA Copy Number"),
                                                  shiny::tags$li("DNA Methylation"),
                                                  shiny::tags$li("mutation"),
                                                  shiny::tags$li("shRNA knock-down data")),
                                   shiny::tags$p("Aberrations in each gene are called for each data type separately and scored as 0 (no significant difference found) or 1 (significant difference found). These scores are summed across all data types giving the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) aberrations and calculates independent scores for both directions, the oncogene score and tumor suppressor score, respectively. Furthermore, a combined score is calculated as oncogene score minus tumor suppressor score."),
                                   shiny::tags$p("OncoScape can be applied to the comparison of any arbitrary groups of samples, such as:"),
                                   shiny::tags$ul(shiny::tags$li("tumors vs. normals"),
                                                  shiny::tags$li("cell lines vs. normals"),
                                                  shiny::tags$li("samples sensitive to treatment vs resistant ones"),
                                                  shiny::tags$li("samples with mutations in gene X vs wild type ones"),
                                                  shiny::tags$li("different cancer subtypes")),
                                   shiny::tags$p("This web page provides access to the results of our comprehensive analysis of tumor samples from 11 cancer types and cell lines from 10 cancer types. Tumor data was obtained from The Cancer Genome Atlas (TCGA) and cell line data from the Cancer Cell Line Encyclopedia (CCLE)."),
                                   shiny::tags$p(shiny::tags$span(style="color:blue", "Summary statistics about OncoScape can be found in the 'Summary Statistics' tab of the app."))
                            ))
                            
                          )  
                        )
             ),
             #),
             
             
             ## comp1
             tabPanel("Top Genes in Cancer types",                      
                      sidebarLayout(
                        column(2,
                          wellPanel(
                            ## cancer type select
                            HTML("<h4>Select Cancer type</h4>"),
                            uiOutput("cancerSelectorC1"),
                            shiny::tags$hr(),
                          
                            ## gene selection criteria
                            HTML("<h4>Gene Selection</h4>"),
                            ## score type
                            uiOutput("scoreSelectInputC1"),                                                    
                            ## number of cutoff genes
                            uiOutput("scoreCutoffSelectorC1"),
                            shiny::tags$h6(" OR "),
                            fileInput('geneListUploadC1', 'Upload Gene List File',accept = c(".tsv")),
                            shiny::tags$hr(),
                            
                            ## tumors or cell-lines
                            HTML("<h4>Select Sample type</h4>"),
                            uiOutput("sampleSelectorC1"),
                            shiny::tags$hr(),
                          
                            ## action button
                            shiny::actionButton("refreshPlot",label="Refresh Results",class='btn btn-primary')
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
                                 shiny::tags$hr(),
                               
                                 ## action button
                                 shiny::actionButton("refreshPlotC2",label="Refresh",class='btn btn-primary')
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
                                   shiny::tags$hr(),
                                   
                                   ## action button
                                   shiny::actionButton("refreshPlotC3",label="Refresh Results",class='btn btn-primary')                         
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
                          shiny::tags$hr(),
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
                                 ## select cancer type
                                 uiOutput("cancerSelectorC4"),
                                 
                                 ## select score type
                                 uiOutput("scoreSelectInputC4"),
                                 
                                 
                                 ## multilpe gene select
                                 selectizeInput('pathwaySelectorChoiceC4', label = "Select Pathway",  choices = NULL, options = list(maxItems = 1,placeholder="Selected Genes")),
                                 
                                 ## tumors or cell lines
                                 uiOutput("sampleSelectorC4"),
                                 shiny::tags$hr(),
                                 
                                 ## action button
                                 shiny:: actionButton("refreshPlotC4",label="Refresh Results",class='btn btn-primary')                         
                               )
                      ),
                        mainPanel(
                          ## plot window
                          fluidRow(
                            column (10, HTML("<h3>Pathway Plot"),
                                    downloadButton('downloadPlotC4', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    imageOutput("pathwayPlot",height="1000"))
                          )
                        )
                      )
             ),
             
             ## comp5
             tabPanel("Summary statistics",
                      sidebarLayout(
                        column(2),
                        mainPanel(
                          ## app about text
                          fluidRow(
                            column(10,
                                   shiny::tags$h3("Number of samples analyzed"),
                                   shiny::tags$p("The following table gives an overview of the number of samples analyzed for each cancer type. Each row contains the number of samples available for the specified data type. Rows 5 and 6 contain the number of samples with available mRNA expression data and copy-number data or DNA methylation data, respectively. The last row contains the number of cell lines available in CCLE for each cancer type."),
                                   dataTableOutput('sampleOverviewC5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 1"),
                                   shiny::tags$p("Number of genes with scores greater of equal to 1. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   dataTableOutput('genesCutoff1C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 2"),
                                   shiny::tags$p("Number of genes with scores greater of equal to 2. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   dataTableOutput('genesCutoff2C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 3"),
                                   shiny::tags$p("Number of genes with scores greater of equal to 3. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   dataTableOutput('genesCutoff3C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 4"),
                                   shiny::tags$p("Number of genes with scores greater of equal to 4. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   dataTableOutput('genesCutoff4C5'),
                                   shiny::tags$br()
                                   
                            )
                          )  
                        )
                      )
                      
             )
             
    )
  
))