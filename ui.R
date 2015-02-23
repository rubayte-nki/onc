
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinysky)

shinyUI(
  bootstrapPage(
  includeCSS("www/styles.css"),
  
  navbarPage("OncoScape", windowTitle ="OncoScape",

             ## comp0
             tabPanel("What is OncoScape?",
                          fluidRow(
                            column(2,wellPanel(
                                   HTML("<h4>App version</h4>"),
                                   HTML("<p>v0.9_beta</p>"),
                                   HTML("<br/>"),
                                   HTML("<h4>License, terms of use, privacy</h4>"),
                                   HTML("<p><a href='http://www.r-project.org/Licenses/GPL-3' target='_blank'> GNU General Public License version 3</a></p>"),
                                   HTML("<br/>"),
                                   HTML("<h4>Contacts</h4>"),
                                   HTML("<p>If you have any questions or suggestions regarding OncoScape or this app, please contact Lodewyk Wessels (l.wessels@nki.nl).</p>"),
                                   HTML("<br/>"),
                                   HTML("<h4>Shiny WebApp</h4>"),
                                   HTML("<p>This shiny web app is developed & hosted by Research - IT department (contact: r.rahman@nki.nl) at the  <a href='http://www.nki.nl/' target='_blank'>Netherlands Cancer Institute - NKI</a>."),
                                   HTML("<br/>"),
                                   HTML("<br/>"),
                                   shiny::tags$hr(),
                                   HTML("<img src='NKIlogo.png' class='img-responsive' />")
                                   )),
                            column(8,
                                   HTML("<h3>What is OncoScape?</h3>"),
                                   HTML("<hr>"),
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
                                   shiny::tags$hr(),
                                   shiny::tags$p(shiny::tags$span(style="color:blue", "Summary statistics about OncoScape can be found in the 'Summary Statistics' tab of the app."))
                            ),
                            column(2)
                            )
             ),
             #),
             
             
             ## comp1
             tabPanel("Top Genes in Cancers", 
                      fluidRow(
                        column(12,
                               HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Display/Refresh Results' to populate table data and generate plots !</strong>
                                    </div>
                                    ")
                               )
                        ),
                      fluidRow(
                        column(2,
                          wellPanel(
                            ## cancer type select
                            HTML("<h4>Select Cancer type</h4>"),
                            uiOutput("cancerSelectorC1"),
                            shiny::tags$hr(),   
                            
                            ## score type select
                            HTML("<h4>Select Score type</h4>"),
                            uiOutput("scoreTypeSelectorC1"),
                            shiny::tags$hr(),
                            
                            ## gene selection criteria
                            HTML("<h4>Gene Selection</h4>"),
                            uiOutput("geneSelectionMethodC1"),
                            uiOutput("geneSelectionPanelC1"),
                            
                            conditionalPanel(
                              condition = "input.geneSelectionMethodC1Value == 'type2'",
                              HTML("
                                    <a class='btn btn-primary' data-toggle='collapse' href='#collapseExample' aria-expanded='false' aria-controls='collapseExample'>
                                   Show help on Input File !
                                   </a>
                                    <div class='collapse' id='collapseExample'>
                                      <div class='well'>
                                          <p><strong>Notice!</strong> Input file must be tab-separated with the first column having the HGNC gene symbols. 
                                          Header line should also be provided.</p>
                                          <p><strong>Another thing!</strong>Once you upload a file, the uploaded data remains active in your current session.
                                          So uploading multiple times is not needed.</p>
                                      </div>
                                    </div>
                                   ")
                              ),
                            
                            conditionalPanel(
                              condition = "input.geneSelectionMethodC1Value == 'type3'",
                              shiny::actionButton("actionAutoFillGeneTextArea", label = "Click to load example !")
                            ),
                            
                            
                            ## gene selection procedure
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type1' && input.selectScoreTypeC1 == 'combined.score'",
#                               selectInput("scoreCutoff", label = "Score cut-off", 
#                                           choices = list(">=2" = 2, ">=3" = 3, ">=4" = 4, "=<-2" = -2, "=<-3" = -3, "=<-4" = -4),
#                                           selected = 3)
#                             ),
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type1' && input.selectScoreTypeC1 != 'combined.score'",
#                               selectInput("scoreCutoff", label = "Score cut-off", 
#                                           choices = list(">=2" = 2, ">=3" = 3, ">=4" = 4),
#                                           selected = 3)
#                             ),
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type2'",
#                               fileInput('geneListUploadC1', 'Upload Gene List File',accept = c(".tsv"))
#                             ),
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type3'",
#                               shiny::tags$textarea(id="geneListValuesC1", rows=10, cols=10, "Copy Paste your genes here separated by comma")
#                             ),
                            shiny::tags$hr(),
                            
                            ## tumors or cell-lines
                            HTML("<h4>Select Sample type</h4>"),
                            uiOutput("sampleSelectorC1"),
                            shiny::tags$hr(),
                          
                            ## action button
                            shiny::actionButton("refreshPlot",label="Display/Refresh Results",class='btn btn-primary')
                            )
                          ),
                        column(10,
                          fluidRow(
                            ## result genes
                            column(10,
                                   HTML("<h3> Result Gene(s)"),
                                   downloadButton('downloadData', 'Download Data', class='btn btn-primary'),
                                   HTML("</h3>"),
                                   shiny::dataTableOutput('genesResTable')
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
             tabPanel("Single Gene in Cancers",
                      fluidRow(
                        column(12,
                               HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Display/Refresh Results' to generate plots !</strong>
                                    </div>
                                    ")
                        )
                      ),
                      fluidRow(
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
                                 shiny::actionButton("refreshPlotC2",label="Display/Refresh Results",class='btn btn-primary')
                               )   
                        ),
                        column(10,
                          ## plot window
                          fluidRow(
                            column (5, HTML("<h3>Gene Expression"),
                                    downloadButton('downloadPlotC2GE', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    plotOutput("geneExpressionPlot")),
                            column (5, HTML("<h3>Copy Numbers"),
                                    downloadButton('downloadPlotC2CNA', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    plotOutput("cnvPlot"))
                            ),
                          fluidRow(
                            column (5, HTML("<h3>Achilles"),
                                    downloadButton('downloadPlotC2A', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    plotOutput("achillesPlot"))
                          )
                        )
                      )
             ),
             
             ## comp3
             tabPanel("Genes across Choromosomes",
                      fluidRow(
                        column(12,
                               HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Display/Refresh Results' to generate plots !</strong>
                                    </div>
                                    ")
                        )
                      ),
                      fluidRow(
                          column(2,
                                 wellPanel(
                                   ## cancer type select
                                   uiOutput("cancerSelectorC3"),
                                 
                                   ## score type select
                                   uiOutput("scoreSelectInputC3"),
                                 
                                   ## sample type select
                                   uiOutput("sampleSelectInputC3"),
                                   
                                   ## chr select
                                   #uiOutput("chrSelector"),
                                   
                                   ## tumors or cell-lines
                                   #uiOutput("sampleSelectorC3"),
                                   shiny::tags$hr(),
                                   
                                   ## action button
                                   shiny::actionButton("refreshPlotC3",label="Display/Refresh Results",class='btn btn-primary')                         
                                 )     
                        ),
                        column(10,
                          ## plot window
                          fluidRow(
                            #column (10, 
                                    HTML("<h3>Result Plot"),
                                    downloadButton('downloadDataView1C3', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    plotOutput("selectedScorePlotNew")#)
                          )#,
                          #shiny::tags$hr(),
                          #fluidRow(
                          #  column (10, HTML("<h3>Affected Samples Plot"),
                          #          downloadButton('downloadDataView2C3', 'Download Plot', class='btn btn-primary'),
                          #          HTML("</h3>"),
                          #          plotOutput("perctAffectedSamplesPlot"))
                          #)
                        )
                      )
             ),
             
             ## comp4
             tabPanel("Pathways",
                      fluidRow(
                        column(12,
                               HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Display/Refresh Results' to generate plots !</strong>
                                    </div>
                                    ")
                        )
                      ),
                      fluidRow(
                        column(2,
                               wellPanel(
                                 ## select cancer type
                                 uiOutput("cancerSelectorC4"),
                                 
                                 ## select score type
                                 uiOutput("scoreSelectInputC4"),
                                 
                                 
                                 ## pathway select
                                 selectizeInput('pathwaySelectorChoiceC4', label = "Select Pathway",  choices = NULL, options = list(maxItems = 1,placeholder="Type and Select a Pathway")),
                                 
                                 ## tumors or cell lines
                                 uiOutput("sampleSelectorC4"),
                                 shiny::tags$hr(),
                                 
                                 ## action button
                                 shiny:: actionButton("refreshPlotC4",label="Display/Refresh Results",class='btn btn-primary')                         
                               )
                      ),
                      column(10,
                             fluidRow(
                          ## plot window
                                    HTML("<h3>Pathway Plot"),
                                    downloadButton('downloadPlotC4', 'Download Plot', class='btn btn-primary'),
                                    HTML("</h3>"),
                                    imageOutput("pathwayPlot")
                            )
                          
                        )
                      )
             ),
             
             ## comp5
             tabPanel("Summary statistics",
                      fluidRow(
                          column(3),
                            column(6,
                                   shiny::tags$h3("Number of samples analyzed"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("The figure shows the number of cell lines, normal samples and tumor samples analyzed for each cancer type and data type. White dots indicate the number of matched normal samples available for DNA methylation (Meth) and copy number altersion (CNA)."),
                                   HTML("<img src='sample_set_sizes.png' class='img-responsive' />"),
                                   shiny::tags$p("The exact numbers of samples are contained in the following table. Each row contains the number of samples available for the specified data type. Rows 5 and 6 contain the number of samples with available mRNA expression data and copy-number data or DNA methylation data, respectively. The last row contains the number of cell lines available in CCLE for each cancer type."),
                                   shiny::dataTableOutput('sampleOverviewC5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Oncogene alterations across cancers"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Overview of the number of genes with activating aberrations for each cancer type and data type. The last row shows the number of genes with aberrations in any of the cancer types."),
                                   HTML("<img src='oncogene_alterations.png' class='img-responsive' />"),
                                   shiny::tags$h3("Tumor suppressor alterations across cancers"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Overview of the number of genes with inactivating aberrations for each cancer type and data type. The last row shows the number of genes with aberrations in any of the cancer types."),
                                   HTML("<img src='tumor_suppressor_alterations.png' class='img-responsive' />"),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 1"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with scores greater of equal to 1. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff1C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 2"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with scores greater of equal to 2. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff2C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 3"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with scores greater of equal to 3. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff3C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 4"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with scores greater of equal to 4. In case of the combined score (calculated as oncogene score (OG) ? tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff4C5'),
                                   shiny::tags$br()
                                   
                            ),
                          column(3)
                          )  
                      )
                      

          )
  
    )
)