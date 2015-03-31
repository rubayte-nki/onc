library(shiny)
library(shinysky)

shinyUI(
  bootstrapPage(
  includeCSS("www/styles.css"),
  
  shinyUI(navbarPage(title="", windowTitle ="OncoScape",header="",

            ## comp0
             tabPanel("OncoScape",
                          fluidRow(
                            column(2),
                            column(8,
                                   fluidRow(class="frontPageBottom1Panel",
                                     HTML("<h3 class=\"frontPageHeader\" ><font color=\"#FFFFFF\">OncoScape</font></h3>"),
                                     HTML("<hr>"),
                                     HTML("<p class='frontPage'><font color=\"#FFFFFF\">OncoScape is a package for cancer gene prioritization for the R statistical programming environment. 
                                          It compares molecular profiling data of two groups of samples in order to identify genes that show significant differences 
                                          between these groups. Currently, OncoScape performs an analysis of the five data types: (A) Gene Expression 
                                          (B) DNA Copy Number (C) DNA Methylation (D) Mutation (E) shRNA knock-down data.</font></p>"),
                                     #shiny::tags$ul(shiny::tags$li("Gene Expression"),
                                      #              shiny::tags$li("DNA Copy Number"),
                                      #              shiny::tags$li("DNA Methylation"),
                                      #              shiny::tags$li("mutation"),
                                      #              shiny::tags$li("shRNA knock-down data")),
                                     HTML("<p class='frontPage'><font color=\"#FFFFFF\">Aberrations in each gene are called for each data type separately and scored as 0 
                                                   (no significant difference found) or 1 (significant difference found). These scores are summed across all 
                                                   data types giving the final score. OncoScape differentiates between activating (oncogene-like) and 
                                                   inactivating (tumor suppressor-like) aberrations and calculates independent scores for both directions, 
                                                   the oncogene score and tumor suppressor score, respectively. Furthermore, a combined score is calculated as 
                                                   oncogene score minus tumor suppressor score. OncoScape can be applied to the comparison of any arbitrary 
                                                   groups of samples, such as: (A) Tumors vs. Normals (B) Cell-lines vs. Normals (C) Samples sensitive
                                                   to treatment vs. resistant ones (D) Samples with mutations in gene X vs. wild type ones 
                                                   (E) Diferent cancer subtypes.</font></p>"),
                                     #shiny::tags$p("OncoScape can be applied to the comparison of any arbitrary groups of samples, such as:"),
                                     #shiny::tags$ul(shiny::tags$li("tumors vs. normals"),
                                    #                shiny::tags$li("cell lines vs. normals"),
                                     #               shiny::tags$li("samples sensitive to treatment vs resistant ones"),
                                    #                shiny::tags$li("samples with mutations in gene X vs wild type ones"),
                                    #                shiny::tags$li("different cancer subtypes")),
                                     ## shiny::tags$p("This web page provides access to the results of our comprehensive analysis of tumor samples from 11 cancer types and cell lines from 10 cancer types. Tumor data was obtained from The Cancer Genome Atlas (TCGA) and cell line data from the Cancer Cell Line Encyclopedia (CCLE)."),
                                     HTML("<p class='frontPage' ><font color=\"#FFFFFF\">This web page provides access to the results of our comprehensive analysis of 
                                          tumor samples from 11 cancer types and cell lines from 10 cancer types. Tumor data were obtained from The Cancer 
                                          Genome Atlas (TCGA) and cell line data from the Cancer Cell Line Encyclopedia (CCLE). The analyses are presented in the 
                                          following publication: <I>Exploring the cancer aberration landscape by genomic data fusion, Andreas Schlicker, 
                                          Magali Michaut, Rubayte Rahman and  Lodewyk FA Wessels (submitted)</I></font></p>"),
                                     
                                     #shiny::tags$p(shiny::tags$span(style="color:blue", "Summary statistics about OncoScape can be found in the 'Summary' tab of the app."))
                                    tags$head(tags$style("
                                                .frontPageBottom1Panel{background-color: #6E6E6E; border-radius: 10px;}"
                                          )   
                                      )
                                     ),
                                   fluidRow(
                                     column(12,
                                            shiny::tags$br()
                                     )
                                   ),
                                   
                                   # fluidRow(shiny::tags$hr()),
                                   fluidRow(class="frontPageBottomPanel",
                                                 column(3, class="frontPageBottomPanelColumn",
                                                            fluidRow(
                                                              HTML("<h4>Publication</h4>"),
                                                              HTML("<p>If you use OncoScape in a publication, please cite:
                                                              Schlicker A, Michaut M, Rahman R and Wessels LFA, Exploring the cancer aberration landscape by genomic data fusion (submitted) </p>"),
                                                              tags$head(tags$style("
                                                                .frontPageBottomPanelColumn{margin: 0.5cm 0.5cm 0.5cm 0.5cm;}"
                                                                )   
                                                              )
                                                            )),
                                                  #column(1),
                                                     column(3,class="frontPageBottomPanelColumn",
                                                            fluidRow(
                                                              HTML("<h4>Contacts</h4>"),
                                                              HTML("<p>If you have any questions or suggestions regarding OncoScape, please contact Lodewyk Wessels (l.wessels@nki.nl).</p>"),
                                                              tags$head(tags$style("
                                                                .frontPageBottomPanelColumn{margin: 0.5cm 0.5cm 0.5cm 0.5cm;}"
                                                              )   
                                                              )
                                                            )),
                                                    #column(1),
                                                     column(3,class="frontPageBottomPanelColumn",
                                                            fluidRow(
                                                              HTML("<h4>Shiny WebApp</h4>"),
                                                              HTML("<p>This shiny web app is developed & hosted by Research - IT department (r.rahman@nki.nl) at the  <a href='http://www.nki.nl/' target='_blank'>Netherlands Cancer Institute - NKI</a>."),
                                                              tags$head(tags$style("
                                                                .frontPageBottomPanelColumn{margin: 0.5cm 0.5cm 0.5cm 0.5cm;}"
                                                              )   
                                                              )
                                                            )),
                                                     #column(3,
                                                    #      )),
                                                tags$head(tags$style("
                                                .frontPageBottomPanel{background-color: #FAFAFA; border-radius: 10px;}"
                                                  )   
                                                ) 
                                                     
                                       #)    
                                    ),
                                   fluidRow(
                                     column(12,
                                            shiny::tags$br()
                                            )
                                     ),
                                   fluidRow(class="frontPageBottomPanel",
                                     column(9,class="frontPageBottomPanelColumn",
                                            fluidRow(
                                              HTML("<h4>License</h4>"),
                                              HTML("<p>Copyright &copy; 2015 Netherlands Cancer Institute</br>
                                                    Licensed under the MIT License.</p>"),
                                                   #<a href='http://www.r-project.org/Licenses/GPL-2' target='_blank'> GNU General Public License version 2</a></p>"),
                                              tags$head(tags$style("
                                                                .frontPageBottomPanelColumn{margin: 0.5cm 0.5cm 0.5cm 0.5cm;}"
                                              )   
                                              )
                                            )
                                          ),
                                     column(2,
                                            fluidRow(
                                              HTML("<img src='NKIlogo.png' class='img-responsive' />")
                                              
                                            )),
                                     tags$head(tags$style("
                                                .frontPageBottomPanel{background-color: #FAFAFA;border-radius: 10px;}"
                                              )   
                                      )
                                     )
                            ),
                            column(2)
                            )
             ),
             #),
             
             
             ## comp1
             tabPanel("Top Candidate Genes",
                      
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
                            
                            ## score cutoff criteria
                            HTML("<h4>Select Cut-off Score</h4>"),
                            uiOutput("scoreCutoffSelectorC1"),
                            shiny::tags$hr(),
                            
#                             ## gene selection criteria
#                             HTML("<h4>Gene Selection</h4>"),
#                             uiOutput("geneSelectionMethodC1"),
#                             
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type3'",
#                               shiny::actionButton("actionAutoFillGeneTextArea", label = "Load an example!", class='btn btn-link')
#                             ),
#                             
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type2'",
#                               HTML("
#                                     <a class='btn btn-link' data-toggle='collapse' href='#collapseExample' aria-expanded='false' aria-controls='collapseExample'>
#                                    Show help on Input File!
#                                    </a>
#                                     <div class='collapse' id='collapseExample'>
#                                       <div class='well'>
#                                           <p><strong>Notice!</strong> Input file must be tab-separated with the first column having the HGNC gene symbols. 
#                                           Header line should also be provided.</p>
#                                           <p><strong>Another thing!</strong>Once you upload a file, the uploaded data remains active in your current session.
#                                           So uploading multiple times is not needed.</p>
#                                       </div>
#                                     </div>
#                                    ")
#                             ),
#                             
#                             uiOutput("geneSelectionPanelC1"),                            
                            
                            
                            
                            ## gene selection procedure
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type1' && input.selectScoreTypeC1 == 'combined.score'",
#                               selectInput("scoreCutoff", label = "Score cut-off", 
#                                           choices = list(">=2" = 2, ">=3" = 3, ">=4" = 4, "=<-2" = -2, "=<-3" = -3, "=<-4" = -4),
#                                           selected = 3)
#                             ),
#                             conditionalPanel(
#                               condition = "input.geneSelectionMethodC1Value == 'type1' && input.selectScoreTypeC1!= 'combined.score'",
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
                            #shiny::tags$hr(),
                            
                            ## tumors or cell-lines
                            HTML("<h4>Select Sample type</h4>"),
                            uiOutput("sampleSelectorC1"),
                            shiny::tags$hr(),
                          
                            ## action button
                            shiny::actionButton("refreshPlot",label="Refresh Results",class='btn btn-primary')
                            )
                          ),
                        column(10,
                              fluidRow(
                                fluidRow(
                                  column(12,
                                         HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Refresh Results' to populate table data and generate plots!</strong>
                                    </div>
                                    ")
                                  )
                                ),
                                 ## plot window 
                                 column(10,                            

                                        tabsetPanel(
                                          ## view type 1
                                          tabPanel("Detailed Aberration Profiles",
                                                   HTML("<h4>Cancer-type-specific top candidate genes"),
                                                   downloadButton('downloadPlotDAPlot', 'Download', class='btn btn-link'),
                                                   HTML("</h4>"),
                                                   plotOutput("distPlot2")),
                                          ## view type 2
                                          tabPanel("Summary Heat-Map",
                                                   HTML("<h4>Cancer-type-specific top candidate genes"),
                                                   downloadButton('downloadPlotHPlot', 'Download', class='btn btn-link'),
                                                   HTML("</h4>"),
                                                   plotOutput("distPlot")),
                                          #tabPanel("Summary Heat-Map",dygraphOutput("distPlot",height="4000px"))
                                          tabPanel("Scores",
                                                   fluidRow(
                                                     column(12,
                                                            HTML("<h4>Cancer-type-specific top candidate genes"),
                                                            downloadButton('downloadData', 'Download', class='btn btn-link'),
                                                            HTML("</h4>")
                                                            )
                                                    # column(6,
                                                    #    HTML("</br><div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                                    #        <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                                    #         <strong>Table header definitions are listed below.</strong>
                                                    #          </div>")
                                                     #)
                                                     ),
                                                   fluidRow(column(12,
                                                                   dataTableOutput('genesResTable')
                                                                   )
                                                    ),
                                                   fluidRow(
                                                     column(12,
                                                            HTML("<hr>"),
                                                            HTML("
                                                          <a class='btn btn-link' data-toggle='collapse' href='#collapseExample' aria-expanded='false' aria-controls='collapseExample'>
                                                          Click to see table header descriptions! 
                                                          </a>
                                                          <div class='collapse' id='collapseExample'>
                                                            <div class='well'>
                                                                <h4>Table column descriptions</h4>
                                                                <ul>
                                                                    <li> <b>Genes</b>: HUGO gene symbol </li>
                                                                    <li> <b>OG Score</b>: Oncogene score for the respective gene</li>
                                                                    <li> <b>TS Score</b>: Tumor suppressor score for the respective gene</li>
                                                                    <li> <b>Combined Score</b>: Combined score, calculated as oncogene score - tumor suppressor score, for the respective gene</li>
                                                                    <li> <b>OG.Meth</b>: Indicates if the gene shows oncogene-like DNA methylation changes</li>
                                                                    <li> <b>OG.CNA</b>: Indicates if the gene shows oncogene-like DNA copy number alterations</li>
                                                                    <li> <b>OG.Mut</b>: Indicates if the gene shows oncogene-like somatic mutation profiles</li>
                                                                    <li> <b>OG.shRNA</b>: Indicates if the gene shows oncogene-like changes after shRNA knock-down</li>
                                                                    <li> <b>OG.Expr</b>: Indicates if the gene shows oncogene-like gene expression changes</li>
                                                                    <li> <b>TS.Meth</b>: Indicates if the gene shows tumor suppressor-like DNA methylation changes</li>
                                                                    <li> <b>TS.CNA</b>: Indicates if the gene shows tumor suppressor-like DNA copy number alterations</li>
                                                                    <li> <b>TS.Mut</b>: Indicates if the gene shows tumor suppressor-like somatic mutation profiles</li>
                                                                    <li> <b>TS.shRNA</b>: Indicates if the gene shows tumor suppressor-like changes after shRNA knock-down</li>
                                                                    <li> <b>TS.Expr</b>: Indicates if the gene shows tumor suppressor-like gene expression changes</li>
                                                                    <li> <b>Cancer</b>: Cancer type showing the association</li>
                                                                    <li> <b>External links</b>: Link to the GeneCards entry for the respective gene</li>
                                                                </ul>
                                                              </div>
                                                            </div>
                                                            ")
                                                      )
                                                     )
                                                    #uiOutput("cancerSelectorC1Sub"),
                                              )
                                        )
                                 )
                              )
                      )
                    )    
             ),
             
        ## comp6
        tabPanel("User-selected Genes", 
         
         fluidRow(
           column(2,
                  wellPanel(                    
                    ## gene selection criteria
                    HTML("<h4>Gene Upload Method</h4>"),
                    uiOutput("geneSelectionMethodC6"),
                    conditionalPanel(
                      condition = "input.geneSelectionMethodC6Value == 'type3'",
                      helpText("Note: Copy paste your genes in the text box below separated by comma(,)")
                    
                    ),
                    conditionalPanel(
                      condition = "input.geneSelectionMethodC6Value == 'type3'",
                      shiny::actionButton("actionAutoFillGeneTextArea", label = "Load an example!", class='btn btn-link')
                    ),
                    
                    conditionalPanel(
                      condition = "input.geneSelectionMethodC6Value == 'type2'",
                      HTML("
                           <a class='btn btn-link' data-toggle='collapse' href='#collapseExampleC6' aria-expanded='false' aria-controls='collapseExampleC6'>
                           Show help on Input File!
                           </a>
                           <div class='collapse' id='collapseExampleC6'>
                           <div class='well'>
                           <p><strong>Notice!</strong> Input file must be tab-separated with the first column having the HGNC gene symbols. 
                           Header line should also be provided.</p>
                           <p><strong>Another thing!</strong>Once your file is uploaded, the meassage 'Upload complete' will be visible. 
                            Please then click the 'Refresh Results' button to generate plots and scores.</p>
                           </div>
                           </div>
                           ")
                    ),
                    uiOutput("geneSelectionPanelC6"),                                                
                    shiny::tags$hr(),
                    
                    ## score type select
                    HTML("<h4>Select Score type</h4>"),
                    uiOutput("scoreTypeSelectorC6"),
                    shiny::tags$hr(),
                    
                    
                    
                    ## tumors or cell-lines
                    HTML("<h4>Select Sample type</h4>"),
                    uiOutput("sampleSelectorC6"),
                    shiny::tags$hr(),
                    
                    ## action button
                    shiny::actionButton("refreshPlotC6",label="Refresh Results",class='btn btn-primary')
                    )
                    ),
           column(10,
                  fluidRow(
                    fluidRow(
                      column(12,
                             HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                       <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                       <strong>Click on the 'Refresh Results' to populate table data and generate plots!</strong>
                       </div>
                       ")
                      )
                    ),
                    
                    ## plot window 
                    column(10,         
                           tabsetPanel(
                             ## view type 1
                             tabPanel("Detailed Aberration Profiles",
                                      HTML("<h4>User-selected genes across cancer-types"),
                                      downloadButton('downloadPlotDAPlotC6', 'Download', class='btn btn-link'),
                                      HTML("</h4>"),
                                      plotOutput("distPlot2C6")),
                             ## view type 2
                             tabPanel("Summary Heat-Map",
                                      HTML("<h4>User-selected genes across cancer-types"),
                                      downloadButton('downloadPlotHPlotC6', 'Download', class='btn btn-link'),
                                      HTML("</h4>"),
                                      plotOutput("distPlotC6")),
                             #tabPanel("Summary Heat-Map",dygraphOutput("distPlot",height="4000px"))
                             tabPanel("Scores",
                                      fluidRow(
                                        column(12,
                                               HTML("<h4>Cancer-type-specific scores for user-selected genes"),
                                               downloadButton('downloadDataC6', 'Download', class='btn btn-link'),
                                               HTML("</h4>"),
                                               ## cancer type select
                                               HTML("<p>Select to view detail scores on scecific cancer type: "),
                                               uiOutput("cancerSelectorC6"),
                                               HTML("</p>"),
                                               shiny::tags$hr()   
                                               
                                        )
                                      ),
                                      fluidRow(
                                        column(12,
                                               #uiOutput("cancerSelectorC1Sub"),
                                               shiny::dataTableOutput('genesResTableC6')
                                               )
                                      ),
                                      fluidRow(
                                        column(12,
                                               HTML("<hr>"),
                                               HTML("
                                                          <a class='btn btn-link' data-toggle='collapse' href='#collapseExampleC6Table' aria-expanded='false' aria-controls='collapseExampleC6Table'>
                                                          Click to see table header descriptions! 
                                                          </a>
                                                          <div class='collapse' id='collapseExampleC6Table'>
                                                            <div class='well'>
                                                                <h4>Table column descriptions</h4>
                                                                <ul>
                                                                    <li> <b>Genes</b>: HUGO gene symbol </li>
                                                                    <li> <b>OG Score</b>: Oncogene score for the respective gene</li>
                                                                    <li> <b>TS Score</b>: Tumor suppressor score for the respective gene</li>
                                                                    <li> <b>Combined Score</b>: Combined score, calculated as oncogene score - tumor suppressor score, for the respective gene</li>
                                                                    <li> <b>OG.Meth</b>: Indicates if the gene shows oncogene-like DNA methylation changes</li>
                                                                    <li> <b>OG.CNA</b>: Indicates if the gene shows oncogene-like DNA copy number alterations</li>
                                                                    <li> <b>OG.Mut</b>: Indicates if the gene shows oncogene-like somatic mutation profiles</li>
                                                                    <li> <b>OG.shRNA</b>: Indicates if the gene shows oncogene-like changes after shRNA knock-down</li>
                                                                    <li> <b>OG.Expr</b>: Indicates if the gene shows oncogene-like gene expression changes</li>
                                                                    <li> <b>TS.Meth</b>: Indicates if the gene shows tumor suppressor-like DNA methylation changes</li>
                                                                    <li> <b>TS.CNA</b>: Indicates if the gene shows tumor suppressor-like DNA copy number alterations</li>
                                                                    <li> <b>TS.Mut</b>: Indicates if the gene shows tumor suppressor-like somatic mutation profiles</li>
                                                                    <li> <b>TS.shRNA</b>: Indicates if the gene shows tumor suppressor-like changes after shRNA knock-down</li>
                                                                    <li> <b>TS.Expr</b>: Indicates if the gene shows tumor suppressor-like gene expression changes</li>
                                                                    <li> <b>Cancer</b>: Cancer type showing the association</li>
                                                                    <li> <b>External links</b>: Link to the GeneCards entry for the respective gene</li>
                                                                  </ul>
                                                              </div>
                                                            </div>
                                                            ")
                                        )
                                     )                                      
                                    )
                           )
                    )
                  )
           )
                  )    
           ),


             ## comp2
             tabPanel("Single Gene Boxplots",
                      
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
                                 shiny::actionButton("refreshPlotC2",label="Refresh Results",class='btn btn-primary')
                               )   
                        ),
                        column(10,
                          ## plot window
                          fluidRow(
                            fluidRow(
                              column(12,
                                     HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Refresh Results' to generate plots!</strong>
                                    </div>
                                    ")
                              )
                            ),
                            column(10,
                                   HTML("<div align='left'>
                                      <p><I>The box plots show the log fold change of the gene expression (resp. copy number data) 
                                        of the tumor and normal samples. Differences were assessed with a Wilcoxon test and p-values were corrected 
                                        with the Benjamin-Hochberg procedure to control the false discovery rate (FDR). 
                                        Differences are considered significant when the adjusted p-value is less than 0.05.
                                      </I></p>
                                    </div>
                                    ")
                            )
                            ),
                          fluidRow(
                            column (5, HTML("<h3>Gene Expression"),
                                    downloadButton('downloadPlotC2GE', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    plotOutput("geneExpressionPlot")),
                            column (5, HTML("<h3>Copy Numbers"),
                                    downloadButton('downloadPlotC2CNA', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    plotOutput("cnvPlot"))
                            ),
                          fluidRow(
                            column (10, HTML("<h3>Achilles"),
                                    downloadButton('downloadPlotC2A', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    plotOutput("achillesPlot"))
                          )
                        )
                      )
             ),
             
             ## comp3
             tabPanel("Genomic Regions",
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
                                   shiny::actionButton("refreshPlotC3",label="Refresh Results",class='btn btn-primary')                         
                                 )     
                        ),
                        column(10,
                               fluidRow(
                                 column(12,
                                        HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Refresh Results' to generate plots!</strong>
                                    </div>
                                    ")
                                 )
                               ),
                               
                               fluidRow(
                                 column(10,
                                        HTML("<div align='left'>
                                      <p><I>Two-dimensional overview of the prioritization scores. Genes are plotted column-wise sorted according to 
                                        chromosomal location, starting with the telomere of chromosome 1p in the upper left corner. Lines delineate the 
                                        different chromosomes and numbers indicate chromosome names.
                                      </I></p>
                                    </div>
                                    ")
                                 )
                               ),
                          ## plot window
                          fluidRow(
                            #column (10, 
                                    HTML("<h3>Result Plot"),
                                    downloadButton('downloadDataView1C3', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    imageOutput("selectedScorePlotNew")#)
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
                                 shiny:: actionButton("refreshPlotC4",label="Refresh Results",class='btn btn-primary')                         
                               )
                      ),
                      column(10,
                             fluidRow(
                               fluidRow(
                                 column(12,
                                        HTML("<div align='center' class='alert alert-info alert-dismissible' role='alert'>
                                      <button type='button' class='close' data-dismiss='alert' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
                                        <strong>Click on the 'Refresh Results' to generate plots!</strong>
                                    </div>
                                    ")
                                 ),
                                 column(10,
                                        HTML("<div align='left'>
                                      <p><I>The pathway diagrams are retrieved from <a href='http://www.kegg.jp' target='_blank'>KEGG</a> using the 
                                      			<a href='http://www.bioconductor.org/packages/release/bioc/html/pathview.html', target='_blank'>pathview</a> package. 
                                            Each square represents a gene family and is coloured according to the maximal absolute value of the selected score 
                                            (oncogene, tumor suppressor or combined) for genes in this family. In the overview with all cancer types, each cancer type is shown as a 
                                            stripe in alphabetical order of the cancer type abbreviations. When you select a specific cancer type in the menu, 
                                              you can see the scores of the tumours, the cell lines, or both at the same time (stripes with tumours on the 
                                              left and cell lines on the right).
                                      </I></p>
                                    </div>
                                    ")
                                 )
                               ),
                               
                          ## plot window
                                    HTML("<h3>Pathway Plot"),
                                    downloadButton('downloadPlotC4', 'Download Plot', class='btn btn-link'),
                                    HTML("</h3>"),
                                    imageOutput("pathwayPlot")
                            )
                          
                        )
                      )
             ),
             
             ## comp5
             tabPanel("Summary",
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
                                   shiny::tags$p("Number of genes with respective scores greater or equal to 1. In case of the combined score (calculated as oncogene score (OG) - tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff1C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 2"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with respective scores greater or equal to 2. In case of the combined score (calculated as oncogene score (OG) - tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff2C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 3"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with respective scores greater or equal to 3. In case of the combined score (calculated as oncogene score (OG) - tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff3C5'),
                                   shiny::tags$br(),
                                   shiny::tags$h3("Number of genes with score >= 4"),
                                   shiny::tags$hr(),
                                   shiny::tags$p("Number of genes with respective scores greater or equal to 4. In case of the combined score (calculated as oncogene score (OG) - tumor suppressor gene score (TS)), the absolute value was taken into account."),
                                   shiny::dataTableOutput('genesCutoff4C5'),
                                   shiny::tags$br()
                                   
                            ),
                          column(3)
                          )  
                      ),

            ## comp7
            tabPanel("FAQ",value='faq',
                     fluidRow(
                       column(3),
                       column(6,
                              shiny::tags$h4("Which cancer types were analyzed and what do the cancer codes mean?"),
                              shiny::tags$hr(),
                              HTML("<ul>
                                    <li>BLCA : Bladder Urothelial Carcinoma</li>
                                    <li>BRCA : Breast invasive carcinoma</li>
                                    <li>COAD : Colon adenocarcinoma</li>
                                    <li>GBM  : Glioblastoma multiforme</li>
                                    <li>HNSC : Head and Neck squamous cell carcinoma</li>
                                    <li>KIRC : Kidney renal clear cell carcinoma</li>
                                    <li>LUAD : Lung adenocarcinoma</li>
                                    <li>LUSC : Lung squamous cell carcinoma</li>
                                    <li>OV  :  Ovarian serous cystadenocarcinoma</li>
                                    <li>READ : Rectum adenocarcinoma</li>
                                    <li>UCEC : Uterine Corpus Endometrial Carcinoma</li>
                                    </ul>
                                   "),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("How are the scores computed / what do the scores mean?"),
                              shiny::tags$hr(),
                              shiny::tags$p("Tumor samples are compared with normal samples to identify differences. If a gene is found to be altered, this gene receives a score of 1 for this data type and else a score of 0. Details on the scoring for each data type are given below. Activating and inactivating alterations are both scored independently for each gene, and the sums of the activating and inactivating aberrations yielded an oncogene score and a tumor suppressor gene score, respectively. Additionally, we calculated the difference between oncogene score and tumor suppressor gene score, referred to as overall score, and the sum between oncogene and tumor suppressor scores, referred to as aberration score. Genes were then ranked based on one of these scores to be classified as potential new oncogene or tumor suppressor gene. Pathway alteration scores were calculated by averaging scores for all genes assigned to the same pathway. Aberrations in cancer cell lines were assessed by comparing the cell lines with normal samples available from TCGA using the same approach as for tumor samples."),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("How are gene expression aberrations scored?"),
                              shiny::tags$hr(),
                              shiny::tags$p("Normalized gene expression data for tumor and normal samples, either from an Illumina sequencing platform or Agilent arrays depending on availability for each cancer type, was obtained from TCGA. We compared expression levels for each gene between tumors and matched normal samples using paired Wilcoxon tests and corrected nominal p-values using the Benjamini-Hochberg procedure. If a gene was significantly (FDR < 0.05) expressed (lower or higher) in tumor samples, it received a +1 towards tumor suppressor or oncogene score, respectively."),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("How are copy-number changes scored?"),
                              shiny::tags$hr(),
                              shiny::tags$p("Segmented DNA copy-number data for tumor and normal samples were obtained using Affymetrix SNP6 arrays by TCGA. We compared log-ratio copy number values between tumor samples and matched normal samples using paired Wilcoxon tests and corrected p-values using Benjamini-Hochberg’s procedure. For genes with a significant difference in copy-number (FDR < 0.05), we calculated the Spearman correlation between copy-number and gene expression. A gene was scored as potential tumor suppressor or oncogene if 1) its copy-number value in tumor samples was significantly lower or higher than in normal samples and 2) the copy-number was significantly (FDR < 0.05) positively correlated with gene expression across the tumor samples."),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("How are DNA methylation changes scored?"),
                              shiny::tags$hr(),
                              shiny::tags$p("We analyzed DNA methylation data in a probe-wiseuser fashion. For each probe, we compared methylation values in tumor and normal samples using unpaired Wilcoxon tests and corrected p-values using the Benjamini-Hochberg procedure. For each probe with a significant (FDR < 0.05) difference between methylation in tumor and normal samples, we computed the Spearman correlation between methylation level and expression of the genes annotated to that probe according to the Illumina annotation across all tumor samples. Correlations with FDR < 0.05 were regarded as significant. A gene was scored as tumor suppressor gene if 1) at least one associated methylation probe in the gene body exhibited significantly lower methylation in tumors and 2) the methylation was positively correlated with gene expression, or if at least one other associated methylation probe showed gain of methylation in tumors and gene expression was negatively correlated. In contrast, a gene was scored as oncogene if 1) at least one probe in the gene body showed gain of methylation and 2) was positively correlated with expression or any other associated probe exhibited lower methylation in tumors and methylation was negatively correlated with gene expression."),
                              shiny::tags$br(),
                              
                              
                              shiny::tags$h4("How are somatic mutations scores?"),
                              shiny::tags$hr(),
                              shiny::tags$p("Mutations were scored according to the 20/20 rule published by Vogelstein and colleagues. We divided mutations according to their classification into oncogene mutations (missense mutations and in frame deletion/insertions) and tumor suppressor mutations (frame-shift deletions/insertions, nonsense mutations and splice site mutations). Then, we calculated an oncogene mutation rate (OGMR) and a tumor suppressor mutation rate (TSMR) for each gene. The OGMR was defined as one minus the number of distinct oncogene mutations divided by the total number of mutations. The TSMR was defined as distinct tumor suppressor mutations divided by the number of total mutations. Genes with an OGMR > 0.2 and TSMR < 0.05 were scored as oncogene. A gene was scored as tumor suppressor if its TSMR > 0.2 or if the OGMR > 0.2 and the TSMR > 0.05. We required a minimum of five oncogene or tumor suppressor mutations in order to score a gene as oncogene or tumor suppressor, respectively."),
                              shiny::tags$br(),
                              
                              
                              shiny::tags$h4("How are data from shRNA knock-down screens scored?"),
                              shiny::tags$hr(),
                              shiny::tags$p("Project Achilles assessed cell viability after knocking down genes using different shRNA hairpins. In order to minimize off-target effects, Project Achilles integrated knock-down results of several hairpins targeting the same gene into so called \"gene solutions\", providing a cell viability score for each gene and cell line combination. We used only genes for which only a single gene solution was provided. For each gene, we derived a distribution of viability values across all cell lines. A gene was scored as potential oncogene for one cancer type if at least 25% of the cell lines for this cancer type had a knock-down viability score that was lower than the 25th percentile across all cell lines. If at least 25% of the cell lines had a knock-down viability score greater than the 75th percentile for that gene, it was scored as potential tumor suppressor gene."),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("What do the '-' indicate in the scores table?"),
                              shiny::tags$hr(),
                              shiny::tags$p("The '-' in the scores table can indicate two things. If a row contains valid oncogene/tumor suppressor/combined score
                                            (e.g. 3,-1,0) but some other column is set to '-' then it indicates a missing value in the database for this column.
                                            On the other hand, if a data row is completely filled with '-' for the score columns then this specific gene is 
                                            not found in the database."),
                              shiny::tags$br(),
                              
                              shiny::tags$h4("Why user uploaded gene not shown in the plots?"),
                              shiny::tags$hr(),
                              shiny::tags$p("If user uploads a gene which is not a valid HGNC/HUGO gene symbol, then it is not shown in the plots. However, the 
                                            scores table will show a data row for this gene but all the scores column will be set to '-'. This incident can 
                                            also occour if some user uploaded gene is not found in the database."),
                              shiny::tags$br()
                              
                              

                              ),
                       column(3)
                       )
              )
                      
          )
  
    ))
)