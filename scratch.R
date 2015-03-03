rm (list=ls())

setwd("C:/Rubayte/Oncoscape/dev3")
source("plotting.R")
source("page1.r")
source("page2.r")
source("page3.r")
source("page4.r")

## page 1

setwd("C:/Rubayte/Oncoscape/data")
# Load the TCGA results
load("prioritize_tcga_pancancer_allgenes_step2.rdata")
tcgaResults = results
# Load the CCLE results
load("prioritize_ccle_pancancer_allgenes_step2.rdata")
ccleResults = results
# Save some memory
rm(results)
setwd("C:/Rubayte/Oncoscape/dev3")



#############################################################################
## prepare data frames
#############################################################################
# TCGA oncogene (OG) scores
tcgaResultsHeatmapOG = heatmapDataframe(tcgaResults, scores=list(Combined="og.score", Meth="og.methylation", CNA="og.cna", Mut="og.mutations",shRNA="og.achilles", Expr="og.exprs"))
# TCGA tumor suppressor (TS) scores
tcgaResultsHeatmapTS = heatmapDataframe(tcgaResults, scores=list(Combined="ts.score", Meth="ts.methylation",CNA="ts.cna", Mut="ts.mutations",shRNA="ts.achilles", Expr="ts.exprs"))
# TCGA combined score
tcgaResultsHeatmapCombined = heatmapDataframe(tcgaResults) 
# CCLE OG scores
ccleResultsHeatmapOG = heatmapDataframe(ccleResults,scores=list(Combined="og.score", Meth="og.methylation",CNA="og.cna", Mut="og.mutations",shRNA="og.achilles", Expr="og.exprs"))
# CCLE OG scores
ccleResultsHeatmapTS = heatmapDataframe(ccleResults,scores=list(Combined="ts.score", Meth="ts.methylation",CNA="ts.cna", Mut="ts.mutations",shRNA="ts.achilles", Expr="ts.exprs"))
# CCLE combined score
ccleResultsHeatmapCombined = heatmapDataframe(ccleResults)

#############################################################################
## update started data frame
#############################################################################
load("starter.RData")
save(cancers,geness,file="starterWidgets.RData")
save(tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,
     ccleResultsHeatmapTS,ccleResultsHeatmapCombined,file="starter.RData")

## add summary statistics data frames
load("starterWidgets.RData")
save(cancers,geness,sampleOverview,genesCutoffOne,genesCutoffTwo,genesCutoffThree,genesCutoffFour,file="starterWidgets.RData")



## test *******************************************************

## subset data frame based on user input
tcgaOG <- page1DataFrame(tcgaResultsHeatmapOG, '2', 'BLCA')
## let's plot
stype=c("combined.score", "ts.score", "og.score")
plotHeatmapPage1(tcgaOG, stype[1])
plotHeatmapPage1(tcgaOG, stype[2])
plotHeatmapPage1(tcgaOG, stype[3])


## page 2
## **********************************************************************************
setwd("C:/Rubayte/Oncoscape/data")

# Load the TCGA data
# Copy number
load("tcga_pancancer4_cna.rdata")
load("tcga_pancancer4_cna_ccle.rdata")
rm(ccle.cna)
# Gene expression
load("tcga_pancancer4_exprs.rdata")
load("tcga_pancancer4_exprs_ccle.rdata")
rm(ccle.exprs)
# DNA methylation
#load("tcga_pancancer4_meth.rdata")
# Methylation annotation data
#load("illumina_infinium450_annotation.rdata")
# Project Achilles data
load("achilles.rdata")

setwd("C:/Rubayte/Oncoscape/dev3")

 

#############################################################################
## prepare data frames
#############################################################################
ccleResultsPrioDetails = NULL
tcgaResultsPrioDetails = NULL
ccleResultsAT = NULL
tcgaResultsAT = NULL

for (n in names(ccleResults))
{
  ccleResultsPrioDetails[[n]] = ccleResults[[n]]$prioritize.details
}

for (n in names(ccleResults))
{
  ccleResultsAT[[n]] = ccleResults[[n]]$achilles.thresholds
}

for (n in names(tcgaResults))
{
  tcgaResultsPrioDetails[[n]] =  tcgaResults[[n]]$prioritize.details
}

for (n in names(tcgaResults))
{
  tcgaResultsAT[[n]] =  tcgaResults[[n]]$achilles.thresholds
}


ccleResultsCLS = NULL
tcgaResultsCLS = NULL

for (n in names(tcgaResults))
{
  tcgaResultsCLS[[n]]$cls =  tcgaResults[[n]]$cls
}

for (n in names(ccleResults))
{
  ccleResultsCLS[[n]]$cls = ccleResults[[n]]$cls
}

save(ccleResultsCLS,tcgaResultsCLS, file="C:/Rubayte/Oncoscape/data/clsdata.RData")



#############################################################################
## update started data frame
#############################################################################
load("starter.RData")

save(tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,
     tcgaResultsPlotTrack,ccleResultsPlotTrack,tcgaScoreMat,ccleScoreMat,achilles,ccle.cna.combat,ccle.exprs.combat,tcga.cna,tcga.cna.combat,tcga.exprs,
     tcga.exprs.combat,tcga.mn.cna,tcga.mn.cna.combat,tcga.mn.exprs,tcga.mn.exprs.combat,ccleResultsPrioDetails,
     tcgaResultsPrioDetails,ccleResultsAT, tcgaResultsAT, file="C:/Rubayte/Oncoscape/data/starterv2.RData")



#############################################################################
## update started data frame
#############################################################################
load("starter.RData")

ccleResultsPlotTrack = NULL
tcgaResultsPlotTrack = NULL

for (n in names(ccleResults))
{
  ccleResultsPlotTrack[[n]] = ccleResults[[n]]$prioritize.combined
}

for (n in names(tcgaResults))
{
  tcgaResultsPlotTrack[[n]] =  tcgaResults[[n]]$prioritize.combined
}

save(tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,
     tcgaResultsPlotTrack,ccleResultsPlotTrack, file="starter.RData")



## test *******************************************************

## let's plot
getPage2Plots("BLCA", "ATM", 1)


## page 3
## **********************************************************************************


#############################################################################
## prepare data frames
#############################################################################

gloc <- sortGenesByLocation(tcgaResults, ccleResults)

## let's plot
res <- getScorePlot(gloc, tcgaResults, cancers$V1[1], '2', 'TS')
  


## page 4
## **********************************************************************************

#############################################################################
## prepare data frames
#############################################################################


cancers = c(names(tcgaResults),"All")

## tcga
tcgaScoreMat = NULL
for(c in cancers)
{
  if (c == "All")
  {
    temp = names(tcgaResults)
    tcgaScoreMat[[c]] = pathviewMat(tcgaResults[intersect(temp, names(tcgaResults))], "combined.score")
  }else{
    tcgaScoreMat[[c]] = pathviewMat(tcgaResults[intersect(c, names(tcgaResults))], "combined.score")
  }
}

cancers = c(names(ccleResults),"All")

## ccle
ccleScoreMat = NULL
for(c in cancers)
{
  if (c == "All")
  {
    temp = names(ccleResults)
    ccleScoreMat[[c]] = pathviewMat(ccleResults[intersect(temp, names(ccleResults))], "combined.score")
  }else{
    ccleScoreMat[[c]] = pathviewMat(ccleResults[intersect(c, names(ccleResults))], "combined.score")
  }
}



#############################################################################
## update started data frame
#############################################################################
load("starter.RData")
save(tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,
     tcgaResultsPlotTrack,ccleResultsPlotTrack,tcgaScoreMat,ccleScoreMat, file="starter.RData")


## test *******************************************************

## let's plot
res <- generatePathview(tcgaResults, ccleResults, '04910', cancers="all",
                            what="tcga",
                            out.dir=".", out.suffix="", kegg.dir=".", 
                            scores="combined.score")




##############################################################################
## update RData object
##############################################################################
load("starterv2.RData")





