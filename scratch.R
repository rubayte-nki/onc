rm (list=ls())

source("plots/plotting.R")
source("page1.r")
source("page2.r")

## page 1
## **********************************************************************************

dir=setwd("path/to")

# Load the TCGA results
load("prioritize_tcga_pancancer_allgenes_step2.rdata")
tcgaResults = results
# Load the CCLE results
load("prioritize_ccle_pancancer_allgenes_step2.rdata")
ccleResults = results
# Save some memory
rm(results)

#############################################################################
## prepare data frames
#############################################################################
# TCGA oncogene (OG) scores
tcgaResultsHeatmapOG = heatmapDataframe(tcgaResults, scores=list(combined="og.score", Meth="og.methylation", CNA="og.cna", Mut="og.mutations",shRNA="og.achilles", Expr="og.exprs"))
# TCGA tumor suppressor (TS) scores
tcgaResultsHeatmapTS = heatmapDataframe(tcgaResults, scores=list(combined="ts.score", Meth="ts.methylation",CNA="ts.cna", Mut="ts.mutations",shRNA="ts.achilles", Expr="ts.exprs"))
# TCGA combined score
tcgaResultsHeatmapCombined = heatmapDataframe(tcgaResults) 
# CCLE OG scores
ccleResultsHeatmapOG = heatmapDataframe(ccleResults,scores=list(combined="og.score", Meth="og.methylation",CNA="og.cna", Mut="og.mutations",shRNA="og.achilles", Expr="og.exprs"))
# CCLE OG scores
ccleResultsHeatmapTS = heatmapDataframe(ccleResults,scores=list(combined="ts.score", Meth="ts.methylation",CNA="ts.cna", Mut="ts.mutations",shRNA="ts.achilles", Expr="ts.exprs"))
# CCLE combined score
ccleResultsHeatmapCombined = heatmapDataframe(ccleResults)

#############################################################################
## update started data frame
#############################################################################
load("www/starter.RData")
save(cancers,geness,file="www/starterWidgets.RData")
save(tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,file="www/starter.RData")


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
load("tcga_pancancer4_meth.rdata")
# Methylation annotation data
load("illumina_infinium450_annotation.rdata")
# Project Achilles data
load("achilles.rdata")


#############################################################################
## prepare data frames
#############################################################################
#achilles,ccle.cna.combat,ccle.exprs.combat,infinium450.probe.ann,tcga.cna,tcga.cna.combat,tcga.exprs,tcga.exprs.combat,tcga.meth,
#tcga.mn.cna,tcga.mn.cna.combat,tcga.mn.exprs,tcga.mn.exprs.combat,tcga.mn.meth

#############################################################################
## update started data frame
#############################################################################
load("www/starter.RData")
save(cancers,geness,tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,
     achilles,ccle.cna.combat,ccle.exprs.combat,infinium450.probe.ann,tcga.cna,tcga.cna.combat,tcga.exprs,tcga.exprs.combat,tcga.meth,
     tcga.mn.cna,tcga.mn.cna.combat,tcga.mn.exprs,tcga.mn.exprs.combat,tcga.mn.meth,tcgaResults,ccleResults, file="www/starter.RData")



## test *******************************************************

## let's plot
getPage2Plots("BLCA", "ATM", 1)



