source("plots/plotting.R")
source("page1.r")

rm (list=ls())

# Load the TCGA results
load("www/prioritize_tcga_pancancer_allgenes_step2.rdata")
tcgaResults = results
# Load the CCLE results
load("www/prioritize_ccle_pancancer_allgenes_step2.rdata")
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
save(cancers,geness,tcgaResultsHeatmapOG,tcgaResultsHeatmapTS,tcgaResultsHeatmapCombined,ccleResultsHeatmapOG,ccleResultsHeatmapTS,ccleResultsHeatmapCombined,file="www/starter.RData")


## test *******************************************************

## subset data frame based on user input
tcgaOG <- page1DataFrame(tcgaResultsHeatmapOG, '2', 'BLCA')
## let's plot
stype=c("combined.score", "ts.score", "og.score")
plotHeatmapPage1(tcgaOG, stype[1])
plotHeatmapPage1(tcgaOG, stype[2])
plotHeatmapPage1(tcgaOG, stype[3])