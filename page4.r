##' Generates pathview plots for all given results. This function can be used to plot 
##' different scores for one cancer type or one score across different cancer types.
##' If cancers == "all", only the first score in scores is plotted. In all other cases,
##' all scores are plotted for the first cancer only. 
##' @param tcgaResults the result list from a TCGA prioritization run
##' @param ccleResults the result list from a CCLE prioritization run
##' @param pathways vector with KEGG pathway IDs to plot; default: NULL (all pathways)
##' @param cancers vector of names of cancer types to plot; default: all
##' @param scores vector of names of the scores to plot; default: combined.score
##' @param out.dir directory for output files; default: "."
##' @param out.suffix suffix to be added to output plots; default: ""
##' @param kegg.dir directory with predownloaded KEGG files; all new downloaded files will be stored there;
##' default: "."
##' @author Andreas Schlicker
generatePathview = function(tcgaResults, ccleResults, pathway, cancers="all",
			    what=c("tcga", "ccle", "both"),
			    out.dir=".", out.suffix="", kegg.dir=".", 
	                    scores="combined.score") {
	# Get the score matrix and combine them if necessary
	if (what == "tcga") {
		if (cancers == "all") {
			cancers = names(tcgaResults)
		}
		scoreMat = pathviewMat(tcgaResults[intersect(cancers, names(tcgaResults))], scores[1])
	} else if (what == "ccle") {
		if (cancers == "all") {
			cancers = names(ccleResults)
		}
		scoreMat = pathviewMat(ccleResults[intersect(cancers, names(ccleResults))], scores[1])
	} else {
		scoreMat = cbind(tcgaResults[[cancers[1]]]$prioritize.combined[, scores],
				 ccleResults[[cancers[1]]]$prioritize.combined[, scores])
	}
	
	if (scores == "combined.score") {
		low = list(gene="#034b87", cpd="blue")
		mid = list(gene="gray98", cpd="gray")
		high = list(gene="#880000", cpd="yellow")
		both.dirs = list(gene=TRUE, cpd=TRUE)
	} else if (scores == "og.score") {
		low = list(gene="gray98", cpd="blue")
		mid = list(gene="gray98", cpd="gray")
		high = list(gene="#880000", cpd="yellow")
		both.dirs = list(gene=FALSE, cpd=FALSE)
	} else {
		low = list(gene="gray98", cpd="blue")
		mid = list(gene="gray98", cpd="gray")
		high = list(gene="#034b87", cpd="yellow")
		both.dirs = list(gene=FALSE, cpd=FALSE)
	}
	
	limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))
	
	# Save working directory and switch to new one
	oldwd = getwd()
	setwd(out.dir)
	
	# Generate plots
	invisible(pathview(gene.data=scoreMat, pathway.id=pathway,
		           kegg.native=TRUE, gene.idtype="SYMBOL",
		           out.suffix=out.suffix, kegg.dir=kegg.dir,
		           limit=list(gene=limit, cpd=1), node.sum="max.abs",
		           multi.state=TRUE, low=low, mid=mid, high=high,
			   both.dirs=both.dirs))
	# Set old working directory
	setwd(oldwd)
	
	paste0(getwd(), "/hsa", pathway, ".", out.suffix, ".multi.png")
}
