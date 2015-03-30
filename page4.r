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
# generatePathview = function(tcgaResults, ccleResults, pathway, cancers="all",
# 			    what=c("tcga", "ccle", "both"),
# 			    out.dir=".", out.suffix="", kegg.dir=".", 
# 	                    scores=c("combined.score","og.score","ts.score")) {
# 	# Get the score matrix and combine them if necessary
# 	if (what == "tcga") {
# 		if (cancers == "all") {
# 			cancers = names(tcgaResults)
# 		}
# 		scoreMat = pathviewMat(tcgaResults[intersect(cancers, names(tcgaResults))], scores[1])
# 	} else if (what == "ccle") {
# 		if (cancers == "all") {
# 			cancers = names(ccleResults)
# 		}
# 		scoreMat = pathviewMat(ccleResults[intersect(cancers, names(ccleResults))], scores[1])
# 	} else {
# 		scoreMat = cbind(tcgaResults[[cancers[1]]]$prioritize.combined[, scores],
# 				 ccleResults[[cancers[1]]]$prioritize.combined[, scores])
# 	}
# 	
# 	if (scores == "combined.score") {
# 		low = list(gene="#034b87", cpd="blue")
# 		mid = list(gene="gray98", cpd="gray")
# 		high = list(gene="#880000", cpd="yellow")
# 		both.dirs = list(gene=TRUE, cpd=TRUE)
# 	} else if (scores == "og.score") {
# 		low = list(gene="gray98", cpd="blue")
# 		mid = list(gene="gray98", cpd="gray")
# 		high = list(gene="#880000", cpd="yellow")
# 		both.dirs = list(gene=FALSE, cpd=FALSE)
# 	} else {
# 		low = list(gene="gray98", cpd="blue")
# 		mid = list(gene="gray98", cpd="gray")
# 		high = list(gene="#034b87", cpd="yellow")
# 		both.dirs = list(gene=FALSE, cpd=FALSE)
# 	}
# 	
# 	limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))
# 	
# 	# Save working directory and switch to new one
# 	oldwd = getwd()
# 	setwd(out.dir)
# 	
# 	# Generate plots
# 	invisible(pathview(gene.data=scoreMat, pathway.id=pathway,
# 		           kegg.native=TRUE, gene.idtype="SYMBOL",
# 		           out.suffix=out.suffix, kegg.dir=kegg.dir,
# 		           limit=list(gene=limit, cpd=1), node.sum="max.abs",
# 		           multi.state=TRUE, low=low, mid=mid, high=high,
# 			   both.dirs=both.dirs))
# 	# Set old working directory
# 	setwd(oldwd)
# 	
# 	paste0(getwd(), "/hsa", pathway, ".", out.suffix, ".multi.png")
# 	
# }


generatePathview2 = function(updateProgress = NULL,pathway, cancer,sample,scores) {
  
  if (substr(pathway,1,3) == "hsa")
  {
    ## get the pathway code out of pathway parameter
    temp <- strsplit(pathway,":")
    code <- temp[[1]][1]
    pathwayCode <- substr(code,4,nchar(code))
    scoreMat = NULL
    what = sample
    low = NULL
    mid = NULL
    high = NULL
    both.dirs = NULL
    limit = NULL
    
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
    
    
    # Get the score matrix and combine them if necessary
    if (what == "tcga") {
      scoreMat <- tcgaScoreMat[[cancer]][[scores]]
      limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))  
    } else if (what == "ccle") {
      scoreMat <- ccleScoreMat[[cancer]][[scores]]
      limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))  
    } else {
      if (nrow(tcgaScoreMat[[cancer]][[scores]]) == nrow(ccleScoreMat[[cancer]][[scores]]))
      {
        scoreMat = cbind(tcgaScoreMat[[cancer]][[scores]],ccleScoreMat[[cancer]][[scores]])
        limit = c(min(scoreMat[,c(2,3)], na.rm=TRUE), max(scoreMat[,c(2,3)], na.rm=TRUE))
      }else{
        tmat <- data.frame(tcgaScoreMat[[cancer]][[scores]],rownames(tcgaScoreMat[[cancer]][[scores]]))
        colnames(tmat) <- c(paste("tumor.",cancer,sep=""),"gene")
        cmat <- data.frame(ccleScoreMat[[cancer]][[scores]],rownames(ccleScoreMat[[cancer]][[scores]]))
        colnames(cmat) <- c(paste("cellline.",cancer,sep=""),"gene")
        scoreMat <- plyr::join(tmat,cmat,type="inner",by="gene")
        rownames(scoreMat) <- scoreMat[,2]
        scoreMat <- scoreMat[,c(1,3)]
        limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))
        rm(tmat)
        rm(cmat)
      }
    }
    
    #limit = c(min(scoreMat, na.rm=TRUE), max(scoreMat, na.rm=TRUE))  
    #limit = c(min(scoreMat[,c(2,3)], na.rm=TRUE), max(scoreMat[,c(2,3)], na.rm=TRUE))
    
    # Save working directory and switch to new one
    #oldwd = getwd()
    #setwd(out.dir)
    
    kdir = paste(getwd(),"/keggdir/",sep="")
    
    # Generate plots
    invisible(pathview(gene.data=scoreMat, pathway.id=pathwayCode,
                       kegg.native=TRUE, gene.idtype="SYMBOL",
                       out.suffix="", kegg.dir=kdir,
                       limit=list(gene=limit, cpd=1), node.sum="max.abs",
                       multi.state=TRUE, low=low, mid=mid, high=high,
                       both.dirs=both.dirs))
    # Set old working directory
    #setwd(oldwd)
    p <- NULL
    if (cancer == "All" || sample == "both"){ 
      p <- paste(getwd(),"/hsa", pathwayCode, "..multi.png",sep="") 
    }else{
      p <- paste(getwd(),"/hsa", pathwayCode, "..png",sep="")
    }
    p
  }else{
    p <- paste(getwd(),"pathwaydefault.png",sep="")
    p
  }

}
