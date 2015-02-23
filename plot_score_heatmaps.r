##' Generates a chromosome plot for the selected score or the number of affected samples.
##' @params geneLocs the data.frame returned by sortGenesByLocation()
##' @params results either tcgaResults or ccleResults
##' @params cancType the selected four letter cancer code
##' @params scoreType which score was selected; one of "OG", "TS", "CO"
##' @return a list with two plotting objects, "score" the actual score and "affected" the percentage of affected samples
##' @author Andreas Schlicker (a.schlicker@nki.nl) 
getScoreHeatmap = function(geneLocs, results, cancType, scoreType) {
	scoreType = match.arg(scoreType, c("TS", "OG", "CO"))
	
	# Get the score vector
	score = results[[cancType]]$prioritize.combined[geneLocs[, "hgnc_symbol"], "og.score"]
	if (scoreType == "TS") {
		score = results[[cancType]]$prioritize.combined[geneLocs[, "hgnc_symbol"], "ts.score"]
		cols = colorpanel(29, low="white", high="#034b87")
		brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
	} else if (scoreType == "OG") {
		score = results[[cancType]]$prioritize.combined[geneLocs[, "hgnc_symbol"], "og.score"]
		cols = colorpanel(29, low="white", high="#880000")
		brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
	} else {
		score = results[[cancType]]$prioritize.combined[geneLocs[, "hgnc_symbol"], "combined.score"]
		cols = colorpanel(29, low="#034b87", mid="white", high="#880000")
		ext = max(abs(min(sortedScores, na.rm=TRUE)), max(sortedScores, na.rm=TRUE))
		#brks=seq(-1*ext, ext, length.out=30)
		brks = 0
	}
	
	matDim = ceil(sqrt(length(score)))
	sortedScores = matrix(c(score, rep(0, times=(matDim^2)-length(sortedGenes))), nrow=matDim, byrow=FALSE)[, 1:(matDim-1)]
	
	aheatmap(sortedScores,
		 color=cols,
		 breaks=brks,
		 scale="none",
		 Rowv=NA, Colv=NA, 
		 cexRow=0, cexCol=0)
	
	#system(paste("composite ", paste(n, "_", sc, "_heatmap.png", sep=""), " -compose Multiply chrom_contours.png ", paste(n, "_", sc, "_heatmap.png", sep=""), sep=""))
}



getScoreHeatmap2 = function(geneLocs, results, cancType, scoreType) {
  scoreType = match.arg(scoreType, c("TS", "OG", "CO"))
  
  # Get the score vector
  score = results[[cancType]][geneLocs[, "hgnc_symbol"], "og.score"]
  if (scoreType == "TS") {
    score = results[[cancType]][geneLocs[, "hgnc_symbol"], "ts.score"]
    cols = colorpanel(29, low="white", high="#034b87")
    brks=seq(0, max(score, na.rm=TRUE), length.out=30)
  } else if (scoreType == "OG") {
    score = results[[cancType]][geneLocs[, "hgnc_symbol"], "og.score"]
    cols = colorpanel(29, low="white", high="#880000")
    brks=seq(0, max(score, na.rm=TRUE), length.out=30)
  } else {
    score = results[[cancType]][geneLocs[, "hgnc_symbol"], "combined.score"]
    cols = colorpanel(29, low="#034b87", mid="white", high="#880000")
    ext = max(abs(min(score, na.rm=TRUE)), max(score, na.rm=TRUE))
    #brks=seq(-1*ext, ext, length.out=30)
    brks = 0
  }
  
  matDim = ceiling(sqrt(length(score)))
  sortedScores = matrix(c(score, rep(0, times=(matDim^2)-length(score))), nrow=matDim, byrow=FALSE)[, 1:(matDim-1)]
  
  aheatmap(sortedScores,
           color=cols,
           breaks=brks,
           scale="none",
           Rowv=NA, Colv=NA, 
           cexRow=0, cexCol=0)
  
  #system(paste("composite ", paste(n, "_", sc, "_heatmap.png", sep=""), " -compose Multiply chrom_contours.png ", paste(n, "_", sc, "_heatmap.png", sep=""), sep=""))
}
