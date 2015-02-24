load("gloc.RData")

##' Fetches the chromosomal locations for all genes and sorts them within each chromosome. 
##' @params tcgaResults the list with TCGA results
##' @params ccleResults the list with CCLE results
##' @return a data.frame with locations for all genes
##' @author Andreas Schlicker (a.schlicker@nki.nl)
sortGenesByLocation = function(tcgaResults, ccleResults) {
	# Genes in all cancer types
	genes = rownames(tcgaResults[[1]]$prioritize.combined)
	for (n in 2:length(tcgaResults)) {
		genes = intersect(genes, rownames(tcgaResults[[n]]$prioritize.combined))
	}
	for (n in 1:length(ccleResults)) {
		genes = intersect(genes, rownames(ccleResults[[n]]$prioritize.combined))
	}

	# Find all gene locations
	mart = useMart(host="ensembl.org", path="/biomart/martservice", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
	geneLoc = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=c("hgnc_symbol"), values=genes, mart=mart)
	# Remove all non-standard chromosome names and make everything numeric
	geneLoc = geneLoc[which(geneLoc[, "chromosome_name"] %in% c("1", "2", "3", "4", "5", "6", "7", "8", 
																															"9", "10", "11", "12", "13", "14", "15", 
																															"16", "17", "18", "19", "20", "21", "22", 
																															"X", "Y")), ]
	geneLoc[which(geneLoc[, "chromosome_name"] == "X"), "chromosome_name"] = "23"
	geneLoc[which(geneLoc[, "chromosome_name"] == "Y"), "chromosome_name"] = "24"
	geneLoc[, "chromosome_name"] = as.integer(geneLoc[, "chromosome_name"])
	
	geneLoc[which(geneLoc[, "strand"] == 1), "strand"] = "+"
	geneLoc[which(geneLoc[, "strand"] == -1), "strand"] = "-"
	
	# Return the sorted gene information
	geneLoc[order(geneLoc$chromosome_name, geneLoc$start_position), ]
}


##' Generates a chromosome plot for the selected score or the number of affected samples.
##' @params geneLocs the data.frame returned by sortGenesByLocation()
##' @params results either tcgaResults or ccleResults
##' @params cancType the selected four letter cancer code
##' @params chrom the selected chromosome; X is coded as 23 and Y as 24
##' @params scoreType which score was selected; one of "OG", "TS", "CO"
##' @return a list with two plotting objects, "score" the actual score and "affected" the percentage of affected samples
##' @author Andreas Schlicker (a.schlicker@nki.nl) 
getScorePlot = function(geneLocs, results, cancType, chrom, scoreType) {
	scoreType = match.arg(scoreType, c("TS", "OG", "CO"))
	
	# Get the genes for the selected chromosome
	genes = subset(geneLocs, chromosome_name==chrom)
	
	# Get the score vector
	score = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "og.score"]
	# And the percentage of affected samples
	affected = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "og.affected.rel"]
	if (scoreType == "TS") {
		score = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "ts.score"] * -1
		affected = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "ts.affected.rel"]
	} else {
		score = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "combined.score"]
		affected = results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "og.affected.rel"] - 
							 results[[cancType]]$prioritize.combined[genes[, "hgnc_symbol"], "ts.affected.rel"]
	}
	
	# Plot two tracks, one for scores and the other one for percent affected samples 
	list(score=plotTracks(DataTrack(GRanges(seqnames=as.character(chrom),
						ranges=IRanges(start=genes[, "start_position"],
						       	       end=genes[, "end_position"]),
						score,
					        strand="*")), 
			      type=c("h", "g"),
			      ylim=c(min(score), max(score)),
			      col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
			      main=""),
	     affected=plotTracks(DataTrack(GRanges(seqnames=as.character(chrom),
		      	      			   ranges=IRanges(start=genes[, "start_position"],
		      	      			   		  end=genes[, "end_position"]),
		      	      			   affected,
						   strand="*")), 
				 type=c("h", "g"),
				 ylim=c(min(affected), max(affected)),
				 col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
				 main=""))
}

## function for plotting the score data track across selected chromosome and cancer
getScorePlot2 = function(geneLocs, results, cancType, chrom, scoreType) {
  scoreType = match.arg(scoreType, c("TS", "OG", "CO"))
  
  # Get the genes for the selected chromosome
  genes = subset(geneLocs, chromosome_name==chrom)
  
  # Get the score vector and the percentage of affected samples
  if (scoreType == "OG"){
    score = results[[cancType]][genes[, "hgnc_symbol"], "og.score"]
    #affected = results[[cancType]][genes[, "hgnc_symbol"], "og.affected.rel"]    
  }else if (scoreType == "TS") {
    score = results[[cancType]][genes[, "hgnc_symbol"], "ts.score"] * -1
    #affected = results[[cancType]][genes[, "hgnc_symbol"], "ts.affected.rel"]
  } else {
    score = results[[cancType]][genes[, "hgnc_symbol"], "combined.score"] 
    #affected = results[[cancType]][genes[, "hgnc_symbol"], "og.affected.rel"] - 
      #results[[cancType]][genes[, "hgnc_symbol"], "ts.affected.rel"]
  }
  
  # Plot two tracks, one for scores and the other one for percent affected samples 
  plotTracks(DataTrack(GRanges(seqnames=as.character(chrom),
                                          ranges=IRanges(start=genes[, "start_position"],
                                                         end=genes[, "end_position"]),
                                          score,
                                          strand="*")), 
                        type=c("h", "g"),
                        ylim=c(min(score), max(score)),
                        col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
                        main="")
}

## function for plotting the percentage of affected samples data track across selected chromosome and cancer
getAffectedPlot2 = function(geneLocs, results, cancType, chrom, scoreType) {
  scoreType = match.arg(scoreType, c("TS", "OG", "CO"))
  
  # Get the genes for the selected chromosome
  genes = subset(geneLocs, chromosome_name==chrom)
  
  # Get the score vector and the percentage of affected samples
  if (scoreType == "OG"){
    #score = results[[cancType]][genes[, "hgnc_symbol"], "og.score"]
    affected = results[[cancType]][genes[, "hgnc_symbol"], "og.affected.rel"]    
  }else if (scoreType == "TS") {
    #score = results[[cancType]][genes[, "hgnc_symbol"], "ts.score"] * -1
    affected = results[[cancType]][genes[, "hgnc_symbol"], "ts.affected.rel"]
  } else {
    #score = results[[cancType]][genes[, "hgnc_symbol"], "combined.score"]
    affected = results[[cancType]][genes[, "hgnc_symbol"], "og.affected.rel"] - 
      results[[cancType]][genes[, "hgnc_symbol"], "ts.affected.rel"]
  }
  
  # Plot two tracks, one for scores and the other one for percent affected samples 
  plotTracks(DataTrack(GRanges(seqnames=as.character(chrom),
                                             ranges=IRanges(start=genes[, "start_position"],
                                                            end=genes[, "end_position"]),
                                             affected,
                                             strand="*")), 
                           type=c("h", "g"),
                           ylim=c(min(affected), max(affected)),
                           col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
                           main="")
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
  
  outputImage = "temp.png"
  aheatmap(sortedScores,
           color=cols,
           breaks=brks,
           scale="none",
           Rowv=NA, Colv=NA,width=13.35,height=10, 
           cexRow=0, cexCol=0,filename=outputImage)

  system('composite "temp.png" -compose Multiply "chrom_contours.png"  "temp.png" ')
  #system('composite -compose Multiply "chrom_contours.png" "temp.png" "temp.png" ')
  outputImage
  
}


## main call to plot function
## view 1
comp3view1Plot <- function (updateProgress = NULL,cancer,scoreType,chr)
{
  if (scoreType == "TS")
  {
    res <- getScorePlot2(gloc, tcgaResultsPlotTrack, cancer, chr, 'TS')
  } else if (scoreType == "OG") {
    res <- getScorePlot2(gloc, ccleResultsPlotTrack, cancer, chr, 'OG')
  } else {
    res <- getScorePlot2(gloc, ccleResultsPlotTrack, cancer, chr, 'CO')
  }  
  res
}

## view 2
comp3view2Plot <- function (updateProgress = NULL,cancer,scoreType,chr)
{
  if (scoreType == "TS")
  {
    res <- getAffectedPlot2(gloc, tcgaResultsPlotTrack, cancer, chr, 'TS')
  } else if (scoreType == "OG") {
    res <- getAffectedPlot2(gloc, ccleResultsPlotTrack, cancer, chr, 'OG')
  } else  {
    res <- getAffectedPlot2(gloc, ccleResultsPlotTrack, cancer, chr, 'CO')
  }
  res
}

## view new
comp3view1PlotNew <- function (updateProgress,cancer,score,sample)
{
  res <- NULL
  if (sample == "tcga")
  {
      res <- getScoreHeatmap2(gloc, tcgaResultsPlotTrack, cancer, score)      

  }else{
    res <- getScoreHeatmap2(gloc, ccleResultsPlotTrack, cancer, score) 
  }
  res  
}

