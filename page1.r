load("starter.RData")

##' Filtering of result frame according to user criteria
##' @param results data.frame with all results
##' @param scoreCutoff threshold that was selected by the user
##' @param cancerType cancer type selected by the user
##' @return a subset of the data.frame that fits the user's selection
##' @author Andreas Schlicker
page1DataFrame = function(results, scoreCutoff, cancerType) {
	# Filter the genes according to user's criteria
	genes = as.character(subset(results, cancer == cancerType & score.type == "combined" & score >= as.integer(scoreCutoff))$gene)
	
	# Sort the genes according to highest sum across all cancer types
	# Get the subset with the selected genes and drop unused levels
	# gene.order = subset(result.df, score.type=="combined" & gene %in% genes)
	gene.order = subset(results, score.type=="combined" & gene %in% genes)
  gene.order$gene = droplevels(gene.order$gene)
	# Do the sorting
	gene.order = names(sort(unlist(lapply(split(gene.order$score, gene.order$gene), sum, na.rm=TRUE))))
	
	# Get the data.frame for plotting
	result.df = subset(results, gene %in% genes)
	result.df$gene = factor(result.df$gene, levels=gene.order)
	result.df$cancer = factor(result.df$cancer, levels=sort(unique(as.character(result.df$cancer))))
	
	result.df
}

##' Get the heatmap for view 1 of page 1.
##' @params results a subsetted data.frame as returned by page1DataFrame()
##' @params colorLow "#034b87" if selected score was TS or combined, else "gray98"
##' @params colorHigh "#880000" if selected score was OG or combined, else "gray98"
##' @return the heatmap object
##' @author Andreas Schlicker
plotHeatmapPage1 = function(results, scoreType=c("combined.score", "ts.score", "og.score")) {
	result.df = results
	colorLow = list(combined.score="#034b87", ts.score="gray98", og.score="gray98") 
	colorMid = list(combined.score="gray98")#,
	colorHigh = list(combined.score="#880000", ts.score="#034b87", og.score="#880000")
	getHeatmap(dataFrame=result.df, yaxis.theme=theme(axis.text.y=element_blank()), 
	   	   color.low=colorLow[[scoreType]], color.mid=colorMid[[scoreType]], color.high=colorHigh[[scoreType]])
}

##' Plots view 2 of page 1
##' @param results a subsetted data.frame as returned by page1DataFrame()
##' @return the ggplot2 object with the plot for view 2 of page 1
##' @author Andreas Schlicker
plotCategoryOverview = function(results) {
	result.df = results
	result.df$score.type = factor(result.df$score.type, levels=c("CNA", "Expr", "Meth", "Mut", "shRNA", "combined"))
	
	# Overwrite the score column with the score type to make it categorical
	# Combined scores are not plotted later
	result.df[, 2] = as.character(result.df[, 2])
	result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 2] = as.character(result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 3])
	result.df[which(is.na(result.df[, 2]) | result.df[, 2] == "0"), 2] = "NONE"
	
	#ggplot(subset(result.df, score.type != "combined" & gene %in% topgenes), aes(x=score.type, y=gene)) + 
	ggplot(subset(result.df, score.type != "combined"), aes(x=score.type, y=gene)) + 
  geom_tile(aes(fill=score), color="white", size=0.7) +
	scale_fill_manual(values=c(NONE="white", CNA="#888888", Expr="#E69F00", Meth="#56B4E9", Mut="#009E73", shRNA="#F0E442"), 
		          breaks=c("CNA", "Expr", "Meth", "Mut", "shRNA")) +
	labs(x="", y="") +
	facet_grid(.~cancer) + 
	theme(panel.background=element_rect(color="white", fill="white"),
	      panel.margin=unit(10, "points"),
	      axis.ticks=element_blank(),
	      axis.text.x=element_blank(),
	      axis.text.y=element_text(color="gray30", size=10, face="bold"),
	      axis.title.x=element_text(color="gray30", size=10, face="bold"),
	      strip.text.x=element_text(color="gray30", size=10, face="bold"),
	      legend.text=element_text(color="gray30", size=10, face="bold"),
	      legend.title=element_blank(),
	      legend.position="bottom")
	#)
}

##' main call to comp1 plots
##' view 1
comp1view1Plot = function(cutoff,cancer,score,sample){
  if (sample == 'tumors'){
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, cutoff, cancer)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }
    }else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }
    }
  }else{
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapOG, cutoff, cancer)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
        
      }
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapTS, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
        
      }
    }else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotHeatmapPage1(resultsSub, score)        
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }
    }
  }
}
##' view 2
comp1view2Plot = function(cutoff,cancer,score,sample){
  if (sample == 'tumors'){
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }
  }else{
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapOG, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapTS, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, cutoff, cancer)      
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }
  }
}

##' main call to page1 gene data frame
geneDataFrameResultSet = function(cutoff,cancer,score,sample){
  if (sample == 'tumors'){
      #rgsog= NULL
      #rgsts= NULL
      #rgscom= NULL
    #if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, cutoff, cancer)
      rgsog <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Oncogene Score","Cancer","External links")
      #dfgenes
    #}else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, cutoff, cancer)
      #rgs <- unique(resultsSub$gene)
      rgsts <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Tumor Suppressor Score","Cancer","External links")
      #dfgenes
    #}else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, cutoff, cancer)
      #rgs <- unique(resultsSub$gene)
      rgscom <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgscom) <- c("Genes","Combined Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Combined Score","Cancer","External links")
      #dfgenes
    #}
    temp <- join(rgsog,rgsts,type="full")
    rgs <- join(temp,rgscom,type="full")
    gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
    temp <- data.frame(rgs[,c(1,2,4,5,3)],gc)
    dfgenes <- replace(temp, is.na(temp), "-")
    rm(temp)
    rm(rgs)
    colnames(dfgenes) <- c("Genes","OG Score","TS Score","Combined Score","Cancer","External links")
    if (nrow(dfgenes)>0){
      dfgenes
    }else{
      dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
      colnames(dfgenes) <- c("")
      dfgenes
    }
  }else{
    #if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapOG, cutoff, cancer)
      #rgs <- unique(resultsSub$gene)
      rgsog <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgsog) <- c("Genes","OG Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Oncogene Score","Cancer","External links")
      #dfgenes
    #}else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapTS, cutoff, cancer)
      #rgs <- unique(resultsSub$gene)
      rgsts <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgsts) <- c("Genes","TS Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Tumor Suppressor Score","Cancer","External links")
      #dfgenes
    #}else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, cutoff, cancer)
      #rgs <- unique(resultsSub$gene)
      rgscom <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
      colnames(rgscom) <- c("Genes","Combined Score","Cancer")
      #rgs <- unique(resultsSub$gene)
      #gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
      #dfgenes <- data.frame(rgs,gc)
      #colnames(dfgenes) <- c("Genes","Combined Score","Cancer","External links")
      #dfgenes
    #}
    temp <- join(rgsog,rgsts,type="full")
    rgs <- join(temp,rgscom,type="full")
    gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
    temp <- data.frame(rgs[,c(1,2,4,5,3)],gc)
    dfgenes <- replace(temp, is.na(temp), "-")
    rm(temp)
    rm(rgs)
    colnames(dfgenes) <- c("Genes","OG Score","TS Score","Combined Score","Cancer","External links")
    if (nrow(dfgenes)>0){
      dfgenes
    }else{
      dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
      colnames(dfgenes) <- c("")
      dfgenes
    }
  }
}