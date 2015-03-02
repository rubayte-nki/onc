load("starter.RData")

##' Filtering of result frame according to user criteria
##' @param results data.frame with all results
##' @param scoreCutoff threshold that was selected by the user
##' @param cancerType cancer type selected by the user
##' @return a subset of the data.frame that fits the user's selection
##' @author Andreas Schlicker
page1DataFrame = function(results, scoreCutoff, cancerType, comstring) {
	# Filter the genes according to user's criteria
	if (scoreCutoff > 0) {
		genes = as.character(subset(results, cancer == cancerType & score.type == comstring & score >= as.integer(scoreCutoff))$gene)
	} else {
		genes = as.character(subset(results, cancer == cancerType & score.type == comstring & score <= as.integer(scoreCutoff))$gene)
	}
	
	# Sort the genes according to highest sum across all cancer types
	# Get the subset with the selected genes and drop unused levels
	# gene.order = subset(result.df, score.type=="combined" & gene %in% genes)
	gene.order = subset(results, score.type == comstring & gene %in% genes)
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
	result.df$gene <- factor(result.df$gene, levels=unique(as.character(result.df$gene)))
	colorLow = list(combined.score="#034b87", ts.score="gray98", og.score="gray98") 
	colorMid = list(combined.score="gray98")
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
	result.df$score.type = factor(result.df$score.type, levels=c("CNA", "Expr", "Meth", "Mut", "shRNA", "Combined"))
	result.df$gene <- factor(result.df$gene, levels=unique(as.character(result.df$gene)))
	
	# Overwrite the score column with the score type to make it categorical
	# Combined scores are not plotted later
	result.df[, 2] = as.character(result.df[, 2])
	result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 2] = as.character(result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 3])
	result.df[which(is.na(result.df[, 2]) | result.df[, 2] == "0"), 2] = "NONE"
	
	#ggplot(subset(result.df, score.type != "combined" & gene %in% topgenes), aes(x=score.type, y=gene)) + 
	ggplot(subset(result.df, score.type != "Combined"), aes(x=score.type, y=gene)) +  #coord_flip()  +
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
comp1view1Plot = function(updateProgress = NULL,cutoff,cancer,score,sample,inputdf = NULL){
  df = NULL
  if (sample == 'tumors'){
    if(score == 'og.score'){
      df = tcgaResultsHeatmapOG
    }else if(score == 'ts.score'){
      df = tcgaResultsHeatmapTS
    }else{
      df = tcgaResultsHeatmapCombined
    }
  }else{
    if(score == 'og.score'){
      df = ccleResultsHeatmapOG
    }else if(score == 'ts.score'){
      df = ccleResultsHeatmapTS
    }else{
      df = ccleResultsHeatmapCombined
    }
  }
  
  ## subset data frame based on user input
  resultsSub <- page1DataFrame(df, cutoff, cancer,"Combined")
  resultsSub <- resultsSub[resultsSub$score.type == "Combined",]
  ## if input dataframe is not null then update the target dataframe with the inputdf genes
  if (!(is.null(inputdf)))
  {
    temp <- as.data.frame(inputdf[,1])
    colnames(temp) <- c("gene")
    resultsSub <- plyr::join(temp,resultsSub,type="inner")          
  }
  ## sort the dataframe to match with results table
  if (cutoff > 0 || cutoff == -10)
  {
    resultsSub <- resultsSub[order(-resultsSub$"score"),]
  }else{
    resultsSub <- resultsSub[order(resultsSub$"score"),]
  }
  temp <- resultsSub[resultsSub$cancer == cancer,]
  temp <- arrange(temp,-row_number())
  resultsSub <- resultsSub[resultsSub$cancer != cancer,]
  resultsSub <- rbind(temp,resultsSub)
  if (nrow(resultsSub) > 0){
    ## call plot function
    plotHeatmapPage1(resultsSub, score)        
  }else{
    plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
    text(1,"Empty result set returned by filter. Nothing to plot.")
  }  

}

## for user file input
comp1view1FilePlot = function(updateProgress = NULL,cancer,inputdf,sample)
{
  if (sample == 'tumors'){
    resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined,-10,cancer,"Combined")
    resultsSub <- subset(resultsSub, cancer == cancer & gene %in% inputdf[,1])
    if (nrow(resultsSub) > 0){
      ## call plot function
      plotHeatmapPage1(resultsSub, score)        
    }else{
      plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
      text(1,"No genes match the app data. Nothing to plot.")
    }  
  }else{
    resultsSub <- page1DataFrame(ccleResultsHeatmapCombined,-10,cancer,"Combined")
    resultsSub <- subset(resultsSub, cancer == cancer & gene %in% inputdf[,1])
    if (nrow(resultsSub) > 0){
      ## call plot function
      plotHeatmapPage1(resultsSub, score)        
    }else{
      plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
      text(1,"No genes match the app data. Nothing to plot.")
    }  
  }
  
}

##' view 2
comp1view2Plot = function(updateProgress = NULL,cutoff,cancer,score,sample,inputdf = NULL){
  if (sample == 'tumors'){
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, cutoff, cancer,"Combined")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("gene")
        resultsSub <- plyr::join(temp,resultsSub,type="inner")          
      }
      ## sort the dataframe to match with results table
      if (cutoff > 0 || cutoff == -10)
      {
        resultsSub <- resultsSub[order(-resultsSub$"score"),]
      }else{
        resultsSub <- resultsSub[order(resultsSub$"score"),]
      }
      temp <- resultsSub[resultsSub$cancer == cancer,]
      temp2 <- temp[temp$score.type == "Combined",]
      temp <- temp[temp$score.type != "Combined",]
      temp2 <- arrange(temp2,-row_number())
      resultsSub <- resultsSub[resultsSub$cancer != cancer,]
      resultsSub <- rbind(temp2,temp,resultsSub)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else if(score == 'ts.score'){
      ## subset data frame based on user inputn 
      resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, cutoff, cancer,"Combined") 
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("gene")
        resultsSub <- plyr::join(temp,resultsSub,type="inner")          
      }
      ## sort the dataframe to match with results table
      if (cutoff > 0  || cutoff == -10)
      {
        resultsSub <- resultsSub[order(-resultsSub$"score"),]
      }else{
        resultsSub <- resultsSub[order(resultsSub$"score"),]
      }
      temp <- resultsSub[resultsSub$cancer == cancer,]
      temp2 <- temp[temp$score.type == "Combined",]
      temp <- temp[temp$score.type != "Combined",]
      temp2 <- arrange(temp2,-row_number())
      resultsSub <- resultsSub[resultsSub$cancer != cancer,]
      resultsSub <- rbind(temp2,temp,resultsSub)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else{
      ## subset data frame based on user input
      com <- page1DataFrame(tcgaResultsHeatmapCombined, cutoff, cancer,"Combined")
      genes <- unique(com$gene)
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- inputdf[,1]
        genes <- intersect(genes,temp)
      }
      og <- subset(tcgaResultsHeatmapOG, cancer == cancer & gene %in% genes)
      colnames(og) <- c('genes','ogs','score.type','cancer')
      ts <- subset(tcgaResultsHeatmapTS, cancer == cancer & gene %in% genes)
      colnames(ts) <- c('genes','tss','score.type','cancer')
      temp <- plyr::join(og,ts,type="inner")
      if (nrow(temp)>0)
      {
        cs <- abs(temp[,2] - temp[,5])
        res <- data.frame(temp[,1],cs,temp[,c(3,4)])
        colnames(res) <- c('gene','score','score.type','cancer')
        ## sort the dataframe to match with results table
        if (cutoff > 0  || cutoff == -10)
        {
          res <- res[order(-res$"score"),]
        }else{
          res <- res[order(res$"score"),]
        }
        temp <- res[res$cancer == cancer,]
        temp2 <- temp[temp$score.type == "Combined",]
        temp <- temp[temp$score.type != "Combined",]
        temp2 <- arrange(temp2,-row_number())
        res <- res[res$cancer != cancer,]
        res <- rbind(temp2,temp,res)
        plotCategoryOverview(res)  
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"No overlapping genes were found using the same cutoff score. Nothing to plot.")        
      }
      
    }
  }else{
    if(score == 'og.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapOG, cutoff, cancer,"Combined")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("gene")
        resultsSub <- plyr::join(temp,resultsSub,type="inner")          
      }
      ## sort the dataframe to match with results table
      if (cutoff > 0)
      {
        resultsSub <- resultsSub[order(-resultsSub$"score"),]
      }else{
        resultsSub <- resultsSub[order(resultsSub$"score"),]
      }
      temp <- resultsSub[resultsSub$cancer == cancer,]
      temp2 <- temp[temp$score.type == "Combined",]
      temp <- temp[temp$score.type != "Combined",]
      temp2 <- arrange(temp2,-row_number())
      resultsSub <- resultsSub[resultsSub$cancer != cancer,]
      resultsSub <- rbind(temp2,temp,resultsSub)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapTS, cutoff, cancer,"Combined")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("gene")
        resultsSub <- plyr::join(temp,resultsSub,type="inner")          
      }
      ## sort the dataframe to match with results table
      if (cutoff > 0  || cutoff == -10)
      {
        resultsSub <- resultsSub[order(-resultsSub$"score"),]
      }else{
        resultsSub <- resultsSub[order(resultsSub$"score"),]
      }
      temp <- resultsSub[resultsSub$cancer == cancer,]
      temp2 <- temp[temp$score.type == "Combined",]
      temp <- temp[temp$score.type != "Combined",]
      temp2 <- arrange(temp2,-row_number())
      resultsSub <- resultsSub[resultsSub$cancer != cancer,]
      resultsSub <- rbind(temp2,temp,resultsSub)
      if (nrow(resultsSub) > 0){
        ## call plot function
        plotCategoryOverview(resultsSub)             
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"Empty result set returned by filter. Nothing to plot.")
      }      
    }else{
      ## subset data frame based on user input
      com <- page1DataFrame(ccleResultsHeatmapCombined, cutoff, cancer,"Combined")
      genes <- unique(com$gene)
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- inputdf[,1]
        genes <- intersect(genes,temp)
      }
      og <- subset(ccleResultsHeatmapOG, cancer == cancer & gene %in% genes)
      colnames(og) <- c('genes','ogs','score.type','cancer')
      ts <- subset(ccleResultsHeatmapTS, cancer == cancer & gene %in% genes)
      colnames(ts) <- c('genes','tss','score.type','cancer')
      temp <- plyr::join(og,ts,type="inner")
      if (nrow(temp)>0)
      {
        cs <- abs(temp[,2] - temp[,5])
        res <- data.frame(temp[,1],cs,temp[,c(3,4)])
        colnames(res) <- c('gene','score','score.type','cancer')
        ## sort the dataframe to match with results table
        if (cutoff > 0  || cutoff == -10)
        {
          res <- res[order(-res$"score"),]
        }else{
          res <- res[order(res$"score"),]
        }
        temp <- res[res$cancer == cancer,]
        temp2 <- temp[temp$score.type == "Combined",]
        temp <- temp[temp$score.type != "Combined",]
        temp2 <- arrange(temp2,-row_number())
        res <- res[res$cancer != cancer,]
        res <- rbind(temp2,temp,res)
        plotCategoryOverview(res)  
      }else{
        plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
        text(1,"No overlapping genes were found using the same cutoff score. Nothing to plot.")        
      }
    }
  }
}
## for user file input
comp1view2FilePlot = function(updateProgress = NULL,cancer,inputdf,sample){

  if (sample == 'tumors'){    
    og <- subset(tcgaResultsHeatmapOG, cancer == cancer & gene %in% inputdf[,1])
    colnames(og) <- c('genes','ogs','score.type','cancer')
    ts <- subset(tcgaResultsHeatmapTS, cancer == cancer & gene %in% inputdf[,1])
    colnames(ts) <- c('genes','tss','score.type','cancer')
    temp <- plyr::join(og,ts,type="inner")
    if (nrow(temp)>0)
    {
      cs <- abs(temp[,2] - temp[,5])
      res <- data.frame(temp[,1],cs,temp[,c(3,4)])
      colnames(res) <- c('gene','score','score.type','cancer')
      plotCategoryOverview(res)  
    }else{
      plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
      text(1,"No overlapping genes were found using the same cutoff score. Nothing to plot.")        
    }
  }else{
    og <- subset(ccleResultsHeatmapOG, cancer == cancer & gene %in% inputdf[,1])
    colnames(og) <- c('genes','ogs','score.type','cancer')
    ts <- subset(ccleResultsHeatmapTS, cancer == cancer & gene %in% inputdf[,1])
    colnames(ts) <- c('genes','tss','score.type','cancer')
    temp <- plyr::join(og,ts,type="inner")
    if (nrow(temp)>0)
    {
      cs <- abs(temp[,2] - temp[,5])
      res <- data.frame(temp[,1],cs,temp[,c(3,4)])
      colnames(res) <- c('gene','score','score.type','cancer')
      plotCategoryOverview(res)  
    }else{
      plot(1,xaxt='n',yaxt='n',ann=FALSE,type="p",col="white")
      text(1,"No overlapping genes were found using the same cutoff score. Nothing to plot.")        
    }
  }
}





##' main call to page1 gene data frame
geneDataFrameResultSet = function(updateProgress = NULL,cutoff,cancer,score,sample,inputdf = NULL){
  
  rgsog= NULL
  rgsts= NULL
  rgscom= NULL
  dfgenes = NULL
  
  if (sample == 'tumors'){
    if(score == 'og.score'){
      
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, cutoff, cancer,"Combined")
      rgsog <- resultsSub[resultsSub[,4]== cancer,]
      rgsog <- reshape(rgsog[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      clist <- NULL
      for (i in 1:nrow(rgsog))
      {
        clist <- c(clist,cancer)
      }
      rgsog <- data.frame(rgsog,clist)
      colnames(rgsog) <- c("Genes","Oncogene Score","Meth","CNA","Mut","shRNA","Expr","Cancer")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
          temp <- as.data.frame(inputdf[,1])
          colnames(temp) <- c("Genes")
          rgsog <- plyr::join(temp,rgsog,type="left")          
      }
      ## handle empty result set
      if (nrow(rgsog)>0){
        ## select others
        resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, -10, cancer,"Combined")
        rgsts <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
        resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, -10, cancer, "Combined")
        rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgscom) <- c("Genes","Combined Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgsog,rgsts,type="left")
        rgs <- plyr::join(temp,rgscom,type="left")
        gc <- paste("<a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=",rgs[,1],"\">","Gene Card","</a>",sep="")
        temp <- data.frame(rgs[,c(1,2,9,10,3,4,5,6,7,8)],gc)
        temp <- temp[order(-temp$"Oncogene.Score"),] 
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgs)
        colnames(dfgenes) <- c("Genes","OG Score","TS Score","Combined Score","OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }      
      
    }else if(score == 'ts.score'){
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, cutoff, cancer,"Combined")
      rgsts <- resultsSub[resultsSub[,4]== cancer,]
      rgsts <- reshape(rgsts[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      clist <- NULL
      for (i in 1:nrow(rgsts))
      {
        clist <- c(clist,cancer)
      }
      rgsts <- data.frame(rgsts,clist)
      colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Meth","CNA","Mut","shRNA","Expr","Cancer")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
          temp <- as.data.frame(inputdf[,1])
          colnames(temp) <- c("Genes")
          rgsts <- plyr::join(temp,rgsts,type="left")          
      }
      ## handle empty result set
      if (nrow(rgsts)>0){
        ## select others
        resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, -10, cancer,"Combined")
        rgsog <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
        resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, -10, cancer, "Combined")
        rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgscom) <- c("Genes","Combined Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgsts,rgsog,type="left")
        rgs <- plyr::join(temp,rgscom,type="left")
        gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
        temp <- data.frame(rgs[,c(1,2,9,10,3,4,5,6,7,8)],gc)
        temp <- temp[order(-temp$"Tumor.Suppressor.Score"),]
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgs)
        colnames(dfgenes) <- c("Genes","TS Score","OG Score","Combined Score","TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      
    }else{
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, cutoff, cancer, "Combined")
      rgscom <- resultsSub[resultsSub[,4]== cancer,]
      rgscom <- reshape(rgscom[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      clist <- NULL
      for (i in 1:nrow(rgscom))
      {
        clist <- c(clist,cancer)
      }
      rgscom <- data.frame(rgscom,clist)
      colnames(rgscom) <- c("Genes","Oncogene Score","Tumor Suppressor Score","Combined Score","OG Score Affected","TS Score Affected","Combined Score Affected","Cancer")
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("Genes")
        rgscom <- plyr::join(temp,rgscom,type="left")            
      }
      ## handle empty result set
      if (nrow(rgscom)>0){
        ## select others
        resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, -10, cancer, "Combined")
        rgsog <- resultsSub[resultsSub[,4]== cancer,]
        rgsog <- reshape(rgsog[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
        colnames(rgsog) <- c("Genes","OG","OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr")
        
        #rgsog <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
        #colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")

        resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, -10, cancer, "Combined")
        rgsts <- resultsSub[resultsSub[,4]== cancer,]
        rgsts <- reshape(rgsts[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
        colnames(rgsts) <- c("Genes","TS","TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr")
        
        #rgsts <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
        #colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgscom,rgsog,type="left")
        rgs <- plyr::join(temp,rgsts,type="left")
        gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
        temp <- data.frame(rgs[,c(1,4,2,3,5,6,7,10,11,12,13,14,15,17,18,19,20,8)],gc)
        if (cutoff > 0 || cutoff == -10)
        {
          temp <- temp[order(-temp$"Combined.Score"),]          
        }else{
          temp <- temp[order(temp$"Combined.Score"),]
        }
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgscom)
        colnames(dfgenes) <- c("Genes","Combined Score","OG Score","TS Score","OG Score Affected","TS Score Affected","Combined Score Affected",
                               "OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr",
                               "TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      
    }
  }else{
    
    if(score == 'og.score'){
      
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapOG, cutoff, cancer, "Combined")
      rgsog <- resultsSub[resultsSub[,4]== cancer,]
      rgsog <- reshape(rgsog[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      if (nrow(rgsog)>0)
      {
        clist <- NULL
        for (i in 1:nrow(rgsog))
        {
          clist <- c(clist,cancer)
        }
        rgsog <- data.frame(rgsog,clist)
        colnames(rgsog) <- c("Genes","Oncogene Score","Meth","CNA","Mut","shRNA","Expr","Cancer")        
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("Genes")
        rgsog <- plyr::join(temp,rgsog,type="left")            
      }
      ## handle empty result set
      if (nrow(rgsog)>0){
        ## select others
        resultsSub <- page1DataFrame(ccleResultsHeatmapTS, -10, cancer, "Combined")
        rgsts <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
        resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, -10, cancer, "Combined")
        rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgscom) <- c("Genes","Combined Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgsog,rgsts,type="left")
        rgs <- plyr::join(temp,rgscom,type="left")
        gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
        temp <- data.frame(rgs[,c(1,2,9,10,3,4,5,6,7,8)],gc)
        temp <- temp[order(-temp$"Oncogene.Score"),]
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgs)
        colnames(dfgenes) <- c("Genes","OG Score","TS Score","Combined Score","OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }      
      
    }else if(score == 'ts.score'){
    
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapTS, cutoff, cancer,"Combined")
      rgsts <- resultsSub[resultsSub[,4]== cancer,]
      rgsts <- reshape(rgsts[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      if (nrow(rgsts)>0)
      {
        clist <- NULL
        for (i in 1:nrow(rgsts))
        {
          clist <- c(clist,cancer)
        }
        rgsts <- data.frame(rgsts,clist)
        colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Meth","CNA","Mut","shRNA","Expr","Cancer")
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("Genes")
        rgsts <- plyr::join(temp,rgsts,type="left")            
      }
      ## handle empty result set
      if (nrow(rgsts)>0){
        ## select others
        resultsSub <- page1DataFrame(ccleResultsHeatmapOG, -10, cancer, "Combined")
        rgsog <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
        resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, -10, cancer, "Combined")
        rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
        colnames(rgscom) <- c("Genes","Combined Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgsts,rgsog,type="left")
        rgs <- plyr::join(temp,rgscom,type="left")
        gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
        temp <- data.frame(rgs[,c(1,2,9,10,3,4,5,6,7,8)],gc)
        temp <- temp[order(-temp$"Tumor.Suppressor.Score"),]
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgs)
        colnames(dfgenes) <- c("Genes","TS Score","OG Score","Combined Score","TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      
    }else{
      
      ## subset data frame based on user input
      resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, cutoff, cancer, "Combined")
      rgscom <- resultsSub[resultsSub[,4]== cancer,]
      rgscom <- reshape(rgscom[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
      if (nrow(rgscom)>0)
      {
        clist <- NULL
        for (i in 1:nrow(rgscom))
        {
          clist <- c(clist,cancer)
        }
        rgscom <- data.frame(rgscom,clist)
        colnames(rgscom) <- c("Genes","Oncogene Score","Tumor Suppressor Score","Combined Score","OG Score Affected","TS Score Affected","Combined Score Affected","Cancer")
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }      
      ## if input dataframe is not null then update the target dataframe with the inputdf genes
      if (!(is.null(inputdf)))
      {
        temp <- as.data.frame(inputdf[,1])
        colnames(temp) <- c("Genes")
        rgscom <- plyr::join(temp,rgscom,type="left")            
      }
      ## handle empty result set
      if (nrow(rgscom)>0){
        ## select others
        resultsSub <- page1DataFrame(ccleResultsHeatmapOG, -10, cancer, "Combined")
        rgsog <- resultsSub[resultsSub[,4]== cancer,]
        rgsog <- reshape(rgsog[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
        colnames(rgsog) <- c("Genes","OG","OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr")
        
        #rgsog <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
        #colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
        
        resultsSub <- page1DataFrame(ccleResultsHeatmapTS, -10, cancer, "Combined")
        rgsts <- resultsSub[resultsSub[,4]== cancer,]
        rgsts <- reshape(rgsts[,c(1,2,3)], direction = "wide", idvar="gene",timevar='score.type')
        colnames(rgsts) <- c("Genes","TS","TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr")
        
        #rgsts <- resultsSub[resultsSub[,3] == 'combined' & resultsSub[,4]== cancer,c(1,2,4)]
        #colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
        ## make final data frame
        temp <- plyr::join(rgscom,rgsog,type="left")
        rgs <- plyr::join(temp,rgsts,type="left")
        gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',rgs[,1],'">','Gene Card','</a>',sep='')
        temp <- data.frame(rgs[,c(1,4,2,3,5,6,7,10,11,12,13,14,15,17,18,19,20,8)],gc)
        if (cutoff > 0 || cutoff == -10)
        {
          temp <- temp[order(-temp$"Combined.Score"),]          
        }else{
          temp <- temp[order(temp$"Combined.Score"),]
        }
        dfgenes <- replace(temp, is.na(temp), "-")
        rm(temp)
        rm(rgscom)
        colnames(dfgenes) <- c("Genes","Combined Score","OG Score","TS Score","OG Score Affected","TS Score Affected","Combined Score Affected",
                               "OG.Meth","OG.CNA","OG.Mut","OG.shRNA","OG.Expr",
                               "TS.Meth","TS.CNA","TS.Mut","TS.shRNA","TS.Expr","Cancer","External links")
        dfgenes
      }else{
        dfgenes <- data.frame(c("Empty result set returned by filter. Nothing to show."))
        colnames(dfgenes) <- c("Empty result set")
        dfgenes
      }
      
    }
    
  }
}

geneFileDataFrameResultSet = function(updateProgress= NULL,cancer,score,inputdf,sample){
  
  if (sample == 'tumors'){
    
    res <- NULL
    ## subset data frame based on user input
    resultsSub <- page1DataFrame(tcgaResultsHeatmapOG, -10, cancer,"Combined")
    rgsog <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
    ## select others
    resultsSub <- page1DataFrame(tcgaResultsHeatmapTS, -10, cancer,"Combined")
    rgsts <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
    resultsSub <- page1DataFrame(tcgaResultsHeatmapCombined, -10, cancer, "Combined")
    rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgscom) <- c("Genes","Combined Score","Cancer")
    ## make final data frame
    temp <- plyr::join(rgsog,rgsts,type="left")
    rgs <- plyr::join(temp,rgscom,type="left")
    temp <- inputdf[,1]
    temp <- as.data.frame(temp)
    colnames(temp) <- c("Genes")
    res <- plyr::join(temp,rgs,type="left")
    cnc <- NULL
    for (i in 1:nrow(res))
    {
      cnc <- c(cnc,cancer)
    }
    gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',res[,1],'">','Gene Card','</a>',sep='')
    res <- data.frame(res[,c(1,2,4,5)],cnc,gc)
    colnames(res) <- c("Genes","OG Score","TS Score","Combined Score","Cancer","External links")
    res <- data.frame(res,inputdf)
    res <- replace(res, is.na(res), "-")
    rm(temp)
    rm(rgs) 
    if (nrow(res)>0){
      res
    }else{
      res <- data.frame(c("Empty result set returned by filter. Nothing to show."))
      colnames(res) <- c("Empty result set")
      res
    }
    
  }else{
    
    res <- NULL
    ## subset data frame based on user input
    resultsSub <- page1DataFrame(ccleResultsHeatmapOG, -10, cancer,"Combined")
    rgsog <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgsog) <- c("Genes","Oncogene Score","Cancer")
    ## select others
    resultsSub <- page1DataFrame(ccleResultsHeatmapTS, -10, cancer,"Combined")
    rgsts <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgsts) <- c("Genes","Tumor Suppressor Score","Cancer")
    resultsSub <- page1DataFrame(ccleResultsHeatmapCombined, -10, cancer, "Combined")
    rgscom <- resultsSub[resultsSub[,3] == 'Combined' & resultsSub[,4]== cancer,c(1,2,4)]
    colnames(rgscom) <- c("Genes","Combined Score","Cancer")
    ## make final data frame
    temp <- plyr::join(rgsog,rgsts,type="left")
    rgs <- plyr::join(temp,rgscom,type="left")
    temp <- inputdf[,1]
    temp <- as.data.frame(temp)
    colnames(temp) <- c("Genes")
    res <- plyr::join(temp,rgs,type="left")
    cnc <- NULL
    for (i in 1:nrow(res))
    {
      cnc <- c(cnc,cancer)
    }
    gc <- paste('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',res[,1],'">','Gene Card','</a>',sep='')
    res <- data.frame(res[,c(1,2,4,5)],cnc,gc)
    colnames(res) <- c("Genes","OG Score","TS Score","Combined Score","Cancer","External links")
    res <- data.frame(res,inputdf)
    res <- replace(res, is.na(res), "-")
    rm(temp)
    rm(rgs)
    if (nrow(res)>0){
      res
    }else{
      res <- data.frame(c("Empty result set returned by filter. Nothing to show."))
      colnames(res) <- c("Empty result set")
      res
    }
    
  }
  
}
