load("clsdata.RData")

getPage2Plots = function(updateProgress = NULL,cancer, gene, sampleSelection) {
  
  if (gene == "")
  {
    return()
  }
  if (is.null(gene))
  {
    return()
  }
  
	if (sampleSelection == 1) {
		#priorDetails = tcgaResults[[cancer]]$prioritize.details
		priorDetails = tcgaResultsPrioDetails[[cancer]]
		exprs.group1 = tcga.exprs[[cancer]][[1]]
		exprs.group2 = tcga.mn.exprs[[cancer]][[1]]
		#meth.group1 = tcga.meth[[cancer]]
		#meth.group2 = tcga.mn.meth[[cancer]]
		acgh.group1 = tcga.cna[[cancer]][[1]]
		acgh.group2 = tcga.mn.cna[[cancer]][[1]]
		achilles.ut = tcgaResultsAT[[cancer]][2]
		achilles.lt = tcgaResultsAT[[cancer]][1]
		color.palette=c("#E69F00", "#56B4E9")
		lab.group1="Tumors" 
		lab.group2="Normals"
		pvalue=TRUE
		cls =  tcgaResultsCLS[[cancer]]$cls # tcgaResults[[cancer]]$cls
	} else if (sampleSelection == 2) {
	  #priorDetails = ccleResults[[cancer]]$prioritize.details
		priorDetails = ccleResultsPrioDetails[[cancer]]
		exprs.group1 = ccle.exprs.combat[[cancer]][[1]]
		exprs.group2 = tcga.mn.exprs.combat[[cancer]][[1]]
		#meth.group1 = tcga.meth[[cancer]]
		#meth.group2 = tcga.mn.meth[[cancer]]
		acgh.group1 = ccle.cna.combat[[cancer]][[1]]
		acgh.group2 = tcga.mn.cna.combat[[cancer]][[1]]
		achilles.ut = ccleResultsAT[[cancer]][2]
		achilles.lt = ccleResultsAT[[cancer]][1]
		color.palette=c("#009E73", "#56B4E9")
		lab.group1="Cell lines" 
		lab.group2="Normals"
		pvalue=TRUE
		cls = ccleResultsCLS[[cancer]]$cls #ccleResults[[cancer]]$cls
	} else {
		priorDetails = tcgaResultsPrioDetails[[cancer]]
		exprs.group1 = tcga.exprs.combat[[cancer]][[1]]
		exprs.group2 = ccle.exprs.combat[[cancer]][[1]]
		#meth.group1 = tcga.meth.combat[[cancer]]
		#meth.group2 = tcga.mn.meth.combat[[cancer]]
		acgh.group1 = tcga.cna.combat[[cancer]][[1]]
		acgh.group2 = ccle.cna.combat[[cancer]][[1]]
		achilles.ut = ccleResultsAT[[cancer]][2]
		achilles.lt = ccleResultsAT[[cancer]][1]
		color.palette=c("#E69F00", "#009E73")
		lab.group1="Tumors" 
		lab.group2="Cell lines"
		pvalue=FALSE
		cls = tcgaResultsCLS[[cancer]]$cls #tcgaResults[[cancer]]$cls
	}
	#meth.anno = infinium450.probe.ann
  if (is.element(gene,rownames(achilles)) || length(cls) == 0)
  {
   achls = achilles[gene, cls, drop=FALSE] # achilles[gene,]    
  }else{
    achls = NULL
  }
	#achls = achilles
	
	plots = plotGene(gene, priorDetails, samples=NULL, sampleSelection,
		 	 exprs.group1, exprs.group2, 
		 	 #meth.group1, meth.group2, meth.anno, 
		 	 acgh.group1, acgh.group2, 
		 	 achls, achilles.ut, achilles.lt, 
		 	 lab.group1=lab.group1, lab.group2=lab.group2, 
		 	 color.palette=color.palette,
		 	 size=4, width=0.2, pvalue=pvalue)
# 	if (sampleSelection != 1) {
# 		plots[["methylation"]] = NULL
# 	}
	if (sampleSelection == 3) {
		plots[["achilles"]] = NULL
	}
	
  rm(priorDetails)
  rm(exprs.group1)
  rm(exprs.group2)
  rm(acgh.group1)
  rm(acgh.group2)
  rm(achilles.ut)
  rm(achilles.lt)
  rm(cls)
  rm(achls)
	plots
}
