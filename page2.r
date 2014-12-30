getPage2Plots = function(cancer, gene, sampleSelection) {
	if (sampleSelection == 1) {
		priorDetails = tcgaResults[[cancer]]$prioritize.details
		exprs.group1 = tcga.exprs[[cancer]][[1]]
		exprs.group2 = tcga.mn.exprs[[cancer]][[1]]
		meth.group1 = tcga.meth[[cancer]]
		meth.group2 = tcga.mn.meth[[cancer]]
		acgh.group1 = tcga.cna[[cancer]][[1]]
		acgh.group2 = tcga.mn.cna[[cancer]][[1]]
		achilles.ut = tcgaResults[[cancer]]$achilles.thresholds[2]
		achilles.lt = tcgaResults[[cancer]]$achilles.thresholds[1]
		color.palette=c("#E69F00", "#56B4E9")
		lab.group1="Tumors" 
		lab.group2="Normals"
		pvalue=TRUE
	} else if (sampleSelection == 2) {
		priorDetails = ccleResults[[cancer]]$prioritize.details
		exprs.group1 = ccle.exprs.combat[[cancer]][[1]]
		exprs.group2 = tcga.mn.exprs.combat[[cancer]][[1]]
		meth.group1 = tcga.meth[[cancer]]
		meth.group2 = tcga.mn.meth[[cancer]]
		acgh.group1 = ccle.cna.combat[[cancer]][[1]]
		acgh.group2 = ccle.mn.cna.combat[[cancer]][[1]]
		achilles.ut = ccleResults[[cancer]]$achilles.thresholds[2]
		achilles.lt = ccleResults[[cancer]]$achilles.thresholds[1]
		color.palette=c("#009E73", "#56B4E9")
		lab.group1="Cell lines" 
		lab.group2="Normals"
		pvalue=TRUE
	} else {
		priorDetails = tcgaResults[[cancer]]$prioritize.details
		exprs.group1 = tcga.exprs.combat[[cancer]][[1]]
		exprs.group2 = ccle.exprs.combat[[cancer]][[1]]
		meth.group1 = tcga.meth.combat[[cancer]]
		meth.group2 = tcga.mn.meth.combat[[cancer]]
		acgh.group1 = tcga.cna.combat[[cancer]][[1]]
		acgh.group2 = ccle.cna.combat[[cancer]][[1]]
		achilles.ut = ccleResults[[cancer]]$achilles.thresholds[2]
		achilles.lt = ccleResults[[cancer]]$achilles.thresholds[1]
		color.palette=c("#E69F00", "#009E73")
		lab.group1="Tumors" 
		lab.group2="Cell lines",
		pvalue=FALSE
	}
	meth.anno = infinium450.probe.ann
	achilles = achilles
	plots = plotGene(gene, priorDetails, samples=NULL, 
		 	 exprs.group1, exprs.group2, 
		 	 meth.group1, meth.group2, meth.anno, 
		 	 acgh.group1, acgh.group2, 
		 	 achilles, achilles.ut, achilles.lt, 
		 	 lab.group1=lab.group1, lab.group2=lab.group2, 
		 	 color.palette=color.palette,
		 	 size=4, width=0.2, pvalue=pvalue)
	if (sampleSelection != 1) {
		plots[["methylation"]] = NULL
	}
	if (sampleSelection == 3) {
		plots[["achilles"]] = NULL
	}
	
	plots
}

