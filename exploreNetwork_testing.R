# Libraries needed:
#	WGCNA
#	RDAVIDWebService
#	biomaRt

# Depends on:
#	.getkMERankAndQuantileForGenes
#	.checkGeneListEnrichmentList
#	.checkTermGeneEnrichmentAllModules
#	.getAllTermsWithGenesAcrossModules
#	.getDAVIDchartForGeneList
#
## Arguments
###
#
## Value
###

.exploreGenes = function(genes, modGenesList, modkMEList, modDAVID, thresh=.05, mfrow=NULL, verbose=T) {
	genesInNet = genes[genes %in% unlist(modGenesList)];
	if (verbose) {cat('Running...\n  .getkMERankAndQuantileForGenes\n')}
	genekMEs = .getkMERankAndQuantileForGenes(genesInNet, modkMEList);
	if (verbose) {cat('  .checkGeneListEnrichmentList\n')}
	modEnrichments = .checkGeneListEnrichmentList(genes, modGenesList, unlist(modGenesList));
	mods = names(table(genekMEs$module));
	if (verbose) {cat('  .checkTermGeneEnrichmentAllModules\n')}
	termEnrichments = .checkTermGeneEnrichmentAllModules(modDAVID=modDAVID, genes=genesInNet, modGenes=modGenesList, thresh=thresh);
	if (verbose) {cat('  .getAllTermsWithGenesAcrossModules\n')}
	terms = .getAllTermsWithGenesAcrossModules(genesInNet, modDAVID);
	if (verbose) {cat('  .getDAVIDchartForGeneList\n')}
	genesDAVID = suppressWarnings(.getDAVIDchartForGeneList(genes=genesInNet, background=unlist(modGenesList)));
	#.verboseBoxplotComparekMEOfGeneSetsWithinAllModules(modkMEList, genesInNet, mfrow)
	#verboseBoxplot(as.numeric(genekMEs), genekMEs$module)
	return(list(kME=genekMEs, modEnrichments=modEnrichments, termEnrichments=termEnrichments, terms=terms, genesDAVID=genesDAVID));
}

# Depends on:
#	.computekME
#	.getModuleGenes
#	.getModGenesRankedBykME
.getNetworkBasics = function (net, data, suffix='.1') {
	assign(paste('dendro', suffix, sep=''), net$dendrograms, envir=.GlobalEnv);
	assign(paste('blockGenes', suffix, sep=''), net$blockGenes, envir=.GlobalEnv);
	assign(paste('colors', suffix, sep=''), net$colors, envir=.GlobalEnv);
	assign(paste('MEs', suffix, sep=''), net$MEs, envir=.GlobalEnv);
	assign(paste('kME', suffix, sep=''), .computekME(data, net$MEs)$all, envir=.GlobalEnv);
	assign(paste('modGenes', suffix, sep=''), .getModuleGenes(data, net$colors), envir=.GlobalEnv);
	assign(paste('modkMEs', suffix, sep=''), 
		   .getModGenesRankedBykME(names(table(net$colors)),
		   						   net$colors,
		   						   .computekME(data, net$MEs)$all
		   						   ), 
		   	envir=.GlobalEnv);
}

#######################################################################################################################################
######################### module genes and IDs
###	.getModuleGenes
###	.convertIDsWithBiomaRt
###	.convertIDsWithBiomaRtList
###	.getAvgNumForGenes
###	.getAvgNumForModules
###	.fundamentalModuleConcepts - deprecate
###	.fundamentalGeneListConcepts
###	.fundamentalGeneListConceptsList
### .buildModuleGenesInfo
#######################################################################################################################################

# Depends on:
#	NA	
#
## Arguments
### DATA: data frame of gene expression data, rows=samples, columns=genes
### colors: character vector of module color assignments for genes in DATA, must be in same order as the columns in DATA
#
## Value
### modGenes: list of character vectors (one for each module) containing colnames from DATA, organized by module color

.getModuleGenes = function(DATA, colors) {
	modNames = names(table(colors));
	modGenes = list();
	for (mod in 1:length(modNames)) {
		modGenes[[mod]] = names(DATA)[colors==modNames[mod]];
		names(modGenes)[mod] = modNames[mod];
	}
	return(modGenes);
}

# Depends on:
#	biomaRt
#
## Arguments
### ids:
#
## Value
### biomart:

.convertIDsWithBiomaRt = function (biomart='ensembl', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', ids) {
	mart = useMart(biomart=biomart, dataset=dataset);
	return(getBM(attributes=c(fromID, toID), filters=fromID, values=ids, mart=mart));
}

# Depends on:
#	.convertIDsWithBiomaRt
#
## Arguments
### modGenes:
#
## Value
### IDs:

.convertIDsWithBiomaRtList = function (biomart='ensembl', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', modGenes) {
	IDs = list();
	for (mod in 1:length(modGenes)) {
		cat(mod, ':', names(modGenes)[mod], ' ', sep='');
		IDs[[mod]] = .convertIDsWithBiomaRt(biomart=biomart, dataset=dataset, fromID=fromID, toID=toID, ids=modGenes[[mod]]);
		names(IDs)[mod] = names(modGenes)[mod];
	}
	cat('\n');
	return(IDs)
}

# Depends on:
#	NA
#
## Arguments
### genes: character vector of gene ids, must be subset of names(nums)
### nums: numeric vector containing some gene-wise measure, must be named
#
## Value
### float:

.getAvgNumForGenes = function (genes, nums) {
	return(mean(nums[names(nums) %in% genes], na.rm=T));
}

# Depends on:
#	.getAvgNumForGenes
#
## Arguments
### modGenes:
### nums: 
#
## Value
### modNum:

.getAvgNumForModules = function (modGenes, nums) {
	modNum = c();
	for ( m in modGenes ) {
		modNum = c(modNum, .getAvgNumForGenes(m, nums))
	}
	names(modNum) = names(modGenes);
	return(modNum)
}

# Depends on:
#	WGCNA
#
## Arguments
### colors:
### DATA:
### type, power:
### verbose:
#
## Value
### modFNC:

.fundamentalModuleConcepts = function(colors, DATA, type='signed', power=18, verbose=T) {
	modNames = names(table(colors));
	modFNC = list();
	if (verbose) {cat('Working on...\n')};
	for (mod in 1:length(modNames)) {
		if (verbose) {cat(modNames[mod], ' ', sep='')};
		modFNC[[mod]] = fundamentalNetworkConcepts(adjacency(DATA[colors==modNames[mod]], type=type, power=power));
		names(modFNC)[mod] = modNames[mod];
	}
	return(modFNC);
}

# Depends on:
#	WGCNA
#
## Arguments
### genes:
### DATA:
### type, power:
#
## Value
### genesFNC:

.fundamentalGeneListConcepts = function(genes, DATA, type='signed', power=18) {
	genesInDATA = intersect(genes, names(DATA));
	gDATA = DATA[, match(genesInDATA, names(DATA))];
	genesFNC = fundamentalNetworkConcepts(adjacency(gDATA, type=type, power=power));
	genesFNC[['Expression']] = gDATA;
	return(genesFNC);
}

# Depends on:
#	.fundamentalGeneListConcepts
#
## Arguments
###			REDUNDANT WITH .fundamentalModuleConcepts?
#
## Value
###
.fundamentalGeneListConceptsList = function(genesList, DATA, type='signed', power=18, verbose=T) {
	genesFNC = list();
	if (verbose) {cat('Working on...\n')};
	for (l in 1:length(genesList)) {
		if (verbose) {cat(names(genesList)[l], ' ', sep='')};
		genesInDATA = intersect(genesList[[l]], names(DATA));
		gDATA = DATA[, match(genesInDATA, names(DATA))];
		genesFNC[[l]] = fundamentalNetworkConcepts(adjacency(gDATA, type=type, power=power));
		genesFNC[[l]][['Expression']] = gDATA;
	}
	names(genesFNC) = names(genesList);
	return(genesFNC);
}

# Depends on:
#	.fundamentalGeneListConcepts

.buildModuleGenesInfo = function (modkMEList, DATA, networkType='signed', power=18) {
	out = modkMEList;
	for (mod in 1:length(out)) {
		mGenes = rownames(out[[mod]]);
		mDATA = DATA[, match(mGenes, names(DATA))];
		mFNC = .fundamentalGeneListConcepts(mGenes, mDATA, type=networkType, power=power);
		out[[mod]] = cbind(out[[mod]], Connectivity=mFNC$Connectivity, ScaledConnectivity=mFNC$ScaledConnectivity, ClusterCoef=mFNC$ClusterCoef, MAR=mFNC$MAR);
	}
	return(out);
}

# Depends on:
#	
.plotModuleFNCrelationships = function (buildModuleGenesInfoOUTPUT, folder='FNC_relationships', width=9, height=6, units='in', quality=100, type='quartz', res=150, mfrow=NULL, verbose=T, ...) {
	fnc = buildModuleGenesInfoOUTPUT;
	dir.create(folder)
	for (mod in 1:length(fnc)) {
		if (verbose) {
			cat(mod, ':', names(fnc)[mod], ' ', sep='');
		}
		tmp = fnc[[mod]][, names(fnc[[mod]]) %in% c(paste('kME', names(fnc)[mod], sep=''), 'ScaledConnectivity', 'ClusterCoef', 'MAR')];
		jpeg(file=paste(folder, '/', names(fnc)[mod], '.jpg', sep=''), width=width, height=height, units=units, quality=quality, type=type, res=res);
		par(mfrow=c(2,3));
		for (j in 1:3) {
			for (i in (j+1):4) {
				verboseScatterplot(tmp[,j], tmp[,i], abline=T, abline.col='red', 
								   frame.plot=F, xlim=c(0,1), ylim=c(0,1), 
								   xlab=colnames(tmp)[j], ylab=colnames(tmp)[i], 
								   col=names(fnc)[mod], ...);
			}
		}
		dev.off();
	}

}

# careful of kME relationships in plots since kME from whole module while others just from this list
.getGenesFNCkMErelationships = function (genes, modkMEs, DATA, type='signed', power=18, plot=T, ...) {
	gkME = .getkMERankAndQuantileForGenes(genes=genes, modkMEList=modkMEs);
	gFNC = .fundamentalGeneListConcepts(genes=genes, DATA=DATA, type=type, power=power);
	gkMEFNC = cbind(gkME, Connectivity=gFNC$Connectivity, ScaledConnectivity=gFNC$ScaledConnectivity, ClusterCoef=gFNC$ClusterCoef, MAR=gFNC$MAR);
	
	if (plot) {
		tmp = gkMEFNC[, names(gkMEFNC) %in% c('kME', 'ScaledConnectivity', 'ClusterCoef', 'MAR')];
		par(mfrow=c(2,3));
		for (j in 1:3) {
			for (i in (j+1):4) {
				verboseScatterplot(tmp[,j], tmp[,i], abline=T, abline.col='red', 
								   frame.plot=F, xlim=c(0,1), ylim=c(0,1), 
								   xlab=colnames(tmp)[j], ylab=colnames(tmp)[i], 
								   col=gkMEFNC$module, 
								   ...);
			}
		}
	}
	return(list(genes=gkMEFNC, density=gFNC$Density, centralization=gFNC$Centralization, heterogeneity=gFNC$Heterogeneity))
}



# mars = c()
# for ( mod in 1:length(x) ) {
	# thisMAR = x[[mod]]$MAR;
	# names(thisMAR) = rownames(x[[mod]]);
	# mars = c(mars, thisMAR)
# }





#######################################################################################################################################
######################### eigengenes and kME
###	.getkMERankAndQuantileForGenes
###	.computekME
###	.getModGenesRankedBykME
###	.getModulekME
###	.getkMEsFromAssignedModules
###	.getTopkMEFromEachModule 
###	.getGenesAvgkMEInMod
###	.plotEigengeneNetworks2
###	.plotMETraitCorsAgainstNums
###	.plotMEsExprWithTraits
#######################################################################################################################################

# Depends on:
#	NA
.getkMERankAndQuantileForGenes = function (genes, modkMEList) {
	mat = as.data.frame(matrix(nrow=length(genes), ncol=4));
	dimnames(mat) = list(genes, c('module','kME','rank','quantile'));
	for (g in 1:length(genes)) {
		gene = genes[g];
		tmp = lapply(modkMEList, function(f) if(any(grepl(gene, rownames(f)))){f});
		tmp = tmp[!sapply(tmp, is.null)][[1]];
		module = gsub('kME', '', names(tmp)[1]);
		row = which(rownames(tmp)==gene);
		kME = tmp[row, 1];
		qqtile = row / nrow(tmp);
		mat[g, ] = c(module, signif(kME,3), row, signif(qqtile,3));
	}
	return(mat);
}

# Depends on:
#	WGCNA
.computekME = function(DATA, MEs, modules=F, colors=NULL) {
	MEs = orderMEs(MEs);
	kME = as.data.frame(cor(DATA, MEs, use='p'));
	names(kME) = paste('k', names(kME), sep='');
	kMEpval = as.data.frame(corPvalueStudent(as.matrix(kME), nrow(DATA)));
	names(kMEpval) = paste('p.', names(kMEpval), sep='');
	
	out = as.data.frame(matrix(nrow=nrow(kME), ncol=(2*ncol(kME))));
	kcols = seq(1, ncol(out), 2); 
	pcols = kcols+1;
	out[, kcols] = kME;
	names(out)[kcols] = names(kME);
	out[, pcols] = kMEpval;
	names(out)[pcols] = names(kMEpval);
	rownames(out) = rownames(kME);
	
	modList = NULL;
	if (modules) {
		modList = list();
		modNames = names(table(colors));
		for (mod in 1:length(modNames)) {
			modList[[mod]] = out[colors==modNames[mod], ];
			names(modList)[mod] = modNames[mod];
		}
	}

	return(list(all=out, mods=modList));
}

# Depends on:
#	.getModulekME
.getModGenesRankedBykME = function(module_names, colors, kME) {
	outList = list();
	for (m in 1:length(module_names)) {
		outList[[m]] = .getModulekME(module_names[m],colors,kME)
	}
	names(outList) = module_names;
	return(outList);
}

# Depends on:
#	NA
.getModulekME = function(module, colors, kMEtable, kMEpvals=NULL, orderByRank=T) {
	modGenes = colors==module;
	if (is.null(kMEpvals)) {
		#modCols = grep(module, names(kMEtable));#######
		modCols = which(gsub('ME|kME|p.kME', '', names(kMEtable)) == module);
		out = kMEtable[modGenes, modCols];
	} else {
		modCol = match(module, gsub('kME', '', names(kMEtable)));
		out = cbind(kMEtable[modGenes, modCol], kMEpvals[modGenes, modCol]);
		dimnames(out) = list(rownames(kMEtable)[modGenes], c(names(kMEtable)[modCol], names(kMEpvals)[modCol]));
	}
	if (orderByRank) {
		out = out[order(out[, 1], decreasing=T), ];
	}
	return(out);
}

# Depends on:
#	NA
.getkMEsFromAssignedModules = function (modkMEList, matchToData=T, data) {
	mat = modkMEList[[1]];
	names(mat) = c('kME','p.kME');
	for (m in 2:length(modkMEList)) {
		names(modkMEList[[m]]) = c('kME','p.kME');
		mat = rbind(mat, modkMEList[[m]])
	}
	if (matchToData) {
		mat = mat[match(names(data), rownames(mat)), ];
	}
	return(mat);
}

# Depends on:
#	NA
.getTopkMEFromEachModule = function(kMEtable, colors) {
	if (nrow(kMEtable)!=length(colors)) {
		stop();
	}
	modNames = names(table(colors));
	hubs = as.data.frame(matrix(nrow=length(modNames), ncol=2));
	rownames(hubs) = modNames;
	names(hubs) = c('id', 'kME');
	for (mod in 1:length(modNames)) {
		ind = match(paste('kME', modNames[mod], sep=''), names(kMEtable));
		tmpkME = kMEtable[colors==modNames[mod], ];
		tmpkME = tmpkME[order(tmpkME[, ind], decreasing=T), ];
		hubs[mod, ] = c(rownames(tmpkME)[1], tmpkME[1, ind]); 
	} 
	return(hubs);
}

# Depends on:
#	NA
.getGenesAvgkMEInMod = function(genes, modkME) {
	return(mean(modkME[rownames(modkME) %in% genes, 1]));
}

# Depends on:
#	WGCNA
.plotEigengeneNetworks2 = function(MEs, heatmapColors=blueWhiteRed(50), returnCors=F) {
	plotEigengeneNetworks(orderMEs(MEs), 
						  setLabels='', marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2), 
						  heatmapColors=heatmapColors);
	if (returnCors) {
		return(cor(MEs));
	}
}

# Depends on:
#	WGCNA
.plotMETraitCorsAgainstNums = function (traitCors, modNums, mfrow, ylim, xlab=NULL, ...) {
	if (is.null(xlab)) {xlab = deparse(substitute(modNums))}
	par(mfrow=mfrow);
	for (tr in 1:ncol(traitCors)) {
		verboseScatterplot(modNums, 
						   traitCors[match(names(modNums), gsub('ME','',rownames(traitCors))), tr],
						   abline=T,
						   pch=19,
						   col=names(modNums),
						   frame.plot=F,
						   xlab=xlab,
						   ylab=colnames(traitCors)[tr],
						   cex=2,
						   ylim=ylim, ...
						   );
	}
}

# Depends on:
#	.checkGeneListEnrichmentList
#	.plotMETraitCorsAgainstNums
.plotMETraitCorsAgainstGeneEnrichmentOddsRatio = function (genes, modGenes, background, traitCors, mfrow=c(2,4), ylim=c(-1,1), return=T) {
	enrich = .checkGeneListEnrichmentList(genes, modGenes, background);
	x = enrich$pvals$ratio;
	names(x) = rownames(enrich$pvals);
	.plotMETraitCorsAgainstNums(traitCors=traitCors,modNums=x,mfrow=mfrow,ylim=ylim,corOptions="method='s'", xlab=paste(deparse(substitute(genes)), ' odds ratio',sep=''));
	if (return) {
		return(enrich);
	}
}

# Depends on:
#	WGCNA
.plotMEsExprWithTraits = function (MEs, factors, mfrow, ...) {
	par(mfrow=mfrow);
	for (i in 1:ncol(MEs)) {
		verboseBoxplot(MEs[,i], as.factor(factors), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i], ...);
	}
}

#######################################################################################################################################
######################### gene list enrichments and analysis
###	.checkGeneListEnrichment
###	.checkGeneListEnrichmentList
#######################################################################################################################################




# Depends on:
#	NA
.checkGeneListEnrichment = function(dataset1, dataset2, refset, alt='two.sided') {
	data1name = as.character(match.call()[2]);
	data2name = as.character(match.call()[3]);
	incommon = sum(dataset1 %in% dataset2);
	outof = length(dataset2);
	innet = sum(dataset1 %in% refset);
	countdimnames = list(c("Y", "N", ""), c("Y", "N", ""));
	names(countdimnames)[[1]] = data2name;
	names(countdimnames)[[2]] = data1name;
	counts0 = matrix(ncol = 3, nrow = 3, dimnames = countdimnames);
	counts0[1,] = c(incommon, outof-incommon, outof);
	counts0[3,] = c(innet, length(refset)-innet, length(refset));
	counts0[2,] = counts0[3,] - counts0[1,];
	counts = counts0[-3,-3];
	test = fisher.test(counts, alternative = alt);
	out = list(counts0, test);
	return(out);
}

# Depends on:
#	.checkGeneListEnrichment
.checkGeneListEnrichmentList = function(dataset1, dataset2list, refset, thresh=.05, alt='two.sided', order=T) {			   	
	data1name = as.character(match.call()[2]);
	out=list(); 
	modpvals = as.data.frame(matrix(nrow=length(dataset2list), ncol=2, dimnames=list(names(dataset2list), c('pval', 'ratio'))));
	for (set in 1:length(dataset2list)) {
		temp = .checkGeneListEnrichment(dataset1, dataset2list[[set]], refset, alt=alt);
		names(dimnames(temp[[1]])) = c(names(dataset2list[set]), data1name);
		out[[set]] = temp;
		names(out)[set] = names(dataset2list[set]);
		modpvals[set, ] = c(signif(temp[[2]]$p.value, 3), signif(temp[[2]]$estimate, 3))
	}
	if(order){modpvals = modpvals[order(as.numeric(modpvals$pval)), ]};
	return(list(details = out, pvals = modpvals));
}

#######################################################################################################################################
######################### explore DAVID results
###	.getTermInfoAcrossMods
###	.checkTermGeneEnrichment
###	.checkTermGeneEnrichmentModule
###	.checkTermGeneEnrichmentAllModules
###	.countTermsWithPatternInModules
###	.modDAVIDByTerm
###	.getAllTermsWithGenesAcrossModules
###	.getAllTermsWithGenes
###	.getAllTermsWithGene
###	.addkMEColToDAVID
###	.addkMEColToDAVIDList
###	.addGeneAvgColToDAVID
###	.getTermGenes
###	.filterModDAVIDList
###	.checkNamesOfModListElements
###	.filterDAVIDChart
#######################################################################################################################################

# Depends on:
#	.modDAVIDByTerm
.getTermInfoAcrossMods = function(term, modDAVID, modDAVIDByTerm.out=NULL) {
	if (is.null(modDAVIDByTerm.out)) {
		cat('...using exn.modDAVIDByTerm() to get distribution of term "', term, '" across modules\n', sep='');
		modDAVIDByTerm.out = .modDAVIDByTerm(modDAVID);
	} else if (any(names(modDAVIDByTerm.out) != c('modsByTerm', 'uniqueToMod'))) {
		stop('modDAVIDByTerm.out must be output from exn.modDAVIDByTerm()');
	}
	termInd = match(term, names(modDAVIDByTerm.out$modsByTerm));
	termMods = modDAVIDByTerm.out$modsByTerm[[termInd]];
	tmpDAVID = modDAVID[match(termMods, names(modDAVID))];
	out = vector(mode='list', length=length(tmpDAVID));
	names(out) = names(tmpDAVID);
	for (m in 1:length(tmpDAVID)) {
		out[[m]] = tmpDAVID[[m]][tmpDAVID[[m]]$Term==term, ];
	}
	return(out);
}

# Depends on:
#	.checkGeneListEnrichment
.checkTermGeneEnrichment = function(modDAVIDmod, row, genes, bg=NULL) {
	termGenes = unlist(strsplit(modDAVIDmod[row,]$Genes, ', '));
	if (is.null(bg)) {
		bg = unique(unlist(strsplit(modDAVIDmod$Genes, ', ')));
	} else if (is.character(bg) & all(termGenes %in% bg)) {
		bg = bg;
	} else {
		warning('Using all term genes in module as background because bg is invalid');
		bg = unique(unlist(strsplit(modDAVIDmod$Genes, ', ')));
	}
	return(.checkGeneListEnrichment(genes, termGenes, bg));
}

# Depends on:
#	.checkTermGeneEnrichment
.checkTermGeneEnrichmentModule = function(modDAVIDmod, genes, bg=NULL, thresh=NULL) {
	outlist = list();
	for (row in 1:nrow(modDAVIDmod)) {
		tmp = .checkTermGeneEnrichment(modDAVIDmod=modDAVIDmod, row=row, genes=genes, bg=bg);
		if (is.numeric(thresh)) {
			if (tmp[[2]]$p.value <= thresh) {
				outlist[[length(outlist)+1]] = tmp;
				names(outlist)[length(outlist)] = modDAVIDmod$Term[row];
			}
		} else {
			outlist[[row]] = tmp;
			names(outlist)[row] = modDAVIDmod$Term[row];
		}
	}
	return(outlist);
}

# Depends on:
#	.checkTermGeneEnrichmentModule
.checkTermGeneEnrichmentAllModules = function(modDAVID, genes, modGenes, bg=NULL, thresh=NULL) {
	if (any(names(modDAVID)!=names(modGenes))) {
		stop('modDAVID and modGenes must have same modules in same order');
	}
	outlist = list();
	for ( mod in 1:length(modDAVID) ) {
		#print(names(modDAVID)[mod])
		if (sum(genes %in% modGenes[[mod]]) == 0) { outlist[[mod]] = 'None in module'; next };
		if (nrow(modDAVID[[mod]]) == 0 | all(is.na(modDAVID[[mod]]))) { outlist[[mod]] = 'No terms'; next };
		outlist[[mod]] = .checkTermGeneEnrichmentModule(modDAVID[[mod]], genes=genes, bg=modGenes[[mod]], thresh=thresh);
	}
	names(outlist) = names(modDAVID);
	return(outlist);
}

# Depends on:
#	.modDAVIDByTerm
.countTermsWithPatternInModules = function (pattern, modDAVID) {
	terms = .modDAVIDByTerm(modDAVID);
	termsU = unlist(terms$uniqueToMod);
	all = sort(table(unlist(terms[[1]][grepl(pattern, names(terms[[1]]))])));
	uni = sort(table(termsU[grepl(pattern, names(termsU))]));
	return(list(all=all,uniqueToMod=uni));
}

# Depends on:
#	.checkNamesOfModListElements
.modDAVIDByTerm = function(modDAVID) {
	.checkNamesOfModListElements(modDAVID);
	tmp = unlist(modDAVID);
	terms = tmp[grepl('Term', names(tmp))];
	terms = terms[!is.na(terms)];
	termList = list();
	for (tt in 1:length(unique(terms))) {
		termMods = terms[terms == unique(terms)[tt]];
		tmp2 = unlist(strsplit(names(termMods), '\\.'));
		mods = tmp2[!grepl('Term', tmp2)];
		termList[[tt]] = mods;
		names(termList)[tt] = unique(terms)[tt];
	}
	out = list(modsByTerm=termList, uniqueToMod=termList[unlist(lapply(termList, length)) == 1]);
	return(out);
}

# Depends on:
#	.getAllTermsWithGenes
.getAllTermsWithGenesAcrossModules = function (genes, modDAVID) {
	outlist = list();
	for (mod in 1:length(modDAVID)) {
		if (nrow(modDAVID[[mod]])==0 | all(is.na(modDAVID[[mod]]))) { 
			outlist[[mod]] = 'No terms'; 
		} else {
			outlist[[mod]] = .getAllTermsWithGenes(genes, modDAVID[[mod]]);
		}
	}
	names(outlist) = names(modDAVID);
	return(outlist);
}

# Depends on:
#	.getAllTermsWithGene
.getAllTermsWithGenes = function (genes, modDAVIDmod) {
	if (sum(genes %in% unlist(strsplit(modDAVIDmod$Genes,', ')))==0) {
		return('No genes in any term');
	}
	mat = .getAllTermsWithGene(genes[1], modDAVIDmod);
	for (g in 2:length(genes)) {
		mat = rbind(mat, .getAllTermsWithGene(genes[g], modDAVIDmod));
	}
	mat = mat[!duplicated(mat), ];
	mat = cbind(mat, qgenes=rep(NA, nrow(mat)));
	orderBy = c();
	for (tt in 1:nrow(mat)) {
		tgenes = strsplit(mat[tt, match('Genes', colnames(mat))], ', ')[[1]];
		g = genes[genes %in% tgenes];
		mat$qgenes[tt] = paste(g, collapse=', ');
		orderBy = c(orderBy, length(g));
	}
	mat = mat[order(orderBy, decreasing=T), ];
	return(mat);
}

# Depends on:
#	NA
.getAllTermsWithGene = function (gene, modDAVIDmod) {
	return(modDAVIDmod[apply(modDAVIDmod, 1, function(f) grepl(gene, f[match('Genes', names(f))])), ]);
}

# Depends on:
#	.getGenesAvgkMEInMod
#	.getTermGenes
.addkMEColToDAVID = function(modDAVIDmod, modkME) {
	modDAVIDmod = cbind(modDAVIDmod, avgkME=1:nrow(modDAVIDmod));
	for (r in 1:nrow(modDAVIDmod)) {
		modDAVIDmod[r, ]$avgkME = .getGenesAvgkMEInMod(.getTermGenes(modDAVIDmod,r), modkME);
	}
	return(modDAVIDmod);
}

# Depends on:
#	.addkMEColToDAVID
.addkMEColsToDAVIDList = function(modDAVID, modkME, order=T) {
	for (m in 1:length(modDAVID)) {
		if (nrow(modDAVID[[m]])==0 | names(modDAVID)[m]=='grey') {
			modDAVID[[m]] = cbind(modDAVID[[m]], avgkME=rep(NA,nrow(modDAVID[[m]])));
		} else {
			modDAVID[[m]] = .addkMEColToDAVID(modDAVID[[m]], modkME[[m]]);
			if (order) {
				modDAVID[[m]] = modDAVID[[m]][order(modDAVID[[m]]$avgkME, decreasing=T), ];
			}
		}
	}
	return(modDAVID);
}

# Depends on:
#	.getTermGenes
#	.getAvgNumForGenes
.addGeneAvgColToDAVID = function(modDAVIDmod, nums, colname='colname') {
	modDAVIDmod = cbind(modDAVIDmod, colname=1:nrow(modDAVIDmod));
	for (row in 1:nrow(modDAVIDmod)) {
		g = .getTermGenes(modDAVIDmod, row);
		modDAVIDmod[row, ncol(modDAVIDmod)] = .getAvgNumForGenes(g, nums);
	}
	names(modDAVIDmod)[ncol(modDAVIDmod)] = colname;
	return(modDAVIDmod)
}

# Depends on:
#	NA
.getTermGenes = function(modDAVIDmod,term) {
	if (is.numeric(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[term], ', ')));
	} else if (is.character(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[modDAVIDmod$Term==term], ', ')));
	} else {
		stop()
	}
}

# Depends on:
#	.checkNamesOfModListElements
#	.filterDAVIDChart
.filterModDAVIDList = function(modDAVID, column, value, lessThan=F, verbose=T) {
	.checkNamesOfModListElements(modDAVID);
	if (is.character(column)) {
		columnName = column;
		column = match(columnName, names(modDAVID[[1]]));
	}
	nrows = unlist(lapply(modDAVID, nrow));
	goodColnames = names(modDAVID[nrows>1][[1]]);
	if (verbose) {
		if (is.character(value)) {
			cat('Filtering $', goodColnames[column], ' to rows that include ', value, '\n', sep='');
		} else if (is.numeric(value)) {
			if (lessThan) {
				cat('Filtering $', goodColnames[column], ' to values < ', value, '\n', sep='');
			} else {
				cat('Filtering $', goodColnames[column], ' to values > ', value, '\n', sep='');
			}
		} else {
			stop('Parameter "value" must be numeric or character');
		}
	} 
	modDAVID_filtered = modDAVID;
	for (m in 1:length(modDAVID_filtered)) {
		if (all(is.na(modDAVID_filtered[[m]]))) {
			modDAVID_filtered[[m]] = modDAVID_filtered[[m]];
			next;
		}
		modDAVID_filtered[[m]] = .filterDAVIDChart(chart=modDAVID_filtered[[m]], column=column, value=value, lessThan=lessThan, verbose=F);
	}
	return(modDAVID_filtered);
}

# Depends on:
#	NA
.checkNamesOfModListElements = function(modList) {
	x = lapply(modList, names);
	x_lengths = unlist(lapply(x, length));
	if ( (!(all(x_lengths[1]==x_lengths))) | length(unique(x_lengths))!=1 ) {
		stop('List elements have different numbers of columns');
	}
	for (element in 2:length(x)) {
		if (!(all(x[[1]]==x[[element]]))) {
			stop('Names in element ', element, ' do not match');
		}
	}
	return(TRUE);
}

# Depends on:
#	NA
.filterDAVIDChart = function(chart, column, value, lessThan=F, verbose=T) {
	if (!is.data.frame(chart)) {
		stop('Input must be a data frame');
	}
	if (is.character(column)) {
		columnName = column;
		column = match(columnName, names(chart));#print(column)
	}
	if (is.character(value)) {
		if (!(is.character(chart[, column]))) {
			stop('Cannot filter ', mode(chart[, column]), ' column by a character value');
		}
		if (verbose) {cat('Filtering $', names(chart)[column], ' to rows that include ', value, '\n', sep='')};
		chart_filtered = chart[grepl(value, chart[, column]), ];
	} else if (is.numeric(value)) {
		if (!(is.numeric(chart[, column]))) {
			stop('Cannot filter ', mode(chart[, column]), ' column by a numeric value');
		}
		if (lessThan) {
			if (verbose) {cat('Filtering $', names(chart)[column], ' to values < ', value, '\n', sep='')};
			chart_filtered = chart[chart[, column]<value, ];
		} else {
			if (verbose) {cat('Filtering $', names(chart)[column], ' to values > ', value, '\n', sep='')};
			chart_filtered = chart[chart[, column]>value, ];
		}
	} else {
		stop('Parameter "value" must be numeric or character');
	}
	return(chart_filtered);
}

#######################################################################################################################################
######################### compute DAVID results
###	.startAndUploadBgAndModListsToDAVID2
###	.getModChartsFromDAVID 
###	.getDAVIDchartForGeneList
###	.getDAVIDchartsForModuleSubsetsInOtherNetwork
#######################################################################################################################################

# Depends on:
#	RDAVIDWebService
.startAndUploadBgAndModListsToDAVID2 = function (colors, login, idType, ids) {
	mods = names(table(colors));
	david = DAVIDWebService$new(email=login);
	addList(david, inputIds=ids, idType=idType, listName='BG', listType='Background');
	for (m in 1:length(mods)) {
		print(mods[m]);
		modGenes = ids[colors==mods[m]];
		addList(david, inputIds=modGenes, idType=idType, listName=mods[m], listType='Gene');
	}
	setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
	setCurrentGeneListPosition(david, 1);
	return(david);
}

# Depends on:
#	RDAVIDWebService
.getModChartsFromDAVID = function(david, setBgPosition=NULL, factorsToChars=T, verbose=T) {
	if (!is.connected(david)) {
		stop('DAVID object not connected to web server');
	}
	if (is.numeric(setBgPosition)) {
		setCurrentBackgroundPosition(david, setBgPosition);
	}
	cat('Using ', getBackgroundListNames(david)[getCurrentBackgroundListPosition(david)], ' as background\n', sep='');
	modDAVID = list();
	geneListNames = getGeneListNames(david);
	cat('Getting charts for module lists... ', sep='');
	for (m in 1:length(geneListNames)) {
		setCurrentGeneListPosition(david, m);
		if (verbose) {cat(getCurrentGeneListPosition(david), ':', geneListNames[getCurrentGeneListPosition(david)], ' ', sep='')};
		modDAVID[[m]] = getFunctionalAnnotationChart(david);
		if (factorsToChars) {
			checkFactor = which(lapply(modDAVID[[m]], class)=='factor');
			if (length(checkFactor)>0) {
				modDAVID[[m]][, checkFactor] = as.character(modDAVID[[m]][, checkFactor]);
			}
		}
		names(modDAVID)[m] = geneListNames[m];
	}
	# check for modules with 0 hits and fix colnames
	nrows = unlist(lapply(modDAVID, nrow));
	goodColnames = names(modDAVID[nrows>0][[1]]);
	if (goodColnames[4]=='X.') {
		goodColnames[4] = 'Percent';
	}
	for (m in 1:length(modDAVID)) {
		if (nrows[m]==0) {
			modDAVID[[m]] = as.data.frame(matrix(nrow=1, ncol=length(goodColnames)));
			names(modDAVID[[m]]) = goodColnames;
		} else {
			names(modDAVID[[m]]) = goodColnames;
		}
	}
	return(modDAVID);
}

# Depends on:
#	RDAVIDWebService
.getDAVIDchartForGeneList = function (genes, background, idType='ENSEMBL_GENE_ID', login='ahilliar@stanford.edu') {
	david = DAVIDWebService$new(email=login);
	addList(david, inputIds=background, idType=idType, listName='BG', listType='Background');
	addList(david, inputIds=genes, idType=idType, listName='genes', listType='Gene');
	setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
	setCurrentGeneListPosition(david, 1);
	chart = getFunctionalAnnotationChart(david);
	checkFactor = which(lapply(chart, class)=='factor');
	if (length(checkFactor)>0) {
		chart[, checkFactor] = as.character(chart[, checkFactor]);
	}
	return(chart);	
}

# Depends on:
#	.getDAVIDchartForGeneList
# right now uploads background everytime, need to fix
.getDAVIDchartsForModuleSubsetsInOtherNetwork = function (moduleRef, modulesTest, modGenesRef, modGenesTest, idType='ENSEMBL_GENE_ID', login='ahilliar@stanford.edu', bg=c('net', 'refModule', 'refAndTestModules'), thresh=10) {
	if (any(sort(unlist(modGenesRef)) != sort(unlist(modGenesTest)))) {stop('Module lists must contain all exact same genes')}
	if (is.vector(bg) | !all(bg %in% c('net', 'refModule', 'refAndTestModules'))) {
		warnings('bg must be either "net", "refModule", or "refAndTestModules". Using "net"');
		bg = 'net';
	} else if (bg == 'net') {
		bg = unlist(modGenesRef);
	} else if (bg == 'refModule') {
		bg = modGenesRef[[moduleRef]];
	} else if (bg == 'refAndTestModules') {
		#bg = unique(c(modGenesRef[[moduleRef]], ))
	} 
	david = DAVIDWebService$new(email=login);
	cat('Uploading background...\n');
	addList(david, inputIds=bg, idType=idType, listName='Background', listType='Background');
	for (m in 1:length(modulesTest)) {
		genes = intersect(modGenesRef[[moduleRef]], modGenesTest[[modulesTest[m]]]);
		print(genes)
		addList(david, inputIds=genes, idType=idType, listName=modulesTest[m], listType='Gene');
	}
	print(david);
	print(getFunctionalAnnotationChart(david))
	#out = list();
	# for (m in 1:length(modulesTest)) {
		# # upload 
		# cat(modulesTest[m], '\n');
		
		# out[[m]] = .getDAVIDchartForGeneList(genes=refgenes[refgenes %in% testgenes], background=bg, idType=idType, login=login);
		# names(out)[m] = modulesTest[m];
	# }
	# return(out);
}


#cl = getClusterReport(d.hs,type='Term', overlap=4L, initialSeed=4L, finalSeed=4L, linkage=.5, kappa=75L)



#######################################################################################################################################
######################### module preservation
###	.computeModulePreservation
#######################################################################################################################################

# Depends on:
#	WGCNA
.computeModulePreservation = function (data1, data2, colors1, colors2, labels=c('1', '2'), networkType='signed', ranksOnly=T, nPermutations=100, verbose=3) {
	multiExpr = list();
	multiExpr[[1]] = list(data=data1);
	multiExpr[[2]] = list(data=data2);
	names(multiExpr) = labels;
	colorList=list();
	colorList[[1]] = colors1;
	colorList[[2]] = colors2;
	names(colorList) = labels;
	
	if (ranksOnly) {
		mp = modulePreservation(multiExpr, colorList, networkType=networkType, nPermutations=0, maxModuleSize=max(table(colors1)), 
								verbose=verbose, savePermutedStatistics=F);
		#tmp = cbind(get(intra,pos=-1), get(inter));
		tmp = cbind(mp[[1]][[1]][[3]][[2]], mp[[1]][[1]][[4]][[2]])
		tmp = tmp[!(rownames(tmp)=='gold'), ];
		tmp2 = apply(-tmp, 2, rank);
		# con - cor.kIM, cor.kME, cor.cor
		medRankCon = apply( tmp2[ , c(8, 9, 11) ], 1, median);
		# density - propVarExplained, meanSignAwareKME, meanSignAwareCorDat, meanAdj
		medRankDen = apply(tmp2[ , c(1, 2, 4, 5) ], 1, median);
		medRankPres = (medRankCon + medRankDen) / 2;
		medRanks = cbind(medRankPres, medRankCon, medRankDen);
		medRanks = medRanks[order(medRanks[,1]), ];
		return(list(medRanks=medRanks, mp=mp));
	} else {
		mp = modulePreservation(multiExpr, colorList, networkType=networkType, maxModuleSize=max(table(colors1)), nPermutations=nPermutations, 
								verbose=verbose, savePermutedStatistics=F);
		return(mp);
	}
}



