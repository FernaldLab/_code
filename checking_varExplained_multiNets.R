setwd('/Volumes/fishstudies/_mammalian_RNAseq/hsa_.33zeros_.5rpkm_log2_nohs.br.M.4_defaultPreProc');

rm(list=ls()); options(stringsAsFactors=F);
source('/Volumes/fishstudies/_code/exploreNetwork.R');
library(WGCNA); allowWGCNAThreads();

toLoad = list.files('_endRuns');
toLoad = gsub('dendro-block1.jpg', '', toLoad);
toLoad = gsub('^_', '', toLoad);

# load nets and data
nets = list();
dats = list();
for (name in toLoad) {
	load(file=paste(name, 'NET.RData', sep=''));
	load(file=paste(name, 'DATA.RData', sep=''));
	subname = gsub('get(stem)_signed_p18_', '', name, fixed=T);
	print(subname)
	nets[[subname]] = net;
	dats[[subname]] = DATA;
	rm(net, DATA);
}; rm(name, subname);

# build traits table and compute ME correlations
# # DATA = dats[[1]];
# # traits = data.frame(br=as.numeric(grepl('br', rownames(DATA))));
# # rownames(traits) = rownames(DATA);
# # traits = cbind(traits, cb=as.numeric(grepl('cb', rownames(DATA))));
# # traits = cbind(traits, br_cb=as.numeric(grepl('br|cb', rownames(DATA))));
# # traits = cbind(traits, ht=as.numeric(grepl('ht', rownames(DATA))));
# # traits = cbind(traits, kd=as.numeric(grepl('kd', rownames(DATA))));
# # traits = cbind(traits, lv=as.numeric(grepl('lv', rownames(DATA))));
# # traits = cbind(traits, ts=as.numeric(grepl('ts', rownames(DATA))));
# # ages = read.table('../Mammalian Samples.txt',header=T,sep='\t');
# # traits = cbind(traits, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATA)), 
												  # # gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);

###################################################################################################
# compute %variance explained by MEs and ME-trait correlations
vexp = list();
#tcors = list();
for (i in 1:length(nets)) {
	tmp = moduleEigengenes(dats[[i]], nets[[i]]$colors);
	tmp2 = tmp$varExplained;
	#names(tmp2) = colnames(tmp$eigengenes);
	vexp[[names(nets)[i]]] = as.numeric(tmp2);
	names(vexp[[names(nets)[i]]]) = colnames(tmp$eigengenes);
	#traitCors.hs = exn.computeAndPlotMETraitCors(traits.hs, MEs.hs, main=stem);
	
	
}; rm(i, tmp, tmp2);

# boxplot of avg ME varExplained
verboseBoxplot(unlist(vexp), gsub('run[0-9]*$', '', names(unlist(vexp))), 
	           cex.axis=.8, cex.lab=1, 
	           ylab='ME varExplained', xlab='', 
	           col='lightgrey', frame.plot=F
	           );
	           
# scatterplot of median varExplained against number of modules in network   
verboseScatterplot(sapply(vexp, median), sapply(vexp, length), 
				   xlab='median varExplained', ylab='number of modules', 
				   frame.plot=F, abline=T, 
				   type='n'); 
text(sapply(vexp, median), sapply(vexp, length), gsub('run[0-9]*$','',names(vexp)), cex=.7);

# histograms of varExplained with median line
par(mfrow=c(3,5));
for (i in 1:length(vexp)) {
	hist(vexp[[i]], main=names(vexp)[i], 
		 xlim=c(.4,1), xlab='ME varExplained', 
		 breaks=length(vexp[[i]]), col='lightgrey', border='grey'
		 );
	abline(v=median(vexp[[i]]), col='blue');
}; rm(i);

###################################################################################################
# compute robustness and separability

mps = list();
for (i in 1:length(nets)) {
	tmp = exn.computeModulePreservation(data1=dats[[i]], data2=dats[[i]], 
										colors1=nets[[i]]$colors, colors2=nets[[i]]$colors, 
										#labels=c('1','1a'), 
										ranksOnly=F, nPermutations=100
										);
	mps[[names(nets)[i]]] = tmp;
}; rm(i, tmp);

seps = list();
quals = list();
zvexp = list();
for (run in 1:length(mps)) {
	s = mps[[run]]$referenceSeparability$Z$ref.1$inColumnsAlsoPresentIn.2;
	seps[[names(mps)[run]]] = s[rownames(s)!='gold', ];
	q = mps[[run]]$quality$Z$ref.1$inColumnsAlsoPresentIn.2
	quals[[names(mps)[run]]] = q[rownames(q)!='gold', ];
	z = mps[[run]]$preservation$Z$ref.1$inColumnsAlsoPresentIn.2[, c(1,5)];
	zvexp[[names(mps)[run]]] = z[rownames(z)!='gold', ];
}; rm(run, s, q, z);

# boxplot of avg separability
verboseBoxplot(unlist(lapply(seps, function(f) f[,2])), 
			   gsub('run[0-9]*$', '', names(unlist(lapply(seps, function(f) f[,2])))), 
	           cex.axis=.8, cex.lab=1, 
	           ylab='Z.separability', xlab='', 
	           col='lightgrey', frame.plot=F
	           );
abline(h=10, col='red', lty='dashed');

# boxplot of avg quality
verboseBoxplot(unlist(lapply(quals, function(f) f[,2])), 
			   gsub('run[0-9]*$', '', names(unlist(lapply(quals, function(f) f[,2])))), 
	           cex.axis=.8, cex.lab=1, 
	           ylab='Z.quality', xlab='', 
	           col='lightgrey', frame.plot=F
	           );
abline(h=10, col='red', lty='dashed');

# scatterplot of median separability against number of modules in network   
verboseScatterplot(sapply(seps, function(f) median(f[,2])), sapply(seps, nrow), 
				   xlab='median separability', ylab='number of modules', 
				   frame.plot=F, abline=T, 
				   type='n'); 
text(sapply(seps, function(f) median(f[,2])), sapply(seps, nrow), gsub('run[0-9]*$','',names(seps)), cex=.7);

# scatterplot of median quality against number of modules in network   
verboseScatterplot(sapply(quals, function(f) median(f[,2])), sapply(quals, nrow), 
				   xlab='median quality', ylab='number of modules', 
				   frame.plot=F, abline=T, 
				   type='n'); 
text(sapply(quals, function(f) median(f[,2])), sapply(quals, nrow), gsub('run[0-9]*$','',names(quals)), cex=.7);

###################################################################################################

# r = list();
# if (length(unique(c(length(seps),length(quals),length(vexp)))) == 1) {
	# for (i in 1:length(seps)) {
		# mat = as.data.frame(matrix(ncol=3, nrow=nrow(seps[[i]]), 
								   # dimnames=list(rownames(seps[[i]]), c('varExp','Z.sep','Z.qual'))
								   # )
							# );
		# mat[, 1] = (vexp[[i]]) #/ max(vexp[[i]]);
		# mat[, 2] = (seps[[i]][,2]) #/ max(seps[[i]][,2]);
		# mat[, 3] = (quals[[i]][,2]) #/ max(quals[[i]][,2]);
		# r[[names(seps)[i]]] = mat;
	# }; rm(i);
# }

r = list();
if (length(unique(c(length(seps),length(quals),length(vexp)))) == 1) {
	for (i in 1:length(seps)) {
		mat = as.data.frame(matrix(ncol=3, nrow=nrow(seps[[i]]), 
								   dimnames=list(rownames(seps[[i]]), c('Z.varExp','Z.sep','Z.qual'))
								   )
							);
		mat[, 1] = (zvexp[[i]][,2]) #/ max(vexp[[i]]);
		mat[, 2] = (seps[[i]][,2]) #/ max(seps[[i]][,2]);
		mat[, 3] = (quals[[i]][,2]) #/ max(quals[[i]][,2]);
		r[[names(seps)[i]]] = mat;
	}; rm(i);
}

rtab = t(sapply(r, function(f) apply(f, 2, median)));
rtab = cbind(rtab, modnum=sapply(r,nrow));
rtabfilt = rtab[!apply(rtab,1,function(f) any(f<10)),];
####################

NEThs = nets$ds2_mm20_mch0.15run9;
DATAhs = dats$ds2_mm20_mch0.15run9;
STATS.hs = r$ds2_mm20_mch0.15run9;
smods = rownames(STATS.hs)[STATS.hs$Z.sep < 10];

.getMostSimilarModule = function (module, MEs) {
	MEcor = cor(MEs);
	tmp = MEcor[gsub('ME','',rownames(MEcor))==module,];
	maxtmp = max(tmp[names(tmp)!=paste('ME',module,sep='')]);
	return( names(tmp)[tmp==maxtmp] );
} 

.getClosestModules = function (modules, MEs) {
	closest = c();
	for ( mod in modules ) {
		closest = c(closest, .getMostSimilarModule(mod, MEs));
	}
	closest = gsub('ME', '', closest);
	closestModulesTable = cbind(modules, closest);
	return(closestModulesTable);
}

.selectiveMergeModules = function (modulesToMerge, colors, DATA) {
	newcolors = colors;
	MEs = moduleEigengenes(DATA, colors)$eigengenes;
	closestModulesTable = .getClosestModules(modulesToMerge, MEs);
	tmp = closestModulesTable;
	for (row in 1:nrow(tmp)) {
		refmod = tmp[row, 1];
		clsmod = tmp[row, 2];
		newcolors[newcolors==refmod] = clsmod;
		tmp[tmp[,2]==refmod, 2] = clsmod;
		print('REFMOD:');
		print(refmod)
		print(cbind(closestModulesTable, tmp));
		print('--------------------------------------')
	}
	return(newcolors);
}











newcolors = NEThs$colors;
for ( mod in smods ) {
	newcolors[newcolors==mod] = gsub('ME', '', .getMostSimilarModule(mod, NEThs$MEs));
}; rm(mod);

# get separability and quality for merged colors
newmp = exn.computeModulePreservation(data1=DATAhs, data2=DATAhs, 
										colors1=newcolors, colors2=newcolors, 
										#labels=c('1','1a'), 
										ranksOnly=F, nPermutations=100
										);
										
tmp = newmp$referenceSeparability$Z$ref.1$inColumnsAlsoPresentIn.2;
newsep = tmp[rownames(tmp)!='gold', ];
tmp = newmp$quality$Z$ref.1$inColumnsAlsoPresentIn.2
newqual = tmp[rownames(tmp)!='gold', ];
rm(tmp);

newMEs = moduleEigengenes(DATAhs, newcolors);
newvexp = as.numeric(newMEs$varExplained);
newSTATS = data.frame(varExp=newvexp, Z.sep=newsep[,2], Z.qual=newqual[,2], row.names=rownames(newsep));

par(mfrow=c(1,3));
fact = c(rep('new',nrow(newSTATS)), rep('old',nrow(STATS.hs)));
for (i in 1:3) {
	verboseBoxplot(c(newSTATS[,i], STATS.hs[,i]), fact, col='grey', frame.plot=F, xlab='', ylab=names(newSTATS)[i]);
}; rm(i)


STATS.hs[rownames(STATS.hs) %in% smods_closest, ];
newSTATS[rownames(newSTATS) %in% smods_closest, ];






newMEs0 = newMEs;
newMEs = orderMEs(newMEs0)$eigengenes;
# # # exn.getNetworkBasics(NEThs, DATAhs, '.hs');
exn.plotDendroAndColors(dendro.hs, colors.hs, block=1, blockGenes=blockGenes.hs, hang=.05);
#exn.plotDendroAndColors(NEThs$dendrograms[[1]], NEThs$colors, block=1, blockGenes=NEThs$blockGenes, hang=.05);
plotDendroAndColors(NEThs$dendrograms[[1]], 
					data.frame(NEThs$colors, newcolors), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs', 'merged'), ylab="Distance (1-TO)", main = '',
					hang=.05, colorHeight=.5, cex.axis=.8
					);



MEnet.hs = exn.plotEigengeneNetworks2(NEThs$MEs, returnCors=T);
newMEnet = exn.plotEigengeneNetworks2(newMEs, returnCors=T);

# build traits table and compute ME correlations
traits.hs = data.frame(br=as.numeric(grepl('br', rownames(DATAhs))));
rownames(traits.hs) = rownames(DATAhs);
traits.hs = cbind(traits.hs, cb=as.numeric(grepl('cb', rownames(DATAhs))));
traits.hs = cbind(traits.hs, br_cb=as.numeric(grepl('br|cb', rownames(DATAhs))));
traits.hs = cbind(traits.hs, ht=as.numeric(grepl('ht', rownames(DATAhs))));
traits.hs = cbind(traits.hs, kd=as.numeric(grepl('kd', rownames(DATAhs))));
traits.hs = cbind(traits.hs, lv=as.numeric(grepl('lv', rownames(DATAhs))));
traits.hs = cbind(traits.hs, ts=as.numeric(grepl('ts', rownames(DATAhs))));
ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits.hs = cbind(traits.hs, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATAhs)), 
												  gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);
												  
traitCor.hs = exn.computeAndPlotMETraitCors(traits.hs, NEThs$MEs, main='old');		
traitCor.hsnew = exn.computeAndPlotMETraitCors(traits.hs, newMEs$eigengenes, main='new');		
kMEnew=exn.computekME(DATAhs,newMEs)

GS.hs = exn.computeGS(traits.hs, DATAhs);	


plotDendroAndColors(NEThs$dendrograms[[1]], 
					data.frame(NEThs$colors, 
							   newcolors, 
							   numbers2colors(GS.hs$GS.br), 
							   numbers2colors(GS.hs$GS.cb), 
							   numbers2colors(GS.hs$GS.br_cb), 
							   numbers2colors(GS.hs$GS.ht),
							   numbers2colors(GS.hs$GS.kd),
							   numbers2colors(GS.hs$GS.lv),
							   numbers2colors(GS.hs$GS.ts)), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('old','new','GS.br', 'GS.cb','GS.br_cb','GS.ht','GS.kd','GS.lv','GS.ts'), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);		

traitCorTable = traitCor.hs$cor;	
GSkMEcors = as.data.frame(matrix(nrow=nrow(traitCorTable), ncol=ncol(traitCorTable), dimnames=dimnames(traitCorTable)));				
for ( tr in 1:ncol(traitCorTable) ) {
	tmp = exn.plotAllModsGSkME(NEThs$colors,colnames(traitCorTable)[tr],GS.hs,kME.hs,c(3,10),plot=F);
	tmp = tmp[match(gsub('ME','',rownames(traitCorTable)), names(tmp))];
	GSkMEcors[, tr] = tmp;
}; rm(tr, tmp)	
diag(cor(traitCor.hs$cor,GSkMEcors))

traitCorTable = traitCor.hsnew$cor;	
newGSkMEcors = as.data.frame(matrix(nrow=nrow(traitCorTable), ncol=ncol(traitCorTable), dimnames=dimnames(traitCorTable)));				
for ( tr in 1:ncol(traitCorTable) ) {
	tmp = exn.plotAllModsGSkME(newcolors,colnames(traitCorTable)[tr],GS.hs,kMEnew$all,c(3,10),plot=F);
	tmp = tmp[match(gsub('ME','',rownames(traitCorTable)), names(tmp))];
	newGSkMEcors[, tr] = tmp;
}; rm(tr, tmp)		  
diag(cor(traitCor.hsnew$cor,newGSkMEcors))