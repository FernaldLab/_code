rm(list=ls()); options(stringsAsFactors=F); 
library(WGCNA); allowWGCNAThreads();
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');

#load('all_primate_brainVar_filtering_DATA.all.higherBrainVar.RData');
load('all_primate_brainVar_filtering_speciesDATA.RData')

#keep = grepl('br|cb',names(DATA.all.higherBrainVar));
#DATA = DATA.all.higherBrainVar[, keep];

#keep = grepl('br', names(DATA.Common)) & grepl('M', names(DATA.Common));
#DATA = DATA.Common[, keep];

DATA = DATA.Common;
remove = c('ggo.cb.M.1');
keepMe = !(names(DATA) %in% remove);
DATA = DATA[, keepMe];

zthresh = floor(ncol(DATA)/3);
out = preProcess(datIN=DATA, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
#####################################################################################
				 
DATA = as.data.frame(t(out$data_Qnorm));
BLOCKSIZE = 2000;
TYPE = 'signed';
sft = pickSoftThreshold(DATA, networkType=TYPE, verbose=3, blockSize=BLOCKSIZE);

POWER = 14;
k = softConnectivity(DATA, type=TYPE, power=POWER, blockSize=BLOCKSIZE);
par(mfrow=c(1,2));
scaleFreePlot(k); hist(k,col='grey',border='darkgrey');

DS = 2;
MM = 10;
MCH = 0.15;
NET = blockwiseModules(datExpr=DATA, maxBlockSize=BLOCKSIZE, networkType=TYPE, power=POWER, deepSplit=DS, minModuleSize=MM, mergeCutHeight=MCH, verbose=3);

dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
library(RDAVIDWebService); library(biomaRt);
exn.plotDendroAndColors(dendro, colors, block=block, blockGenes=blockGenes);
exn.plotEigengeneNetworks2(MEs)

kME = exn.computekME(DATA, MEs)$all;
mod.genes=exn.getModuleGenes(DATA, colors);

.getModGenesRankedBykME = function(module_names,colors,kME) {
	outList = list();
	for (m in 1:length(module_names)) {
		outList[[m]] = exn.getModulekME(module_names[m],colors,kME)
	}
	names(outList) = module_names;
	return(outList);
}
modkMEs = .getModGenesRankedBykME(names(table(colors)),colors,kME);





















length(table(colors))
MFROW = c(3,7);
factors = unlist(strsplit(rownames(DATA),'\\.'))[seq(2,length(unlist(strsplit(rownames(DATA),'\\.'))),4)];
par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(factors), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=1);
}

factors2 = factors; factors2[grepl('br|cb', factors2)] = 'br/cb';
par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(factors2), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=1);
}

factors3 = factors2; factors3[grepl('ht|kd|lv', factors3)] = 'ht/kd/lv';
par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(factors3), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=1);
}





tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
primates_hiv1_interactions = primates_hiv1_interactions[primates_hiv1_interactions==1]
checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA));
hiv = checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA))$pvals;hiv