##### assumes raw fpkm files were processed using processCufflinksFPKM.py #####

##############################################################
###### read in data and get genes common to all samples ######
##############################################################
options('stringsAsFactors'=F);
# read in fpkm
rm(list=ls());
#setwd('~/Documents/_Fernald_lab/_hybridRNAseq/_newFPKM');
setwd('~/Documents/_Fernald_lab/_hybridRNAseq/_masterFPKM');
dirs = list.files();
#files = paste(dirs,'/genes.fpkm_tracking.noZeros',sep='');
files = paste(dirs,'/isoforms.fpkm_tracking.noZeros',sep='');

DAT0 = list();
for (file in 1:length(files)) {
	DAT0[[file]] = read.table(files[file], header=T, row.names=1);
	cat(files[file], '\n');
	print(head(DAT0[[file]]));
	cat('\n');
}; rm(file);

#names(DAT0) = gsub('Cufflinkswithannotation', '', dirs);
names(DAT0) = dirs;
lapply(DAT0,dim);

###
# insert FPKM_status filter
for (set in 1:length(DAT0)) {
	cat(names(DAT0)[set], ':\n', sep='');
	print(dim(DAT0[[set]]));
	print(sum(DAT0[[set]]$FPKM_status!='OK'));
	DAT0[[set]] = DAT0[[set]][DAT0[[set]]$FPKM_status=='OK',]
	print(dim(DAT0[[set]]));
}; rm(set);


# possibly some kind of confidence interval filter
###

# get genes common to all datasets
genes = rownames(DAT0[[1]]);
for (set in 2:length(DAT0)) {
	genes = genes[genes %in% rownames(DAT0[[set]])];
}; rm(set);


DAT = data.frame(matrix(nrow=length(genes), ncol=length(DAT0)));
rownames(DAT) = genes;
for (set in 1:length(DAT0)) {
	thisSet = DAT0[[set]];
	thisSet = thisSet[match(genes, rownames(thisSet)), ];
	DAT[, set] = thisSet$FPKM;
	names(DAT)[set] = names(DAT0)[set];
}; rm(set, thisSet)

##############################################################
###### remove outliers and normalize #########################
##############################################################

source('../../_code/preProcATH-for_web_noVSN.R');
x.probes=removeOutlierProbesIterate(DAT,2);
x.NA=removeTooManyNAs(x.probes$dataClean);
x.samples=outlierSamplesIterate(x.NA$dataClean);
x.qnorm=normalize.quantiles(as.matrix(x.samples$dataClean));
DATA = as.data.frame(x.qnorm);
dimnames(DATA) = dimnames(x.samples$dataClean);
DATA = as.data.frame(t(DATA));
##############################################################
###### test parameters for network construction ##############
##############################################################
library(WGCNA); allowWGCNAThreads();
#sft = pickSoftThreshold(DATA,networkType='signed',blockSize=3500);
sft = pickSoftThreshold(DATA,networkType='signed',blockSize=5000,verbose=3);
plot(sft$fitIndices$Power,sft$fitIndices$SFT.R.sq,type='n');text(sft$fitIndices$Power,sft$fitIndices$SFT.R.sq, sft$fitIndices$Power);

#k15 = softConnectivity(DATA,type='signed',power=15,blockSize=3500);
#names(k15) = names(DATA);
k14 = softConnectivity(DATA,type='signed',power=14,blockSize=5000,verbose=3);
names(k14) = names(DATA);
k16 = softConnectivity(DATA,type='signed',power=16,blockSize=5000,verbose=3);
names(k16) = names(DATA);
toKeep = k16 > quantile(k16, .2);
DATA2=DATA[,toKeep];
k16.2 = softConnectivity(DATA2,type='signed',power=16,blockSize=5000,verbose=3);
names(k16.2) = names(DATA2);
# toKeep = names(k15)[k15>quantile(k15,.25)];
# DATAk15 = DATA[, names(DATA) %in% toKeep];
# sft15 = pickSoftThreshold(DATAk15,networkType='signed',blockSize=3500);
# k15.25=softConnectivity(DATAk15,type='signed',power=15,blockSize=3500);
# names(k15.25)=names(DATAk15);
# par(mfrow=c(2,2));hist(k15);scaleFreePlot(k15);hist(k15.25);scaleFreePlot(k15.25)

##############################################################
###### construct network and relate to traits ################
##############################################################

net=blockwiseModules(DATA,maxBlockSize=6000,power=16,networkType='signed',minModuleSize=10,mergeCutHeight=.2,verbose=3);
table(net$colors);

b=1;
plotDendroAndColors(net$dendrograms[[b]], net$colors[net$blockGenes[[b]]],"Module colors",rowText=net$colors[net$blockGenes[[b]]],
 					main = "",dendroLabels = FALSE, hang = 0.03, 
 					addGuide = TRUE, guideHang = 0.05);


plotDendroAndColors(net$dendrograms[[1]], net$colors,"Module colors",rowText=net$colors,
					main = "",dendroLabels = FALSE, hang = 0.03, 
					addGuide = TRUE, guideHang = 0.05);
MET=orderMEs(net$MEs);
plotEigengeneNetworks(MET,setLabels='',marDendro=c(0,4,1,2),marHeatmap=c(3,4,1,2),heatmapColors=blueWhiteRed(50));

# adj=adjacency(DATAk15,power=15,type='signed');
# TOMdist=TOMdist(adj,TOMType='signed');
# rm(adj);
# rownames(TOMdist)=rownames(adj); colnames(TOMdist)=colnames(adj);
# clust=flashClust(as.dist(TOMdist),method='average');

# modLabels=cutreeDynamic(clust,distM=TOMdist,minClusterSize=10);
# MEs=moduleEigengenes(DATAk15,modLabels)$eigengenes;
# MEsMerge=mergeCloseModules(DATAk15,colors=modLabels,cutHeight=0.2,MEs=MEs,relabel=T);
# mergedColors=labels2colors(MEsMerge$colors);

# plotDendroAndColors(clust, mergedColors,"Module colors",rowText=mergedColors,
					# main = "",dendroLabels = FALSE, hang = 0.03, 
					# addGuide = TRUE, guideHang = 0.05);
# table(mergedColors);

# moduleColors=mergedColors;
# nSamples = nrow(DATAk15);
# MEs0=moduleEigengenes(DATAk15,moduleColors)$eigengenes
# MEs=orderMEs(MEs0);
# # MEs should be ordered
# plotEigengeneNetworks(MEs,setLabels='',marDendro=c(0,4,1,2),marHeatmap=c(3,4,1,2));

# modNames=substring(names(MEs),3);


kME = as.data.frame(cor(DATAk15, MEs, use='p'));
names(kME)=paste('k',names(kME),sep='');
kMEpval = as.data.frame(corPvalueFisher(as.matrix(kME), nrow(DATAk15)));
names(kMEpval)=paste('p.',names(kMEpval),sep='');


# purplekME = kME[moduleColors=='purple',match('purple',gsub('kME','',names(kME)))];
# names(purplekME)=rownames(kME)[moduleColors=='purple'];


##############################################################
###### get genes from each module ############################
##############################################################
modGenes=list();
for (mod in 1:length(modNames)) {
	modGenes[[mod]] = names(DATAk15)[moduleColors==modNames[mod]];
	names(modGenes)[mod]=modNames[mod];
}; rm(mod);

modGeneskMEsorted=list();
for (mod in 1:length(modGenes)) {
	moduleGenes=modGenes[[mod]];
	moduleColumn=match(names(modGenes)[mod], gsub('kME','',names(kME)));
	modulekME=kME[match(moduleGenes,rownames(kME)),moduleColumn];
	modGeneskMEsorted[[mod]]=modulekME;
	names(modGeneskMEsorted[[mod]])=moduleGenes;
	names(modGeneskMEsorted)[mod]=names(modGenes)[mod];
	modGeneskMEsorted[[mod]]=sort(modGeneskMEsorted[[mod]],decreasing=T);
}; rm(mod,moduleGenes,modulekME,moduleColumn);

annos=read.table('../../_broadftp/cichlid_geneNames_MZnoNONE.txt',header=F,sep='\t',stringsAsFactors=F,quote="");
#annos[,1]=as.character(annos[,1]);
annos=annos[grepl('^mz',annos[,1]),];
#annos=annos[,1:3];

#annos[gsub('mz.gene.','',annos[,1]) %in% gsub('.[0-9]$', '', gsub('mz.mrna.','',names(DATA))), ]


allGenes=names(DATA);
allGenesAnnos=annos[gsub('mz.gene.','',annos[,1]) %in% gsub('.[0-9]$', '', gsub('mz.mrna.','',allGenes)), ];
allGenesSym=toupper(as.character(allGenesAnnos[,3]));
allGenesSym=strsplit(gsub('\\s','',allGenesSym), '\\(');
temp=c();
for(gene in 1:length(allGenesSym)) {
	if(length(allGenesSym[[gene]])==0) {
		next;
	} else {
		temp = c(temp, allGenesSym[[gene]][1]);
	}
}; rm(gene);
allGenesSym=temp; rm(temp);
allGenesIDs=IDs[names(IDs) %in% allGenesSym];
IDvec = c();
for (gene in 1:length(allGenesIDs)) {
	IDvec=c(IDvec, allGenesIDs[[gene]]);
}; allGenesIDs=IDvec; rm(gene,IDvec);
write.table(allGenesIDs, file='isoNet.allGenesEntrez.txt', quote=F, row.names=F, col.names=F)

modGenesAnnos=list();
for (mod in 1:length(modGenes)) {
	#modGenesAnnos[[mod]] = annos[annos[,1] %in% modGenes[[mod]],];
	modGenesAnnos[[mod]] = annos[gsub('mz.gene.','',annos[,1]) %in% gsub('.[0-9]$', '', gsub('mz.mrna.','',modGenes[[mod]])), ];
	names(modGenesAnnos)[mod] = names(modGenes)[mod];
}; rm(mod);

modGenesSym=list();
for (mod in 1:length(modGenesAnnos)) {
	x=toupper(as.character(modGenesAnnos[[mod]][,3]));
	x=strsplit(gsub('\\s','',x), '\\(');
	xclean=c();
	for (gene in 1:length(x)) {
		if (length(x[[gene]]) == 0) {
			next;
		} else {
			xclean=c(xclean,x[[gene]][1]);
		}
	}
	modGenesSym[[mod]]=xclean;
	names(modGenesSym)[mod]=names(modGenesAnnos)[mod];
}; rm(mod, x, gene, xclean);


library(org.Hs.eg.db);
IDs=as.list(org.Hs.egSYMBOL2EG);

# convert for each module and store in list
modIDs=list();
for (mod in 1:length(modGenes)) {
	temp=IDs[names(IDs) %in% modGenesSym[[mod]]];
	IDvec=c();
	for (gene in 1:length(temp)) {
		IDvec=c(IDvec,temp[[gene]]);
	}
	modIDs[[mod]]=IDvec;
	names(modIDs)[mod]=names(modGenesSym)[mod];
}
rm(mod,temp,IDvec,gene);
# set filenames for outputing each module to own file
filenames=c();
for (mod in 1:length(modIDs)){
	filenames=c(filenames,paste('isoNet.', names(modIDs)[mod], 'Entrez.txt', sep=''));
}
rm(mod);
# write files
for (mod in 1:length(modIDs)) {
	write.table(modIDs[[mod]], file=filenames[mod], quote=F, row.names=F, col.names=F)
}
rm(mod);

#or combine mods into 1 file
tmplengths = unlist(lapply(modIDs,length));
tmpMat = matrix(nrow=max(tmplengths), ncol=length(tmplengths));
colnames(tmpMat)=names(tmplengths); rm(tmplengths);
for (mod in 1:ncol(tmpMat)) {
	thismod = modIDs[[mod]];
	diff = nrow(tmpMat) - length(thismod);
	if (diff>0) {
		thismod = c(thismod, rep('',diff));
	}
	tmpMat[, mod] = thismod;
}; rm(mod, thismod, diff);
write.table(tmpMat, file='isoNet.allModsEntrez.txt', quote=F, row.names=F, sep='\t');
#####################
snp=read.table('../CichlidBAMs/Metriaclima_zebra.BROADMZ2cp.V2.gtfSNPs',header=F,sep='\t');
snps=c();
for (row in 1:nrow(snp)) {
	#snps = c(snps, gsub('gene_id ', '', strsplit(as.character(snp[row,9]), ';')[[1]][1]));
	snps = c(snps, gsub(' transcript_id ', '', strsplit(as.character(snp[row,9]), ';')[[1]][2]));
}; rm(row);

for (mod in 1:length(modGenes)) {
	print(sum(modGenes[[mod]] %in% snps));
}; rm(mod);

snp2=read.table('../CichlidBAMs/Metriaclima_zebra.BROADMZ2cp.V2.gtfSNPsASE',header=F,sep='\t');

source('~/Documents/_analysis_compare/_code/checkGeneListEnrichment.R');

checkGeneListEnrichmentList(unique(snps),modGenes,names(DATAk15))


snpGeneCounts=sort(table(snps),decreasing=T);
snpGeneCounts=snpGeneCounts[names(snpGeneCounts)%in%names(DATAk15)];
snpkME=c();
for (gene in 1:length(snpGeneCounts)) {
	thisgene=names(snpGeneCounts)[gene];
	mod = moduleColors[names(DATAk15)==thisgene];
	#kMErank=which(names(modGeneskMEsorted[names(modGeneskMEsorted)==mod][[1]])==thisgene);
	#cat(thisgene,', ',mod,', ',kMErank,'\n',sep='');
	kMEvals=modGeneskMEsorted[names(modGeneskMEsorted)==mod][[1]];
	snpkME=c(snpkME, kMEvals[names(kMEvals)==thisgene]);
}; rm(gene,mod,thisgene,kMEvals);


####################
lnc=read.table('/Users/ath/Documents/_Fernald_lab/_broadftp/Annotation/lncRNA/mzeb.lnc.final.gtf',header=F,sep='\t');
lncs=c();
for (row in 1:nrow(lnc)) {
	lncs = c(lncs, gsub(' transcript_id ', '', strsplit(as.character(lnc[row,9]), ';')[[1]][2]));
}; rm(row);
lncs=unique(lncs);

names(DATA)[gsub('mz.mrna.','',names(DATA)) %in% gsub('mz.lncrna.scaffold_','s',lncs) ];

for (mod in 1:length(modGenes)) {
	print(names(modGenes)[mod]);
	print(length(modGenes[[mod]]))
	print(sum(gsub('mz.mrna.','', modGenes[[mod]]) %in% gsub('mz.lncrna.scaffold_','s',lncs)));
}; rm(mod);

for (mod in 1:length(modGeneskMEsorted)) {
	inmod = gsub('mz.lncrna.scaffold_','s',lncs) %in% gsub('mz.mrna.','', names(modGeneskMEsorted[[mod]]));
	
}
