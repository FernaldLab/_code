##############################################################################
###########################################################
##########################################
### load data 
rm(list=ls());
options(stringsAsFactors=F);
setwd('/Volumes/fishstudies-1/_Elim_RNAseq/4_expression/NCBI/htseq-count');

# get filenames to read in
dirs = c('Jan2015/','Mar2015/');
files = c();
for (d in dirs) {
	files = c(files, paste(d, list.files(d), sep=''));
}; rm(d);

# get data
rawData = list();
for (f in files) {
	rawData[[f]] = read.table(f);
}; rm(f);

counts = rawData[[1]];
for (r in rawData) {
	counts = cbind(counts, r[,2]);
}; rm(r);
rownames(counts) = counts[,1];
counts = counts[, -c(1,2)];

# get read groups info file
# make sure rows are in same order as "files"
rg0 = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);
rg = rg0; 
rg$BatchSeq = gsub('2015', '', rg$BatchSeq);

names(counts) = gsub(' ','',apply(rg[,c(1,2,4,5,6)], 1, function(f) paste(f, collapse='-')));

# restrict to only rows with "LOCxxxxxxxxx" as gene name
counts = counts[grepl('^LOC', rownames(counts)), ];
write.table(counts, file='countsLOC.txt', sep='\t', quote=F, row.names=T, col.names=T);
##############################################################################
###########################################################
##########################################
### normalize with cqn
# match to loci from ncbi gff file
counts = read.table('countsLOC.txt',header=T,sep='\t');
names(counts) = gsub('X','',names(counts));

gff = read.table('~/Documents/_Burtoni_annotations/ref_AstBur1.0_scaffolds.clean.translate.final.combo.gff3_bedtools_nucLOC',sep='\t',header=F);
gffSplit = unlist(strsplit(gff$V9, ';'));
gffLOCs = grep('gene=', gffSplit);
gffLOCs = unlist(strsplit(gffSplit[gffLOCs], 'gene='));
gffLOCs = gffLOCs[seq(2,length(gffLOCs),2)];
rownames(gff) = gffLOCs;

gff = gff[rownames(gff) %in% rownames(counts), ];
counts = counts[match(rownames(gff), rownames(counts)), ];
counts0 = counts;

library(cqn);
cqnOut = cqn(counts=counts, x=gff$V10, lengths=gff$V11, verbose=T);
cqnRPKM = cqnOut$y + cqnOut$offset;
RPKM = 2^cqnRPKM;

#Back calculate counts from RPKM if needed (Normalized counts)
countsFromRPKM = t(t(RPKM)*(colSums(counts)/1000000))*(gff$V11/1000);

##############################################################################
###########################################################
##########################################
### remove noise and outliers
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

datIN = as.data.frame(RPKM);

# remove bottom x% of values
centile = 0.25;
noyze = quantile(WGCNA::vectorizeMatrix(as.matrix(datIN), T), centile);
datIN[datIN < noyze] = NA;

# remove genes with >NAthresh noise values
sumNAs = apply(datIN, 1, function(f) sum(is.na(f)));
summary(sumNAs); table(sumNAs);
NAthresh = 0;
datIN = datIN[-which(sumNAs > NAthresh), ];

# log2 transform before more processing
datIN = log2(datIN);

# define cv helper function
.cv = function(x){return(mean(x,na.rm=T)/sd(x,na.rm=T))}

# remove genes in bottom 10% variance
cv = apply(datIN, 1, .cv);
datIN = datIN[cv > quantile(cv, .1), ];

# remove outlier measurements and samples, check for batch effects
out = preProcess(datIN=datIN, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
## if some measurements ("probes") removed, need to impute missing values for ComBat
library(impute);
datImpute = impute.knn(as.matrix(out$data_Qnorm))$data;
	 
# didn't remove samples but major batch effects of library prep and sequencing run

##############################################################################
###########################################################
##########################################
### remove batch effects
source('/Volumes/fishstudies/_code/ComBat.R');

# format expression data and write to file
datToWrite = as.data.frame(datImpute);
datToWrite = as.data.frame(cbind(gene=rownames(datToWrite), datToWrite));
write.table(datToWrite, file='datForComBat.txt', quote=F, sep='\t', row.names=F);

# format sample info table and write to file
# assumes samples in rows of rg match samples in cols of datToWrite
rg0 = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);

db0 = read.csv('/Volumes/fishstudies/dbsATH.csv');
db0$Dissection.id[db0$Dissection.id == 'RS102414'] = 'RS102414.1';
db = db0[match(rg0$Sample, gsub('RS','',db0$Dissection.id)), ];

rg = rg0; 
rg$BatchSeq = gsub('2015', '', rg$BatchSeq);
rg = cbind(names(datToWrite)[2:ncol(datToWrite)], rg);
rg = as.data.frame(cbind(rg[, c(1,6,2,3,7)], db$Tank));
names(rg) = c('Array name', 'Sample name', 'Batch', 'Covariate 1', 'Covariate 2', 'Covariate 3');
write.table(rg, file='rgForComBat.txt', quote=F, sep='\t', row.names=F);

# run ComBat to correct for library batch effect
combatout = ComBat(expression_xls='datForComBat.txt', 
				   sample_info_file='rgForComBat.txt', 
				   filter=F, write=F, skip=1
				   );
				   
# save output and setup to run ComBat again for sequencing run batch effect
write.table(combatout, file='datForComBat_run2.txt', quote=F, sep='\t', row.names=F);

rg = rg[, c(1,2,4,3,5,6)];
names(rg)[c(3,4)] = c('Batch', 'Covariate 1');
rg$Batch[rg$Batch=='Jan'] = 1;
rg$Batch[rg$Batch=='Mar'] = 2;
write.table(rg, file='rgForComBat_run2.txt', quote=F, sep='\t', row.names=F);

combatout2 = ComBat(expression_xls='datForComBat_run2.txt', 
				    sample_info_file='rgForComBat_run2.txt', 
				    filter=F, write=F, skip=1
				    );
				    
				    
# save output and setup to run ComBat again for dissection order batch effect
write.table(combatout2, file='datForComBat_run3.txt', quote=F, sep='\t', row.names=F);

rg = rg[, c(1,2,6,4,3,5)];
names(rg)[3:6] = c('Batch', 'Covariate 1', 'Covariate 2', 'Covariate 3');
write.table(rg, file='rgForComBat_run3.txt', quote=F, sep='\t', row.names=F);

combatout3 = ComBat(expression_xls='datForComBat_run3.txt', 
				    sample_info_file='rgForComBat_run3.txt', 
				    filter=F, write=F, skip=1
				    );
				    
				    
				    
				    
				    
newDAT = combatout3;
rownames(newDAT) = combatout3[,1];
newDAT = newDAT[, -1];

outCheck = preProcess(datIN=newDAT, 
				 removeOutlierProbes=F, 
				 removeTooManyNAs=F, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=F);
				 
##############################################################################
###########################################################
##########################################
### test scale-freeness
library(WGCNA);allowWGCNAThreads();
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');

DAT = as.data.frame(t(newDAT));
exn.plotPowers(DAT, blockSize=5000);

k20 = exn.computeConnectivityAndPlotScaleFreeness(DAT, power=20, blockSize=5000);
toRemove = which(k20 < quantile(k20, .25));

DAT = DAT[, -toRemove];
collectGarbage();
net = blockwiseModulesEnriched(DAT, power=20, minModuleSize=100, maxBlockSize=ncol(DAT)+1, densityPermTest=F);

# took 41 runs to remove grey genes, load that data
load('run41DATA.RData');
#load('run3DATA.RData');

# extract relevant info from net
exn.getNetworkBasics(net, DATA, '');

exn.plotEigengeneNetworks2(MEs);


### build trait table
db0 = read.csv('/Volumes/fishstudies/dbsATH.csv');
db0$Dissection.id[db0$Dissection.id == 'RS102414'] = 'RS102414.1';
db = db0[match(rg0$Sample, gsub('RS','',db0$Dissection.id)), ];

traits0 = as.data.frame(cbind(rg0[,-3], db[,c(1,4,5,10:15,23:25,31,38)]));
rownames(traits0) = traits0$Sample;
traits0 = traits0[, -c(3,4,7,8,10)];

traits0$BatchSeq[traits0$BatchSeq == 'Jan2015'] = 1;
traits0$BatchSeq[traits0$BatchSeq == 'Mar2015'] = 2;

traits = traits0;
traits = as.data.frame(cbind(traits, as.data.frame(model.matrix(~0+Condition, traits))));
traits = as.data.frame(cbind(traits, as.data.frame(model.matrix(~0+Tank, traits))));
traits = traits[, -c(3,4,5)];		# ignoring time of sacrifice for now, need to figure out daylight savings issue
names(traits) = gsub('Condition|Batch','',names(traits));

###
traitCors = exn.computeAndPlotMETraitCors(traits, MEs, main='');
GS = exn.computeGS(traits, DATA);

colorFrame = data.frame(colors, 
						numbers2colors(GS$GS.Fish.length.cm),  
						numbers2colors(GS$GS.GSI), 
						numbers2colors(GS$GS.ASC),
						numbers2colors(GS$GS.D),
						numbers2colors(GS$GS.ND),
						numbers2colors(GS$GS.F),
						numbers2colors(GS$GS.Lib));
plotDendroAndColors(dendro[[1]], colorFrame,
					dendroLabels=F, addGuide=F,
					groupLabels=gsub('numbers2colors.GS.GS.','',names(colorFrame),fixed=T), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);
					
					
exn.plotAllModsGSkME(colors=colors, trait='F', GStable=GS, kMEtable=kME, mfrow=c(5,7), order=T, returnCors=F);
exn.plotModuleAllGSkME('orange', colors, GS, kME, c(3,6), horiz=T);

par(mfrow=c(2,2)); 
for (i in c(6:9)) {
	verboseBoxplot(traits0[,i], traits0$Condition, col='grey', ylab=names(traits0)[i], xlab='', frame.plot=F);
}; rm(i);



traits2 = cbind(traits, GSI.F=traits$GSI, GSI.M=traits$GSI);
traits2$GSI.F[traits2$F==0] = NA;
traits2$GSI.M[traits2$F==1] = NA;
traits2 = traits2[, c(1,2,7:11,16:18,3,4,13,15,12,14,6,19,20)];

traits2 = cbind(traits2, Fish.length.cm.D=traits2$Fish.length.cm, Fish.length.cm.REST=traits2$Fish.length.cm);
traits2$Fish.length.cm.D[traits2$D==0] = NA;
traits2$Fish.length.cm.REST[traits2$D==1] = NA;
###############
###########

library(biomaRt);
library(RDAVIDWebService);

anno0 = read.table('~/Documents/_Burtoni_annotations/H_burtoni_rna_blastx_FISH_ENS_top1_reciprocalBackFrom_Drer_Olat_Onil_Trub_ENS_pep_noComments_passRecipHsENS',header=T,sep='\t');

anno = anno0[match(names(DATA),anno0$gene), ];

modGenesHs = modGenes;
for (m in 1:length(modGenesHs)) {
	modGenesHs[[m]] = anno$hsENS[anno$gene %in% modGenesHs[[m]]];
	modGenesHs[[m]] = as.vector(na.omit(unlist(strsplit(modGenesHs[[m]], ','))));
}; rm(m);

dHs = .startAndUploadBgAndModGeneListToDAVID(modGenes=modGenesHs, login='ahilliar@stanford.edu', idType='ENSEMBL_GENE_ID');
modDAVID0 = exn.getModChartsFromDAVID(dHs);

modDAVID = exn.filterModDAVIDList(modDAVID0, 13, 10, T);
modDAVIDterms = exn.modDAVIDByTerm(modDAVID);

modDAVIDconvertedGenes = .convertTermGenesForDAVIDChartList(modDAVID, anno, 7, 1)

modDAVIDkME = .addkMEColsToDAVIDList(modDAVIDconvertedGenes,modkMEs);

modDAVIDclust = .getTermClustersFromDAVID(dHs, EASEthresh=1.2)



modGenesHsEntrez0 = .convertIDsWithBiomaRtList(modGenes=modGenesHs);
modGenesHsEntrez = modGenesHsEntrez0;
for (m in 1:length(modGenesHsEntrez)) {
	modGenesHsEntrez[[m]] = modGenesHsEntrez[[m]]$entrezgene;
}; rm(m);

library(GOFunction);
modGOFunction = .GOFunctionModules(modGenesHsEntrez, unlist(modGenesHsEntrez));