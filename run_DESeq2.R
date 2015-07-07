rm(list=ls());
options(stringsAsFactors=F);
library(DESeq2);library(cqn);library(WGCNA);allowWGCNAThreads();
setwd('/Volumes/fishstudies-1/_Elim_RNAseq/4_expression/NCBI/htseq-count');

# get filenames to read in
dirs = c('Jan2015/','Mar2015/');
files = c();
for (d in dirs) {
	files = c(files, paste(d, list.files(d), sep=''));
}; rm(d);

# get read groups info file
rg = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);
# make sampleTable for DESeqDataSetFromHTSeqCount
rg = cbind(SeqInd=apply(rg[,2:3],1,function(f) paste(f, collapse='/')), rg);
rg = cbind(SeqInd=rg[,1], File=files, rg[,2:ncol(rg)], Sex=rep('M',nrow(rg)));
rg$Sex[rg$Condition=='F'] = 'F';
#rg = rg[,-4];

# RG = subset(rg, Condition %in% c('D','F'));
# DESIGN = formula(~ Sex);
# dds = DESeqDataSetFromHTSeqCount(sampleTable=RG, directory='./', design=DESIGN);
# ddsDE = DESeq(dds);
# res = results(ddsDE);

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
names(counts) = gsub(' ','',apply(rg[,c(3,4,5,7,8)], 1, function(f) paste(f, collapse='-')));

# restrict to only rows with "LOCxxxxxxxxx" as gene name
counts = counts[grepl('^LOC', rownames(counts)), ];
write.table(counts, file='countsLOC.txt', sep='\t', quote=F, row.names=T, col.names=T);
##############################################################################
###########################################################

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


cqnOut = cqn(counts=counts, x=gff$V10, lengths=gff$V11, verbose=T);
cqnRPKM = cqnOut$y + cqnOut$offset;
RPKM = 2^cqnRPKM

#Back calculate counts from RPKM (Normalized counts)
countsFromRPKM = t(t(RPKM)*(colSums(counts)/1000000))*(gff$V11/1000);

#############

RPKMlib2 = RPKM[, grep('^2', colnames(RPKM))];


# gMeans = rowSums(counts)/ncol(counts);
# counts = counts[-which(gMeans<median(gMeans)), ];

# countsNorm = counts2;
# for (s in 1:ncol(countsNorm)) {
	# countsNorm[,s] = countsNorm[,s]/max(countsNorm[,s])
# }; rm(s)


########################################
#################################
########################












SAMPLES = which(rg$Condition %in% c('D','F'));
COUNTS = counts[, SAMPLES];
RG = rg[SAMPLES, ];
DESIGN = formula(~ Sex);
countsDE = DESeqDataSetFromMatrix(COUNTS, RG, DESIGN);
countsDE = DESeq(countsDE, test='LRT');
res = results(countsDE);



###################################

source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');

DAT = as.data.frame(t(log2(RPKMlib2)));

# try on rpkm values from cqn, no filtering
exn.plotPowers(DAT, blockSize=5000);

k14 = exn.computeConnectivityAndPlotScaleFreeness(DAT, power=14, blockSize=5000);
k16 = exn.computeConnectivityAndPlotScaleFreeness(DAT, power=16, blockSize=5000);
k18 = exn.computeConnectivityAndPlotScaleFreeness(DAT, power=18, blockSize=5000);

net = blockwiseModulesEnriched(DAT, power=18, minModuleSize=20, maxBlockSize=5000, densityPermTest=F);

# try filtering 
thresh.rpkm = 0.5;
thresh.num = floor(nrow(DAT)/3);

numlow = apply(DAT, 2, function(f) sum(f < thresh.rpkm));
remove = which(numlow > thresh.num);
DAT = DAT[, -remove];

exn.plotPowers(DAT, blockSize=5000);
k18 = exn.computeConnectivityAndPlotScaleFreeness(DAT, power=18, blockSize=5000);
DAT = DAT[,k18 >= median(k18)];
net = blockwiseModulesEnriched(DAT, power=18, minModuleSize=20, maxBlockSize=5000, densityPermTest=F);

exn.getNetworkBasics(net, DAT, '');

exn.plotEigengeneNetworks2(MEs);


MEs.2tmp = mergeCloseModules(DAT, colors, MEs=MEs, cutHeight=.75)
colors2 = MEs.2tmp$colors
MEs2 = MEs.2tmp$newMEs

rg = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);
# traits = rg[,c(1,2,6)];
# traits = traits[grep('^2', colnames(RPKM)), ];		# if using only lib2 batch
# traits$BatchLib = traits$BatchLib - 1;
# traits$BatchSeq = c(0,0,0,0,0,1,1,1,1,1,1,1,1,1);
# traits$Condition[traits$Condition=='F'] = 0;
# traits$Condition[traits$Condition!=0] = 1;

traits = rg[, c(4,6)];
traits = traits[grep('^2', colnames(RPKM)), ];		# if using only lib2 batch
rownames(traits) = rownames(DAT);
traits = as.data.frame(model.matrix(~0+Condition+Lane, data=traits));
names(traits) = gsub('Condition','',names(traits));

traitCors = exn.computeAndPlotMETraitCors(traits, MEs, main='');

factors = rg$Condition[grep('^2', colnames(RPKM))];
exn.plotMEsExprWithTraits(MEs, factors, c(5,8), cex.axis=.8, cex.main=1.1);
GS = exn.computeGS(traits, DAT);

block = 2;
colorFrame = data.frame(colors[blockGenes[[block]]], 
					    numbers2colors(GS$GS.ASC[blockGenes[[block]]]), 
					    numbers2colors(GS$GS.D[blockGenes[[block]]]), 
						numbers2colors(GS$GS.F[blockGenes[[block]]]), 
						numbers2colors(GS$GS.Lane[blockGenes[[block]]])
						);
plotDendroAndColors(dendro[[block]], colorFrame,
					dendroLabels=F, addGuide=F,
					groupLabels=gsub('numbers2colors.GS.','',names(colorFrame),fixed=T), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);

#####################################
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

datIN = as.data.frame(RPKM);

## no zeros post cqn
# zeros = apply(datIN, 1, function(f) sum(f==0));
# zthresh = floor(ncol(datIN)/3);
# remove = which(zeros > zthresh);
# datIN = datIN[-remove, ];


noyze = quantile(vectorizeMatrix(as.matrix(datIN),T), .25)
datIN[datIN < noyze] = NA;

sumNAs = apply(datIN, 1, function(f) sum(is.na(f)));summary(sumNAs);
NAthresh = 0;
datIN = datIN[-which(sumNAs > NAthresh), ];

datIN = log2(datIN);
out = preProcess(datIN=datIN, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
				 
###########################
# source('~/Downloads/RemoveOutliers/RemoveOutliers_0.9');
# library(cluster)
# library(impute)
# library(preprocessCore)

# DAT = countsFromRPKM;
# DAT = cbind(1:nrow(DAT), DAT);
# sampleInfo = as.data.frame(matrix(nrow=(ncol(DAT)-1), ncol=7));
# names(sampleInfo) = c('sampleID','batchSeq','batchLib','lane','condition','sex','dummyGroup');
# sampleInfo$sampleID = colnames(DAT)[2:ncol(DAT)];
# sampleInfo$batchSeq = rg$BatchSeq;
# sampleInfo$batchLib = rg$BatchLib;
# sampleInfo$lane = rg$Lane;
# sampleInfo$condition = rg$Condition;
# sampleInfo$sex = rg$Sex;
# sampleInfo$dummyGroup = rep('A',nrow(sampleInfo));
# indices = list(c(2:ncol(DAT)));

# INFO = sampleInfo;
# RemoveOutliers(datexpr1=DAT,
			   # subset1=NULL,
			   # impute1=FALSE,
			   # skip1=1,
			   # indices1=indices,
			   # sampleinfo1=INFO,
			   # samplelabels1=1,
			   # grouplabels1=7,
			   # trait1=c(2,3,5,6),
			   # btrait1=c(2,3,5,6),
			   # asfactors1=c(2,3,5,6),
			   # projectname1='combatTest',
			   # cexlabels=.9,
			   # normalize1=TRUE,
			   # replacenegs1=TRUE,
			   # fitmodels1=TRUE,
			   # whichfit1='mean',
			   # exportfigures1=FALSE
			   # )
			   
			   
############################
#######################
library(sva);

# DAT = datIN;
# BATCH = as.numeric(substr(names(datIN),1,1));
# MOD = model.matrix(~as.factor(BatchLib)+as.factor(BatchSeq)+as.factor(Condition), data=rg);

datToWrite = datIN;
datToWrite = as.data.frame(cbind(rownames(datToWrite), datToWrite));
write.table(datToWrite, file='datForComBat.txt', quote=F, sep='\t', row.names=F);

tmp = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);
tmp = cbind(gsub(' ','',apply(tmp[,c(1,2,3,5,6)], 1, function(f) paste(f, collapse='-'))), tmp);
tmp = tmp[, c(1,6,2,3,7)];
tmp$Condition[tmp$Condition != 'F'] = 'M';
names(tmp) = c('Array name', 'Sample name', 'Batch', 'Covariate 1', 'Covariate 2');
write.table(tmp, file='rgForComBat.txt', quote=F, sep='\t', row.names=F);

combatout = ComBat(expression_xls='datForComBat.txt', sample_info_file='rgForComBat.txt', filter=F, write=F, skip=1);

write.table(combatout, file='datForComBat_run2.txt', quote=F, sep='\t', row.names=F);

tmp = tmp[, c(1,2,4,3,5)];
names(tmp)[c(3,4)] = c('Batch', 'Covariate 1');
write.table(tmp, file='rgForComBat_run2.txt', quote=F, sep='\t', row.names=F);

combatout2 = ComBat(expression_xls='datForComBat_run2.txt', sample_info_file='rgForComBat_run2.txt', filter=F, write=F, skip=1);

#newDAT = as.data.frame(cbind(datIN[,which(!(names(datIN) %in% names(combatout)))], combatout));
newDAT = combatout2;
rownames(newDAT) = combatout2[,1];
newDAT = newDAT[, -1];

out = preProcess(datIN=newDAT, 
				 removeOutlierProbes=F, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=F);




#################
############
#######
###
#

# # cqn_data<- cqn(Peromyscus_Expression_Raw_Counts[,1:8], lengths = Peromyscus_Expression_Raw_Counts$LENGTH,x = Peromyscus_Expression_Raw_Counts$PERC_GC, sizeFactors = colSums(Peromyscus_Expression_Raw_Counts[,1:8]), verbose = TRUE)
# png(file="DASE_Data_CQN_Norm_Output_140819.png",width=2000,height=1000)
# par(mfrow=c(1,2))
# cqnplot(cqn_data, n = 1, xlab = "GC content", lwd=5, lty = 1, ylim = c(1,7))
# cqnplot(cqn_data, n = 2, xlab = "length", lwd=5, lty = 1, ylim = c(1,7))
# dev.off()
# RPKM.cqn <- cqn_data$y + cqn_data$offset

# #RPKM
# Peromyscus_HTSEQ_RPKM<-2^RPKM.cqn

# #Back calculate counts from RPKM (Normalized counts)
# Peromyscus_HTSEQ_Back_Counts<-t(t(Peromyscus_HTSEQ_RPKM)*(colSums(Peromyscus_Expression_Raw_Counts[,1:8])/1000000))*(Peromyscus_Expression_Raw_Counts$LENGTH/1000)

# #Normalized counts
# write.table(Peromyscus_HTSEQ_Back_Counts,"~/Desktop/HTSEQ/Peromyscus_HTSEQ_Back_Counts.txt", sep="\t", quote=FALSE)

# #Normalized RPKM
# write.table(Peromyscus_HTSEQ_RPKM,"~/Desktop/HTSEQ/Peromyscus_HTSEQ_RPKM.txt", sep="\t", quote=FALSE)

