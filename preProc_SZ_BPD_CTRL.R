rm(list=ls());
setwd('~/Documents/_Fernald_lab');
source('_code/preProcATH-for_web_noVSN.R');

load('SZ_DATA_LOOKUPs.RData');

lapply(DATAlist, dim); lapply(LOOKUPlist, dim);
summary(DATAlist); summary(LOOKUPlist);

#################################
### outlier transcript values ###
#################################

outTranscripts = list();
for (set in 2:4)	#skip dat0
{
	outTranscripts[[set-1]] = removeOutlierProbesIterate(as.data.frame(DATAlist[[set]]), deviate=3, rowORcol=1);
	names(outTranscripts)[set-1] = names(DATAlist)[set];
}
rm(set);
save(outTranscripts, file = 'SZ_DATASETS_outTranscripts_deviate3.RData');

#######################################################
### remove transcripts with too many missing values ###
#######################################################

outNAs = list();
for (set in 1:3)
{
	DAT = outTranscripts[[set]]$dataClean;
	transcript_thresh = floor(ncol(DAT)/2);
	sample_thresh = floor(nrow(DAT)/2);
	outNAs[[set]] = removeTooManyNAs(DAT, probe_thresh=transcript_thresh, sample_thresh=sample_thresh);
	names(outNAs)[set] = names(outTranscripts)[set];
}
rm(set, DAT, transcript_thresh, sample_thresh);
save(outNAs, file = 'SZ_DATASETS_outTranscripts_deviate3_outNAs_moreThanHalf.RData');

##############################
### remove outlier samples ###
##############################

outSamples = list();

outSamples[[1]] = outlierSamplesIterate(outNAs[[1]]$dataClean);
outSamples[[2]] = outlierSamplesIterate(outNAs[[2]]$dataClean);
outSamples[[3]] = outlierSamplesIterate(outNAs[[3]]$dataClean);

names(outSamples) = names(outNAs);

save(outSamples, file = 'SZ_DATASETS_outTranscripts_deviate3_outNAs_moreThanHalf_outSamples_sd2.RData');

##############################
### quantile normalization ###
##############################

outQnorm = list();
for (set in 1:3)
{
	DAT = outSamples[[set]]$dataClean;
	outQnorm[[set]] = as.data.frame(normalize.quantiles(as.matrix(DAT)));
	dimnames(outQnorm[[set]]) = dimnames(DAT);
}
rm(set, DAT);
names(outQnorm) = names(outSamples);

lapply(outQnorm, dim);

save(outQnorm, file = 'SZ_DATASETS_outTranscripts_deviate3_outNAs_moreThanHalf_outSamples_sd2_Qnorm.RData');

###########################
### build test networks ###
###########################
rm(list=ls());
setwd('~/Documents/_Fernald_lab');
load('SZ_DATASETS_outTranscripts_deviate3_outNAs_moreThanHalf_outSamples_sd2_Qnorm.RData');
library(WGCNA);
options(stringsAsFactors = F);
allowWGCNAThreads();
collectGarbage();

pickSoft = list();
for (set in 1:3)
{
	DAT = as.data.frame(t(outQnorm[[set]]));
	pickSoft[[set]] = pickSoftThreshold(DAT, networkType='signed', blockSize=6000, verbose=2);
	collectGarbage();
}
rm(set, DAT);
names(pickSoft) = names(outQnorm);

save(pickSoft, file = 'SZ_DATASETS_outTranscripts_deviate3_outNAs_moreThanHalf_outSamples_sd2_Qnorm_pickSoft.RData');

# plot scale-free fit, mean k, and slope of fit as function of power

#jpeg(file = '.jpg', width = 8, height = 10, units = 'in', quality = 100, type = 'quartz', res = 150);

par( mfrow = c(3, 3) );
for( j in 1:3 ) {
	set = pickSoft[[j]]$fitIndices;
	x = set$Power;
	y = set$SFT.R.sq;
	plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = 0.8, col = 'red', lty = 'dashed' );
	
	y = set$mean.k.;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = 50, col = 'red', lty = 'dashed' );
	
	y = set$slope;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = -2, col = 'red', lty = 'dashed' );
	abline( h = -1, col = 'red', lty = 'dashed' );
	}
	rm(j, set, x, y);
	
#dev.off();

# or to do one at a time
set0 = 'dat.SZ';

set0 = match(set0, names(pickSoft));
par( mfrow = c(1, 3) );
set = pickSoft[[set0]]$fitIndices;
x = set$Power;
y = set$SFT.R.sq;
plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = 0.8, col = 'red', lty = 'dashed' );

y = set$mean.k.;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = 50, col = 'red', lty = 'dashed' );
	
y = set$slope;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

rm(set0, set, x, y);

####
k14.SZ = softConnectivity(as.data.frame(t(outQnorm$dat.SZ)), type='signed', power=14, blockSize=6000, verbose=3);
names(k14.SZ) = names(as.data.frame(t(outQnorm$dat.SZ)));
par(mfrow=c(1,2));hist(k14.SZ);scaleFreePlot(k14.SZ);
collectGarbage();

net.SZ0 = blockwiseModules(as.data.frame(t(outQnorm$dat.SZ)),
 						   maxBlockSize = 6000,
 						   power = 14,
 						   networkType = 'signed',
 						   deepSplit = 2,
 						   minModuleSize = 10,
 						   verbose = 3,
 						   saveTOMs = F
 						   );
collectGarbage(); 						   

block = 1;
plotDendroAndColors(net.SZ0$dendrograms[[block]],
					net.SZ0$colors[ net.SZ0$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = net.SZ0$colors[ net.SZ0$blockGenes[[block]] ],
					main = paste('block', block),
					dendroLabels = F,
					#hang = 0.03,
					addGuide = T,
					guideHang = 0.05);
rm(block); 

bgGenes = net.SZ0$colors=='grey';
DAT = as.data.frame(t(outQnorm$dat.SZ));

######
# net=blockwiseModulesEnriched(DAT, maxBlockSize=6000, skipThresh=500, skipGrey=T, onlyGreyThresh=2, verbose=2, saveFileBase='SZ_dev3_sd2_qnorm_permTest1000Skip500')
######

rm(list=ls());
setwd('~/Documents/_Fernald_lab');
library(WGCNA);
options(stringsAsFactors = F);
allowWGCNAThreads();

load('SZ_dev3_sd2_qnormrun25DATA.RData');
source('_code/blockwiseModulesEnriched-Feb2013.R');
net=blockwiseModulesEnriched(DAT,maxBlockSize=6000,densityPermTest=F,verbose=2,saveFileBase='SZ_dev3_sd2_qnorm',setRun=25);