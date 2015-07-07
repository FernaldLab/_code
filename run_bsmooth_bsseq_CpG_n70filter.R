rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');

##########################################################################################################################
##########################################################################################################################

load(file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50.RData');
rm(bsdalln50);
bsdall;

##########################################################################################################################
##########################################################################################################################

scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
l70 = scaffoldCounts[scaffoldCounts$num < 70,];
l70n = na.omit(as.numeric(unlist(strsplit(as.character(l70$scaffold), '_'))));
bsdalln70 = bsdall[-which(seqnames(bsdall) %in% l70$scaffold)];
bsdalln70;
save(bsdalln70, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70.RData');

##########################################################################################################################
##########################################################################################################################

bsdalln70.fit = BSmooth(bsdalln70,
					 ns=70, 				# min number of CpGs in window
					 h=1000,					# half min window size
					 maxGap=10^8,			# longest distance between 2 CpGs before cluster is broken into 2
					 mc.cores=4, 
					 parallelBy='sample');
gc();
save(bsdalln70.fit, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.RData');


##########################################################################################################################
##########################################################################################################################

bsdalln70.fitCov = getCoverage(bsdalln70.fit);
# keep only sites that had 5x or higher coverage in all samples
reqCov = 5;
reqSamples = 4;
keep = which(apply(bsdalln70.fitCov, 1, function(f) sum(f>=reqCov)==reqSamples));
bsdalln70.fit5x = bsdalln70.fit[keep,];
scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdalln70.fit5x))$values, 
							num=attributes(seqnames(bsdalln70.fit5x))$lengths
							);
l70 = scaffoldCounts[scaffoldCounts$num < 70,];
l70n = na.omit(as.numeric(unlist(strsplit(as.character(l70$scaffold), '_'))));
bsdalln70.fit5x.n70 = bsdalln70.fit5x[-which(seqnames(bsdalln70.fit5x) %in% l70$scaffold)];
save(bsdalln70.fit5x, bsdalln70.fit5x.n70, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed_5x_n70.RData');


##########################################################################################################################
##########################################################################################################################


bsdalln70.fit5x.n70.tstat = BSmooth.tstat(bsdalln70.fit5x.n70, 
								 		  group1=c('3165_BRISCOE','3581_LYNLEY'),
								 		  group2=c('3157_TENNISON','3677_MONK'),
								 		  estimate.var='same',
										  local.correct=T,
										  mc.cores=4,
										  verbose=T
										  );gc();
save(bsdalln70.fit5x.n70.tstat, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed_5x_n70_tstat.RData');


##########################################################################################################################
##########################################################################################################################

plot(bsdalln70.fit5x.n70.tstat);
apply(attributes(bsdalln70.fit5x.n70.tstat)$stats, 2, summary);

for (pct in c(.01, .001, .0001)) {
	cat(pct, '\n');
	int = c(pct/2, 1-(pct/2));
	tvals = quantile(attributes(bsdalln70.fit5x.n70.tstat)$stats[,5], int, na.rm=T);
	assign(paste('dmrs0', pct, sep=''), dmrFinder(bsdalln70.fit5x.n70.tstat, cutoff=tvals));gc();
}; rm(pct,int,tvals);

dmrs01e04 = get('dmrs01e-04');
save(dmrs00.001, dmrs00.01, dmrs01e04, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed_5x_n70_tstat_dmrs0.RData');

nThresh = 3;
meanDiffThresh = 0.1;
dmrs.01 = subset(dmrs00.01, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.001 = subset(dmrs00.001, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.1e04 = subset(dmrs01e04, n>=nThresh & abs(meanDiff)>=meanDiffThresh);

par(mfrow=c(3,6));
for (i in c('dmrs.01', 'dmrs.001', 'dmrs.1e04')) {
	for (j in 7:12) {
		hist(get(i)[,j], xlab='', main=paste(i, names(get(i))[j]), breaks=nrow(get(i)))
	}
}; rm(i,j);


##########################################################################################################################
##########################################################################################################################

DMRS = dmrs.001;
cols = which(names(DMRS) %in% c('n','width','invdensity','areaStat','maxStat','meanDiff','tstat.sd'));
par(mfrow=c(2,4));
for (i in cols) {
	if (names(DMRS)[i] %in% c('areaStat','maxStat','meanDiff')) {
		WGCNA::verboseBoxplot(abs(DMRS[,i]), DMRS$direction, frame.plot=F, col='grey', ylab=paste('abs(',names(DMRS)[i],')',sep=''), xlab='');
	} else {
		WGCNA::verboseBoxplot(DMRS[,i], DMRS$direction, frame.plot=F, col='grey', ylab=names(DMRS)[i], xlab='');
	}
}; rm(cols,i);

BSD = bsdalln70.fit5x.n70;
pData = pData(BSD);
pData$col = c('blue','red','red','blue');
pData(BSD) = pData;

pdf(file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed_5x_n70_tstat_dmrs.001.pdf', width=10, height=5);
plotManyRegions(BSD, DMRS, extend=5000, addRegions=DMRS);
dev.off();

##########################################################################################################################
##########################################################################################################################

#exonDb = makeTranscriptDbFromGFF('../_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix_CDS.gtf', format='gtf');gc();
exonDb = read.table('../_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix_CDS.gtf', sep='\t', header=F);
exonDbGR = GRanges(seqnames=Rle(exonDb$V1), ranges=IRanges(start=exonDb$V4, end=exonDb$V5), strand=Rle(exonDb$V7), mcols=exonDb$V9);
utrDb = read.table('../_Burtoni_annotations/Astatotilapia_burtoni.BROADAB1.UTRs.gff3', sep='\t', header=F);
utrDbGR = GRanges(seqnames=Rle(utrDb$V1), ranges=IRanges(start=utrDb$V4, end=utrDb$V5), strand=Rle(utrDb$V7), mcols=utrDb$V9);

### make mcols into data frame with column for color to plot that interval, edit bsseq:::plotAnnoTrack

plotRegion(BSD, DMRS[1,], extend=5000, addRegions=DMRS, annoTrack=list(UBXN4=c(utrDbGR, exonDbGR)));


#plotRegion(BSD, DMRS[1,], extend=5000, addRegions=DMRS, annoTrack=list(UBXN4=exonsBy(txdb.ab,'gene')$ab.gene.s1314.1));



db = read.table('../_Burtoni_annotations/ref_AstBur1.0_scaffolds.clean.translate.final.combo.gff3', header=F, sep='\t');
dbGR = GRanges(seqnames=Rle(db$V1), ranges=IRanges(start=db$V4, end=db$V5), strand=Rle(db$V7), mcols=db$V9);

pdf(file='test.pdf', width=10, height=5);
plotManyRegions(BSD, DMRS[1:10,], extend=5000, addRegions=DMRS, annoTrack=list(annos=dbGR));
dev.off();