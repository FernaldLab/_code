rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');
subjects = list.files()[grep('^3', list.files())];

##########################################################################################################################
##########################################################################################################################

# ### analyze CpGs only on scaffolds0-1494
# ## load data, ignore some columns
# d0 = list();
# for (s in subjects) {
	# cat(s,'\n');
	# d0[[s]] = read.table(paste(s, '/aligned.adapters.q30.m0.methratio.CG.clean.upto1494', sep=''), 
						 # header=F, 
						 # colClasses=c('character','numeric',rep('NULL',3),'numeric','numeric',rep('NULL',5))
						 # );
# }; rm(s); gc();
# save(d0, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494_d0.RData');

### analyze CpGs only 
## load data, ignore some columns
d0 = list();
for (s in subjects) {
	cat(s,'\n');
	d0[[s]] = read.table(paste(s, '/aligned.adapters.q30.m0.methratio.CG.clean', sep=''), 
						 header=F, 
						 colClasses=c('character','numeric',rep('NULL',3),'numeric','numeric',rep('NULL',5))
						 );
}; rm(s); gc();
save(d0, file='aligned.adapters.q30.m0.methratio.CG.clean_d0.RData');

##########################################################################################################################
##########################################################################################################################

bsd = list();
# add colnames and create BSseq objects
for (s in 1:length(d0)) {
	cat(names(d0)[s], '\n');
	names(d0[[s]]) = c('chr','pos','Cov','M');
	if (names(d0)[s] %in% c('3157_TENNISON','3677_MONK')) {
		group = 'ND';
	} else if (names(d0)[s] %in% c('3165_BRISCOE','3581_LYNLEY')) {
		group = 'D'
	} else {
		stop('Check samples');
	}
	bsd[[names(d0)[s]]] = BSseq(M=as.matrix(d0[[s]]$M, ncol=1), 
							    Cov=as.matrix(d0[[s]]$Cov, ncol=1), 
							    pos=d0[[s]]$pos, 
							    chr=d0[[s]]$chr,
							    sampleNames=names(d0)[s],
							    pData=data.frame(group=as.character(group), row.names=names(d0)[s])
							    );
}; rm(s,group);
save(bsd, file='aligned.adapters.q30.m0.methratio.CG.clean_bsd.RData');
gc();

##########################################################################################################################
##########################################################################################################################


bsdall = combine(bsd[[1]], bsd[[2]]);gc();
bsdall = combine(bsdall, bsd[[3]]);gc();
bsdall = combine(bsdall, bsd[[4]]);gc();


scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
l50 = scaffoldCounts[scaffoldCounts$num < 50,];		
l50n = na.omit(as.numeric(unlist(strsplit(as.character(l50$scaffold), '_'))));
bsdalln50 = bsdall[-which(seqnames(bsdall) %in% l50$scaffold)];

save(bsdall, bsdalln50, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50.RData');

##########################################################################################################################
##########################################################################################################################

bsdalln50.fit = BSmooth(bsdalln50,
					 ns=50, 				# min number of CpGs in window
					 h=1000,					# half min window size
					 maxGap=10^8,			# longest distance between 2 CpGs before cluster is broken into 2
					 mc.cores=4, 
					 parallelBy='sample');
gc();
save(bsdalln50.fit, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50Smoothed.RData');

##########################################################################################################################
##########################################################################################################################

bsdalln50.fitCov = getCoverage(bsdalln50.fit);

# keep only sites that had 5x or higher coverage in all samples
reqCov = 5;
reqSamples = 4;
keep = which(apply(bsdalln50.fitCov, 1, function(f) sum(f>=reqCov)==reqSamples));
bsdalln50.fit5x = bsdalln50.fit[keep,];

scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdalln50.fit5x))$values, 
							num=attributes(seqnames(bsdalln50.fit5x))$lengths
							);
l50 = scaffoldCounts[scaffoldCounts$num < 50,];
l50n = na.omit(as.numeric(unlist(strsplit(as.character(l50$scaffold), '_'))));
bsdalln50.fit5x.n50 = bsdalln50.fit5x[-which(seqnames(bsdalln50.fit5x) %in% l50$scaffold)];

save(bsdalln50.fit5x, bsdalln50.fit5x.n50, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50Smoothed_5x_n50.RData');

##########################################################################################################################
##########################################################################################################################


bsdalln50.fit5x.n50.tstat = BSmooth.tstat(bsdalln50.fit5x.n50, 
								 		  group1=c('3165_BRISCOE','3581_LYNLEY'),
								 		  group2=c('3157_TENNISON','3677_MONK'),
								 		  estimate.var='same',
										  local.correct=T,
										  mc.cores=4,
										  verbose=T
										  );gc();
save(bsdalln50.fit5x.n50.tstat, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50Smoothed_5x_n50_tstat.RData');

##########################################################################################################################
##########################################################################################################################

plot(bsdalln50.fit5x.n50.tstat);
apply(attributes(bsdalln50.fit5x.n50.tstat)$stats, 2, summary);

for (pct in c(.01, .001, .0001)) {
	cat(pct, '\n');
	int = c(pct/2, 1-(pct/2));
	tvals = quantile(attributes(bsdalln50.fit5x.n50.tstat)$stats[,5], int, na.rm=T);
	assign(paste('dmrs0', pct, sep=''), dmrFinder(bsdalln50.fit5x.n50.tstat, cutoff=tvals));gc();
}; rm(pct,int,tvals);

dmrs01e04 = get('dmrs01e-04');
save(dmrs00.001, dmrs00.01, dmrs01e04, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50Smoothed_5x_n50_tstat_dmrs0.RData');

nThresh = 3;
meanDiffThresh = 0.1;
dmrs.01 = subset(dmrs00.01, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.001 = subset(dmrs00.001, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.1e04 = subset(dmrs01e04, n>=nThresh & abs(meanDiff)>=meanDiffThresh);

##########################################################################################################################
##########################################################################################################################

par(mfrow=c(1,3)); 
hist(abs(dmrs.01$meanDiff), breaks=nrow(dmrs.01)); 
hist(abs(dmrs.001$meanDiff), breaks=nrow(dmrs.001)); 
hist(abs(dmrs.1e04$meanDiff), breaks=nrow(dmrs.1e04));

pData = pData(bsdalln50.fit5x.n50);
pData$col = c('blue','red','red','blue');
pData(bsdalln50.fit5x.n50) = pData;				 

pdf(file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n50Smoothed_5x_n50_tstat_dmrs.001.pdf', width=10, height=5);
plotManyRegions(bsdalln50.fit5x.n50, dmrs.001, extend=5000, addRegions=dmrs.001);
dev.off();





