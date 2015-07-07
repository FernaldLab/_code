rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');
subjects = list.files()[grep('^3', list.files())];

### analyze CpG data that was combined across strands when symmetric
## load data, ignore some columns
d0 = list();
for (s in subjects) {
	cat(s,'\n');
	d0[[s]] = read.table(paste(s, '/aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean', sep=''), 
						 header=F, 
						 colClasses=c('character','numeric',rep('NULL',3),'numeric','numeric',rep('NULL',5))
						 );
}; rm(s); gc();
save(d0, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean_d0.RData');

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
save(bsd, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean_bsd.RData');
gc();

# smooth
bsdall = combine(bsd[[1]], bsd[[2]]);gc();
bsdall = combine(bsdall, bsd[[3]]);gc();
bsdall = combine(bsdall, bsd[[4]]);gc();
save(bsdall, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombined.RData')
bsdall.fit = BSmooth(bsdall,
					 ns=10, 				# min number of CpGs in window
					 h=250,					# half min window size
					 maxGap=1000,			# longest distance between 2 CpGs before cluster is broken into 2
					 mc.cores=4, 
					 parallelBy='sample');
gc();
save(bsdall, bsdall.fit, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombinedAndSmoothed.RData');
bsdall.fit2 = BSmooth(bsdall,
 					 ns=20, 				# min number of CpGs in window
 					 h=500,					# half min window size
 					 maxGap=1000,			# longest distance between 2 CpGs before cluster is broken into 2
 					 mc.cores=4, 
 					 parallelBy='sample');
save(bsdall.fit2, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanSmoothed2.RData');

# compute t-stats
bsdall.fitCov = getCoverage(bsdall.fit);
bsdall.fitCov2 = getCoverage(bsdall.fit2);
# keep only sites that had 5x or higher coverage in all samples
reqCov = 5;
reqSamples = 4;
keep = which(apply(bsdall.fitCov, 1, function(f) sum(f>=reqCov)==reqSamples));
keep2 = which(apply(bsdall.fitCov2, 1, function(f) sum(f>=reqCov)==reqSamples));
gc();

bsdall.fit = bsdall.fit[keep,];			 
bsdall.fit.tstat = BSmooth.tstat(bsdall.fit, 
								 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 group2=c('3157_TENNISON','3677_MONK'),
								 estimate.var='same',
								 local.correct=T,
								 mc.cores=4,
								 verbose=T
								 );
save(bsdall.fit, bsdall.fit.tstat, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanFilteredAndTstat.RData');

bsdall.fit2 = bsdall.fit2[keep2,];






allNAs = c();
for (j in 1:ncol(attributes(bsdall.fit.tstat)$stats)) {
	allNAs = c(allNAs, which(is.na( attributes(bsdall.fit.tstat)$stats[,j] )));
}; rm(j);
allNAs = unique(allNAs);

tvals = quantile(attributes(bsdall.fit.tstat)$stats[,5],na.rm=T,c(.005,.995));		 
dmrs0 = dmrFinder(bsdall.fit.tstat, cutoff=tvals);	
gc();

nthresh = 3;
diffThresh = 0.1;
dmrs = subset(dmrs0, n>=nthresh & abs(meanDiff)>=diffThresh);
dmrs = dmrs[order(abs(dmrs$areaStat),decreasing=T), ];
save(dmrs0, dmrs, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanDMRs.RData');

pData = pData(bsdall.fit);
pData$col = c('blue','red','red','blue');
pData(bsdall.fit) = pData;
								 
plotRegion(bsdall.fit, dmrs[4,], extend=2500, addRegions=dmrs, addPoints=F);

pdf(file='dmrs_top10.pdf', width=10, height=5);
plotManyRegions(bsdall.fit, dmrs[1:10,], extend=3000, addRegions=dmrs);
dev.off();


dmrs100bp = dmrs[dmrs$width>=100, ];

pdf(file='dmrs100bp_top100.pdf', width=10, height=5);
plotManyRegions(bsdall.fit, dmrs100bp[1:100,], extend=5000, addRegions=dmrs);
dev.off();


pdf(file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombinedSmoothed_ns20_h250_maxGap1000_2x.tstat_estvarSame_dmrs0_n20.pdf', width=10, height=5);
plotManyRegions(bsdall.fit, dmrs0[dmrs0$n>=20,], extend=5000, addRegions=dmrs);
dev.off();