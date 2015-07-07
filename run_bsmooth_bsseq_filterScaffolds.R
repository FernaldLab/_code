rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');

load('aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombined.RData');

scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
						
l70 = scaffoldCounts[scaffoldCounts$num < 70,];
l70n = na.omit(as.numeric(unlist(strsplit(as.character(l70$scaffold), '_'))));

bsdall2 = bsdall[-which(seqnames(bsdall) %in% l70$scaffold)];

scaffoldCounts2 = data.frame(scaffold=attributes(seqnames(bsdall2))$values, 
							num=attributes(seqnames(bsdall2))$lengths
							);
							
bsdall2.fit = BSmooth(bsdall2,
					 ns=70, 				# min number of CpGs in window
					 h=1000,					# half min window size
					 maxGap=10^8,			# longest distance between 2 CpGs before cluster is broken into 2
					 mc.cores=4, 
					 parallelBy='sample');
					 
reqCov = 5;
reqSamples = 4;
keep = which(apply(getCoverage(bsdall2.fit), 1, function(f) sum(f>=reqCov)==reqSamples));

bsdall2.fit5x = bsdall2.fit[keep,];		

tmp = data.frame(scaffold=attributes(seqnames(bsdall2.fit5x))$values, 
							num=attributes(seqnames(bsdall2.fit5x))$lengths
							);
bsdall2.fit5x2 = bsdall2.fit5x[-which(seqnames(bsdall2.fit5x) %in% as.character(tmp$scaffold[tmp$num<70]))]
	 
	 
	 
	 
bsdall2.fit5x2.tstat = BSmooth.tstat(bsdall2.fit5x2, 
								 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 group2=c('3157_TENNISON','3677_MONK'),
								 estimate.var='same',
								 local.correct=T,
								 mc.cores=4,
								 verbose=T
								 );
								 
dmrs0 = dmrFinder(bsdall2.fit5x2.tstat, cutoff=quantile(attributes(bsdall2.fit5x2.tstat)$stats[,6],na.rm=T,c(.001,.999)));
dmrs = dmrs0[dmrs0$n>=3 & abs(dmrs0$meanDiff)>=.1, ];						 
								
pData = pData(bsdall2.fit5x2);
pData$col = c('blue','red','red','blue');
pData(bsdall2.fit5x2) = pData;				 

pdf(file='run_bsmooth_bsseq_filterScaffolds_dmrs_n3_md.1.pdf', width=10, height=5);
plotManyRegions(bsdall2.fit5x2, dmrs, extend=5000, addRegions=dmrs);
dev.off();

#################################


rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');
load('aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined.RData');

scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
l70 = scaffoldCounts[scaffoldCounts$num < 70,];
bsdall2 = bsdall[-which(seqnames(bsdall) %in% l70$scaffold)];
bsdall2.fit = BSmooth(bsdall2,
					 ns=70, 				# min number of CpGs in window
					 h=1000,					# half min window size
					 maxGap=10^8,			# longest distance between 2 CpGs before cluster is broken into 2
					 mc.cores=4, 
					 parallelBy='sample');
save(bsdall2, bsdall2.fit, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults.RData');

######################		 
reqCov = 5;
reqSamples = 4;
keep = which(apply(getCoverage(bsdall2.fit), 1, function(f) sum(f>=reqCov)==reqSamples));
bsdall2.fit5x = bsdall2.fit[keep,];

tmp = data.frame(scaffold=attributes(seqnames(bsdall2.fit5x))$values, 
							num=attributes(seqnames(bsdall2.fit5x))$lengths
							);
bsdall2.fit5x2 = bsdall2.fit5x[-which(seqnames(bsdall2.fit5x) %in% as.character(tmp$scaffold[tmp$num<70]))];
bsdall2.fit5x2.tstat = BSmooth.tstat(bsdall2.fit5x2, 
								 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 group2=c('3157_TENNISON','3677_MONK'),
								 estimate.var='same',
								 local.correct=T,
								 mc.cores=4,
								 verbose=T
								 );gc()
save(bsdall2.fit5x2.tstat, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_5x.tstat.RData');
dmrs0_5x = dmrFinder(bsdall2.fit5x2.tstat, cutoff=quantile(attributes(bsdall2.fit5x2.tstat)$stats[,6],na.rm=T,c(.001,.999)));
save(dmrs0_5x, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_5x.tstat_dmrs0.RData');



######################
reqCov = 2;
reqSamples = 4;
keep = which(apply(getCoverage(bsdall2.fit), 1, function(f) sum(f>=reqCov)==reqSamples));
bsdall2.fit2x = bsdall2.fit[keep,];

tmp = data.frame(scaffold=attributes(seqnames(bsdall2.fit2x))$values, 
							num=attributes(seqnames(bsdall2.fit2x))$lengths
							);
bsdall2.fit2x2 = bsdall2.fit2x[-which(seqnames(bsdall2.fit2x) %in% as.character(tmp$scaffold[tmp$num<70]))];
bsdall2.fit2x2.tstat = BSmooth.tstat(bsdall2.fit2x2, 
								 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 group2=c('3157_TENNISON','3677_MONK'),
								 estimate.var='same',
								 local.correct=T,
								 mc.cores=4,
								 verbose=T
								 );gc()
save(bsdall2.fit2x2.tstat, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_2x.tstat.RData');
dmrs0_2x = dmrFinder(bsdall2.fit2x2.tstat, cutoff=quantile(attributes(bsdall2.fit2x2.tstat)$stats[,6],na.rm=T,c(.001,.999)));
save(dmrs0_2x, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_2x.tstat_dmrs0.RData');

######################


nThresh = 5;
meanDiffThresh = 0.15;
dmrs5x = dmrs0_5x[dmrs0_5x$n>=nThresh & abs(dmrs0_5x$meanDiff)>=meanDiffThresh, ];
write.table(dmrs5x, file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_5x.tstat_dmrs_n5_md.15', quote=F, row.names=F, col.names=T, sep='\t');

pData = pData(bsdall2.fit5x2);
pData$col = c('blue','red','red','blue');
pData(bsdall2.fit5x2) = pData;

pdf(file='aligned.adapters.q30.m0.methratio.CG.clean.upto1494Combined_filt_BSmooth_defaults_5x.tstat_dmrs_n5_md.15.pdf', width=10, height=5);
plotManyRegions(bsdall2.fit5x2, dmrs5x, extend=5000, addRegions=dmrs5x);
dev.off();