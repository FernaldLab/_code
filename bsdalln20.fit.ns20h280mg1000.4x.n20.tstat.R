rm(list=ls());


load('aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean_bsdCombined_n70.RData')
rm(bsdalln70);


NUMLOCI = 20;
### get number of methylation loci per scaffold
scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
							
### get scaffolds to keep
scaffoldsCountsKeep = scaffoldCounts[scaffoldCounts$num > NUMLOCI,];
NAME = paste('bsdalln', NUMLOCI, sep='');
assign(NAME, bsdall[which(seqnames(bsdall) %in% scaffoldsCountsKeep$scaffold)]);

save(bsdalln20, file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean_bsdCombined_n20.RData')

###########################################################################################

### with NUMLOCI=20, leaves 4549 scaffolds covering 98.8% of genome, shortest scaffold is 1002bp

### read in scaffold lengths
sl = read.table('/Volumes/fishstudies/_Burtoni_genome_files/scaffold_lengths');
### filter to scaffolds from above
slFilt = sl[match(scaffoldsCountsKeep$scaffold, sl$V1), ];
slFilt = slFilt[order(as.numeric(rownames(slFilt))), ];
### compute percentage of genome covered by remaining scaffolds
sum(slFilt[,2]) / sum(sl[,2]);
###########################################################################################

bsdalln20.fit.ns20h280mg1000 = BSmooth(bsdalln20, ns=20, h=280, maxGap=1000, mc.cores=4, parallelBy='sample');

gc();
save(bsdalln20.fit.ns20h280mg1000, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h280mg1000.RData'
	 );
	 
###########################################################################################
	 
	 
FIT = bsdalln20.fit.ns20h280mg1000;

FIT.Cov = getCoverage(FIT);
# keep only sites that had 4x or higher coverage in all samples
reqCov = 4;
reqSamples = 4;
keep = which(apply(FIT.Cov, 1, function(f) sum(f>=reqCov)==reqSamples));
FIT.FiltX = FIT[keep, ];

NUMLOCI = 20;
### get number of methylation loci per scaffold
scaffoldCountsPostFit = data.frame(scaffold=attributes(seqnames(FIT.FiltX))$values, 
								   num=attributes(seqnames(FIT.FiltX))$lengths
							       );
### get scaffolds to keep
scaffoldsCountsPostFitKeep = scaffoldCountsPostFit[scaffoldCountsPostFit$num > NUMLOCI,];
FIT.FiltXn = FIT.FiltX[which(seqnames(FIT.FiltX) %in% scaffoldsCountsPostFitKeep$scaffold)];

### compute percentage of genome covered by remaining scaffolds
sl = read.table('/Volumes/fishstudies/_Burtoni_genome_files/scaffold_lengths');
slFilt = sl[match(scaffoldsCountsPostFitKeep$scaffold, sl$V1), ];
slFilt = slFilt[order(as.numeric(rownames(slFilt))), ];
sum(slFilt[,2]) / sum(sl[,2]);


bsdalln20.fit.ns20h280mg1000.4x.n20 = FIT.FiltXn;
save(bsdalln20.fit.ns20h280mg1000.4x.n20, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h280mg1000_4x_n20.RData'
	 );
	 
	 
	 
FIT.FiltXn = bsdalln20.fit.ns20h280mg1000.4x.n20;

groupD = c('3165_BRISCOE','3581_LYNLEY');
groupND = c('3157_TENNISON','3677_MONK');
Sys.time()
FIT.FiltXn.T = BSmooth.tstat(FIT.FiltXn, 
							 group1=groupD,
						     group2=groupND,
			       	 		 estimate.var='same',
							 local.correct=T,
							 mc.cores=4,
							 verbose=T
							 );gc();
							 
							 
							 
bsdalln20.fit.ns20h280mg1000.4x.n20.tstat = FIT.FiltXn.T;
save(bsdalln20.fit.ns20h280mg1000.4x.n20.tstat, 
     file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h280mg1000_4x_n20_tstat.RData'
     );
     
     
FIT.FiltXn.T = bsdalln20.fit.ns20h280mg1000.4x.n20.tstat;

plot(FIT.FiltXn.T);
apply(attributes(FIT.FiltXn.T)$stats, 2, summary);

for (pct in c(.01, .001, .0001)) {
	cat(pct, '\n');
	int = c(pct/2, 1-(pct/2));
	tvals = quantile(attributes(FIT.FiltXn.T)$stats[,5], int, na.rm=T);
	assign(paste('dmrs0', pct, sep=''), dmrFinder(FIT.FiltXn.T, cutoff=tvals));gc();
}; rm(pct,int,tvals);

dmrs01e04 = get('dmrs01e-04');
save(dmrs00.001, dmrs00.01, dmrs01e04, 
     file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h280mg1000_4x_n20_tstat_dmrs0.RData'
     );
     
     
nThresh = 10;
meanDiffThresh = 0.1;
dmrs.01 = subset(dmrs00.01, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.001 = subset(dmrs00.001, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.1e04 = subset(dmrs01e04, n>=nThresh & abs(meanDiff)>=meanDiffThresh);