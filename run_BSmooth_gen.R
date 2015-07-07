###################################################################################
##### load data and setup workspace
###################################################################################
rm(list=ls());
options(stringsAsFactors=F);
library('bsseq');
print(getwd())
load('data/aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombined.RData');
bsdall;

###################################################################################
##### set parameters
###################################################################################
##### window settings for BSmooth
ns = ;
h = ;
maxGap = ;

##### filtering before BSmooth.tstat
reqCov = ;

##### BSmooth.tstat variance estimation
estimate.var = ;

##### filtering tstats before dmrFinder
tFilter = ;

#####

###################################################################################
##### run BSmooth
###################################################################################
bsdall.fit = BSmooth(bsdall, ns=ns, h=h, maxGap=maxGap, mc.cores=4, parallelBy='sample');				 
bsdall.fit;
filename = paste('aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombinedSmoothed_ns', ns,
				 '_h', h,
				 '_maxGap', maxGap,
				 '.RData',
				 sep=''
				 );
cat('writing smoothed data to file: ', filename, '\n', sep='');
save(bsdall.fit, file=filename);

###################################################################################
##### run BSmooth.tstat
###################################################################################
##### filter smoothed data to only sites with required coverage
bsdall.fitCov = getCoverage(bsdall.fit);
reqSamples = 4;
keep = which(apply(bsdall.fitCov, 1, function(f) sum(f>=reqCov)==reqSamples));
bsdall.fit = bsdall.fit[keep, ];
bsdall.fit;

bsdall.fit.tstat = BSmooth.tstat(bsdall.fit, 
								 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 group2=c('3157_TENNISON','3677_MONK'),
				 				 estimate.var=estimate.var,
								 local.correct=T,
							     mc.cores=4,
								 verbose=T
								 );
bsdall.fit.tstat;
summary(attributes(bsdall.fit.tstat)$stats[,6]);

filename = gsub('.RData', '', filename, fixed=T);
filename = paste(filename, '_', reqCov, 'x.RData', sep='');
cat('writing tstats to file: ', filename, '\n', sep='');
save(bsdall.fit.tstat, file=filename);

filename = gsub('.RData', '', filename, fixed=T);
filename = paste(filename, '_tDist.pdf', sep='');
cat('saving plot of tstat distribution to file: ', filename, '\n', sep='');
pdf(file=filename, width=10, height=5);
plot(bsdall.fit.tstat[!is.na(attributes(bsdall.fit.tstat)$stats[,6]), ]);
dev.off();

###################################################################################
##### run dmrFinder
###################################################################################
tFilterInt = c(abs(1 - 1 - (tFilter/2)), 
			   1 - (tFilter/2) 
			   );
tvals = quantile(attributes(bsdall.fit.tstat)$stats[,6], tFilterInt, na.rm=T);
dmrs0 = dmrFinder(bsdall.fit.tstat, cutoff=tvals);
