rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');

load(file='aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombined.RData');


for (ns in c(20,50)) {
	for (h in c(250,500)) {
		for (maxGap in c(1000,5000)) {
			cat('=============================\n', sep='');
			cat('==================\n', sep='');
			cat('===== Smoothing...\n', sep='');
			cat('        ns=', ns, ', h=', h, ', maxGap=', maxGap, '\n', sep='');
			bsdall.fit = BSmooth(bsdall, ns=ns, h=h, maxGap=maxGap, mc.cores=4, parallelBy='sample');gc();
			filename = paste('aligned.adapters.q30.m0.methratio_CpGcombined.CG.cleanCombinedSmoothed_ns', ns,
							 '_h', h,
							 '_maxGap', maxGap,
							 '.RData',
							 sep=''
							 );
			cat('        writing file: ', filename, '\n', sep='');
			save(bsdall.fit, file=filename);gc();
		}
	}
}
rm(ns,h,maxGap,filename);


#####################################################################################
rm(list=ls());
filesSmoothed = list.files()[grep('maxGap',list.files())];

for (fs in filesSmoothed) {
	cat('=============================\n', sep='');
	print(fs)
	load(fs);
	print(bsdall.fit);
	bsdall.fitCov = getCoverage(bsdall.fit);
	for (reqCov in c(2,5)) {
		cat('==================\n', sep='');
		cat('Using only sites with ', reqCov, 'x coverage\n', sep='');
		keep = which(apply(bsdall.fitCov, 1, function(f) sum(f>=reqCov)==4));
		bsdall.fit = bsdall.fit[keep,];gc();
		bsdall.fit.tstat = BSmooth.tstat(bsdall.fit, 
								 		 group1=c('3165_BRISCOE','3581_LYNLEY'),
								 		 group2=c('3157_TENNISON','3677_MONK'),
									     estimate.var='same',
								 		 local.correct=T,
		     					 		 mc.cores=4,
								 		 verbose=T
								 		 );
		filename = gsub('.RData', '', fs);
		filename = paste(filename, '_', reqCov, 'x.tstat_estvarSame.RData', sep='');
		cat('        writing file: ', filename, '\n', sep='');
		save(bsdall.fit.tstat, file=filename);gc();
	}
}; rm(fs,reqCov,keep,bsdall.fitCov, filename);
