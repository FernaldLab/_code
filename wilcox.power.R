.wilcox.power = function(group1, group2, replace=T,change, trials = 10000, cex.main=1.2, verbose=T,...) {
	if (is.null(group2)) {
		if (is.numeric(change)) {group2 = change * group1}
	} 
	W = wilcox.test(group1, group2)$statistic;
	
	#--- null ---
	boxNULL = c(group1, group2);
	statsNULL = c();
	if (verbose) {cat('Working on null distribution...\n')}
	for (trial in 1:trials)
	{
		#if (verbose) {printRunNumber(trials = trials, trial = trial)}
		pseudo_gp1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_gp2 = sample(boxNULL, length(group2), replace = replace);
		pseudoW.NULL = wilcox.test(pseudo_gp1, pseudo_gp2)$statistic;
		statsNULL = c(statsNULL, pseudoW.NULL);
	}
	confintNULL = quantile(statsNULL, c(0.025, 0.975));	

	#--- alt ---
	statsALT = c();
	if (verbose) {cat('Working on alternate distribution...\n')}
	for (trial in 1:trials)
	{
		#if (verbose) {printRunNumber(trials = trials, trial = trial)}
		pseudo_gp1 = sample(group1, length(group1), replace = replace);
		pseudo_gp2 = sample(group2, length(group2), replace = replace);
		pseudoW.ALT = wilcox.test(pseudo_gp1, pseudo_gp2)$statistic;
		statsALT = c(statsALT, pseudoW.ALT);
	}
	confintALT = quantile(statsALT, c(0.025, 0.975));
	
	right_count = sum(statsALT > confintNULL[2]);
	left_count = sum(statsALT < confintNULL[1]);
	power = (right_count + left_count) / trials;
	
	histNULL = hist(statsNULL,breaks=length(statsNULL)/2,plot=F);
	histALT = hist(statsALT,breaks=length(statsALT)/2,plot=F);
	ymax = max(c(histNULL$counts, histALT$counts));
	xmax = max(c(histNULL$breaks, histALT$breaks));
	xmin = min(c(histNULL$breaks, histALT$breaks));
	par(mfrow=c(2,1))
	plot(statsNULL,main='NULL',col='grey',border='grey', ylim = c(0, ymax), xlim = c(xmin, xmax));
	abline(v=confintNULL[1],col='red');
	abline(v=confintNULL[2],col='red');
	plot(statsALT,main='ALT',col='grey',border='grey', ylim = c(0, ymax), xlim = c(xmin, xmax));
	abline(v=confintALT[1],col='red');
	abline(v=confintALT[2],col='red');
	return(power)

}
