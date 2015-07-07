dat=read.csv('behavioraldatacontrolcort-t-kt-e.csv');

for (h in 1:4) {
	box = as.numeric(c(dat[,h], dat[,h+4]));
	
	print(wilcox.test(dat[,h],dat[,h+4]))
}

wilcoxRankSumPower = function(dat1, dat2, trials=10000, replace=T) {
	box.null = as.numeric(dat1, dat2);
	stats.null = c();
	stats.alt = c();
	for (trial in 1:trials) {
		if(trial%%10000==0) {print(trial)}
		pdat1 = sample(box.null, length(dat1), replace=replace);
		pdat2 = sample(box.null, length(dat2), replace=replace);
		stats.null = c(stats.null, suppressWarnings(wilcox.test(pdat1,pdat2)$statistic[[1]]));
		
		pdat1.alt = sample(dat1, length(dat1), replace=replace);
		if (length(unique(pdat1.alt))==1) {pdat1.alt = sample(dat1, length(dat1), replace=replace)}
		pdat2.alt = sample(dat2, length(dat2), replace=replace);
		if (length(unique(pdat2.alt))==1) {pdat2.alt = sample(dat2, length(dat2), replace=replace)}
		stats.alt = c(stats.alt, suppressWarnings(wilcox.test(pdat1.alt,pdat2.alt)$statistic[[1]]));
	}
	confint.null = quantile(stats.null, c(.025, .975));
	
	right_count = sum(stats.alt > confint.null[2]);
	left_count = sum(stats.alt < confint.null[1]);
	power = (right_count + left_count) / trials;	
	
	par(mfrow = c(2, 1), oma = c(0,0,3,0), mar = c(2,4,3,2));
	hist.null = hist(stats.null, plot=F);
	hist.alt = hist(stats.alt, plot=F);
	
	ymax = max(c(hist.null$counts, hist.alt$counts));
	xmax = max(c(hist.null$breaks, hist.alt$breaks));
	xmin = min(c(hist.null$breaks, hist.alt$breaks));
	plot(hist.null, main='null', col='grey', border='darkgrey', xlab='', ylim=c(0, ymax), xlim=c(xmin, xmax));
	abline(v=confint.null[1], col='red');
	abline(v=confint.null[2], col='red');
	plot(hist.alt, main='alt', col='grey', border='darkgrey', xlab='', ylim=c(0, ymax), xlim=c(xmin, xmax));
	abline(v=confint.null[1], col='red');
	abline(v=confint.null[2], col='red');
	
	return(power);
	
}


# > wilcoxRankSumPower(dat[,1],dat[,5],trials=100000)
# [1] 10000
# [1] 20000
# [1] 30000
# [1] 40000
# [1] 50000
# [1] 60000
# [1] 70000
# [1] 80000
# [1] 90000
# [1] 100000
# [1] 0.08575
# > wilcoxRankSumPower(dat[,2],dat[,6],trials=100000)
# [1] 10000
# [1] 20000
# [1] 30000
# [1] 40000
# [1] 50000
# [1] 60000
# [1] 70000
# [1] 80000
# [1] 90000
# [1] 100000
# [1] 0.1111
# > wilcoxRankSumPower(dat[,3],dat[,7],trials=100000)
# [1] 10000
# [1] 20000
# [1] 30000
# [1] 40000
# [1] 50000
# [1] 60000
# [1] 70000
# [1] 80000
# [1] 90000
# [1] 100000
# [1] 0.3873
# > wilcoxRankSumPower(dat[,4],dat[,8],trials=100000)
# [1] 10000
# [1] 20000
# [1] 30000
# [1] 40000
# [1] 50000
# [1] 60000
# [1] 70000
# [1] 80000
# [1] 90000
# [1] 100000
# [1] 0.22089