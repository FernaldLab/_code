runWilcoxAndBootstrap = function(inputMatrix, groupingVec, trials=10000, replace=T, sigThresh=.05, plot=F, warnings=F)
{
	if (is.numeric(groupingVec))
	{
		groupingVec = as.logical(groupingVec);
	}
	
	group1 = inputMatrix[groupingVec, ];
	group0 = inputMatrix[!groupingVec, ];
	
	outputMatrix = as.data.frame(matrix(nrow=ncol(inputMatrix), ncol=9));
	rownames(outputMatrix) = colnames(inputMatrix);
	colnames(outputMatrix) = c('wilcoxon', 'stat', 'sig', 'bootstrap_means', 'stat', 'sig', 'bootstrap_medians', 'stat', 'sig');
	for (test in 1:ncol(inputMatrix))
	{
		cat(rownames(outputMatrix)[test], '...\n', sep='');
		#cat('Wilcoxon...\n');
		if (warnings)
		{
			wilcoxOut = wilcox.test(group1[, test], group0[, test]);
		} 
		else
		{
			wilcoxOut = suppressWarnings(wilcox.test(group1[, test], group0[, test]));
		}
		wilcoxPval = wilcoxOut$p.value
		wilcoxStat = wilcoxOut$statistic;
		
		#cat('bootstrap means...\n');
		bootstrapMeanOut = bootstrap_test(group1[, test], group0[, test], trials=trials, Func='mean', replace=replace, plot=plot);
		bootstrapMeanPval = bootstrapMeanOut$p;
		bootstrapMeanStat = signif(bootstrapMeanOut$stat, 3);
		
		#cat('bootstrap medians...\n');
		bootstrapMedianOut = bootstrap_test(group1[, test], group0[, test], trials=trials, Func='median', replace=replace, plot=plot);
		bootstrapMedianPval = bootstrapMedianOut$p;
		bootstrapMedianStat = signif(bootstrapMedianOut$stat, 3);
		
		outputMatrix[test, 1] = wilcoxPval;
		outputMatrix[test, 2] = wilcoxStat;
		outputMatrix[test, 3] = wilcoxPval < sigThresh;
		outputMatrix[test, 4] = bootstrapMeanPval;
		outputMatrix[test, 5] = bootstrapMeanStat;
		outputMatrix[test, 6] = bootstrapMeanPval < sigThresh;
		outputMatrix[test, 7] = bootstrapMedianPval;
		outputMatrix[test, 8] = bootstrapMedianStat;
		outputMatrix[test, 9] = bootstrapMedianPval < sigThresh;
	}
	return(outputMatrix);
} 


#################

bootstrap_test = function(ctrl, exp, trials, Func, replace = F, plot=T)
{
	boxNULL = c(ctrl, exp);
	statsNULL = c();
	
	stat = eval(call(Func, ctrl)) - eval(call(Func, exp));
	
	for (trial in 1:trials)
	{
		pseudo_ctrl = sample(boxNULL, length(ctrl), replace = replace);
		pseudo_exp = sample(boxNULL, length(exp), replace = replace);
		pseudo_stat_null = eval(call(Func, pseudo_ctrl)) - eval(call(Func, pseudo_exp));
		statsNULL = c(statsNULL, pseudo_stat_null);
	}

	if (stat < 0)
	{ 
		p = sum(statsNULL < stat) / trials
	} else {
		p = sum(statsNULL > stat) / trials
	}
	
	if (plot)
	{
		hist(statsNULL);
		abline(v = stat, col = 'red');
	}
	return(list(stat = stat, p = p));
}