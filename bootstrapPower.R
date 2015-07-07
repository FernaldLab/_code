bootstrapPower = function(ctrl, 
						 exp = NULL, 
						 change, 
						 Func,
						 trials = 10000,
						 replace = T,
						 abs = F,
						 plot = T,
						 breaks = NULL,
						 cex.main = 1.2,
						 ...
						 )
{
	if (is.null(exp)) {exp = change * ctrl}
	old = eval(call(Func, ctrl));
	new = eval(call(Func, exp));
	stat = eval(call(Func, ctrl)) - eval(call(Func, exp));
	if (abs) {stat = abs(stat)}
	delta = (new - old) / old;
	
	#--- null ---
	boxNULL = c(ctrl, exp);
	statsNULL = c();
	
	for (trial in 1:trials)
	{
		pseudo_ctrl = sample(boxNULL, length(ctrl), replace = replace);
		pseudo_exp = sample(boxNULL, length(exp), replace = replace);
		pseudo_stat_null = eval(call(Func, pseudo_ctrl)) - eval(call(Func, pseudo_exp));
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	
	if (abs) 
	{
		statsNULL = abs(statsNULL);
		NULLconfint = quantile(statsNULL, 0.95);
	} 
	else 
	{
		NULLconfint = quantile(statsNULL, c(0.025, 0.975));	
	}
	
	#--- alt ---
	statsALT = c();
	
	for (trial in 1:trials)
	{
		pseudo_ctrlALT = sample(ctrl, length(ctrl), replace = replace);
		pseudo_expALT = sample(exp, length(exp), replace = replace);
		pseudo_stat_alt = eval(call(Func, pseudo_ctrlALT)) - eval(call(Func, pseudo_expALT));
		statsALT = c(statsALT, pseudo_stat_alt);
	}
	
	if (abs) 
	{
		statsALT = abs(statsALT);
		right_count = sum(statsALT > NULLconfint);
		left_count = NULL;
		power = right_count / trials;
	} 
	else 
	{
		right_count = sum(statsALT > NULLconfint[2]);
		left_count = sum(statsALT < NULLconfint[1]);
		power = (right_count + left_count) / trials;
	}
	
	if (plot)
	{
		if (is.null(breaks)) {breaks = trials / 2}
		par(mfrow = c(2, 1), oma = c(0,0,3,0), mar = c(2,4,3,2));
		histNULL = hist(statsNULL, breaks = breaks, plot = F);
		histALT = hist(statsALT, breaks = breaks, plot = F);
		ymax = max(c(histNULL$counts, histALT$counts));
		xmax = max(c(histNULL$breaks, histALT$breaks));
		xmin = min(c(histNULL$breaks, histALT$breaks));
		
		#hist(statsNULL, breaks = breaks, main = 'NULL', col = 'grey', border = 'grey', xlab = '', ...);
		plot(histNULL, main = 'null distribution', col = 'grey', border = 'darkgrey', xlab = '', ylim = c(0, ymax), xlim = c(xmin, xmax), cex.main = cex.main*.9, ...);
		if (abs)
		{
			abline(v = NULLconfint, col = 'red', lty = 'dashed');
		} else {
			abline(v = NULLconfint[1], col = 'red', lty = 'dashed');
			abline(v = NULLconfint[2], col = 'red', lty = 'dashed');
		}
		#hist(statsALT, breaks = breaks, main = 'ALT', col = 'grey', border = 'grey', xlab = '', ...);
		plot(histALT, main = 'alternate distribution', col = 'grey', border = 'darkgrey', xlab = '', ylim = c(0, ymax), xlim = c(xmin, xmax), cex.main = cex.main*.9, ...);
		if (abs)
		{
			abline(v = NULLconfint, col = 'red', lty = 'dashed');
		} else {
			abline(v = NULLconfint[1], col = 'red', lty = 'dashed');
			abline(v = NULLconfint[2], col = 'red', lty = 'dashed');
		}
		title(main = paste('Test statistic: difference of group ', Func, 's\npower = ', signif(power, 2), sep = ''), outer = T, cex.main = cex.main);
	}
	
	return(list(ctrl = ctrl,
				exp = exp,
				stat = stat,
				delta = delta,
				confint = NULLconfint,
				rightcount = right_count,
				leftcount = left_count,
				power = power,
				null.dist = statsNULL,
				alt.dist = statsALT
				)
		  )
}
##
#
##
bootstrapPower.wilcox = function(group1, group2 = NULL, change, trials = 10000, replace = T, verbose = T, plot = T, cex.main = 1.2, ...)
{
	if (is.null(group2))
	{
		if (is.numeric(change)) {group2 = change * group1}
	}  

	W = wilcoxStat(group1, group2);
	
	#--- null ---
	boxNULL = c(group1, group2);
	statsNULL = c();
	if (verbose) {cat('Working on null distribution...\n')}
	for (trial in 1:trials)
	{
		if (verbose) {printRunNumber(trials = trials, trial = trial)}
		pseudo_gp1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_gp2 = sample(boxNULL, length(group2), replace = replace);
		pseudoW.NULL = wilcoxStat(pseudo_gp1, pseudo_gp2);
		statsNULL = c(statsNULL, pseudoW.NULL);
	}
	confintNULL = quantile(statsNULL, c(0.025, 0.975));	
	
	#--- alt ---
	statsALT = c();
	if (verbose) {cat('Working on alternate distribution...\n')}
	for (trial in 1:trials)
	{
		if (verbose) {printRunNumber(trials = trials, trial = trial)}
		pseudo_gp1 = sample(group1, length(group1), replace = replace);
		pseudo_gp2 = sample(group2, length(group2), replace = replace);
		pseudoW.ALT = wilcoxStat(pseudo_gp1, pseudo_gp2);
		statsALT = c(statsALT, pseudoW.ALT);
	}
	confintALT = quantile(statsALT, c(0.025, 0.975));
	
	right_count = sum(statsALT > confintNULL[2]);
	left_count = sum(statsALT < confintNULL[1]);
	power = (right_count + left_count) / trials;
	
	if (plot)
	{
		bootstrapPlot.powerHistograms(statsNULL = statsNULL, statsALT = statsALT, 
									  breaks = NULL, trials = trials, 
									  NULLconfint = confintNULL, Func = 'Wilcoxon', 
									  power = power, cex.main, abs = F, col, border, col.line, lty, 
									  title.main = 'Wilcoxon W'
									  );
	}
	
	return(list(group1 = group1, group2 = group2, stat = W, conf.int = confintNULL, rightcount = right_count, leftcount = left_count, power = power, null.dist = statsNULL, alt.dist = statsALT));
}