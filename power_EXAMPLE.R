power_EXAMPLE = function(ctrl, 
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
		NULLconfint = quantile(statsNULL, 0.95, na.rm=T);
	} else {
		NULLconfint = quantile(statsNULL, c(0.025, 0.975), na.rm=T);	
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
		right_count = sum(statsALT > NULLconfint,na.rm=T);
		left_count = NULL;
		power = right_count / trials;
	} else {
		right_count = sum(statsALT > NULLconfint[2],na.rm=T);
		left_count = sum(statsALT < NULLconfint[1],na.rm=T);
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