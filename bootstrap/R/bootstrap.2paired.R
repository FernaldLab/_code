bootstrap.2paired <-
function(condition1, condition2, 
						  		 trials = 10000, 
						  		 Func = 'mean',
						  		 plots = T, 
						  		 savePlot = F,
						  		 abs.diffs = F,
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, 
						  		 conditionNames = c('condition1', 'condition2'), 
						  		 pch = 21,
						  		 verbose = T,
						  		 ...
						  		 )
{
	allParam = as.list(sys.frame(sys.nframe()));

	tmp = bootstrapUtil.checkNumericAndNAs(condition1, condition2, paired = T);
	condition1 = tmp$condition1;
	condition2 = tmp$condition2;
	
	# compute statistic
	diffs = condition1 - condition2; 
	stat = eval(call(Func, diffs));
	if (abs.diffs) {stat = abs(stat)}
	
	# store data for output
	data = data.frame(condition1, condition2, diffs);
	names(data) = c(conditionNames[1], conditionNames[2], 'diffs');
	
	# resample to create null distribution of statistic
	if (verbose) {
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: ', Func, ' of individual differences = ', stat, '\n\n', sep = '');
	}
	statsNULL = c();
	for (trial in 1:trials) {
		if (verbose) {bootstrapUtil.printRunNumber(trials = trials, trial = trial)}
		signs = sample(c(-1, 1), length(diffs), replace = T);
		pseudo_diffs = signs * diffs;
		pseudo_stat_null = eval(call(Func, pseudo_diffs));
		if (abs.diffs) {pseudo_stat_null = abs(pseudo_stat_null)}
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	if (any(is.na(statsNULL))) {stop('NA IN NULL DISTRIBUTION, CHECK YOUR DATA')}
	
	# compute p-values
	if (abs.diffs) {p = sum(statsNULL > stat) / trials; ptemp = list(p = p)}
	else {ptemp = bootstrapUtil.compute2sidePval(statsNULL = statsNULL, stat = stat, trials = trials)}
	
	if (plots) {
		par(mfrow = c(1, 3));
		
		outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2pairedGroupDiffHist))];
		data.lab = do.call(bootstrapPlot.2pairedGroupDiffHist, c(outerArgs, diffs = list(diffs), 
																 returnDataLabel = T, ...
																 ));
		
		# plot histogram of null distribution with lines at stat and reflected stat
		if (abs.diffs) {xlabNULL = paste('statistic: abs(', Func, ' of individual differences)', sep = '')}
		else {xlabNULL = paste('statistic: ', Func, ' of individual differences', sep = '')}
		
		hist(statsNULL, 
		 	 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 	 xlab = xlabNULL, col = col, border = border, breaks = floor(trials / 100),
		 	 #cex.lab = cex.lab,
			 #cex.axis = cex.axis,
			 ...);
			 
		if (abs.diffs) {abline(v = stat, col = col.line)}
		else {
			abline(v = stat, col = col.line);
			abline(v = ptemp$reflect, col = col.line, lty = 'dashed');
		}

		# boxplot of conditions with individual data points
		toPlot = data[, 1:2];
		
		if (ptemp$p == 0) {toPaste = paste('p < 1e-5 (n=', length(diffs), ')', sep = '')}
		else {toPaste = paste('p = ', signif(ptemp$p, 2), ' (n=', length(diffs), ')', sep = '')}
		
		# draw boxes
		boxplot(toPlot, ylab = data.lab,
				#medlty = medlty,
				main = toPaste, frame.plot = F,
				#cex.lab = cex.lab,
				#cex.axis = cex.axis,
				...);
				
		# add data points
		stripchart(toPlot, vertical = T, add = T, pch = pch, bg = col, cex = 1.5);
				   
		# connect paired points across conditions
		for (row in 1:nrow(data)) {segments(1, data[row, 1], 2, data[row, 2], col = border)}
	}
	
	# build output list
	if (abs.diffs) {
		output = list(stat = stat,
					  p = ptemp$p,
					  null.dist = statsNULL,
					  data = data,
					  call.param = match.call(),
					  all.param = allParam);
	}
	else {
		output = list(stat = stat, 
					  stat.reflect = ptemp$reflect,
					  p = ptemp$p, 
					  p.left = ptemp$p_left,
				  	  p.right = ptemp$p_right,
				  	  null.dist = statsNULL,
				  	  null.dist.mean = ptemp$midNULL,
				  	  data = data,
					  call.param = match.call(),
					  all.param = allParam);
	}
	if (verbose) {cat('\np = ', ptemp$p, '\n', sep = '')}
	return(output);
}
