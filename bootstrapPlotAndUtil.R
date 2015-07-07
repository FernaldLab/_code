##########################
### Plotting functions ###
##########################
##
#
##
bootstrapPlot.nullDist = function(stat, abs.diffs, reflect, statsNULL, trials, Func, replace, col, border, col.line)
{
	hist_null = hist(statsNULL, breaks = floor(trials / 100), plot = F);
	plot(hist_null,
		 xlim = c(min(hist_null$breaks), max(hist_null$breaks)), 
		 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 xlab = paste('statistic (', Func, 's, replace = ', substr(replace, 1, 1), ')', sep = ''),
		 col = col, border = border
		 );

	if (abs.diffs)
	{
		abline(v = stat, col = col.line);
	}
	else if (!abs.diffs)
	{
		abline(v = stat, col = col.line);
		abline(v = reflect, col = col.line, lty = 'dashed');
	}
}
##
#
##
bootstrapPlot.confIntDist = function(pseudo, Func, trials, value, conf.int, col, line.col, border)
{
	hist(pseudo, 
	     breaks = length(pseudo)/100, 
	     col = col, border = border, 
	     main = paste('bootstrapped ', Func, '\n(n=', trials, ')', sep = ''), 
	     xlab = ''
	     );
	abline(v = value, col = line.col);
	abline(v = conf.int[1], col = line.col, lty = 'dashed');
	abline(v = conf.int[2], col = line.col, lty = 'dashed');
}
##
#
##
bootstrapPlot.groupHist = function(group, groupName, Func, dataDescriptor = NULL, lims = NULL, col, border, col.line)
{
	# add n to group name
	groupName = paste(groupName, ' (n=', length(group), ')', sep = '');
	
	# compute histogram
	hist_group = hist(group, breaks = length(group), plot = F);
	
	if (is.null(lims))
	{
		# compute axis limits
		ymax = max(hist_group$counts);
		xmin = min(hist_group$breaks);
		xmax = max(hist_group$breaks);
	}
	else if (!is.null(lims))
	{
		if (is.list(lims) & is.numeric(unlist(lims)))
		{
			ymax = lims$ymax;
			xmin = lims$xmin;
			xmax = lims$xmax;
		}
	}
	
	# check for dataDescriptor 
	data.lab = checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# plot group histogram with line at mean/median
	plot(hist_group, 
		 main = groupName, 
		 ylim = c(0, ymax), xlim = c(xmin, xmax), 
		 col = col, border = border, 
		 xlab = data.lab
		 );
	if (Func %in% c('mean', 't')) {v = mean(group)}
	else if (Func %in% c('median', 'cv', 'wilcox', 'ks')) {v = median(group)}
	#abline(v = eval(call(Func, group)), col = col.line);	
	abline(v = v, col = col.line);
}
##
#
##
bootstrapPlot.groupStripchart = function(data, interval, value, conf.int, Func, line.col, pch, col)
{
	stripchart(data, vertical = T, 
			   pch = pch, bg = col, 
			   frame.plot = F, 
			   main = paste(interval*100, '% conf.int around ', round(value, 1), 
			   				'\n[', round(conf.int[1], 1), ', ', 
			   				round(conf.int[2], 1), ']', 
			   				sep = ''
			   				)
			   );
	segments(.8, value, 1.2, value, col = line.col);
	text(1.26, value, Func);
	segments(.9, conf.int[1], 1.1, conf.int[1], col = line.col, lty = 'dashed');
	segments(.9, conf.int[2], 1.1, conf.int[2], col = line.col, lty = 'dashed');
}
##
#
##
bootstrapPlot.2pairedGroupDiffHist = function(dataDescriptor = NULL, diffs, Func, conditionNames, col, border, cex.lab, cex.axis, col.line, returnDataLabel = T)
{
	# check for dataDescriptor 
	if (!is.null(dataDescriptor) & is.character(dataDescriptor))
	{
		data.lab = dataDescriptor;
	}
	else
	{
		data.lab = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
	}
	
	# compute histogram of individual differences
	hist.diffs = hist(diffs, breaks = length(diffs), plot = F);
	xmax = max(abs(hist.diffs$breaks));

	# plot histogram
	plot(hist.diffs,
		 xlim = c(-xmax, xmax),
		 main = paste('Individual differences with ', Func, '\n(n=', length(diffs), ')', sep = ''),
		 xlab = paste(data.lab, 
		 			  ' (', conditionNames[1], ' - ', conditionNames[2], ')',
		 			   sep = ''
		 			   ),
		 col = col, border = border,
		 cex.lab = cex.lab,
		 cex.axis = cex.axis
		 );
	abline(v = mean(diffs), col = col.line);
	if (returnDataLabel) {return(data.lab)}
}
##
#
##
bootstrapPlot.2groupBoxplot = function(group1, group2, groupNames, paired, p, dataDescriptor=NULL, col, border, pch, jitter, confidenceIntervals, Func)
{
	toPlot = c(group1, group2);
	grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], length(group2)));
	boxLineType = Func;
	if (boxLineType == 'mean') {medlty = 'blank'}
	else {medlty = 'solid'}
	toPaste = getPvalToPaste(p);
	data.lab = checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# draw boxes
	boxplot(toPlot ~ grp,
			ylab = data.lab,
			medlty = medlty,
			main = toPaste
			);
	
	# draw mean lines if specified
	if (boxLineType == 'mean')
	{
		segments(0.6, mean(group1), 1.4, mean(group1), lwd = 3);
		segments(1.6, mean(group2), 2.4, mean(group2), lwd = 3);
	}
	
	# add data points
	stripchart(toPlot ~ grp,
			   vertical = T,
			   add = T,
			   method = 'jitter',
			   jitter = jitter,
			   pch = pch,
			   bg = col,
			   cex = 1.5
			   );
			   
	if (paired)
	{
		data = data.frame(group1, group2);
		# connect paired points across conditions
		for (row in 1:nrow(data))
		{
			segments(1, data[row, 1], 2, data[row, 2], col = border);
		}
	}
	
	# if (confidenceIntervals)
	# {
		# if (verbose)
		# {
			# cat('  Computing 95% confidence intervals\n');
		# }
		# values = c();
		# ints = matrix(nrow = ncol(data), ncol = 2);
		# for (gp in 1:ncol(data))
		# {
			# temp = confidenceInterval(data[, gp], plots = F, Func = Func);###NEED TO DEFINE data
			# values = c(values, temp$value);
			# ints[gp, ] = temp$conf.int;
		# }
		# coMat = matrix(nrow = ncol(data), ncol = 2);
		# coMat[1, ] = c(.75, 1.25);
		# for (gp in 2:nrow(coMat))
		# {
			# coMat[gp, ] = coMat[gp - 1, ] + 1;
		# }
		# for (gp in 1:length(values))
		# {
			# segments(coMat[gp, 1], values[gp], 
				     # coMat[gp, 2], values[gp], 
				     # lwd = 1, col = line.col
				     # );
			# segments(coMat[gp, 1], ints[gp, 1], 
					 # coMat[gp, 2], ints[gp, 1], 
					 # lwd = 1, lty = 'dashed', col = line.col
					 # );
			# segments(coMat[gp, 1], ints[gp, 2], 
				     # coMat[gp, 2], ints[gp, 2], 
				     # lwd = 1, lty = 'dashed', col = line.col
				     # );
		# }
	# }
}
##
#
##
bootstrapPlot.NgroupBoxplot = function(data, groupNames = NULL, points = T, confidenceIntervals = T, Func, col, pch, verbose)
{
	if (is.null(groupNames))
	{
		groupNames = colnames(data);
	}
	else if (!is.null(groupNames) & length(groupNames)==ncol(data))
	{
		groupNames = groupNames;
	}
	boxplot(data, 
		    names = groupNames, 
		    main = paste('p = ', p, sep = ''), 
		    border = 'lightgrey', 
		    medlty = 'blank', 
		    boxwex = .5
		    );
		
	if (points)
	{
		stripchart(as.data.frame(data), vertical = T, add = T, pch = pch, bg = col);
	}
	
	if (confidenceIntervals)
	{
		if (verbose)
		{
			cat('  Computing 95% confidence intervals\n');
		}
		values = c();
		ints = matrix(nrow = ncol(data), ncol = 2);
		for (gp in 1:ncol(data))
		{
			temp = confidenceInterval(data[, gp], plots = F, Func = Func);
			values = c(values, temp$value);
			ints[gp, ] = temp$conf.int;
		}
		coMat = matrix(nrow = ncol(data), ncol = 2);
		coMat[1, ] = c(.75, 1.25);
		for (gp in 2:nrow(coMat))
		{
			coMat[gp, ] = coMat[gp - 1, ] + 1;
		}
		for (gp in 1:length(values))
		{
			segments(coMat[gp, 1], values[gp], 
				     coMat[gp, 2], values[gp], 
				     lwd = 1, col = line.col
				     );
			segments(coMat[gp, 1], ints[gp, 1], 
					 coMat[gp, 2], ints[gp, 1], 
					 lwd = 1, lty = 'dashed', col = line.col
					 );
			segments(coMat[gp, 1], ints[gp, 2], 
				     coMat[gp, 2], ints[gp, 2], 
				     lwd = 1, lty = 'dashed', col = line.col
				     );
		}
	}
}
##
#
##
bootstrapPlot.powerHistograms = function(statsNULL, statsALT, breaks = NULL, trials, NULLconfint, Func, power, cex.main, abs = F, col, border, col.line, lty, title.main = NULL, ...)
{
	if (is.null(breaks)) {breaks = trials / 2}
	par(mfrow = c(2, 1), oma = c(0,0,3,0), mar = c(2,4,3,2));
	histNULL = hist(statsNULL, breaks = breaks, plot = F);
	histALT = hist(statsALT, breaks = breaks, plot = F);
	ymax = max(c(histNULL$counts, histALT$counts));
	xmax = max(c(histNULL$breaks, histALT$breaks));
	xmin = min(c(histNULL$breaks, histALT$breaks));
	
	#hist(statsNULL, breaks = breaks, main = 'NULL', col = 'grey', border = 'grey', xlab = '', ...);
	plot(histNULL, 
		 main = 'NULL', 
		 col = 'grey', border = 'darkgrey', 
		 xlab = '', 
		 ylim = c(0, ymax), xlim = c(xmin, xmax)
		 );
	if (abs)
	{
		abline(v = NULLconfint, col = 'red', lty = 'dashed');
	} else {
		abline(v = NULLconfint[1], col = 'red', lty = 'dashed');
		abline(v = NULLconfint[2], col = 'red', lty = 'dashed');
	}
	#hist(statsALT, breaks = breaks, main = 'ALT', col = 'grey', border = 'grey', xlab = '', ...);
	plot(histALT,
	     main = 'ALT', 
	     col = 'grey', border = 'darkgrey', 
	     xlab = '', ylim = c(0, ymax), xlim = c(xmin, xmax)
	     );
	if (abs)
	{
		abline(v = NULLconfint, col = 'red', lty = 'dashed');
	} else {
		abline(v = NULLconfint[1], col = 'red', lty = 'dashed');
		abline(v = NULLconfint[2], col = 'red', lty = 'dashed');
	}
	if (is.null(title.main))
	{
		toPaste = paste('Test statistic: difference of group ', Func, 's\npower = ', signif(power, 2), sep = '');
	}
	else 
	{
		toPaste = paste(title.main, '\npower = ', signif(power, 2), sep  = '');
	}
	title(main = toPaste, outer = T, cex.main = cex.main);
}
##
#
##
#################
### Utilities ###
#################
##
#
##
compute2sidePval = function(statsNULL, stat, trials)
{
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	#if (stat < 0)
	if (stat < midNULL)
	{ 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} 
	#else if (stat > 0)
	else if (stat > midNULL)
	{
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	}
	else
	{
		stop('Either statistic==0 or midNULL or something else is wrong...\n  Figure it out or find Austin!');
	}
	return(list(midNULL = midNULL, reflect = reflect, p_left = p_left, p_right = p_right, p = p))
}
##
#
##
printRunNumber = function(trials, trial)
{
	if (trials <= 20000)
	{
		if (trial %% 1000 == 0) {cat('  Run ', trial, '\n', sep = '')}
	}
	else if (trials > 20000)
	{
		if (trial %% 10000 == 0) {cat('  Run ', trial, '\n', sep = '')}
	}
}
##
#
##
checkYesOrNo = function(answer)
{
	check = tolower(answer) %in% c('y', 'yes', 'n', 'no');
	if (!check) 
	{
		while (!check) 
		{
			answer = tolower(readline(prompt = '  Please answer yes or no'));
			check = tolower(answer) %in% c('y', 'yes', 'n', 'no');
		}	
	}
	if (tolower(answer) %in% c('y', 'yes')) 
	{
		answer = 'y';
	} 
	else if (tolower(answer) %in% c('n', 'no'))
	{
		answer = 'n';
	}
	return(answer);
}
##
#
##
match2histAxisLims = function (group1, group2)
{
	hist_group1 = hist(group1, breaks = length(group1), plot = F);
	hist_group2 = hist(group2, breaks = length(group2), plot = F);
	
	# compute axis limits
	ymax = max(c(hist_group1$counts, hist_group2$counts));
	xmin = min(c(hist_group1$breaks, hist_group2$breaks));
	xmax = max(c(hist_group1$breaks, hist_group2$breaks));
	
	return(list(ymax = ymax, xmin = xmin, xmax = xmax));
}
##
#
##
checkDataDescriptorToGetLabel = function (dataDescriptor)
{
	if (!is.null(dataDescriptor) & is.character(dataDescriptor))
	{
		data.lab = dataDescriptor;
	}
	else
	{
		data.lab = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
	}
	return(data.lab);
}
##
#
##
getPvalToPaste = function(p)
{
	if (p == 0) {toPaste = paste('p < 1e-5', sep = '')}
	else {toPaste = paste('p = ', signif(p, 2), sep = '')}
	return(toPaste);
}
##
#
##
checkNumericAndNAs = function(group1, group2)
{
	not_numeric = !is.numeric(c(group1, group2));
	NA_check = sum(is.na(c(group1, group2))) > 0;
	if (not_numeric)
	{
		stop('DATA CONTAINS NON-NUMERIC VALUES...\n   I ONLY EAT NUMBERS!!!\n    GIVE ME NUMBERS!!!!!');
	}
	if (NA_check)
	{
		group1 = group1[!is.na(group1)];
		group2 = group2[!is.na(group2)];
		boxNULL = c(group1, group2);
		warning('NAs removed from one or both groups, check your data');
	}
	else
	{
		boxNULL = c(group1, group2);
	}
	
	return(list(group1 = group1, group2 = group2, boxNULL = boxNULL));
}
##
#
##
tStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(t.test(x = group1, y = group2,
											 alternative = 'two.sided',
											 paired = F,
											 var.equal = F
											)$statistic
									  )
						   );
	return(stat);
}
##
#
##
wilcoxStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(wilcox.test(x = group1, y = group2, 
											      exact = T, 
											      alternative = 'two.sided'
											     )$statistic
									  )
							);
	return(stat);
}
##
#
##
ksStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(ks.test(x = group1, y = group2, alternative= 'two.sided'
											 )$statistic
									 )
						   );
	return(stat);
}
##
#
##
cv = function(data)
{
	value = sd(as.numeric(data), na.rm = T) / mean(as.numeric(data), na.rm = T);
	return(value);	
}
##
#
##
checkForNegativeDataIfCV = function(group1, group2, Func)
{
	if (Func == 'cv')
	{
		if (any(c(group1, group2) < 0))
		{
			stop('CV IS NONSENSICAL FOR DATA WITH NEGATIVE VALUES\n ...USE MEAN OR MEDIAN');
		}
	}
}
##
#
##
computeStat2independent = function(group1, group2, Func = c('mean.diff', 'median.diff', 'cv.diff', 't.test', 'wilcox.test', 'ks.test'), verbose, dataDescriptor, replace)
{
	choiceVec = c('mean.diff', 'median.diff', 'cv.diff', 't.test', 'wilcox.test', 'ks.test');
	title = 'Please select one of the following statistics to use:';
	err = 'SOMETHING IS WRONG IN computeStatOther2independent(), GET AUSTIN';
	if (length(Func) > 1)
	{
		checkArg = menu(choices = choiceVec, title = title);
		Func = choiceVec[checkArg];
		cat('Using ', Func, ' statistic...\n', sep = '');
	}
	else if (length(Func) == 1) 
	{
		checkArg = pmatch(Func, choiceVec);
		Func = choiceVec[checkArg];
	}
	else {stop(err)}
	
	if (is.na(checkArg)) 
	{
		checkArg = menu(choices = choiceVec, title = title);print(checkArg)
		Func = choiceVec[checkArg];
	}
	
	if (checkArg %in% 1:3)
	{
		Func = gsub('.diff', '', Func);
		checkForNegativeDataIfCV(group1, group2, Func);
		stat = eval(call(Func, group1)) - eval(call(Func, group2));
		if (verbose)
		{
			cat('...........................................\n');
			cat('Testing: ', dataDescriptor, '\n\n', sep = '');
			cat('Building null distribution...\n');
			cat(' Test statistic: difference of group ', Func, 's\n', sep = '');
			cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
		}
	}
	else if (checkArg %in% 4:6)
	{
		Func = gsub('.test', '', Func);
		if (verbose)
		{
			cat('...........................................\n');
			cat('Testing: ', dataDescriptor, '\n\n', sep = '');
			cat('Building null distribution...\n');
			cat(' Test statistic: ', Func, '\n', sep = '');
			cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
		}
		if (checkArg == 4) {stat = tStat(group1, group2)}
		else if (checkArg == 5) {stat = wilcoxStat(group1, group2)}
		else if (checkArg == 6) {stat = ksStat(group1, group2)}
	}
	else {stop(err)}
	
	return(list(stat = stat, Func = Func));
}
##
#
##
