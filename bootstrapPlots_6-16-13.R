###########################################################################
###### plot bootstrap distribution for computing confidence interval ######
###########################################################################
##
###

bootstrapPlot.confIntDist = function(pseudo, Func, trials, value, conf.int, plotParams, ...)
{
	col = plotParams$col;
	border = plotParams$border;
	line.col = plotParams$line.col;
	hist(pseudo, breaks = length(pseudo)/100, 
	     col = col, border = border, 
	     main = paste('bootstrapped ', Func, '\n(n=', trials, ')', sep = ''), 
	     xlab = '', ...);
	abline(v = value, col = line.col);
	abline(v = conf.int[1], col = line.col, lty = 'dashed');
	abline(v = conf.int[2], col = line.col, lty = 'dashed');
}

########################################################################
###### plot stripchart with confidence interval around descriptor ######
########################################################################
##
###

bootstrapPlot.groupStripchart = function(data, interval, value, conf.int, Func, plotParams, ...)
{
	line.col = plotParams$line.col;
	pch = plotParams$pch;
	col = plotParams$col;
	stripchart(data, vertical = T, 
			   pch = pch, bg = col, 
			   frame.plot = F, 
			   main = paste(interval*100, '% conf.int around ', round(value, 1), 
			   				'\n[', round(conf.int[1], 1), ', ', 
			   				round(conf.int[2], 1), ']', 
			   				sep = ''
			   				), ...
			   );
	segments(.8, value, 1.2, value, col = line.col);
	text(1.26, value, Func, ...);
	segments(.9, conf.int[1], 1.1, conf.int[1], col = line.col, lty = 'dashed');
	segments(.9, conf.int[2], 1.1, conf.int[2], col = line.col, lty = 'dashed');
}

###########################################################################
###### plot  ######
###########################################################################
##
###

bootstrapPlot.groupHist = function(group, groupName, Func, dataDescriptor = NULL, lims = NULL, col, border, col.line, ...)
{
	# add n to group name
	groupName = paste(groupName, ' (n=', length(group), ')', sep = '');
	
	# compute histogram
	hist_group = hist(group, breaks = length(group), plot = F);
	
	if (is.null(lims)) {
		# compute axis limits
		ymax = max(hist_group$counts);
		xmin = min(hist_group$breaks);
		xmax = max(hist_group$breaks);
	}
	else if (!is.null(lims)) {
		if (is.list(lims) & is.numeric(unlist(lims))) {
			ymax = lims$ymax;
			xmin = lims$xmin;
			xmax = lims$xmax;
		}
	}
	
	# check for dataDescriptor 
	data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# plot group histogram with line at mean/median
	plot(hist_group, 
		 main = groupName, 
		 ylim = c(0, ymax), xlim = c(xmin, xmax), 
		 col = col, border = border, 
		 xlab = data.lab, ...);
	if (Func %in% c('mean', 't')) {v = mean(group)}
	else if (Func %in% c('median', 'cv', 'wilcox', 'ks')) {v = median(group)}
	#abline(v = eval(call(Func, group)), col = col.line);	
	abline(v = v, col = col.line);
}

###########################################################################
###### plot  ######
###########################################################################
##
###

bootstrapPlot.nullDist = function(stat, abs.diffs, reflect, statsNULL, trials, Func, replace, col, border, col.line, ...)
{
	hist_null = hist(statsNULL, breaks = floor(trials / 100), plot = F);
	
	xmax = max(max(hist_null$breaks), abs(min(hist_null$breaks)));
	if (all(statsNULL > 0)) {xlim = c(0, xmax)}
	else if (all(statsNULL < 0)) {xlim = c(-xmax, 0)}
	else {xlim = c(-xmax, xmax)}
	
	if (Func %in% c('mean', 'median', 'cv')) {
		xlab = paste('statistic (', Func, '.diff, replace = ', substr(replace, 1, 1), ')', sep = '');
	}
	else {
		xlab = paste('statistic (', Func, ', replace = ', substr(replace, 1, 1), ')', sep = '');
	}
	 
	plot(hist_null, xlim = xlim, xlab = xlab,
		 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
		 col = col, border = border, ...);

	if (abs.diffs) {abline(v = stat, col = col.line)}
	else {
		abline(v = stat, col = col.line);
		abline(v = reflect, col = col.line, lty = 'dashed');
	}
}

###########################################################################
###### plot  ######
###########################################################################
##
###

bootstrapPlot.2groupBoxplot = function(group1, group2, groupNames, paired, p, dataDescriptor=NULL, col, border, pch, jitter, Func, ...)
{
	# add n to group names
	groupNames[1] = paste(groupNames[1], ' (n=', ( length(group1)-sum(is.na(group1)) ), ')', sep = '');
	groupNames[2] = paste(groupNames[2], ' (n=', ( length(group2)-sum(is.na(group2)) ), ')', sep = '');
	
	toPlot = c(group1, group2);
	grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], length(group2)));
	boxLineType = Func;
	if (boxLineType == 'mean') {medlty = 'blank'}
	else {medlty = 'solid'}
	toPaste = bootstrapUtil.getPvalToPaste(p);
	data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# draw boxes
	boxplot(toPlot ~ grp,
			ylab = data.lab,
			medlty = medlty,
			main = toPaste,...);
	# draw mean lines if specified
	if (boxLineType == 'mean') {
		segments(0.6, mean(group1,na.rm=T), 1.4, mean(group1,na.rm=T), lwd = 3, col = 'black');
		segments(1.6, mean(group2,na.rm=T), 2.4, mean(group2,na.rm=T), lwd = 3);
	}
	
	# add data points
	stripchart(toPlot ~ grp,
			   vertical = T,
			   add = T,
			   method = 'jitter',
			   jitter = jitter,
			   pch = pch,
			   bg = col,
			   #cex = 1.5, 
			   ...);
			   
	if (paired) {
		data = data.frame(group1, group2);
		# connect paired points across conditions
		for (row in 1:nrow(data)) {segments(1, data[row, 1], 2, data[row, 2], col = border)}
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

###########################################################################
###### plot  ######
###########################################################################
##
###

bootstrapPlot.2pairedGroupDiffHist = function(dataDescriptor = NULL, diffs, Func, conditionNames, col, border, col.line, returnDataLabel = T, ...)
{
	# check for dataDescriptor 
	if (!is.null(dataDescriptor) & is.character(dataDescriptor)) {data.lab = dataDescriptor}
	else {
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
		 #cex.lab = cex.lab,
		 #cex.axis = cex.axis,
		 ...);
	abline(v = mean(diffs), col = col.line);
	if (returnDataLabel) {return(data.lab)}
}
