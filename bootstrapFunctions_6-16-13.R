source('/Volumes/fishstudies/_code/bootstrapUtils_6-16-13.R');
source('/Volumes/fishstudies/_code/bootstrapPlots_6-16-13.R');

###########################################################
###### compute confidence interval ########################
###########################################################
# bootstrap.confidenceInterval()
# Depends on:	bootstrapUtil.printRunNumber()
#				bootstrapPlot.confIntDist()
#				bootstrapPlot.groupStripchart()
## Arguments
### data:		numeric vector
### interval:	number defining confidence interval
### trials:		number of resampling runs
### Func:		char string of descriptor around which confidence interval will be computed
### plots:		logical indicating whether plot(s) is desired
### showDist:	logical indicating whether to plot bootstrap distribution or only stripchart
### savePlot:	logical indicating whether to save plot as jpeg
### plotParams:	list made character strings (col, border, line.col) and a number (pch) defining plot parameters
### verbose:	logical indicating whether info printed to screen during computation

bootstrap.confidenceInterval = function(data, interval = 0.95, trials = 10000, Func = 'mean', 
									    plots = T, showDist = T, savePlot = F,
									    plotParams = list(col = 'grey', border = 'darkgrey', line.col = 'red', pch = 21),
									    verbose = T,
									    ...
									    )
{
	# check data is a numeric vector
	if (!is.numeric(data) | !is.vector(data)) {stop('CHECK THAT DATA IS A NUMERIC VECTOR')}
	
	# check for and remove NAs
	if (any(is.na(data))) {
		data0 = data;
		data = data[!(is.na(data))];
		warning(paste(sum(is.na(data0)), ' NAs detected in data', sep = ''));
	}
	
	# compute actual value
	value = eval(call(Func, data));
	if (verbose) {
		cat('....................................................................\n');
		cat('Computing ', interval*100, '% confidence interval around the ', Func, '...\n', sep = '');
		cat(' Resampling...\n');
	}
	
	# bootstrap
	pseudo = c();
	for (trial in 1:trials) {
		# print run number to console
		if (verbose) {bootstrapUtil.printRunNumber(trials = trials, trial = trial)}
		
		# resample data with replacement to create pseudo data
		pseudoData = sample(data, length(data), replace = T);
		
		# compute pseudo value
		pseudo = c(pseudo, eval(call(Func, pseudoData)));
	}
	
	# compute confidence interval of pseudo value distribution
	half = (1 - interval) / 2;
	int = c(half, interval + half);
	conf.int = quantile(pseudo, int);

	# plot stripchart and pseudo value distribution
	if (plots) {
		if (showDist) {
			if (savePlot) {
				jpeg(file = paste(deparse(match.call()$data), '.jpg', sep = ''), 
				 	 width = 8, height = 4, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			par(mfrow = c(1, 2));
			bootstrapPlot.confIntDist(pseudo = pseudo, Func = Func, trials = trials, 
									  value = value, conf.int = conf.int, plotParams = plotParams, ...);
			bootstrapPlot.groupStripchart(data = data, interval = interval, value = value, 
										  conf.int = conf.int, Func = Func, plotParams = plotParams, ...);
		}
		else {
			if (savePlot) {
				jpeg(file = paste(deparse(match.call()$data), '.jpg', sep = ''), 
				 	 width = 4, height = 4, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			bootstrapPlot.groupStripchart(data = data, interval = interval, value = value, 
										  conf.int = conf.int, Func = Func, plotParams = plotParams, ...);
		}
		if (savePlot) {dev.off()}
	}
	output = list(value = value, conf.int = conf.int, dist = pseudo, call.param = match.call());
	return(output);
}

bootstrap.multiConfidenceIntervals = function(data, Func, groupNames=NULL)
{
	if (is.null(groupNames)) {groupNames = colnames(data)}
	values = c();
	ints = matrix(nrow = ncol(data), ncol = 2, dimnames = list(groupNames, c('lower','upper')));
	
	for (gp in 1:ncol(data)) {
		cat('  ', groupNames[gp], '\n', sep = '');
		temp = bootstrap.confidenceInterval(data[, gp], plots = F, Func = Func, verbose=F);
		values = c(values, temp$value);
		ints[gp, ] = temp$conf.int;
	}
	return(list(ints=ints, values=values));
}
	
	
####################################################################################
###### test two numeric data vectors representing independent groups ###############
####################################################################################
# bootstrap.2independent()
# Depends on:	bootstrapUtil.checkNumericAndNAs()
#				bootstrapUtil.computeStat2independent()
#				bootstrapUtil.printRunNumber()
#				bootstrapUtil.compute2sidePval()
#				bootstrapUtil.match2histAxisLims()
#				bootstrapPlot.groupHist()
#				bootstrapPlot.2groupBoxplot()
#
## Arguments
### group1, group2:			numeric vectors representing independent datasets
### Func:					char string of function name that statistic will be based on
### trials:					number of resampling runs
### replace:				logical indicating whether to resample with/without replacement
### groupNames:				2 element character vector holding names for group1 and group2
### plots:					logical indicating whether plot(s) are desired
### plotNullDist:			logical indicating whether to plot bootstrapped null distribution
### savePlot:				logical indicating whether to save plot as jpeg
### dataDescriptor:			character string describing the data, used for axis labels
### col, border, col.line:	character strings setting plot parameters
### jitter, pch:			numbers setting plot parameters
### verbose: 				logical indicating whether info printed to screen during computation

bootstrap.2independent = function(group1, group2, 
						  		 Func = 'mean', trials = 10000, replace = T, 
						  		 groupNames = c('group1', 'group2'),
						  		 plots = T, plotNullDist = T, savePlot = F, dataDescriptor = NULL, 
						  		 col = 'grey', border = 'darkgrey', col.line = 'red',
						  		 jitter = .15, pch = 21,
						  		 verbose = T, ...
						  		 )
{
	# get all input values for later
	allParam = as.list(sys.frame(sys.nframe()));#print(allParam)
	
	# check data
	tmp = bootstrapUtil.checkNumericAndNAs(group1, group2, paired = F);
	group1 = tmp$group1;
	group2 = tmp$group2;
	boxNULL = tmp$boxNULL;
	
	# compute actual stat
	tmp = bootstrapUtil.computeStat2independent(group1=group1, group2=group2, 
												Func=Func, replace=replace, 
												verbose=verbose, dataDescriptor=dataDescriptor);
	stat = tmp$stat; 
	Func = tmp$Func;
	
	# bootstrap
	statsNULL = c();
	for (trial in 1:trials) {
		if (verbose) {bootstrapUtil.printRunNumber(trial = trial, trials = trials)}
		pseudo_group1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_group2 = sample(boxNULL, length(group2), replace = replace);
		pseudo_stat_null = bootstrapUtil.computeStat2independent(group1 = pseudo_group1, 
												   				 group2 = pseudo_group2, 
												   			 	 Func = Func, 
												   				 verbose = F, 
												   				 dataDescriptor = dataDescriptor, 
												   				 replace = replace
												   				 )$stat;
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	if (any(is.na(statsNULL))) {stop('NA IN NULL DISTRIBUTION, CHECK YOUR DATA')}
	
	# compute p-values
	ptemp = bootstrapUtil.compute2sidePval(statsNULL = statsNULL, stat = stat, trials = trials);
	
	if (plots) {
		lims = bootstrapUtil.match2histAxisLims(group1, group2);
		if (plotNullDist) {
			if (savePlot) {
				tmp = paste(deparse(match.call()$group1), '-', deparse(match.call()$group2), '_', Func, sep = '');
				jpeg(file = paste(tmp, '.jpg', sep = ''), 
				 	 width = 8, height = 8, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			par(mfrow = c(2, 2));		
			
			# plot histograms of data for each group
			bootstrapPlot.groupHist(group=group1, groupName=groupNames[1], 
									Func=Func, dataDescriptor=dataDescriptor, 
									lims=lims, col=col, border=border, col.line=col.line);		
			bootstrapPlot.groupHist(group=group2, groupName=groupNames[2], 
									Func=Func, dataDescriptor=dataDescriptor, 
									lims=lims, col=col, border=border, col.line=col.line);	
									
			# plot histogram of null distribution with lines at stat and reflected stat
			bootstrapPlot.nullDist(stat=stat, abs.diffs=F, reflect=ptemp$reflect, 
								   statsNULL=statsNULL, trials=trials, Func=Func, replace=replace, 
								   col=col, border=border, col.line=col.line);
								   
			# boxplot of groups with individual data points
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2groupBoxplot))];#print(outerArgs)
			do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, p = ptemp$p, paired = F, ...));
		}
		else {
			if (savePlot) {
				tmp = paste(deparse(match.call()$group1), '-', deparse(match.call()$group2), '_', Func, sep = '');
				jpeg(file = paste(tmp, '.jpg', sep = ''), 
				 	 width = 12, height = 4, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			par(mfrow = c(1, 3));
			
			# plot histograms of data for each group
			bootstrapPlot.groupHist(group=group1, groupName=groupNames[1], 
									Func=Func, dataDescriptor=dataDescriptor, 
									lims=lims, col=col, border=border, col.line=col.line);
			bootstrapPlot.groupHist(group=group2, groupName=groupNames[2], 
									Func=Func, dataDescriptor=dataDescriptor, 
									lims=lims, col=col, border=border, col.line=col.line);

			# boxplot of groups with individual data points
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2groupBoxplot))];
			do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, p = ptemp$p, paired = F, ...));
		}
		if (savePlot) {dev.off()}
	}
	
	if (verbose) {cat('\np = ', ptemp$p, '\n', sep = '')}
	
	# build output list
	output = list(stat = stat, 
				stat.reflect = ptemp$reflect,
				p = ptemp$p, 
				p.left = ptemp$p_left,
				p.right = ptemp$p_right,
				null.dist = statsNULL,
				null.dist.mean = ptemp$midNULL, 
				data = setNames(list(group1, group2), groupNames),
				param = list(call.param = match.call(), all.param = allParam)
				);
	
	return(output);
}

####################################################################################
###### test two numeric data vectors representing dependent conditions #############
####################################################################################
# bootstrap.2paired()
# Depends on:	bootstrapUtil.checkNumericAndNAs()
#				bootstrapUtil.printRunNumber()
#				bootstrapUtil.compute2sidePval()
#
## Arguments
### condition1, condition2:
### trials:
### Func:
### plots:
### savePlot:
### abs.diffs:
### col, border, col.line:
### pch
### dataDescriptor:
### conditionNames:
### verbose:

bootstrap.2paired = function(condition1, condition2, 
						  		 trials = 10000, 
						  		 Func = 'mean',
						  		 plots = T, noDist = F, 
						  		 savePlot = F,
						  		 abs.diffs = F,
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, 
						  		 conditionNames = c('condition1', 'condition2'), 
						  		 pch = 21,
						  		 verbose = T, main=NULL,
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
		if (noDist) {
			data.lab = dataDescriptor;
			# boxplot of conditions with individual data points
			toPlot = data[, 1:2];
			
			if (ptemp$p == 0) {toPaste = paste('p < 1e-5 (n=', length(diffs), ')', sep = '')}
			else {toPaste = paste('p = ', signif(ptemp$p, 2), ' (n=', length(diffs), ')', sep = '')}
			toPaste=paste(main, ': ', toPaste, sep='')
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
		} else {
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

#######################################
####  ####
#######################################
# bootstrap.ANOVA()
# Depends on:	bootstrapUtil.calculateF()
#				bootstrapUtil.printRunNumber()
#				bootstrapUtil.2sidePval()
#				bootstrap.confidenceInterval()
#
## Arguments
### data:
### Func:
### absDiffs:
### trials:
### replace:
### plots:
### dataDescriptor:
### groupNames:
### col, border, col.line:
### pch
### verbose:

bootstrap.ANOVA = function(data = 'matrix or dataframe with groups in columns', 
					  	  Func = 'median',
					  	  absDiffs = T,
					  	  trials = 10000, 
					  	  replace = F,
					  	  Ftype = 'raw',
					  	  plots = T,
					  	  dataDescriptor = NULL,
					  	  groupNames = NULL,
					  	  pch = 21,
					  	  col = 'grey', border = 'darkgrey',
					  	  line.col = 'red', 
					  	  verbose = T, ...
					  	  )
{
	if (! (Ftype %in% c('raw','SS')) ) {
		choices = c('raw', 'SS');
		checkF = menu(choices=choices, title='Pick a method for computing "F":');
		Ftype = choices[checkF];#print(Ftype)
	}
	
	# calculate actual statistic and get group info
	if (Ftype=='SS') {
		calcF = bootstrapUtil.calculateF_SS(data = data);
		Fstat = calcF$Fstat; groupLengths = calcF$groupLengths; nGroups = calcF$nGroups;
	} else if (Ftype=='raw') {
		calcF = bootstrapUtil.calculateF(data = data, Func = Func, absDiffs = absDiffs);
		Fstat = calcF$FstatInfo$Fstat; groupLengths = calcF$groupInfo$groupLengths; nGroups = calcF$groupInfo$nGroups;
	} else {
		stop('SOMETHING WRONG WITH F CALCULATION SELECTION');
	}
	
	# combine groups
	#boxNULL = unlist(calcF$groupInfo$data);
	boxNULL = data[!(is.na(data))];#print(boxNULL)
	statsNULL = c();
	statNAs=c();
	
	# bootstrap
	if (verbose) {
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		if (Ftype=='raw') {
			cat(' \"F\" statistic: based on group ', Func, 's = ', Fstat, '\n', sep = '');
		} else if (Ftype=='SS') {
			cat(' \"F\" statistic: based on SS = ', Fstat, '\n', sep = '');
		}
		cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
	}
	
	for (trial in 1:trials) {
		if (verbose) {bootstrapUtil.printRunNumber(trials = trials, trial = trial)}
		pseudoNULL = matrix(nrow = max(groupLengths), ncol = nGroups);
		# add NAs to pseudo data before resampling if needed		
		for (gp in 1:ncol(pseudoNULL)) {
			if (groupLengths[gp] < nrow(pseudoNULL)) {
				numNA = nrow(pseudoNULL) - groupLengths[gp];
				pseudoNULL[1:numNA, gp] = NA;
				pseudoNULL[(numNA+1):nrow(pseudoNULL), gp] = sample(boxNULL, groupLengths[gp], replace = replace);
			}
			else {pseudoNULL[, gp] = sample(boxNULL, groupLengths[gp], replace = replace)}	
		}
		
		if (Ftype=='SS') {
			pseudoF = bootstrapUtil.calculateF_SS(data = pseudoNULL)$Fstat;
		} else if (Ftype=='raw') {
			pseudoF = bootstrapUtil.calculateF(data = pseudoNULL, Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat;
		} else {
			stop('SOMETHING WRONG WITH F CALCULATION SELECTION');
		}
			
		if (is.infinite(pseudoF) | is.nan(pseudoF) | is.na(pseudoF)) {statNAs = c(statNAs, trial); next}
		statsNULL = c(statsNULL, pseudoF);#print(mean(statsNULL))
	}
	if (length(statNAs) > 0) {warning('Check $statNAs for run #s where pseudoF was NA, consider setting absDiffs=T')}

	# compute p-values
	if (absDiffs) {p = sum(statsNULL > Fstat) / trials; p_left = NULL; p_right = NULL; ptemp = list(p = p)}
	else {ptemp = bootstrapUtil.compute2sidePval(statsNULL=statsNULL, stat=Fstat, trials=trials)}
	
	if (verbose) {cat('\nComputing 95% confidence intervals...\n')}
	citmp = bootstrap.multiConfidenceIntervals(data=data, Func=Func, groupNames=groupNames);
	ints = citmp$ints; 
	values = citmp$values;
	
	if (plots) {
		par(mfrow = c(1, 2));
		hist(statsNULL, 
			 xlim = c(min(statsNULL), max(statsNULL, Fstat)), xlab = '',
			 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
			 col = col, border = border, ...);
			 
		if (absDiffs) {abline(v = Fstat, col = line.col)}
		else {
			abline(v = Fstat, col = line.col);
			abline(v = ptemp$reflect, col = line.col, lty = 'dashed');
		}
		
		if (is.null(groupNames)) {groupNames = colnames(data)}
		else if (!is.null(groupNames) & length(groupNames)==ncol(data)) {groupNames = groupNames}
		boxplot(data, names = groupNames, 
			    main = paste('p = ', ptemp$p, sep = ''), 
			    border = 'lightgrey', medlty = 'blank', boxwex = .5, ...);
			
		stripchart(as.data.frame(data), vertical = T, add = T, pch = pch, bg = col, ...);
		
		coMat = matrix(nrow = ncol(data), ncol = 2);
		coMat[1, ] = c(.75, 1.25);
		for (gp in 2:nrow(coMat)) {coMat[gp, ] = coMat[gp - 1, ] + 1}
		for (gp in 1:length(values)) {
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
	
	output = list(stat = Fstat,
			      data = data,
			      null.dist = statsNULL,
			      p = ptemp$p, p.left = ptemp$p_left, p.right = ptemp$p_right,
			      conf.ints = as.data.frame(ints), values=values, calcF = calcF,
			      parameters = match.call(), statNAs = statNAs
			      );
	return(output);
}