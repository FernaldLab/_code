bootstrap.2independent <-
function(group1, group2, 
						  		 Func = 'mean', trials = 10000, replace = T, 
						  		 groupNames = c('group1', 'group2'),
						  		 plots = T, plotNullDist = T, savePlot = F, dataDescriptor = NULL, 
						  		 col = 'grey', border = 'darkgrey', col.line = 'red',
						  		 jitter = .15, pch = 21,
						  		 verbose = T, ...
						  		 )
{
	# get all input values for later
	allParam = as.list(sys.frame(sys.nframe()));
	
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
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2groupBoxplot))];
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
