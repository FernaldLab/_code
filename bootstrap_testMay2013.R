################################################################################################
#### bootstrap2independent() tests two numeric data vectors representing independent groups ####
################################################################################################

bootstrap2independent = function(group1, group2, 
						  		 trials = 10000, 
						  		 Func = 'median', 
						  		 replace = F, 
						  		 plots = T, 
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, 
						  		 groupNames = c('group1', 'group2'),
						  		 jitter = .15, 
						  		 pch = 21,
						  		 boxLineMedian = T,
						  		 printResults = T,
						  		 verbose = T
						  		 )
{
	# combine datasets and check that data is numeric
	boxNULL = c(group1, group2);
	not_numeric = !is.numeric(boxNULL);
	NA_check = sum(is.na(boxNULL)) > 0;
	if (not_numeric)
	{
		stop('DATA CONTAINS NON-NUMERIC VALUES...\n   I ONLY EAT NUMBERS!!!\n    GIVE ME NUMBERS!!!!!');
	}
	if (NA_check)
	{
		group1 = group1[!is.na(group1)];
		group2 = group2[!is.na(group2)];
		boxNULL = c(group1, group2);
		warning('NAs in one or both groups, check your data');
	}
	data = list(group1, group2);
	names(data) = c(groupNames[1], groupNames[2]);
	
	# compute statistic
	stat = eval(call(Func, group1)) - eval(call(Func, group2));
	
	# resample to create null distribution of statistic
	if (verbose)
	{
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: difference of group ', Func, 's\n', sep = '');
		cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
	}
	statsNULL = c();
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trials <= 20000)
			{
				if (trial %% 1000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
			else if (trials > 20000)
			{
				if (trial %% 10000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
		}
		pseudo_group1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_group2 = sample(boxNULL, length(group2), replace = replace);
		pseudo_stat_null = eval(call(Func, pseudo_group1)) - eval(call(Func, pseudo_group2));
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	
	# compute p-values
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	if (stat < 0)
	{ 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} 
	else if (stat > 0)
	{
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	}
	else
	{
		stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
	}
	
	# histograms of datasets and null distribution, and boxplot with p-value
	if (plots)
	{
		# add n to group names
		groupNames[1] = paste(groupNames[1], ' (n=', length(group1), ')', sep = '');
		groupNames[2] = paste(groupNames[2], ' (n=', length(group2), ')', sep = '');
		
		par(mfrow = c(2, 2));
		
		# compute histograms for each group
		hist_group1 = hist(group1, breaks = length(group1), plot = F);
		hist_group2 = hist(group2, breaks = length(group2), plot = F);
		
		# compute axis limits
		ymax = max(c(hist_group1$counts, hist_group2$counts));
		xmin = min(c(hist_group1$breaks, hist_group2$breaks));
		xmax = max(c(hist_group1$breaks, hist_group2$breaks));
		
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
		
		# plot group histograms with lines at means/medians
		plot(hist_group1, 
			 main = groupNames[1], 
			 ylim = c(0, ymax), xlim = c(xmin, xmax), 
			 col = col, border = border, 
			 xlab = data.lab
			 );
		abline(v = eval(call(Func, group1)), col = col.line);
		plot(hist_group2, 
			 main = groupNames[2], 
			 ylim = c(0, ymax), xlim = c(xmin, xmax), 
			 col = col, border = border, 
			 xlab = data.lab
			 );
		abline(v = eval(call(Func, group2)), col = col.line);
		
		# plot histogram of null distribution with lines at stat and reflected stat
		hist(statsNULL, 
		 	 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 	 xlab = paste('statistic (', Func, 's, replace = ', substr(replace, 1, 1), ')', sep = ''),
		 	 col = col, border = border,
		 	 breaks = floor(trials / 100)
		 	 );
		abline(v = stat, col = col.line);
		abline(v = reflect, col = col.line, lty = 'dashed');
		
		# boxplot of groups with individual data points
		toPlot = boxNULL;
		grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], length(group2)));
		
		# check if line in box should be mean or median (default)
		if (boxLineMedian)
		{
			medlty = 'solid';
		}
		else
		{
			medlty = 'blank';
		}
		
		# draw boxes
		boxplot(toPlot ~ grp,
				ylab = data.lab,
				medlty = medlty,
				main = paste('p = ', signif(p, 2), sep = '')
				);
		
		# draw mean lines if specified
		if (!boxLineMedian)
		{
			segments(0.6, mean(group1), 1.4, mean(group1), lwd = 2);
			segments(1.6, mean(group2), 2.4, mean(group2), lwd = 2);
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
	}
	
	# build output list
	output = list(stat = stat, 
				stat.reflect = reflect,
				p = p, 
				p.left = p_left,
				p.right = p_right,
				null.dist = statsNULL,
				null.dist.mean = midNULL, 
				replacement = replace,
				data = data,
				parameters = match.call()
				);
	
	# print summary of results to screen if specified			
	if (printResults)
	{
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist and $parameters\n');
	}
	
	return(output);
}