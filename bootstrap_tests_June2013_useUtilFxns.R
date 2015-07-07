#############
####  ####
#############
confidenceInterval = function(data, 
							  interval = 0.95, 
							  trials = 10000, 
							  Func = 'mean', 
							  plots = T, 
							  showDist = T,
							  col = 'grey', 
							  border = 'darkgrey', 
							  line.col = 'red',
							  pch = 21,
							  verbose = T
							  )
{
	# get input parameter values
	allParam = as.list(sys.frame(sys.nframe()));
	
	# check data is a numeric vector
	if (!is.numeric(data) | !is.vector(data)) {stop('CHECK THAT DATA IS A NUMERIC VECTOR')}
	
	# compute actual value
	value = eval(call(Func, data));
	
	# bootstrap
	pseudo = c();
	for (trial in 1:trials)
	{
		# print run number to console
		if (verbose) {printRunNumber(trials = trials, trial = trial)}
		
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
	if (plots)
	{
		if (showDist)
		{
			par(mfrow = c(1, 2));
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.confIntDist))];
			do.call(bootstrapPlot.confIntDist, c(outerArgs, 
												 pseudo = list(pseudo), 
												 value = value, conf.int = list(conf.int)
												 ));
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.groupStripchart))];
			do.call(bootstrapPlot.groupStripchart, c(outerArgs, value = value, conf.int = list(conf.int)));
		}
		else 
		{
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.groupStripchart))];
			do.call(bootstrapPlot.groupStripchart, c(outerArgs, value = value, conf.int = list(conf.int)));
		}

	}
	
	output = list(value = value, conf.int = conf.int, dist = pseudo, parameters = match.call());
	return(output);
}

################################################################################################
#### bootstrap2independent() tests two numeric data vectors representing independent groups ####
################################################################################################

bootstrap2independent = function(group1, group2, 
						  		 Func = 'mean', trials = 10000, replace = T, 
						  		 groupNames = c('group1', 'group2'),
						  		 plots = T, plotNullDist = T, dataDescriptor = NULL, 
						  		 col = 'grey', border = 'darkgrey', col.line = 'red',
						  		 jitter = .15, pch = 21,
						  		 verbose = T
						  		 )
{
	# get all input values for later
	allParam = as.list(sys.frame(sys.nframe()));
	
	# check data
	tmp = checkNumericAndNAs(group1, group2);
	group1 = tmp$group1;
	group2 = tmp$group2;
	boxNULL = tmp$boxNULL;
	checkForNegativeDataIfCV(group1, group2, Func);

	# compute statistic
	# stat = eval(call(Func, group1)) - eval(call(Func, group2));
	
	# # resample to create null distribution of statistic
	# if (verbose)
	# {
		# cat('...........................................\n');
		# cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		# cat('Building null distribution...\n');
		# cat(' Test statistic: difference of group ', Func, 's\n', sep = '');
		# cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
	# }
	
	args = allParam[names(allParam) %in% names(formals(computeStat2independent))];
	tmp = do.call(computeStat2independent, args);
	stat = tmp$stat; Func = tmp$Func;
	
	statsNULL = c();
	for (trial in 1:trials)
	{
		if (verbose) {printRunNumber(trial = trial, trials = trials)}
		pseudo_group1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_group2 = sample(boxNULL, length(group2), replace = replace);
		pseudo_stat_null = computeStat2independent(group1 = pseudo_group1, 
												   group2 = pseudo_group2, 
												   Func = Func, 
												   verbose = F, 
												   dataDescriptor = dataDescriptor, 
												   replace = replace
												   )$stat;
		#pseudo_stat_null = eval(call(Func, pseudo_group1)) - eval(call(Func, pseudo_group2));
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	if (any(is.na(statsNULL))) {stop('NA IN NULL DISTRIBUTION, CHECK YOUR DATA')}
	
	# compute p-values
	ptemp = compute2sidePval(statsNULL = statsNULL, stat = stat, trials = trials);
	
	if (plots)
	{
		lims = match2histAxisLims(group1, group2);
		if (plotNullDist)
		{
			par(mfrow = c(2, 2));
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.groupHist))];
			do.call(bootstrapPlot.groupHist, c(outerArgs[!names(outerArgs)=='Func'], 
											   group = list(group1), 
											   groupName = groupNames[1], 
											   lims = list(lims),
											   Func = Func
											   ));		
			do.call(bootstrapPlot.groupHist, c(outerArgs[!names(outerArgs)=='Func'], 
											   group = list(group2), 
											   groupName = groupNames[2], 
											   lims = list(lims),
											   Func = Func
											   ));				
			# plot histogram of null distribution with lines at stat and reflected stat
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.nullDist))];
			do.call(bootstrapPlot.nullDist, c(outerArgs[!names(outerArgs)=='Func'],
											  stat = stat, abs.diffs = F,
											  reflect = ptemp$reflect, 
							   				  statsNULL = list(statsNULL),
											   Func = Func
							   				  ));
			# boxplot of groups with individual data points
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2groupBoxplot))];
			do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, p = ptemp$p, paired = F));
		}
		else
		{
			par(mfrow = c(1, 3));
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.groupHist))];
			do.call(bootstrapPlot.groupHist, c(outerArgs[!names(outerArgs)=='Func'], 
											   group = list(group1), 
											   groupName = groupNames[1], 
											   Func = Func));		
			do.call(bootstrapPlot.groupHist, c(outerArgs[!names(outerArgs)=='Func'], 
											   group = list(group2), 
											   groupName = groupNames[2], 
											   Func = Func));	
			# boxplot of groups with individual data points
			outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2groupBoxplot))];
			do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, p = ptemp$p, paired = F));
		}
	}
	
	if (verbose) {cat('\nstatistic = ', stat, ', p = ', ptemp$p, '\n', sep = '')}
	
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

################################################
#### bootstrap2paired  ####
################################################

bootstrap2paired = function(condition1, condition2, 
						  		 trials = 10000, 
						  		 Func = 'mean',
						  		 plots = T, 
						  		 abs.diffs = F,
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, 
						  		 conditionNames = c('condition1', 'condition2'), 
						  		 pch = 21,
						  		 cex.lab = 1.2, cex.axis = 1.2,
						  		 printResults = T,
						  		 verbose = T,
						  		 ...
						  		 )
{
	allParam = as.list(sys.frame(sys.nframe()));
	# check that data is numeric
	numCheck1 = is.numeric(condition1);
	numCheck2 = is.numeric(condition2);
	if (!numCheck1 | !numCheck2)
	{
		stop('DATA CONTAINS NON-NUMERIC VALUES...\n   I ONLY EAT NUMBERS!!!\n    GIVE ME NUMBERS!!!!!');
	}
	
	# check that conditions have same number of data points
	if (length(condition1) != length(condition2))
	{
		stop('GROUPS ARE DIFFERENT SIZES...\n   THIS DATA SHOULD BE PAIRED!!!');
	}
	
	# check for missing data
	NAcheck1 = is.na(condition1);
	NAcheck2 = is.na(condition2);
	removeMe = NAcheck1 | NAcheck2;
	if (sum(removeMe) > 0)
	{
		condition1 = condition1[!removeMe];
		condition2 = condition2[!removeMe];
		warning('NAs in one/both conditions, corresponding data points removed from BOTH');
	}
	
	# compute statistic
	diffs = condition1 - condition2; 
	stat = eval(call(Func, diffs));
	if (abs.diffs) {stat = abs(stat)}
	
	# store data for output
	data = data.frame(condition1, condition2, diffs);
	names(data) = c(conditionNames[1], conditionNames[2], 'diffs');
	
	# resample to create null distribution of statistic
	if (verbose)
	{
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: ', Func, ' of individual differences\n\n', sep = '');
	}
	statsNULL = c();
	for (trial in 1:trials)
	{
		if (verbose) {printRunNumber(trials = trials, trial = trial)}
		signs = sample(c(-1, 1), length(diffs), replace = T);
		pseudo_diffs = signs * diffs;
		pseudo_stat_null = eval(call(Func, pseudo_diffs));
		if (abs.diffs) {pseudo_stat_null = abs(pseudo_stat_null)}
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	
	# compute p-values
	if (abs.diffs) {p = sum(statsNULL > stat) / trials}
	else if (!abs.diffs)
	{
		ptemp = compute2sidePval(statsNULL = statsNULL, stat = stat, trials = trials);
	}
	
	if (plots)
	{
		par(mfrow = c(1, 3));
		
		outerArgs = allParam[names(allParam) %in% names(formals(bootstrapPlot.2pairedGroupDiffHist))];
		data.lab = do.call(bootstrapPlot.2pairedGroupDiffHist, c(outerArgs, 
																 diffs = list(diffs), 
																 returnDataLabel = T
																 ));
		
		# plot histogram of null distribution with lines at stat and reflected stat
		if (abs.diffs)
		{
			xlabNULL = paste('statistic: abs(', Func, ' of individual differences)', sep = '');
		}
		else if (!abs.diffs)
		{
			xlabNULL = paste('statistic: ', Func, ' of individual differences', sep = '');
		}
		hist(statsNULL, 
		 	 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 	 xlab = xlabNULL,
		 	 col = col, border = border,
		 	 breaks = floor(trials / 100),
		 	 cex.lab = cex.lab,
			 cex.axis = cex.axis
		 	 );
		if (abs.diffs)
		{
			abline(v = stat, col = col.line);
		}
		else if (!abs.diffs)
		{
			abline(v = stat, col = col.line);
			abline(v = ptemp$reflect, col = col.line, lty = 'dashed');
		}
	
		# boxplot of conditions with individual data points
		toPlot = data[, 1:2];
		
		if (ptemp$p == 0)
		{
			toPaste = paste('p < 1e-5 (n=', length(diffs), ')', sep = '');
		}
		else
		{
			toPaste = paste('p = ', signif(ptemp$p, 2), ' (n=', length(diffs), ')', sep = '');
		}
		
		# draw boxes
		boxplot(toPlot,
				ylab = data.lab,
				#medlty = medlty,
				main = toPaste,
				frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,
				...
				);
				
		# add data points
		stripchart(toPlot,
				   vertical = T,
				   add = T,
				   pch = pch,
				   bg = col,
				   cex = 1.5
				   );
				   
		# connect paired points across conditions
		for (row in 1:nrow(data))
		{
			segments(1, data[row, 1], 2, data[row, 2], col = border);
		}
		
	}
	
	# build output list
	if (abs.diffs)
	{
		output = list(stat = stat,
					  p = ptemp$p,
					  null.dist = statsNULL,
					  data = data,
					  call.param = match.call(),
					  all.param = allParam
					  );
	}
	else if (!abs.diffs)
	{
		output = list(stat = stat, 
					  stat.reflect = ptemp$reflect,
					  p = ptemp$p, 
					  p.left = ptemp$p_left,
				  	  p.right = ptemp$p_right,
				  	  null.dist = statsNULL,
				  	  null.dist.mean = ptemp$midNULL,
				  	  data = data,
					  call.param = match.call(),
					  all.param = allParam
				  	  );
	}
				  
	# print summary of results to screen if specified			
	if (printResults)
	{
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist, $data, and $parameters\n');
	}

	return(output);
}

################################################################################################
####  ####
################################################################################################

calculateF = function(data = 'matrix or dataframe with groups in columns', 
					  Func = 'mean',
					  absDiffs = T
					  )
{
	# check data
	if (!is.data.frame(data) & !is.matrix(data))
	{
		stop('Data is not stored in matrix or dataframe...\n   FIX IT!');
	}
	if (mode(data) != 'numeric')
	{
		stop('DATA ISN\'T NUMERIC...\n   I ONLY EAT NUMBERS!!!');
	}
	if (nrow(data) < ncol(data))
	{
		warning('More groups than data points per group...\n   ARE YOU SURE?');
	}
	data = as.data.frame(data);
	
	# group info
	nGroups = ncol(data);
	groupNums = apply(data, 2, get(Func), na.rm = T);
	groupNAs = apply(apply(data, 2, is.na), 2, sum);
	groupLengths = nrow(data) - groupNAs;	
	# store for output
	groupInfo = list(data = data, 
					 nGroups = nGroups, 
					 groupNums = groupNums, 
					 groupNAs = groupNAs, 
					 groupLengths = groupLengths
					 );		
	
	# grand mean/median/Func
	grandTop = sum(groupLengths * groupNums);
	grandBot = sum(groupLengths);
	grandNum = grandTop / grandBot;
	# store for output
	grandNumInfo = list(grandTop = grandTop, 
						grandBot = grandBot, 
						grandNum = grandNum
						);
	
	# test statistic numerator
	diffs = grandNum - groupNums;
	if (absDiffs)
	{
		diffs = abs(diffs);
	}
	numerator = sum(groupLengths * diffs);
	# store for output
	Fnum = list(diffs = diffs,
				numerator = numerator
				);
	
	# test statistic denominator
	toSum = list();
	sums = c();
	for (gp in 1:nGroups)
	{
		temp = data[, gp] - groupNums[gp];
		temp = temp[!is.na(temp)];
		if (absDiffs)
		{
			temp = abs(temp);
		}
		toSum[[gp]] = temp;
		sums[gp] = sum(toSum[[gp]]);
	}
	denominator = sum(sums);
	# store for output
	Fden = list(toSum = toSum,
				sums = sums,
				denominator = denominator
				);
	
	Fstat = numerator / denominator; 
	# store for output
	FstatInfo = list(Fnum = Fnum,
					 Fden = Fden,
					 Fstat = Fstat
					 );
	
	out = list(groupInfo = groupInfo, 
			   grandNumInfo = grandNumInfo,
			   FstatInfo = FstatInfo,
			   parameters = match.call()
			   );
	
	return(out);
}

#######################################
####  ####
#######################################

bootstrapANOVA = function(data = 'matrix or dataframe with groups in columns', 
					  	  Func = 'mean',
					  	  absDiffs = T,
					  	  trials = 10000, 
					  	  replace = F,
					  	  plots = T,
					  	  groupNames = NULL,
					  	  pch = 21,
					  	  col = 'grey', border = 'darkgrey',
					  	  line.col = 'red', 
					  	  verbose = T
					  	  )
{
	# calculate actual statistic and get group info
	calcF = calculateF(data = data, Func = Func, absDiffs = absDiffs);
	Fstat = calcF$FstatInfo$Fstat;
	groupLengths = calcF$groupInfo$groupLengths;
	nGroups = calcF$groupInfo$nGroups;
	
	# combine groups
	boxNULL = unlist(calcF$groupInfo$data);#print(boxNULL)
	statsNULL = c();
	
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trial %% 1000 == 0)
			{
				cat('  Run ', trial, '\n', sep = '');
			}
			
		}
		pseudoNULL = matrix(nrow = max(groupLengths), ncol = nGroups);		
		for (gp in 1:ncol(pseudoNULL))
		{
			pseudoNULL[, gp] = sample(boxNULL, groupLengths[gp], replace = replace);
		}
		pseudoF = calculateF(data = pseudoNULL, Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat;
		statsNULL = c(statsNULL, pseudoF);
	}
	
	# compute p-values
	if (absDiffs)
	{
		p = sum(statsNULL > Fstat) / trials;
		p_left = NULL;
		p_right = NULL;
	}
	else
	{
		midNULL = mean(statsNULL);
		reflect = midNULL - (Fstat - midNULL);
			if (Fstat < 0)
			{
				p_left = sum(statsNULL < Fstat) / trials;
				p_right = sum(statsNULL > reflect) / trials;
				p = p_left + p_right;
			}
			else if (Fstat > 0)
			{
				p_right = sum(statsNULL > Fstat) / trials;
				p_left = sum(statsNULL < reflect) / trials;
				p = p_right + p_left;
			}
			else
			{
				stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
			}
	}
	

	
	if (plots)
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
		print(Fstat);print(summary(statsNULL))
		par(mfrow = c(1, 2));
		hist(statsNULL, 
			 xlim = c(min(statsNULL), max(statsNULL, Fstat)), 
			 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
			 col = col, border = border, 
			 xlab = ''
			 );
		if (absDiffs)
		{
			abline(v = Fstat, col = line.col);
		}
		else 
		{
			abline(v = Fstat, col = line.col);
			abline(v = reflect, col = line.col, lty = 'dashed');
		}
		
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
			    main = paste('p = ', signif(p,2), sep = ''), 
			    border = 'lightgrey', 
			    medlty = 'blank', 
			    boxwex = .5
			    );
			
		stripchart(as.data.frame(data), vertical = T, add = T, pch = pch, bg = col);
		
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
	
	output = list(stat = Fstat,
			      data = data,
			      null.dist = statsNULL,
			      p = p,
			      p.left = p_left,
			      p.right = p_right,
			      calcF = calcF,
			      conf.ints = ints,
			      parameters = match.call()
			      );
	
	return(output);
}
#######################################
####  ####
#######################################
bootstrapStats = function(file)
{
	if (!is.character(file)) {stop('WRONG KIND OF INPUT\n...MUST BE FILE NAME IN QUOTES')}

	this_ext = substring(file, first = (nchar(file)-3));
	if (this_ext != '.csv') {stop('WRONG FILE TYPE DETECTED: ', file, '\n...DATA MUST BE IN .csv FILE!')}
	
	data = data.frame(read.csv(file));
	cat('..............................................\n');
	if (!is.null(nrow(data))) {cat('Input data has ', nrow(data), ' rows and ', ncol(data), ' columns\n', sep = '')}
	
	if (ncol(data)==2 & nrow(data)>2)
	{
		answer = checkYesOrNo(readline(prompt = ' Do the 2 columns represent different groups or conditions? y/n'));
		if (answer == 'y')
		{
			answer = checkYesOrNo(readline(prompt = '   Is this paired data? y/n'));
			if (answer == 'y')
			{
				
			}
			else if (answer == 'n')
			{
				checkCor = cor.test(data[, 1], data[, 2]);
				if (checkCor$p.value < 0.05)
				{
					cat('    Correlation between columns is ', signif(checkCor$statistic, 2), ', p = ', checkCor$p.value, '\n', sep = '');
					answer = checkYesOrNo(readline(prompt = '    Are you sure the data isn\'t paired? y/n'));
				}
			}
		}
	}
	return(T);
}	