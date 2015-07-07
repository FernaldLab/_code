   ###############################################################
###### print run number to console every 1000/10000 runs ######
###############################################################
##
###

bootstrapUtil.printRunNumber = function (trials, trial)
{
	if (trials <= 20000) {if (trial %% 1000 == 0) {cat('  Run ', trial, '\n', sep = '')}}
	else if (trials > 20000) {if (trial %% 10000 == 0) {cat('  Run ', trial, '\n', sep = '')}}
}

###############################################################
###### check that groups are numeric and remove any NAs #######
###############################################################
##
###

bootstrapUtil.checkNumericAndNAs = function (group1, group2, paired = F)
{
	not_numeric = !is.numeric(c(group1, group2));
	if (not_numeric) {stop('DATA CONTAINS NON-NUMERIC VALUES... FIX IT!\n')}
	
	if (paired) {
		nCheck = length(group1) != length(group2);
		if (nCheck) {stop('GROUPS ARE DIFFERENT SIZES...\n   DATA SHOULD BE PAIRED!!!')}
		
		# check for missing data
		NAcheck1 = is.na(group1);
		NAcheck2 = is.na(group2);
		removeMe = NAcheck1 | NAcheck2;
		if (sum(removeMe) > 0) {
			group1 = group1[!removeMe];
			group2 = group2[!removeMe];
			warning('NAs in one/both conditions, corresponding data points removed from BOTH');
		}
		output = list(condition1 = group1, condition2 = group2);
	}
	else {
		NA_check = sum(is.na(c(group1, group2))) > 0;
		if (NA_check) {
			group1 = group1[!is.na(group1)];
			group2 = group2[!is.na(group2)];
			boxNULL = c(group1, group2);
			warning('NAs removed from one or both groups, check your data');
		}
		else {boxNULL = c(group1, group2)}
		output = list(group1 = group1, group2 = group2, boxNULL = boxNULL);
	}
	return(output);
}

##########################################################################
###### check that data is non-negative if CV is used as descriptor #######
##########################################################################
##
###

bootstrapUtil.checkForNegativeDataIfCV = function(group1, group2, Func)
{
	if (Func == 'cv') {
		if (any(c(group1, group2) < 0)) {
			stop('CV IS NONSENSICAL FOR DATA WITH NEGATIVE VALUES\n ...USE MEAN OR MEDIAN');
		}
	}
}

##########################################################################
###### compute value of actual statistic for 2 independent groups ########
##########################################################################
##
###

bootstrapUtil.computeStat2independent = function(group1, group2, Func = c('mean.diff', 'median.diff', 'cv.diff', 't.test', 'wilcox.test', 'ks.test'), verbose, dataDescriptor, replace)
{
	choiceVec = c('mean.diff', 'median.diff', 'cv.diff _DO_NOT_USE_', 't.test', 'wilcox.test', 'ks.test');
	title = 'Please select one of the following statistics to use:';
	err = 'SOMETHING IS WRONG IN bootstrapUtil.computeStat2independent(), GET AUSTIN';
	
	# if Func not explicitly defined, ask user to select one
	if (length(Func) > 1) {
		checkArg = menu(choices = choiceVec, title = title);
		Func = choiceVec[checkArg];
		cat('Using ', Func, ' statistic...\n', sep = '');
	}
	# if Func defined, make sure it's a valid choice 
	else if (length(Func) == 1) {
		checkArg = pmatch(Func, choiceVec);
		if (is.na(checkArg)) {checkArg = menu(choices = choiceVec, title = title)}
		Func = choiceVec[checkArg];
	}
	else {stop(err)}
	
	# if Func is mean, median, or cv
	if (checkArg %in% 1:3) {
		Func = gsub('.diff', '', Func);
		bootstrapUtil.checkForNegativeDataIfCV(group1, group2, Func);
		stat = eval(call(Func, group1)) - eval(call(Func, group2));
		if (verbose) {
			cat('...........................................\n');
			cat('Testing: ', dataDescriptor, '\n\n', sep = '');
			cat('Building null distribution...\n');
			cat(' Test statistic: difference of group ', Func, 's = ', stat, '\n', sep = '');
			cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
		}
	}
	# if Func is t.test, wilcox.test, or ks.test
	else if (checkArg %in% 4:6) {
		Func = gsub('.test', '', Func);
		if (checkArg == 4) {stat = bootstrapUtil.tStat(group1, group2)}
		else if (checkArg == 5) {stat = bootstrapUtil.wilcoxStat(group1, group2)}
		else if (checkArg == 6) {stat = bootstrapUtil.ksStat(group1, group2)}
		if (verbose) {
			cat('...........................................\n');
			cat('Testing: ', dataDescriptor, '\n\n', sep = '');
			cat('Building null distribution...\n');
			cat(' Test statistic: ', Func, ' = ', stat, '\n', sep = '');
			cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
		}
	}
	else {stop(err)}
	
	return(list(stat = stat, Func = Func));
}

##########################################################################
######  ########
##########################################################################
##
###

bootstrapUtil.tStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(t.test(x = group1, y = group2,
											 alternative = 'two.sided',
											 paired = F,
											 var.equal = F
											)$statistic));
	return(stat);
}
##
bootstrapUtil.wilcoxStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(wilcox.test(x = group1, y = group2, 
											      exact = T, 
											      alternative = 'two.sided'
											     )$statistic));
	return(stat);
}
##
bootstrapUtil.ksStat = function(group1, group2)
{
	stat = suppressWarnings(as.double(ks.test(x = group1, y = group2, alternative= 'two.sided'
											 )$statistic));
	return(stat);
}
##
cv = function(data)
{
	value = sd(as.numeric(data), na.rm = T) / mean(as.numeric(data), na.rm = T);
	return(value);	
}

##########################################################################
######  ########
##########################################################################
##
###

bootstrapUtil.compute2sidePval = function(statsNULL, stat, trials)
{
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	#if (stat < 0)
	if (stat < midNULL) { 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} 
	#else if (stat > 0)
	else if (stat > midNULL) {
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	}
	else {
		stop('Either statistic==0 or midNULL or something else is wrong...\n');
	}
	return(list(midNULL = midNULL, reflect = reflect, p_left = p_left, p_right = p_right, p = p));
}

##########################################################################
######  ########
##########################################################################
##
###

bootstrapUtil.match2histAxisLims = function (group1, group2)
{
	hist_group1 = hist(group1, breaks = length(group1), plot = F);
	hist_group2 = hist(group2, breaks = length(group2), plot = F);
	
	# compute axis limits
	ymax = max(c(hist_group1$counts, hist_group2$counts));
	xmin = min(c(hist_group1$breaks, hist_group2$breaks));
	xmax = max(c(hist_group1$breaks, hist_group2$breaks));
	
	return(list(ymax = ymax, xmin = xmin, xmax = xmax));
}

bootstrapUtil.checkDataDescriptorToGetLabel = function (dataDescriptor)
{
	if (!is.null(dataDescriptor) & is.character(dataDescriptor)) {data.lab = dataDescriptor}
	else {
		data.lab = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
	}
	return(data.lab);
}

bootstrapUtil.getPvalToPaste = function(p)
{
	if (p == 0) {toPaste = paste('p < 1e-5', sep = '')}
	else {toPaste = paste('p = ', signif(p, 2), sep = '')}
	return(toPaste);
}
######################
#######
#######################
##
###
bootstrapUtil.calculateF = function(data, Func = 'mean', absDiffs = T)
{
	# check data
	if (!is.data.frame(data) & !is.matrix(data)) {stop('DATA MUST BE IN MATRIX OR DATAFRAME')}
	if (mode(data) != 'numeric') {stop('DATA MUST BE NUMERIC')}
	if (nrow(data) < ncol(data)) {warning('More groups than data points per group...\n   ARE YOU SURE?')}
	data = as.data.frame(data);
	
	# group info
	nGroups = ncol(data);
	groupNums = apply(data, 2, get(Func), na.rm = T);
	groupNAs = apply(apply(data, 2, is.na), 2, sum);
	groupLengths = nrow(data) - groupNAs;	
	# store for output
	groupInfo = list(data=data, nGroups=nGroups, groupNums=groupNums, groupNAs=groupNAs, groupLengths=groupLengths);		
	
	# grand mean/median/Func
	grandTop = sum(groupLengths * groupNums);
	grandBot = sum(groupLengths);
	grandNum = grandTop / grandBot;
	# store for output
	grandNumInfo = list(grandTop = grandTop, grandBot = grandBot, grandNum = grandNum);
	
	# test statistic numerator
	diffs = grandNum - groupNums;
	if (absDiffs) {diffs = abs(diffs)}
	numerator = sum(groupLengths * diffs);
	# store for output
	Fnum = list(diffs = diffs, numerator = numerator);
	
	# test statistic denominator
	toSum = list();
	sums = c();
	for (gp in 1:nGroups) {
		temp = data[, gp] - groupNums[gp];
		temp = temp[!is.na(temp)];
		if (absDiffs) {temp = abs(temp)}
		toSum[[gp]] = temp;
		sums[gp] = sum(toSum[[gp]]);
	}
	denominator = sum(sums);
	# store for output
	Fden = list(toSum = toSum, sums = sums, denominator = denominator);
	
	Fstat = numerator / denominator; 
	# store for output
	FstatInfo = list(Fnum = Fnum, Fden = Fden, Fstat = Fstat);
	
	out = list(groupInfo = groupInfo, 
			   grandNumInfo = grandNumInfo,
			   FstatInfo = FstatInfo,
			   parameters = match.call()
			   );
	
	return(out);
}
######################
#######
#######################
##
###

bootstrapUtil.calculateF_SS = function (data)
{	
	dat = data;
	
	#get group lengths
	gp_lengths = apply(dat,2,length) - apply(apply(dat,2,is.na), 2, sum);
	
	#compute mean of each group
	gp_means = apply(dat,2,mean,na.rm=T);
	
	#compute overall mean
	#don't just take mean of group means since group sizes may be unbalanced
	overall_mean = mean(dat, na.rm=T);
	
	#compute between-group sum of squares
	df_b = ncol(dat)-1;
	ss_b = sum( gp_lengths * ((gp_means - overall_mean)^2) ) / df_b;
	
	#compute within-group sum of squares
	df_w = (length(as.vector(dat)) - sum(is.na(as.vector(dat)))) - ncol(dat);
	ss_w = c();
	for (col in 1:ncol(dat)) {ss_w = c(ss_w, dat[, col] - gp_means[col])}; rm(col);
	ss_w = ss_w[!is.na(ss_w)]^2;
	ss_w = sum(ss_w) / df_w;
	
	#compute F-ratio
	f = ss_b / ss_w;
	
	return(list(Fstat=f, groupLengths=gp_lengths, nGroups=ncol(dat)));
}


#############
#####
#############
bootstrapUtil.checkConfIntOverlap = function(ints)
{
	mat = matrix(nrow=nrow(ints), ncol=nrow(ints), 
				 dimnames=list(rownames(ints), rownames(ints))
				 );
	for (gp in 1:nrow(ints)) {
		bigger = ints$lower[gp] > ints$upper;
		smaller = ints$upper[gp] < ints$lower;
		mat[gp, ] = bigger | smaller;
	}
	diag(mat) = NA;
	return(mat);
}

