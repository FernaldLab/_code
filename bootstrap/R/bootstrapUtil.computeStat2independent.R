bootstrapUtil.computeStat2independent <-
function(group1, group2, Func = c('mean.diff', 'median.diff', 'cv.diff', 't.test', 'wilcox.test', 'ks.test'), verbose, dataDescriptor, replace)
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
