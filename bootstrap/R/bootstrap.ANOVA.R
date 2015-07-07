bootstrap.ANOVA <-
function(data = 'matrix or dataframe with groups in columns', 
					  	  Func = 'mean',
					  	  absDiffs = T,
					  	  trials = 10000, 
					  	  replace = F,
					  	  plots = T,
					  	  dataDescriptor = NULL,
					  	  groupNames = NULL,
					  	  pch = 21,
					  	  col = 'grey', border = 'darkgrey',
					  	  line.col = 'red', 
					  	  verbose = T, ...
					  	  )
{
	# calculate actual statistic and get group info
	calcF = bootstrapUtil.calculateF(data = data, Func = Func, absDiffs = absDiffs);
	Fstat = calcF$FstatInfo$Fstat;
	groupLengths = calcF$groupInfo$groupLengths;#print(groupLengths)
	nGroups = calcF$groupInfo$nGroups;
	
	# combine groups
	boxNULL = unlist(calcF$groupInfo$data);
	boxNULL = boxNULL[!(is.na(boxNULL))];#print(boxNULL)
	statsNULL = c();
	statNAs=c();
	
	# bootstrap
	if (verbose) {
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' \"F\" statistic: based on group ', Func, 's = ', Fstat, '\n', sep = '');
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
		pseudoF = bootstrapUtil.calculateF(data = pseudoNULL, Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat;
		
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
			      conf.ints = as.data.frame(ints), calcF = calcF,
			      parameters = match.call(), statNAs = statNAs
			      );
	return(output);
}
