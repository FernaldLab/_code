bootstrapUtil.calculateF <-
function(data, Func = 'mean', absDiffs = T)
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
