bootstrapUtil.checkNumericAndNAs <-
function (group1, group2, paired = F)
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
