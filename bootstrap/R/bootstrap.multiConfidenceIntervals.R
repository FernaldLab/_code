bootstrap.multiConfidenceIntervals <-
function(data, Func, groupNames=NULL)
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
