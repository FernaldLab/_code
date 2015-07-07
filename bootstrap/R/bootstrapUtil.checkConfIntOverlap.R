bootstrapUtil.checkConfIntOverlap <-
function(ints)
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
