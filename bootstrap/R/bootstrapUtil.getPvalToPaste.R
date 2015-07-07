bootstrapUtil.getPvalToPaste <-
function(p)
{
	if (p == 0) {toPaste = paste('p < 1e-5', sep = '')}
	else {toPaste = paste('p = ', signif(p, 2), sep = '')}
	return(toPaste);
}
