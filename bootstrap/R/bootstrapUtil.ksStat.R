bootstrapUtil.ksStat <-
function(group1, group2)
{
	stat = suppressWarnings(as.double(ks.test(x = group1, y = group2, alternative= 'two.sided'
											 )$statistic));
	return(stat);
}
