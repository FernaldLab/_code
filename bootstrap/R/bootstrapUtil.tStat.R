bootstrapUtil.tStat <-
function(group1, group2)
{
	stat = suppressWarnings(as.double(t.test(x = group1, y = group2,
											 alternative = 'two.sided',
											 paired = F,
											 var.equal = F
											)$statistic));
	return(stat);
}
