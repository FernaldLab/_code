bootstrapUtil.wilcoxStat <-
function(group1, group2)
{
	stat = suppressWarnings(as.double(wilcox.test(x = group1, y = group2, 
											      exact = T, 
											      alternative = 'two.sided'
											     )$statistic));
	return(stat);
}
