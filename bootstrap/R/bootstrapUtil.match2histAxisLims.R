bootstrapUtil.match2histAxisLims <-
function (group1, group2)
{
	hist_group1 = hist(group1, breaks = length(group1), plot = F);
	hist_group2 = hist(group2, breaks = length(group2), plot = F);
	
	# compute axis limits
	ymax = max(c(hist_group1$counts, hist_group2$counts));
	xmin = min(c(hist_group1$breaks, hist_group2$breaks));
	xmax = max(c(hist_group1$breaks, hist_group2$breaks));
	
	return(list(ymax = ymax, xmin = xmin, xmax = xmax));
}
