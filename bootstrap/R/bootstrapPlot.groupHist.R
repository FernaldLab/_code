bootstrapPlot.groupHist <-
function(group, groupName, Func, dataDescriptor = NULL, lims = NULL, col, border, col.line, ...)
{
	# add n to group name
	groupName = paste(groupName, ' (n=', length(group), ')', sep = '');
	
	# compute histogram
	hist_group = hist(group, breaks = length(group), plot = F);
	
	if (is.null(lims)) {
		# compute axis limits
		ymax = max(hist_group$counts);
		xmin = min(hist_group$breaks);
		xmax = max(hist_group$breaks);
	}
	else if (!is.null(lims)) {
		if (is.list(lims) & is.numeric(unlist(lims))) {
			ymax = lims$ymax;
			xmin = lims$xmin;
			xmax = lims$xmax;
		}
	}
	
	# check for dataDescriptor 
	data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# plot group histogram with line at mean/median
	plot(hist_group, 
		 main = groupName, 
		 ylim = c(0, ymax), xlim = c(xmin, xmax), 
		 col = col, border = border, 
		 xlab = data.lab, ...);
	if (Func %in% c('mean', 't')) {v = mean(group)}
	else if (Func %in% c('median', 'cv', 'wilcox', 'ks')) {v = median(group)}
	#abline(v = eval(call(Func, group)), col = col.line);	
	abline(v = v, col = col.line);
}
