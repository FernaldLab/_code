bootstrapPlot.2pairedGroupDiffHist <-
function(dataDescriptor = NULL, diffs, Func, conditionNames, col, border, col.line, returnDataLabel = T, ...)
{
	# check for dataDescriptor 
	if (!is.null(dataDescriptor) & is.character(dataDescriptor)) {data.lab = dataDescriptor}
	else {
		data.lab = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
	}
	
	# compute histogram of individual differences
	hist.diffs = hist(diffs, breaks = length(diffs), plot = F);
	xmax = max(abs(hist.diffs$breaks));

	# plot histogram
	plot(hist.diffs,
		 xlim = c(-xmax, xmax),
		 main = paste('Individual differences with ', Func, '\n(n=', length(diffs), ')', sep = ''),
		 xlab = paste(data.lab, 
		 			  ' (', conditionNames[1], ' - ', conditionNames[2], ')',
		 			   sep = ''
		 			   ),
		 col = col, border = border,
		 #cex.lab = cex.lab,
		 #cex.axis = cex.axis,
		 ...);
	abline(v = mean(diffs), col = col.line);
	if (returnDataLabel) {return(data.lab)}
}
