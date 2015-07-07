bootstrapPlot.nullDist <-
function(stat, abs.diffs, reflect, statsNULL, trials, Func, replace, col, border, col.line, ...)
{
	hist_null = hist(statsNULL, breaks = floor(trials / 100), plot = F);
	
	xmax = max(max(hist_null$breaks), abs(min(hist_null$breaks)));
	if (all(statsNULL > 0)) {xlim = c(0, xmax)}
	else if (all(statsNULL < 0)) {xlim = c(-xmax, 0)}
	else {xlim = c(-xmax, xmax)}
	
	if (Func %in% c('mean', 'median', 'cv')) {
		xlab = paste('statistic (', Func, '.diff, replace = ', substr(replace, 1, 1), ')', sep = '');
	}
	else {
		xlab = paste('statistic (', Func, ', replace = ', substr(replace, 1, 1), ')', sep = '');
	}
	 
	plot(hist_null, xlim = xlim, xlab = xlab,
		 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
		 col = col, border = border, ...);

	if (abs.diffs) {abline(v = stat, col = col.line)}
	else {
		abline(v = stat, col = col.line);
		abline(v = reflect, col = col.line, lty = 'dashed');
	}
}
