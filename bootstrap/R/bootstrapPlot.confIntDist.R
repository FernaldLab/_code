bootstrapPlot.confIntDist <-
function(pseudo, Func, trials, value, conf.int, plotParams, ...)
{
	col = plotParams$col;
	border = plotParams$border;
	line.col = plotParams$line.col;
	hist(pseudo, breaks = length(pseudo)/100, 
	     col = col, border = border, 
	     main = paste('bootstrapped ', Func, '\n(n=', trials, ')', sep = ''), 
	     xlab = '', ...);
	abline(v = value, col = line.col);
	abline(v = conf.int[1], col = line.col, lty = 'dashed');
	abline(v = conf.int[2], col = line.col, lty = 'dashed');
}
