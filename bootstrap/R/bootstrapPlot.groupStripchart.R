bootstrapPlot.groupStripchart <-
function(data, interval, value, conf.int, Func, plotParams, ...)
{
	line.col = plotParams$line.col;
	pch = plotParams$pch;
	col = plotParams$col;
	stripchart(data, vertical = T, 
			   pch = pch, bg = col, 
			   frame.plot = F, 
			   main = paste(interval*100, '% conf.int around ', round(value, 1), 
			   				'\n[', round(conf.int[1], 1), ', ', 
			   				round(conf.int[2], 1), ']', 
			   				sep = ''
			   				), ...
			   );
	segments(.8, value, 1.2, value, col = line.col);
	text(1.26, value, Func, ...);
	segments(.9, conf.int[1], 1.1, conf.int[1], col = line.col, lty = 'dashed');
	segments(.9, conf.int[2], 1.1, conf.int[2], col = line.col, lty = 'dashed');
}
