bootstrap.confidenceInterval <-
function(data, interval = 0.95, trials = 10000, Func = 'mean', 
									    plots = T, showDist = T, savePlot = F,
									    plotParams = list(col = 'grey', border = 'darkgrey', line.col = 'red', pch = 21),
									    verbose = T,
									    ...
									    )
{
	# check data is a numeric vector
	if (!is.numeric(data) | !is.vector(data)) {stop('CHECK THAT DATA IS A NUMERIC VECTOR')}
	
	# check for and remove NAs
	if (any(is.na(data))) {
		data0 = data;
		data = data[!(is.na(data))];
		warning(paste(sum(is.na(data0)), ' NAs detected in data', sep = ''));
	}
	
	# compute actual value
	value = eval(call(Func, data));
	if (verbose) {
		cat('....................................................................\n');
		cat('Computing ', interval*100, '% confidence interval around the ', Func, '...\n', sep = '');
		cat(' Resampling...\n');
	}
	
	# bootstrap
	pseudo = c();
	for (trial in 1:trials) {
		# print run number to console
		if (verbose) {bootstrapUtil.printRunNumber(trials = trials, trial = trial)}
		
		# resample data with replacement to create pseudo data
		pseudoData = sample(data, length(data), replace = T);
		
		# compute pseudo value
		pseudo = c(pseudo, eval(call(Func, pseudoData)));
	}
	
	# compute confidence interval of pseudo value distribution
	half = (1 - interval) / 2;
	int = c(half, interval + half);
	conf.int = quantile(pseudo, int);

	# plot stripchart and pseudo value distribution
	if (plots) {
		if (showDist) {
			if (savePlot) {
				jpeg(file = paste(deparse(match.call()$data), '.jpg', sep = ''), 
				 	 width = 8, height = 4, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			par(mfrow = c(1, 2));
			bootstrapPlot.confIntDist(pseudo = pseudo, Func = Func, trials = trials, 
									  value = value, conf.int = conf.int, plotParams = plotParams, ...);
			bootstrapPlot.groupStripchart(data = data, interval = interval, value = value, 
										  conf.int = conf.int, Func = Func, plotParams = plotParams, ...);
		}
		else {
			if (savePlot) {
				jpeg(file = paste(deparse(match.call()$data), '.jpg', sep = ''), 
				 	 width = 4, height = 4, units = 'in', quality = 100, type = 'quartz', res = 150);
			}
			bootstrapPlot.groupStripchart(data = data, interval = interval, value = value, 
										  conf.int = conf.int, Func = Func, plotParams = plotParams, ...);
		}
		if (savePlot) {dev.off()}
	}
	output = list(value = value, conf.int = conf.int, dist = pseudo, call.param = match.call());
	return(output);
}
