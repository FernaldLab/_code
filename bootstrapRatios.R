bootstrapRatios = function(group1, group2, ratios, trials)
{
	# count ratios > 1
	bigger = sum(ratios > 1);
	#smaller_equal = sum(ratios <= 1);
	
	# compute proportion of ratios > 1
	actualRatio = bigger / length(ratios);
	
	data = cbind(group1, group2);
	pseudo = c();
	for(trial in 1:trials)
	{
		# resample to randomize winner
		resampled = sample(c(1, 2), nrow(data), replace = T, prob = c(.2857143, 1-.2857143));
		p_ratios = c();
		for (r in 1:nrow(data))
		{
			if (resampled[r] == 1)
			{
				p_ratios = c(p_ratios, data[r, 1] / data[r, 2]);
			}
			else if (resampled[r] == 2)
			{
				p_ratios = c(p_ratios, data[r, 2] / data[r, 1]);
			}
		}
		p_bigger = sum(p_ratios > 1) / length(p_ratios);
		pseudo = c(pseudo, p_bigger);
	}
	p = sum(pseudo > actualRatio) / trials;
	hist(pseudo, col = 'grey', border = 'darkgrey', main = paste('p = ', p, sep = ''));
	abline(v = actualRatio, col = 'red');
	output = list(actual.ratio = actualRatio, data = data, null = pseudo, p = p);
	return(output);
}