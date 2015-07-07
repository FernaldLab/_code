bootstrapUtil.compute2sidePval <-
function(statsNULL, stat, trials)
{
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	#if (stat < 0)
	if (stat < midNULL) { 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} 
	#else if (stat > 0)
	else if (stat > midNULL) {
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	}
	else {
		stop('Either statistic==0 or midNULL or something else is wrong...\n');
	}
	return(list(midNULL = midNULL, reflect = reflect, p_left = p_left, p_right = p_right, p = p));
}
