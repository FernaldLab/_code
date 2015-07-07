bootstrapUtil.printRunNumber <-
function (trials, trial)
{
	if (trials <= 20000) {if (trial %% 1000 == 0) {cat('  Run ', trial, '\n', sep = '')}}
	else if (trials > 20000) {if (trial %% 10000 == 0) {cat('  Run ', trial, '\n', sep = '')}}
}
