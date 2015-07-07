bootstrapUtil.checkForNegativeDataIfCV <-
function(group1, group2, Func)
{
	if (Func == 'cv') {
		if (any(c(group1, group2) < 0)) {
			stop('CV IS NONSENSICAL FOR DATA WITH NEGATIVE VALUES\n ...USE MEAN OR MEDIAN');
		}
	}
}
