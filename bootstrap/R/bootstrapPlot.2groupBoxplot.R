bootstrapPlot.2groupBoxplot <-
function(group1, group2, groupNames, paired, p, dataDescriptor=NULL, col, border, pch, jitter, Func, ...)
{
	# add n to group names
	groupNames[1] = paste(groupNames[1], ' (n=', ( length(group1)-sum(is.na(group1)) ), ')', sep = '');
	groupNames[2] = paste(groupNames[2], ' (n=', ( length(group2)-sum(is.na(group2)) ), ')', sep = '');
	
	toPlot = c(group1, group2);
	grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], length(group2)));
	boxLineType = Func;
	if (boxLineType == 'mean') {medlty = 'blank'}
	else {medlty = 'solid'}
	toPaste = bootstrapUtil.getPvalToPaste(p);
	data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor);
	
	# draw boxes
	boxplot(toPlot ~ grp,
			ylab = data.lab,
			medlty = medlty,
			main = toPaste,...);
	# draw mean lines if specified
	if (boxLineType == 'mean') {
		segments(0.6, mean(group1,na.rm=T), 1.4, mean(group1,na.rm=T), lwd = 3, col = 'black');
		segments(1.6, mean(group2,na.rm=T), 2.4, mean(group2,na.rm=T), lwd = 3);
	}
	
	# add data points
	stripchart(toPlot ~ grp,
			   vertical = T,
			   add = T,
			   method = 'jitter',
			   jitter = jitter,
			   pch = pch,
			   bg = col,
			   #cex = 1.5, 
			   ...);
			   
	if (paired) {
		data = data.frame(group1, group2);
		# connect paired points across conditions
		for (row in 1:nrow(data)) {segments(1, data[row, 1], 2, data[row, 2], col = border)}
	}
	
	# if (confidenceIntervals)
	# {
		# if (verbose)
		# {
			# cat('  Computing 95% confidence intervals\n');
		# }
		# values = c();
		# ints = matrix(nrow = ncol(data), ncol = 2);
		# for (gp in 1:ncol(data))
		# {
			# temp = confidenceInterval(data[, gp], plots = F, Func = Func);###NEED TO DEFINE data
			# values = c(values, temp$value);
			# ints[gp, ] = temp$conf.int;
		# }
		# coMat = matrix(nrow = ncol(data), ncol = 2);
		# coMat[1, ] = c(.75, 1.25);
		# for (gp in 2:nrow(coMat))
		# {
			# coMat[gp, ] = coMat[gp - 1, ] + 1;
		# }
		# for (gp in 1:length(values))
		# {
			# segments(coMat[gp, 1], values[gp], 
				     # coMat[gp, 2], values[gp], 
				     # lwd = 1, col = line.col
				     # );
			# segments(coMat[gp, 1], ints[gp, 1], 
					 # coMat[gp, 2], ints[gp, 1], 
					 # lwd = 1, lty = 'dashed', col = line.col
					 # );
			# segments(coMat[gp, 1], ints[gp, 2], 
				     # coMat[gp, 2], ints[gp, 2], 
				     # lwd = 1, lty = 'dashed', col = line.col
				     # );
		# }
	# }
}
