plotCors = function(dataMat, test1, test2, abline=T, frame.plot=F, xlab=NULL, ylab=NULL, main='')
{
	if (is.null(xlab))
	{
		xlab = test1;
	}
	else
	{
		xlab = xlab;
	}
	if (is.null(ylab))
	{
		ylab = test2;
	}
	else
	{
		ylab = ylab;
	}
	
	test1col = match(test1, colnames(dataMat));
	test2col = match(test2, colnames(dataMat));
	
	# par(mfrow=c(1,3));
	
	# verboseScatterplot(datamat[, test1col], 
						# dataMat[, test2col], 
						# abline=abline, 
						# main='Pwinners',
						# xlab=xlab,
						# ylab=ylab,
						# frame.plot=F, 
						# type='n'
						# );
	# text(datP$Length.diff, 
		 # datP$'11.KT', 
		 # gsub('Dyad ', '', rownames(datP))
		 # );
	
	# verboseScatterplot(datNPmatch$Length.diff, 
						# datNPmatch$'11.KT', 
						# abline=T, 
						# main='NP match',
						# xlab=xlab,
						# ylab=ylab,
						# frame.plot=F, 
						# type='n'
						# );
	# text(datNPmatch$Length.diff, 
		 # datNPmatch$'11.KT', 
		 # gsub('Dyad ', '', rownames(datNPmatch))
		 # );
	
	# verboseScatterplot(datNPnomatch$Length.diff, 
						# datNPnomatch$'11.KT', 
						# abline=T, 
						# main='NP older',
						# xlab=xlab,
						# ylab=ylab,
						# frame.plot=F, 
						# type='n'
						# );
	# text(datNPnomatch$Length.diff, 
		 # datNPnomatch$'11.KT', 
		 # gsub('Dyad ', '', rownames(datNPnomatch))
		 # );
		 
}