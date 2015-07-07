
bootstrapANOVA = function(data = 'matrix or dataframe with groups in columns', 
					  	  Func = 'mean',
					  	  absDiffs = T,
					  	  trials = 10000, 
					  	  replace = F,
					  	  plots = T,
					  	  groupNames = NULL,
					  	  pch = 21,
					  	  col = 'grey', border = 'darkgrey',
					  	  line.col = 'red', 
					  	  verbose = T
					  	  )
{
	# calculate actual statistic and get group info
	calcF = calculateF(data = data, Func = Func, absDiffs = absDiffs);
	Fstat = calcF$FstatInfo$Fstat;
	groupLengths = calcF$groupInfo$groupLengths;print(groupLengths)
	nGroups = calcF$groupInfo$nGroups;
	
	# combine groups
	boxNULL = unlist(calcF$groupInfo$data);print(boxNULL)
	statsNULL = c();
	
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trial %% 1000 == 0)
			{
				cat('  Run ', trial, '\n', sep = '');
			}
			
		}
		pseudoNULL = matrix(nrow = max(groupLengths), ncol = nGroups);		
		for (gp in 1:ncol(pseudoNULL))
		{
			pseudoNULL[, gp] = sample(boxNULL, groupLengths[gp], replace = replace);
		}
		pseudoF = calculateF(data = pseudoNULL, Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat;
		statsNULL = c(statsNULL, pseudoF);
	}
	
	# compute p-values
	if (absDiffs)
	{
		p = sum(statsNULL > Fstat) / trials;
		p_left = NULL;
		p_right = NULL;
	}
	else
	{
		midNULL = mean(statsNULL);
		reflect = midNULL - (Fstat - midNULL);
			if (Fstat < 0)
			{
				p_left = sum(statsNULL < Fstat) / trials;
				p_right = sum(statsNULL > reflect) / trials;
				p = p_left + p_right;
			}
			else if (Fstat > 0)
			{
				p_right = sum(statsNULL > Fstat) / trials;
				p_left = sum(statsNULL < reflect) / trials;
				p = p_right + p_left;
			}
			else
			{
				stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
			}
	}
	

	
	if (plots)
	{
		
		if (verbose)
		{
			cat('  Computing 95% confidence intervals\n');
		}
		values = c();
		ints = matrix(nrow = ncol(data), ncol = 2);
		for (gp in 1:ncol(data))
		{
			temp = confidenceInterval(data[, gp], plots = F, Func = Func);
			values = c(values, temp$value);
			ints[gp, ] = temp$conf.int;
		}
		
		par(mfrow = c(1, 2));
		hist(statsNULL, 
			 xlim = c(min(statsNULL), max(statsNULL, Fstat)), 
			 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
			 col = col, border = border, 
			 xlab = ''
			 );
		if (absDiffs)
		{
			abline(v = Fstat, col = line.col);
		}
		else 
		{
			abline(v = Fstat, col = line.col);
			abline(v = reflect, col = line.col, lty = 'dashed');
		}
		
		if (is.null(groupNames))
		{
			groupNames = colnames(data);
		}
		else if (!is.null(groupNames) & length(groupNames)==ncol(data))
		{
			groupNames = groupNames;
		}
		boxplot(data, 
			    names = groupNames, 
			    main = paste('p = ', p, sep = ''), 
			    border = 'lightgrey', 
			    medlty = 'blank', 
			    boxwex = .5
			    );
			
		stripchart(as.data.frame(data), vertical = T, add = T, pch = pch, bg = col);
		
		coMat = matrix(nrow = ncol(data), ncol = 2);
		coMat[1, ] = c(.75, 1.25);
		for (gp in 2:nrow(coMat))
		{
			coMat[gp, ] = coMat[gp - 1, ] + 1;
		}
		for (gp in 1:length(values))
		{
			segments(coMat[gp, 1], values[gp], 
				     coMat[gp, 2], values[gp], 
				     lwd = 1, col = line.col
				     );
			segments(coMat[gp, 1], ints[gp, 1], 
					 coMat[gp, 2], ints[gp, 1], 
					 lwd = 1, lty = 'dashed', col = line.col
					 );
			segments(coMat[gp, 1], ints[gp, 2], 
				     coMat[gp, 2], ints[gp, 2], 
				     lwd = 1, lty = 'dashed', col = line.col
				     );
		}
		
	}
	
	output = list(stat = Fstat,
			      data = data,
			      null.dist = statsNULL,
			      p = p,
			      p.left = p_left,
			      p.right = p_right,
			      calcF = calcF,
			      conf.ints = ints,
			      parameters = match.call()
			      );
	
	return(output);
}
