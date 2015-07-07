# assumes files organized as in "logBoop_cleaned_CH.txt"
getData = function(filename)
{
	dat = read.table(filename);
	dat = dat[, 1:2];
	dat[,1] = gsub('[', '', dat[,1], fixed=T);
	dat[,1] = gsub(',', '', dat[,1], fixed=T); #print(head(dat))
	
	temp = table(dat[, 2]); #print(temp)
	lookup = as.data.frame(matrix(nrow=length(temp), ncol=2));
	lookup[, 1] = names(temp);
	lookup[, 2] = 1:length(temp); #print(lookup)
	
	dat_num = dat;
	
	for (row in 1:nrow(dat_num))
	{
		lookup_ind = match(dat_num[row, 2], lookup[,1]);
		dat_num[row, 2] = lookup[lookup_ind, 2];
	}
	
	return(list(rawData=dat, numData=dat_num, freqs=temp));
}
#################

buildTransitionProbsFromFile = function(filename)
{
	data = getData(filename);
	raw_data = data$rawData[,2]; 
	num_data = data$numData[,2];
	freqs = data$freqs;
	
	temp = syntaxEntropy(num_data);
	prob_matrix = temp$transProbs;
	rownames(prob_matrix) = names(freqs);
	colnames(prob_matrix) = names(freqs);
	
	return(list(probMatrix=prob_matrix));
}
#################

buildGraphFromFile = function()

###################
buildDotFile = function(probMatrix, originalDataVec, file='', title='untitled', fontsize=24)
{
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# get behavior frequencies
	freqs = table(originalDataVec);
	
	# check that behaviors are in same order in freqs and probMatrix
	if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
	# compute proportions of behaviors for relative node size
	for (beh in 1:length(freqs))
	{
		prop = freqs[row] / sum(freqs) * 10; #beh not row?
		cat('		', names(freqs)[beh], ' [width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
	}
	
	# loop through probMatrix to get probabilities
	probMatrix=probMatrix*10;
	for (row in 1:nrow(probMatrix))
	{
		for (col in 1:ncol(probMatrix))
		{
			val = probMatrix[row,col];
			if (val > 0)
			{
				leader = rownames(probMatrix)[row];#print(leader);
				follower = colnames(probMatrix)[col];#print(follower);
		
				cat('		', leader, ' -> ', follower, ' [label="", style="setlinewidth(', val, ')", arrowsize=1];','\n' ,sep='', file=file, append=T);	
		 	}
		}
	}
	
	# write last line of file
	cat('	}', file=file, append=T);
	return();
}

#################
#Computes transition probability for string-based analysis.
#################
transProbString = function(data, leader, follower) #Used by syntaxSimilarityDiffs()
{
	#set initial values for # of transitions
	transitions = 0
	terminations = 0
	length=length(data)
	syl = c(1:length(data))
	
	for (pos in syl)
	{
		#check if leader is in current position
		check = data[pos] == leader;
		
		#if leader is in current position
		if (check)
		{
			#check whether current position is the end of the string
			if(is.na(data[pos+1]==follower))
			{
				transitions = transitions+0		# what is this?
			}
			
			else 
			if (data[pos+1]==follower)
			{
				transitions = transitions+1
			}
		}
	}
	
	#compute probability
	total_leader = sum(data == leader)
	if (data[length(data)]==leader)
	{
		total_leader = total_leader-1
	}
	prob = transitions/total_leader
	
	output = list(transitions, total_leader, prob);
	names(output) = c('transitions', 'total_leader', 'trans_prob');
	return(output);
}

#################



#################
#Creates transition probability table and computes weighted and unweighted syntax entropy scores as in Miller et al. (2010).
#################
syntaxEntropy <- function(syntaxIn, verbose=0)
{
	#find the unique syllables and count them
	uniqueSyls <- sort(unique(as.numeric(syntaxIn)))
	nSyls = length(uniqueSyls)
	
	# intialize empty matrix, each row is a syllable
	# columns denote transition probabilities to other syllables and end of motif
	probMatrix = matrix(nrow = nSyls, ncol=nSyls);
	rownames(probMatrix) = paste('syl', uniqueSyls, sep = '');
	colnames(probMatrix) = paste('syl', uniqueSyls, sep = '');
	
	# loop over syllables, treat each as leader in turn
	for (leader in uniqueSyls)
	{
		if (verbose > 0)
		{
			cat('leader:',leader,'\n');
		}		
		
		# find row for current leader in probability matrix
		leaderRow = match(paste('syl', leader, sep = ''), rownames(probMatrix));
		#cat('leaderRow:',leaderRow,'\n')
		
		# loop over syllables, treat each as follower in turn
		for (follower in uniqueSyls)
		{
			if (verbose > 0)
			{
				cat('follower:',follower,', ');	
			}
			
			# find column for current follower in probability matrix
			followerCol = match(paste('syl', follower, sep = ''), colnames(probMatrix));
			#cat('followerCol:',followerCol,'\n')
			
			# compute transition probability and store in matrix
			tmp = transProbString(syntaxIn, leader, follower);
			probMatrix[leaderRow, followerCol] = tmp$trans_prob;	
		}
		
		# add termination probability for current leader to vector
		#termProbs[leaderRow] = tmp$terminations / tmp$total_leader;
		if (verbose > 0)
		{
			cat('termProbs: ',termProbs,'\n\n');
		}		
	}
	
	#for each syllable type, calculate transition entropy with normalization by motif max possible entropy
	transEntropy <- vector()
	transEntropyNoNorm <- vector()
	for (row in 1:nrow(probMatrix))
	{
		h <- vector()
		
		for (col in 1:ncol(probMatrix))
		{
			if (probMatrix[row,col]==0)
			{
				h <- c(h,0)
			}
			
			if (!probMatrix[row,col]==0)
			{
				tempEntropy <- (-probMatrix[row,col])*log2(probMatrix[row,col])
				h <- c(h,tempEntropy)
			}	
		}
		cols <- ncol(probMatrix)
		h.norm <- sum(h)/(ncol(probMatrix)*((-1/ncol(probMatrix))*log2(1/ncol(probMatrix))))
		transEntropy <- c(transEntropy,h.norm)
		transEntropyNoNorm <- c(transEntropyNoNorm,sum(h))
	}
	
	#weighted syllable entropy calculation step
	count.totals <- vector()
	for (value in unique(syntaxIn))
	{
		count.temp <- sum(syntaxIn==value)
		count.totals <- c(count.totals,count.temp)
	}
	norm <- count.totals/max(count.totals)
	transEntropyWeighted <- transEntropy*norm
	
	out <- list(transProbs = probMatrix,
				entropy = mean(transEntropy),
				transEntropy = transEntropy,
				stereotypy = 1-mean(transEntropy),
				entropy.weighted = mean(transEntropyWeighted),
				stereotypy.weighted = 1-mean(transEntropyWeighted),
				entropyNoNorm = mean(transEntropyNoNorm),
				stereotypyNoNorm = 1-mean(transEntropyNoNorm)
				)
	return(out)	
}

####################

