# assumes files organized as in "logBoop_cleaned_CH.txt"
# e.g. ...
# 	[22.9845, 'f']
# 	[23.0245, 'p']
# 	[23.0395, 'p']
# 	...
########

getData = function(filename, sigFigs=6)
{
	dat = read.table(filename, header=F, sep='');
	dat = dat[, 1:2];
	colnames(dat) = c('timestamp', 'behavior');
	dat[,1] = gsub('[', '', dat[,1], fixed=T);
	dat[,1] = gsub(',', '', dat[,1], fixed=T); 
	dat[,1] = signif(as.numeric(dat[,1]), digits=sigFigs);
	dat[,2] = as.character(dat[,2]);					
	dat_table = table(dat[, 2]);						
	
	dat_num = dat;
	for (row in 1:nrow(dat_num))
	{													
		dat_num[row, 2] = as.numeric(match(dat_num[row, 2], names(dat_table)));	
	}
	dat_num_table = table(dat_num[, 2]);				
	
	##### check output for errors #####
	if ( length(dat_table) != length(dat_num_table) ) 
	{ 
		stop('RAW AND NUMERIC DATASETS DONT MATCH');
	} 
	ok = c();
	for (ind in 1:length(dat_table))
	{
		check_dat = which(dat[, 2] == names(dat_table)[ind]);
		check_dat_num = which(dat_num[, 2] == names(dat_num_table)[ind]);
		check = sum(check_dat == check_dat_num);
		if (check != dat_table[ind] | check != dat_num_table[ind])
		{
			stop('RAW AND NUMERIC DATASETS DONT MATCH');
		}
		else
		{
			ok = c(ok, T);
		}
	}
	##### done checking output #####
	
	return(list(rawData = dat, 
				numData = dat_num, 
				freqs = dat_table,
				check = ok)
				);			
}
#################

getDataRA = function (filename) {
	dat = as.character(read.table(filename, fill=T, sep=' ', colClasses=(c('character','character','character'))));
	freqs = table(dat);
	dat_num = match(dat, names(freqs));
	return(list(rawData=dat, numData=dat_num, freqs=freqs));
}
#################
buildTransitionProbsFromFile = function(filename, datVersion=c('CH','RA'), ...)
{
	if (datVersion=='CH') {
		data = getData(filename);
		raw_data = data$rawData[,2]; 
		num_data = data$numData[,2];
		if (sum(data$check) != length(freqs)){stop('PROBLEM WITH getData() FUNCTION\n\t...FIND AUSTIN OR FIGURE IT OUT!')}
	} else if (datVersion=='RA') {
		data = getDataRA(filename);
		raw_data = data$rawData;
		num_data = data$numData;
	} else {
		stop('INVALID ARG datVersion, must be "CH" or "RA"');
	}
	freqs = data$freqs;
	
	
	temp = syntaxEntropy(num_data, ...);
	prob_matrix = temp$transProbs;
	rownames(prob_matrix) = names(freqs);
	colnames(prob_matrix) = names(freqs);
	
	return(list(probMatrix=prob_matrix, data=data));
}

###################
buildDotFile = function(probMatrix, originalDataVec, file='', title='untitled', fontsize=24)
{
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# # get behavior frequencies
	freqs = table(originalDataVec);
	#print(freqs);print(probMatrix)
	# # check that behaviors are in same order in freqs and probMatrix
	if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
	# # compute proportions of behaviors for relative node size
	for (beh in 1:length(freqs))
	{
		prop = freqs[beh] / sum(freqs) * 10; #print(file)
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
	
	return(NULL)
}

##################
buildGraphFromFile = function(filename, dotFileBase=NULL, datVersion=c('CH','RA'), ...)
{
	out = buildTransitionProbsFromFile(filename, datVersion=datVersion, ...);
	if ( !(is.null(dotFileBase)) )
	{
		outfilename = paste(dotFileBase, '.dot', sep='');
	}
	else
	{
		outfilename = paste(gsub('.txt', '', filename, fixed=T), '.dot', sep='');
	}
	if (datVersion=='CH') {
		rawData = out$data$rawData[,2];
	} else if (datVersion=='RA') {
		rawData = out$data$rawData;
	} else {
		stop('INVALID ARG datVersion, must be "CH" or "RA"');
	}
	buildDotFile(out$probMatrix, 
				 rawData, 
				 file=outfilename
				 );
	return(list(file=outfilename, matrix=out$probMatrix, data=out$data$rawData));
}

#######################################################################################################################################
########################################################################################################################################
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
	#print(data)
	#compute probability
	total_leader = sum(data == leader)
	if (data[length(data)]==leader) {
		if (total_leader==1) {
			total_leader=total_leader;
		} else {
			total_leader = total_leader-1;
		}
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
syntaxEntropy <- function(syntaxIn=NULL, probabilityMatrix=NULL, verbose=0)
{
	
	if (is.null(syntaxIn) & !is.null(probabilityMatrix) & is.matrix(probabilityMatrix)) 
	{
		probMatrix = probabilityMatrix;		###ADD ERROR CHECKING, E.G. CHECK ROW SUMS###
	}
	else
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
		#print(uniqueSyls)
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
				tmp = transProbString(syntaxIn, leader, follower);#print(tmp)
				probMatrix[leaderRow, followerCol] = tmp$trans_prob;	
			}
		
			# add termination probability for current leader to vector
			#termProbs[leaderRow] = tmp$terminations / tmp$total_leader;
			#if (verbose > 0)
			#{
			#	cat('termProbs: ',termProbs,'\n\n');
			#}		
		}
	}
	
	#for each syllable type, calculate transition entropy with normalization by motif max possible entropy
	transEntropy <- vector()
	#transEntropyNoNorm <- vector()
	#print(probMatrix)
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
		cols <- ncol(probMatrix);
		maxH = ncol(probMatrix) * ( (-1/ncol(probMatrix))*log2(1/ncol(probMatrix)) );
		h.norm <- sum(h) / maxH;
		transEntropy <- c(transEntropy,h.norm)
		#transEntropyNoNorm <- c(transEntropyNoNorm,sum(h))
	}
	
	#weighted syllable entropy calculation step
	# count.totals <- vector()
	# for (value in unique(syntaxIn)) ###NEED TO FIX REFERENCE TO syntaxIn IN CASE NULL###
	# {
		# count.temp <- sum(syntaxIn==value)
		# count.totals <- c(count.totals,count.temp)
	# }
	# norm <- count.totals/max(count.totals)
	# transEntropyWeighted <- transEntropy*norm
	
	out <- list(transProbs = probMatrix,
				entropy = mean(transEntropy),
				transEntropy = transEntropy,
				stereotypy = 1-mean(transEntropy), 
				maxH = maxH
				#entropy.weighted = mean(transEntropyWeighted),
				#stereotypy.weighted = 1-mean(transEntropyWeighted),
				#entropyNoNorm = mean(transEntropyNoNorm),
				#stereotypyNoNorm = 1-mean(transEntropyNoNorm)
				)
	return(out)	
}

####################

