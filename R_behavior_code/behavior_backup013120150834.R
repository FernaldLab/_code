library(stringr)
options(stringsAsFactors = FALSE)
source("~/Desktop/Katrina/behavior_code/bootstrap_tests_June2013_STABLE.R")
#use color=blue in .dot output script to make separate sets of lines for 1st follower, 2nd, etc.
# Find a way to represent AB -> C, ABC -> D instead of only A -> B probabilities
# collapse statistics across animals. Arrows only allowed from subj to different subj.
# frame counts!!!
# Ask scott - what defines start of spawning???
# maybe paths leading up to spawning - size of circle reps # of paths. Forks and such.
#heatmaps with k!
# hella basic stats to convince Scott this is legitimite
# HEATMAP OF COUNTS FOR INDIVIDUAL SAMPLES. 6x6. DO THE ONES CLUSTER.
#bouts

# behaviors <- names(table(sjan_data[[6]]$dat$beh))
# behaviors <- c(behaviors[1:2], 'c', behaviors[3:length(behaviors)]) #hacky and not generally correct
# mat = matrix(nrow = length(behaviors), ncol = 6)
# rownames(mat)<- behaviors
# colnames(mat)<- names(sjan_data)
# for (i in 1:length(sjan_data)) {tab = table(sjan_data[[i]]$dat$beh);
# for (j in 1:length(behaviors)) { if (behaviors[j] %in% names(tab)) {mat[behaviors[j],i] <- tab[behaviors[j]]} else {mat[behaviors[j],i] <- 0}
# }}

# > eval(parse(text="tmp[order(tmp[,4], tmp[,3], tmp[,2], tmp[,6], tmp[,7]),]"))



# heatmap(cor(t(probma)), symm = TRUE)










# Source: ethograms_from_scorevideo.R
# Calls .getData() on every file in <folderPath> and returns the results in a list.
# If groups = TRUE, will recurse through directories. Use with one folder for control animals,
# one folder for experimental condition 1, etc.
.getDataBatch = function (folderPath, groups = FALSE) {
	filenames = if (groups) {paste(folderPath,list.files(folderPath, pattern = "txt$", recursive = TRUE),sep = "");}
	            else {paste(folderPath,list.files(folderPath, pattern = "txt$"),sep = "");}
	data = list();
	for (f in 1:length(filenames)) {
		print(filenames[f])
		data[[f]] = .getData(filenames[f]);
		names(data)[f] = gsub(folderPath, "", filenames[f]);
	}
	return(data);
}


# Reads in a score log (.txt file) and returns a data frame with columns
#   time    behavior    subject    type    pair_time    duration
.getData = function (filename) {
	data0 = read.table(filename, fill=T, colClasses='character', sep='\t', header=F, quote='', blank.lines.skip=T, strip.white=T);
	desc_table = .getDescriptionTableFromRawData(data0);
	
	fpsRow = which(grepl('^Frames/sec of files:', data0[,1]));
	framesPerSecond = as.numeric(gsub('[^0-9]', '', data0[fpsRow,1]));
	
	df = .parseFullLog(data0, desc_table, framesPerSecond);
	
	
	df = df[order(as.numeric(df$time)),];
	dimnames(df)[[1]] <- 1:(dim(df)[1]);
	
	if (nrow(df)<1) {
		warning(paste('No data in ', filename, sep=''));
		return(list(data=df, codes=NULL));
	}
#	desc_table = desc_table[as.character(desc_table[, 1]) %in% as.character(df$behavior), ];

#   could return desc_table here if desired
	return(df);
}

# Source: ethograms_from_scorevideo.R
# Helper function for .getData
# Gets description table
.getDescriptionTableFromRawData = function (data0) {
	desc_row = which(grepl('^start.*description$', data0[,1])) + 2;
	desc_table = matrix(nrow=1, ncol=3);
	while (TRUE) {
		if (grepl('^-.*-$', data0[desc_row, ])) {
			break;
		} else {
			tmp = unlist(strsplit(data0[desc_row, ], '[0123]'));
			code = unlist(strsplit(gsub(' ','',tmp[1]),''));
			if (length(code) > 1) {
				code <- code [1];
			}
			description = gsub('  ','',tmp[length(tmp)]);
			subject = unlist(strsplit(data0[desc_row, ], '[^0123]'));
			desc_table = rbind(desc_table, c(code, description, subject[15]));
			desc_row = desc_row + 1;
		}
	}
	desc_table = desc_table[-1, ];
	dimnames(desc_table)[[2]] <- c("code", "description", "subj")
	return(desc_table);
}

# Helper function for .getData
# Parses full log to populate data frame
.parseFullLog = function (data0, desc_table, framesPerSecond) {
	start_full = which(data0=='FULL LOG') + 4;
	end_full = which(data0=='NOTES') - 2;

	full = strsplit(data0[start_full:end_full, ], '    ');
	full = lapply(full, function(f) gsub('^ *','', gsub(' *$','', f)));
	full = lapply(full, function(f) f[f!=""]);
	
	# Get time, behavior, subject, type
	df = data.frame(matrix(nrow=1, ncol=6));
	for (entry in full) {
		isStart = if (length(entry)==5) {entry[5]} else {NA}
		newRow = c(as.numeric(entry[1]), as.character(entry[3]), as.character(entry[4]), isStart, NA, NA)
		df = rbind(df, newRow);
	}
	names(df) = c('time', 'behavior', 'subject', 'type', 'pair_time', 'duration');
	df$type <- as.factor(df$type);
	df$time <- as.numeric(df$time) / framesPerSecond;
	df$pair_time <- as.numeric(df$pair_time);
	df = df[-1, ];
	
	# Match up starts and stops; calculate durations
	for (i in 1:(dim(df)[1])) {
		entry = df[i,];
		if (!is.na(entry$type) && entry$type == 'stop') {
			this_time = entry$time;
			startIndex = 0;
			for (j in i:1) {
				if (df[j,]$behavior == entry$behavior && df[j,]$type=='start') {
					startIndex = j;
					break;
				}
			}
			df$pair_time[startIndex] <- entry$time;
			df$duration[startIndex] <- df$pair_time[startIndex] - df$time[startIndex]
			df$pair_time[i] <- df$time[startIndex]
			df$duration[i] <- df$duration[startIndex]
		}
	}
	return(df);
}

# Renames behaviors to differentiate between starts and stops
.renameStartStop = function(data) {
	for (i in 1:dim(data)[1]) {
		if (!is.na(data$type[i])) data$behavior[i] <- paste(data$behavior[i], " ", as.character(data$type[i]), sep = "");  
	}
	return(data);
}

.groupLevelProbMats = function(data, byTotal = FALSE) {
	probMats <- lapply(data, function(d) {.getProbabilityMatrix(d$behavior, byTotal=byTotal)});
	counts <- lapply(data, function(d) {table(d$behavior)})
	behnames = names(table(unlist(lapply(probMats, function(d){dimnames(d)[[1]]}))));
	groupNames = names(table(gsub("/.+", "", names(data))));
	probMatsByGroup = list();
	for (i in 1:length(groupNames)) {
		groupData = list(probMats = probMats[grepl(paste("^", groupNames[i], sep = ""), names(probMats))],
		                 counts = counts[grepl(paste("^", groupNames[i], sep = ""), names(counts))]);
		probMatsByGroup[[i]] = .combineProbabilityMatrices(groupData, behnames, byTotal);
		names(probMatsByGroup)[i] <- groupNames[i];
	}
	return(probMatsByGroup);
}

.combineProbabilityMatrices = function(data, behnames = NULL, byTotal) {
	if (is.null(behnames)) {behnames = names(table(unlist(lapply(probmas, function(d){dimnames(d)[[1]]}))));}
	probmas = data$probMats;
	counts = data$counts;
	probMat = matrix(data = numeric(length(behnames) * length(behnames)),
					 nrow = length(behnames), ncol = length(behnames), dimnames = list(behnames, behnames));
	countVec = numeric(length(behnames));
	names(countVec) <- behnames;
	nfish <- countVec;
	for (row in rownames(probMat)) {
		for (subject in 1:length(probmas)) {
			if (row %in% rownames(probmas[[subject]])) {
				countVec[row] <- countVec[row] + counts[[subject]][row];
				nfish[row] <- nfish[row] + 1;
				for (col in colnames(probmas[[subject]])) {
					probMat[row,col] <- probMat[row,col] + probmas[[subject]][row,col];
				}
			}
		}
	}
	
	if (byTotal) {
		probMat = probMat / length(probmas);
	} else {
		for (col in colnames(probMat)) {
			probMat[,col] <- probMat[,col] / countVec;
		}
	}
	return(list(probMat = probMat, counts = countVec, nfish = nfish));
}

# Source: ethograms_from_scorevideo.R
# Reads in a vector giving a sequence of behaviors and returns a matrix giving the transitional
#  probability for each pair of behaviors.
# Usage: .getProbabilityMatrix(.renameStartStop(dataFrame)$behavior)      to include behavior starts and ends
.getProbabilityMatrix = function (data, removeZeroCol=F, ...) {
	beh = names(table(data[-length(data)]));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeTransitionProbability(data=data, leader=leader, follower=follower, ...);
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$probability;
		}
	}
	if (removeZeroCol) {
		colSums = apply(probMat, 2, sum);
		if (any(colSums==0)) {
			probMat = probMat[-which(colSums==0), -which(colSums==0)];
		}
	}
	return(probMat);
}

# Source: ethograms_from_scorevideo.R
# Helper function for .getProbabilityMatrix
# Computes the transition probability between leader and follower
.computeTransitionProbability = function (data, leader, follower, byTotal = FALSE) {
	count = 0;
	termination = 0;
	for (i in 1:length(data)) {
		if (data[i] == leader) {
			if (i == length(data)) {
				termination = 1;
			} else if (data[i+1] == follower) {
				count = count + 1;
			}
		}
	}
	total_leader = if(termination) {sum(data == leader) - 1} else {sum(data == leader)};
	prob = if (byTotal) {count / length(data)} else if(total_leader) {count / total_leader} else {0} ;
	return(list(probability=prob, termination=termination, count_transitions=count, count_leader=total_leader));
}


# Reads in a data frame giving a sequence of behaviors with time numbers and returns a matrix giving a list containing:
#   1. timeDistances, the average amount of time between two behaviors,
#   2. behDistances, the average number of behaviors between two behaviors (minimum 1)
#   3. probFollowed, the probability that the leading behavior was followed by the follower at some point. 
# Usage: .getIntervalMatrix(.renameStartStop(dataFrame))
.getIntervalMatrix = function (data) {
	beh = names(table(data$behavior));
	timeMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	behMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeAverageInterval(data=data, leader=leader, follower=follower);
			timeMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$average_times;
			behMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$average_behaviors;
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] =
				(tmp$count_with / (tmp$count_with + tmp$count_without + tmp$count_end));
		}
	}
	return(list(timeDistances = timeMat, behDistances = behMat, probFollowed = probMat));
}


# Helper function for .getIntervalMatrix
# Computes average time interval, count interval, and returns counts of different situations.
# count_without is the number of times leader was followed by leader before it was followed by follower
# count_end is 1 if there is no follower after the last leader, and 0 otherwise
# count_with is the number of times leader was followed by follower before the end of behavior and before
#       another leader occurred
.computeAverageInterval = function (data, leader, follower) {
	count = 0;
	countWithout = 0; 
	totaltimes = 0;
	totalBehaviors = 0;
	for (i in 1:dim(data)[1]) {
		if (data$behavior[i] == leader) {
			lastStartTime = as.numeric(data$time[i]);
			j = i+1;
			while (j <= dim(data)[1] && data$behavior[j] != follower && data$behavior[j] != leader) {j = j + 1;}
			if (j > dim(data)[1]) {
				countEnd = countEnd + 1;
			} else if (data$behavior[j] == follower) {
				count = count + 1;
				totalTimes = totalTimes + (data$time[j] - lastStartTime);
				totalBehaviors = totalBehaviors + (j-i);
			} else if (data$behavior[j] == leader) {
				countWithout = countWithout + 1;
			} else {
				stop('Something is wrong in ".computeAverageInterval()". GO FIND KATRINA!!!')
			}
		}
	}
	if (count > 0) {
		avgF = totalTimes / count;
		avgB = totalBehaviors / count;
	} else {
		avgF = NA;
		avgB = NA;
	}
	return(list(average_times=avgF, average_behaviors = avgB, count_with = count, count_without=countWithout, count_end=countEnd));
}

# Returns a vector giving the distances between occurances of <centerBeh> and <varBeh> in data frame <data>.
# If <noRepCenterBeh> is false, the vector has length count(centerBeh)*count(varBeh), and it has distances between every occurance of each behavior.
# If <noRepCenterBeh> is true, the vector has length roughly 2*count(varBeh), as it only measures distances between a varBeh and its two neighboring centerBehs.
# Optionally you can provide a buffer, specifying not to count centerBehs that occur too close to the beginning or end of data. This buffer
# can be given in times or behavior counts.
# This function returns a list containing:
#     1. <timeDists>, a vector of the distances in times
#     2. <behDists>, a vector of the distances in behaviors
#     3. <centerCount>, the number of occurances of behavior centerBeh
#     4. <varCount>, the number of occurances of behavior varBeh
.getAllIntervals = function (data, centerBeh, varBeh, startBuffer = 0, endBuffer = 0, bufferInTimes = TRUE, noRepCenterBeh = TRUE) {
	varBehLocs = which(data$behavior == varBeh);
	varBehTimes = data$time[varBehLocs];
	centerBehLocs = which(data$behavior == centerBeh);
	if (bufferInTimes) {
		max_time = max(data$time);
		centerBehTimes = data$time[centerBehLocs];
		centerBehLocs = centerBehLocs[centerBehTimes > startBuffer & centerBehTimes <= max_time - endBuffer];
	} else {
		centerBehLocs = centerBehLocs[as.numeric(centerBehLocs) > startBuffer & as.numeric(centerBehLocs) <= length(data$behavior - endBuffer)];
	}
	centerBehTimes = data$time[centerBehLocs];
	
	
	behDists = numeric();
	timeDists = numeric();
	
	for (i in 1:length(centerBehLocs)) {
		if (noRepCenterBeh) {
			if (i == 1) {
				behDists = c(behDists, varBehLocs[varBehLocs <= centerBehLocs[i+1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes <= centerBehTimes[i+1]] - centerBehTimes[i]);				
			} else if (i == length(centerBehLocs)) {
				behDists = c(behDists, varBehLocs[varBehLocs >= centerBehLocs[i-1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes >= centerBehTimes[i-1]] - centerBehTimes[i]);				
			} else {
				behDists = c(behDists, varBehLocs[varBehLocs >= centerBehLocs[i-1] & varBehLocs <= centerBehLocs[i+1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes >= centerBehTimes[i-1] & varBehTimes <= centerBehTimes[i+1]] - centerBehTimes[i]);
			}
		} else {
			behDists = c(behDists, varBehLocs - centerBehLocs[i]);
			timeDists = c(timeDists, varBehTimes - centerBehTimes[i]);
		}
	}
	
	return(list(behDists=behDists, timeDists=timeDists, centerCount = length(centerBehLocs), varCount = length(varBehLocs)));
}
#TODO raster plot - ethograms.dataSummary()

# Calls .getIntervalsAcrossFish on every data frame in the list <data>, and combines the results into one
# list. behDists and timeDists contain all the distances from all the fish, and varCount and centerCount
# contain the total counts across all fish.
.getIntervalsAcrossFish = function (data, ...) {
	alldata <- list(behDists = numeric(), timeDists = numeric(), centerCount = 0, varCount = 0);
	for (i in 1:length(data)) {
		tmp <- .getAllIntervals(data[[i]], ...);
		alldata$behDists <- c(alldata$behDists, tmp$behDists);
		alldata$timeDists <- c(alldata$timeDists, tmp$timeDists);
		alldata$centerCount <- alldata$centerCount + tmp$centerCount;
		alldata$varCount <- alldata$varCount + tmp$varCount;
	}
	return(alldata);
}

# behaviorsToPlotAndColors should have column 1 as a list of the behaviors to plot and column 2 as the color each behavior should be.
# Usage: .plotGraphs(.renameStartStop(mariana_dat[[1]]), bc, "Male In Pot start")
.plotGraphs = function(data, behaviorsToPlotAndColors, centerBeh, lim, noRepCenterBeh = TRUE, multifish = FALSE, timesPerBin = 30, ymax = 1) {
	histBreaks = ((-ceiling(lim/timesPerBin) - 1):(ceiling(lim/timesPerBin)) + 0.5) * timesPerBin;
	for (i in 1:(dim(behaviorsToPlotAndColors)[1])) {
		x <- if (multifish) {.getIntervalsAcrossFish(data, centerBeh = centerBeh, varBeh = behaviorsToPlotAndColors[i, 1], noRepCenterBeh = noRepCenterBeh);}
			 else {.getAllIntervals(data, centerBeh, behaviorsToPlotAndColors[i, 1], noRepCenterBeh = noRepCenterBeh);}
		h <- hist(x$timeDists[x$timeDists > -lim & x$timeDists < lim], breaks = histBreaks, plot = FALSE);
		count <- x$varCount;
		#print(paste("mids length =", length(h$mids), ", counts length =", length(h$counts), ", count =", count, ",extra length =", length(h$counts / count)))
		plot(x = h$mids, y = h$counts / count, type = "l", col = behaviorsToPlotAndColors[i, 2], xlim = c(-lim, lim), ylim = c(0, ymax));
		par(new = TRUE);
	}
	par(new = FALSE);
}




# Source: ethograms_from_scorevideo.R
# Reads in a probability matrix and returns a list containing:
#    1. $probMat, the input probability matrix
#    2. $hMat, a matrix of entropy
#    3. $h, a vector of the raw entropy for each behavioral code
#    4. $h_norm, a vector of the normalized entropy for each behavioral code
#    5. $h_max, a value representing the maximum possible entropy (which was used to normalize
#        the entropy values)
.computeEntropyProbMatrix = function (probMat) {
	hMat = probMat;
	for (row in 1:nrow(probMat)) {
		hMat[row, ] = .computeEntropyOneState(probMat[row, ]);
	}
	num_states = ncol(probMat);
	h_max = num_states * ( (-1/num_states)*log2(1/num_states) );
	h = apply(hMat, 1, sum);
	h_norm = h / h_max;
	return(list(probMat=probMat, hMat=hMat, h=h, h_norm=h_norm, h_max=h_max));
}

# Source: ethograms_from_scorevideo.R
# Helper function for .computeEntropyProbMatrix
.computeEntropyOneState = function (probVec) {
	h = c();
	for (follower in 1:length(probVec)) {
		prob = probVec[follower];
		if (prob == 0) {
			h = c(h, 0);
		} else {
			h = c(h, (-prob)*log2(prob));
		}
	}
	names(h) = names(probVec);
	return(h);
}




# Source: behavior_syntax.R
# Given a probability matrix and a data vector, writes a dot file to <file> (default is to console)
.buildDotFile = function (probMatrix, originalDataVec, file='', title='untitled', fontsize=24, minValForLine = 0, weird = FALSE, singleCharLabels = TRUE, byTotal = FALSE) {
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# # get behavior frequencies
	if (is.numeric(originalDataVec)) {freqs = originalDataVec;}
	else if (is.character(originalDataVec)) {freqs = table(originalDataVec);}
	else {stop("DATA VEC IS WRONG IN .buildDotFile!!!!!!!");}
	#print(freqs);print(probMatrix)
	# # check that behaviors are in same order in freqs and probMatrix
	if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
	# # compute proportions of behaviors for relative node size
	for (beh in 1:length(freqs))
	{
#		if (names(freqs)[beh] %in% c("[","]")) {next;}
		prop = freqs[beh] / sum(freqs) * 10; #print(file)
		if (!singleCharLabels || prop < 0.7) { # 0.7 is the magic size for single characters in 24pt font
			cat('		', gsub(' ', '', names(freqs)[beh]), ' [label="", xlabel="', gsub(' ', '', names(freqs)[beh]),'", width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
		} else {
			cat('		', gsub(' ', '', names(freqs)[beh]), ' [width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
		}
	}
	
	# loop through probMatrix to get probabilities
	# probMatrix=probMatrix*10;
	for (row in 1:nrow(probMatrix))
	{
		for (col in 1:ncol(probMatrix))
		{
			val = if (byTotal) {(probMatrix[row,col] / sum(probMatrix)) * 100} else {(probMatrix[row,col] / sum(probMatrix[row,])) * 10};
			if (weird) {val = 10 * (1 - (probMatrix[row,col] / max(probMatrix[row,])));} #FIX THIS !!!TODO
			if (val > minValForLine)
			{
				leader = rownames(probMatrix)[row];
				follower = colnames(probMatrix)[col];#print(paste(leader, follower, sep = ','));
				if (leader == "STOP") {next;}
#				if (follower %in% c("[","]")) {next;}
		
				cat('		', gsub(' ', '', leader), ' -> ', gsub(' ', '', follower),
				    ' [label="", style="setlinewidth(', val, ')", arrowsize=1];','\n' ,sep='', file=file, append=T);	
		 	}
		}
	}
	
	# write last line of file
	cat('	}', file=file, append=T);
	
	return(NULL)
}

# Takes a data frame, converts the descriptions to single letters using key (or using
# a smart algorithm to assign letters if no key is provided), and returns a list of the
# data frame and the key used. A key should be a list whose names are the behavior names
# and whose entries are single-letter codes.
.descriptionsToLetters = function (data, key = NULL) {
	frequencies = table(data$behavior);
	if (is.null(key)){key = .assignLetters(frequencies);}
	else {.validateKey(key, names(frequencies))}
	for (i in 1:dim(data)[1]) {
		data$behavior[i] <- key[[data$behavior[i]]];
	}
	descTable = data.frame(code = as.character(key), description = names(key));
	#dimnames(descTable)[[1]] <- 1:dim(descTable)[1];
	return(list(data, descTable));
}

.validateKey = function (key, behaviors) {
	if (!is.list(key)) stop("Error: Key provided to .descriptionsToLetters is not a list.\nKey should be a list of single-letter codes, with the names of key set to the behavior names.");
	if (sum(!(behaviors %in% names(key))) != 0) stop("Error: Key provided to .descriptionsToLetters does not contain a code for all behaviors in the data provided.");
	asvector = as.character(key);
	if(sum(nchar(asvector) != 1) != 0) {warning("Not all codes in key are single characters");}
}

.assignLetters = function (freqtable) {
	behList = names(freqtable[order(freqtable, decreasing = TRUE)]);
	key = list();
	orphans = character();
	for (behaviorName in behList) {
		preferredCodes = toupper(substr(unlist(strsplit(behaviorName, " ")), 1,1));
		for (code in preferredCodes) {
			if ( !(code %in% as.character(key)) ) {
				key[[behaviorName]] <- code;
				break;
			}
		}
		if (! behaviorName %in% names(key)) {orphans <- c(orphans, behaviorName);}
	}
	
	for (behaviorName in orphans) {
		for (code in c(toupper(letters), letters)) {
			if ( !(code %in% as.character(key)) ) {
				key[[behaviorName]] <- code;
				break;
			}
		}
		if (! behaviorName %in% names(key)) {stop("More than 52 behaviors are present. .assignLetters() needs modification.")}
	}
	return(key);
}

.sepSubject = function(data, toSeparate, separationTable) {
	for ( i in 1:length(data$behavior)) {
		if (data$behavior[i] == toSeparate)
			data$behavior[i] <- separationTable[data$subject[i]]
	}
	return(data);
}


# OPTIONS:
# startTime, endTime - only get behavior from times [startTime, endTime]
# subject - only consider this subject's behavior
# startOnly - only consider starts of behaviors (ignore ends)
# boutInterval - interval to separate bouts, in times. Default is null (no bout separation)
# minNumBehaviors - remove behaviors that do not occur at least this many times.
# toExclude - remove all behaviors in this character vector.
.filterData = function(data, startTime = NA, endTime = NA, subject = NULL, startOnly = FALSE,
					   boutInterval = NULL, minNumBehaviors = NULL, toExclude = NULL, splitPot = FALSE) {
	if (!is.na(startTime)) {data <- data[data$time >= startTime,];}
	if (!is.na(endTime)) {data <- data[data$time <= endTime,];}
	if (!is.null(subject)) {
		data <- data[data$subject == subject,]
		if(length(data$behavior) == 0) {stop(paste("Error: Incorrect Subject Name:", subject))}
	}
	if (startOnly) {data <- data[is.na(data$type) | data$type != "stop",];}
	if (!is.null(boutInterval)) {
		data = .separateBouts(data, boutInterval);
	}
	if (!is.null(minNumBehaviors)) {
		freqtable = table(data$behavior);
		behaviorsToTrash = names(freqtable)[freqtable < minNumBehaviors]
		data <- data[!(data$behavior %in% behaviorsToTrash),];
		if(length(data$behavior) == 0) {stop(paste("Error: minNumBehaviors is set too high. No behavior occurs", minNumBehaviors, "times."))}
	}
	if (!is.null(toExclude)) {
		data <- data[!(data$behavior %in% toExclude),];
	}
	if (splitPot) {
		v <- c("male IN POT", "female IN POT");
		names(v) <- c("male", "female");
		data <- .sepSubject(data, "inside POT", v)
	}
	return(data);
}

.filterDataList = function(data, ...) {
	for (i in 1:length(data)) {
		data[[i]] <- .filterData(data[[i]], ...);
	}
	return(data);
}


#----------------------
# OPTIONS:
# k - number of spots that are ok to count for transitional probabilities. weighted makes behaviors closer to leader count more.
# minValForLine - minimum probability*10 value to justify drawing an arrow
.makeDotPlot = function (data, filename, k = 1, weighted = FALSE, minValForLine = 0) {
	data = .filterData(data, ...)
	
	cleanerDataForPlot <- .descriptionsToLetters(.renameStartStop(data));
	behaviorVec <- cleanerDataForPlot[[1]]$behavior;
	probMatrix <- .getProbabilityMatrixK(behaviorVec, k = k, weighted = weighted);
	.buildDotFile(probMatrix, behaviorVec, file = filename, minValForLine = minValForLine);
	
	
	return(list(probMatrix = probMatrix, codes = cleanerDataForPlot[[2]], entropyData = .computeEntropyProbMatrix(probMatrix)));
}

.makeDotPlots = function (listOfDataLists, filename_ext, ...) {
	filenames = paste(names(listOfDataLists), filename_ext, sep = '_')
	for (i in 1:length(filenames)) {
		.makeDotPlot(listOfDataLists[[i]], filenames[i], ...)
	}
	#TODO make a list of the lists returned by .makeDotPlots. so entropy, probmatrix, codes, etc are not lost.
}

.makeDotPlotsFromProbMas = function(probMatDat, fileprefix, ...) {
	for (i in 1:length(probMatDat)) {
		.buildDotFile(probMatDat[[i]]$probMat, probMatDat[[i]]$counts, file = paste(fileprefix, names(probMatDat)[[i]], ".dot", sep = ""), ...);
	}
}


.separateBouts = function (data, intervalToSeparate) {
	# 	names(df) = c('time', 'behavior', 'subject', 'type', 'pair_time', 'duration');
	newData = data.frame(time = data$time[1], behavior = "START", subject = NA, type = NA, pair_time = NA, duration = NA); #TODO maybe change type to start and pair with stop
	newData = rbind(newData, data[1,]);
	for (i in 2:length(data[,1])) {
		if (as.numeric(data$time[i]) - as.numeric(data$time[i-1]) >= intervalToSeparate) {
			stopRow = data.frame(time = data$time[i-1], behavior = "STOP", subject = NA, type = NA, pair_time = NA, duration = NA);
			newData = rbind(newData, stopRow);
			startRow = data.frame(time = data$time[i], behavior = "START", subject = NA, type = NA, pair_time = NA, duration = NA);
			newData = rbind(newData, startRow);
		}
		newData = rbind(newData, data[i,]);
	}
	stopRow = data.frame(time = data$time[length(data$time)], behavior = "STOP", subject = NA, type = NA, pair_time = NA, duration = NA);
	newData = rbind(newData, stopRow);
	dimnames(newData)[[1]] <- 1:length(dimnames(newData)[[1]])
	return(newData);
}


# Reads in a vector giving a sequence of behaviors and returns a matrix giving the probability that the
# behavior in a given column follows the behavior in a fiven row within k behaviors.
.getProbabilityMatrixK = function (data, removeZeroCol=F, k = 2, weighted = FALSE) {
	beh = names(table(data));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeTransitionProbabilityK(data=data, leader=leader, follower=follower, k=k, weighted = weighted);
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$probability;
		}
	}
	if (removeZeroCol) {
		colSums = apply(probMat, 2, sum);
		if (any(colSums==0)) {
			probMat = probMat[-which(colSums==0), -which(colSums==0)];
		}
	}
	return(probMat);
}

# Helper function for .getProbabilityMatrixK
.computeTransitionProbabilityK = function (data, leader, follower, k, weighted) {
	count = 0;
	termination = 0;
	for (i in 1:length(data)) {
		if (data[i] == leader) {
			if (i == length(data)) {
				termination = 1;
			} else {
				increment = 1;
				while (increment <= k && i + increment <= length(data)) {
					if (data[i + increment] == follower) {
						count = if(weighted) {count + (k + 1) - increment} else {count + 1}
					}
					increment = increment + 1;
				}
			}
		}
	}
	total_leader = sum(data == leader);
	prob = count / total_leader;
	return(list(probability=prob, termination=termination, count_transitions=count, count_leader=total_leader));
}




# TODO figure out how to preserve time numbers, etc in 3D data structure. list beh=current, time=, subj=, type=, .......
# TODO sort.
.getAllContexts = function(behaviorName, data, k=NULL, kbefore=0, kafter=0) {
	if (!is.null(k)) {
		kbefore = k;
		kafter = k;
	}
	indicesVector = which(data$behavior == behaviorName);
	indicesVector = indicesVector[indicesVector > kbefore & indicesVector <= (length(data$behavior) - kafter)];
	if (length(indicesVector) == 0) stop(paste("Error: Zero occurances of", behaviorName, "within margins in data provided to .getAllContexts"));
	if (length(indicesVector) < 3) warning(paste("Only", length(indicesVector), "occurances of", behaviorName, "in data provided to .getAllContexts"));
	df = data.frame(centerTimeNum = data$time[indicesVector]);
	for (i in (-kbefore):kafter) {
		df = cbind(df, data$behavior[indicesVector + i]);
		dimnames(df)[[2]][length(dimnames(df)[[2]])] <- paste("n", if(i>=0){"+"}else{""}, as.character(i), sep="");
	}
	print(df);
	centerCol = which(dimnames(df)[[2]] == "n+0");
	# df = df[order(df[,centerCol + 1], df[,centerCol + 2], df[,centerCol + 3]),]      but tailored to actual num of cols.
	return(df);
}



# Replaces behavior description <toReplace> with description <replacement> in data
# frame data. An example use of this function would be to replace "Male in pot" with
# "Male In Pot"
.replaceBeh = function(data, toReplace, replacement) {
	for (i in 1:length(data$behavior)) {
		if (data$behavior[i] == toReplace) data$behavior[i] <- replacement;
	}
	return(data);
}

# Calls .replaceBeh on every data frame of a data list.
.replaceBehAll = function(data, toReplace, replacement) {
	for (i in 1:length(data)) {
		data[[i]] <- .replaceBeh(data[[i]], toReplace, replacement);
	}
	return(data);
}

# This is an example of how to use .replaceBehAll to clean up data such that all logs
# have the same descriptions for the same behavior. This accomplishes the task for the
# logs provided by Mariana.
# To find behaviors in your data that may be named differently in different logs, call
# table(unlist(lapply(<your_data>, function(f) {names(table(f$behavior))}))), which will give
# a vector with all the behaviors you scored as the names and the number of logs that behavior
# appears in as the value. You can visually inspect this output to look for duplicates.
.cleanUpMariana = function(mariana_dat) {
	mariana_clean <- .replaceBehAll(mariana_dat, "Male darts", "Male Darts");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male Darting", "Male Darts");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male in Pot", "Male In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male in pot", "Male In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male to Pot after Female", "Male In Pot")
	mariana_clean <- .replaceBehAll(mariana_clean, "Female flees", "Female Flees");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female follows", "Female Follows");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female in pot", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female in Pot", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female to Pot w/o Male", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male bites/rams", "Male Bites/Rams");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male chases", "Male Chases");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female Acceptance Still, Faces Male", "Female Follows");
	return(mariana_clean);
}

.cleanUpPGF2A = function(dat) {
	clean <- .replaceBehAll(dat, "approach", "APPROACH");
	clean <- .replaceBehAll(clean, "bite", "BITE");
	clean <- .replaceBehAll(clean, "male CHASE", "CHASE");
	clean <- .replaceBehAll(clean, "chase", "CHASE");
	clean <- .replaceBehAll(clean, "flee", "FLEE");
	clean <- .replaceBehAll(clean, "female FOLLOW", "FOLLOW");
	clean <- .replaceBehAll(clean, "inside pot", "inside POT");
	clean <- .replaceBehAll(clean, "lead", "LEAD");
	clean <- .replaceBehAll(clean, "male LEAD", "LEAD");
	clean <- .replaceBehAll(clean, "male BITES", "BITE");
	clean <- .replaceBehAll(clean, "male QUIVER", "QUIVER");
	clean <- .replaceBehAll(clean, "quiver", "QUIVER");
	clean <- .replaceBehAll(clean, "spawning", "SPAWNING");
	return(clean);
}


.getMatMatNames = function(behnames) {
	vec = character();
	for (i in 1:length(behnames)) {
			vec = c(vec, paste(behnames[i], ": Latency"), paste(behnames[i], ": Count"));
	}
	return(vec);
}


.extractBasicStats = function(data, behnames) {
	tables = lapply(data, function(f) {table(f$behavior)});
	countsMat = matrix(nrow = length(names(tables)), ncol = length(behnames), dimnames = list(names(tables), behnames));
	latMat = matrix(nrow = length(names(tables)), ncol = length(behnames), dimnames = list(names(tables), behnames));
	matmat = matrix(nrow = length(names(tables)), ncol = length(behnames) * 2, dimnames = list(names(tables), .getMatMatNames(behnames)));
	for (i in 1:length(names(tables))) {
		for (j in 1:length(behnames)) {
			if (behnames[j] %in% names(tables[[i]])) {
				countsMat[i,j] <- tables[[i]][[behnames[j]]];
				latMat[i,j] <- data[[i]][data[[i]]$behavior == behnames[j],][1,]$time;
				matmat[i, 2*j-1] <- data[[i]][data[[i]]$behavior == behnames[j],][1,]$time;
				matmat[i, 2*j] <- tables[[i]][[behnames[j]]];
			}
			else {
				countsMat[i,j] <- 0;
				matmat[i,2*j] <- 0;
			}
		}
	}
	return(list(counts = t(countsMat), latencies = t(latMat), total = t(matmat)));
	# write.csv(matmat, file = "~/Desktop/Katrina/MARIANA_spawning_analysis/basic_stats.csv")
}


#maybe flag groups in other fxn, have user give group1Names, group2Names: colname %in% group#names
# and then separate here.
.calcBasicStats = function(data, bootstrapTrials = 10000) {
	groupNames = names(table(gsub("/.+", "", names(data))));
	behnames = names(table(unlist(lapply(data, function(f) {names(table(f$behavior))}))));
	dataByGroup = list();
	for (i in 1:length(groupNames)) {
		groupData = data[grepl(paste("^", groupNames[i], sep = ""), names(data))];
		dataByGroup[[i]] = .extractBasicStats(groupData, behnames);
		names(dataByGroup)[i] <- groupNames[i];
	}
	avgCounts = lapply(dataByGroup, function(d) {apply(d$counts, 1, mean)});
	stdDevCounts = lapply(dataByGroup, function(d) {apply(d$counts, 1, sd)});
	
	wilcoxCounts = numeric();
	studentCounts = numeric();
	wilcoxLat = numeric();
	bootstrapCounts = numeric();
	# studentLat = numeric();
	for (i in 1:length(behnames)) {
		print(behnames[i])
		# print(dataByGroup[[1]]$counts[i,]);
		# print(dataByGroup[[2]]$counts[i,]);
		# print(wilcox.test(x=dataByGroup[[1]]$counts[i,], y=dataByGroup$counts[[2]][i,]));
		# print(wilcox.test(x=dataByGroup[[1]]$counts[i,], y=dataByGroup$counts[[2]][i,])$p.value);
		wilcoxCounts = c(wilcoxCounts, wilcox.test(x=dataByGroup[[1]]$counts[i,], y=dataByGroup$counts[[2]][i,])$p.value);
		studentCounts = c(studentCounts, t.test(x=dataByGroup[[1]]$counts[i,], y=dataByGroup$counts[[2]][i,])$p.value);
		wilcoxLat = c(wilcoxLat, wilcox.test(x=dataByGroup[[1]]$latencies[i,], y=dataByGroup$latencies[[2]][i,])$p.value);
		bootstrapCounts = c(bootstrapCounts,  bootstrap2independent(group1 = dataByGroup[[1]]$counts[i,],
																	group2 = dataByGroup[[2]]$counts[i,],
																	dataDescriptor = behnames[i],
																	groupNames = names(dataByGroup)[1:2],
																	printResults = F, verbose = F,
																	trials = bootstrapTrials)$p );
		# studentLat = c(studentLat, t.test(x=dataByGroup[[1]]$latencies[i,], y=dataByGroup$latencies[[2]][i,])$p.value);
		names(wilcoxCounts)[i] <- behnames[i];
		names(studentCounts)[i] <- behnames[i];
	}
	return(data.frame(averageCount = avgCounts,
					  stddevCounts = stdDevCounts,
					  wilcox_p_counts = wilcoxCounts,
					  studentt_p_counts = studentCounts,
					  bootstrap_p_counts = bootstrapCounts,
					  wilcox_p_lat = wilcoxLat));
}

#TODO double check seconds/time is OK



.makeRasterPlot = function (dataFrame, plots=T, ...) {
	data = dataFrame$behavior;
	# if (is.null(.checkInputDataVecOK(data))) {
		# return(NULL);
	# } 
	codes = names(table(dataFrame$behavior));
	frames = as.numeric(dataFrame$time);
	counts = table(data);
	num_beh = length(counts);
	diffs = diff(frames);
	
	avg_diffs = counts;
	for (beh in 1:num_beh) {
		avg_diffs[beh] = mean(diffs[data==names(avg_diffs)[beh]], na.rm=T);
	}
	
	if (plots) {
		par(mfrow=c(2, 1));
		par(oma=c(0,5,0,0));
		plot(frames, as.numeric(as.factor(data)), frame.plot=F, 
			 axes=F, xlab='time (seconds)', ylab='', 
			 col='blue', cex=3, pch=3, 
			 ...);
		axis(2, at=1:num_beh, labels=codes, tick=F, las=2);
		axis(1, yaxp=c(0, max(frames), 10), col='white', col.ticks='black');
		for (i in 1:num_beh) { 
			abline(h=i, col='darkgrey');
		}
		
		boxplot(frames ~ as.factor(data), frame.plot=F,
				col='grey', notch=F, width=table(data)/length(data), 
				horizontal=T, names=codes, las=1
				);
	}
	probMat = .getProbabilityMatrix(data);
	return(list(beh_counts=counts, 
				frame_diffs=diffs, 
				frame_diffsAvg=avg_diffs, 
				frames_per_beh=max(frames)/sum(counts),
				ethogram=.computeEntropyProbMatrix(probMat),
				data=dataFrame
				)
			);
}


