# Source: ethograms_from_scorevideo.R
# Reads in a score log and returns a list with:
#    1. $data, a data frame with columns behavior, start frame, (end frame), (duration)
#    2. $codes, a matrix with column 1 containing the codes, column 2 containing the meanings, and
#          column 3 containing the subjects.
.getData = function (filename) {
	data0 = read.table(filename, fill=T, colClasses='character', sep='\t', header=F, quote='', blank.lines.skip=T, strip.white=T);
	desc_table = .getDescriptionTableFromRawData(data0);
	
	df = .parseRawLog(data0, desc_table);
	df = .parseFullLog(data0, desc_table, df);
	.checkStartEndMatch(df);
	
	
	df = df[order(as.numeric(df[,2])),];
	dimnames(df)[[1]] <- 1:(dim(df)[1]);
	
	if (nrow(df)<1) {
		warning(paste('No data in ', filename, sep=''));
		return(list(data=df, codes=NULL));
	}
	check1 = which(as.numeric(df$end_frame)<as.numeric(df$start_frame));
	check2 = which(as.numeric(df$end_frame)==as.numeric(df$duration));
	df[unique(c(check1, check2)), 3:4] = NA;
#	desc_table = desc_table[as.character(desc_table[, 1]) %in% as.character(df$behavior), ];
	return(list(data=df, codes=desc_table));
}



# Helper function for .getData
# Makes preliminary data frame using raw log
.parseRawLog = function (data0, desc_table) {
	startcodes = desc_table[,1]
	
	start_raw = which(data0=='RAW LOG') + 4;
	end_raw = which(data0=='FULL LOG') - 2;
	
	raw = strsplit(data0[start_raw:end_raw, ], ' ');
	raw = lapply(raw, function(f) f[c(1, length(f))]);
	df = data.frame(matrix(nrow=1, ncol=4));
	for (entry in raw) {
		if (nchar(entry)[2]==1 && entry[2] %in% startcodes) {
			df = rbind(df, c(as.character(entry[2]), as.numeric(entry[1]), NA, NA));
		}
	}
	df = df[-1, ];
	names(df) = c('behavior', 'start_frame', 'end_frame', 'duration');
	return(df);
}

# Helper function for .getData
# Fills in end times of behaviors using full log.
.parseFullLog = function (data0, desc_table, df) {
	start_full = which(data0=='FULL LOG') + 4;
	end_full = which(data0=='NOTES') - 2;

	full = strsplit(data0[start_full:end_full, ], '    ');
	full = lapply(full, function(f) c(f[1], gsub('^ *','', gsub(' *$','', f[3])), gsub(' ', '', f[length(f)])));
	for (i in 1:length(full)) {
		entry = full[[i]];
		if (entry[3]=='stop') {
			this_frame = as.numeric(entry[1]);
			startIndex = 0;
			for (j in 1:i) {
				if (full[[j]][2] == entry[2] && full[[j]][3]=='start') {  ####[[3]][315:365,]
					startIndex = j;
				}
			}
			last_frame = as.numeric(full[[startIndex]][1]);
			rrow = which(as.numeric(df[,2])==last_frame & df[,1] == desc_table[desc_table[,2] == entry[2],1]);
			if(length(rrow) > 1) {
				rrow = rrow[1];
				warning(paste("There are two of the same behavior (", desc_table[desc_table[,2] == entry[2],1], ") in the same frame (", last_frame,")", sep = ''), immediate. = TRUE)
			}
			df[rrow, 3:4] = c(as.numeric(this_frame), as.numeric(this_frame)-as.numeric(last_frame));
			df = rbind(df, c(paste("end_", df[rrow, 1], sep = ""), this_frame, NA, NA));
		}
	}
	return(df);
}

# Helper function for .getData
# Ensures there are no ends of behaviors that haven't started yet
# Throws a warning if input was malformed
.checkStartEndMatch = function (df) {
	dimnames(df)[[1]] <- 1:(dim(df)[1]);
	behvec <- df$behavior
	names(behvec) <- dimnames(df)[[1]]
	freqtable <- table(behvec);
	for (end in 1:length(names(freqtable))) {
		if (grepl('^end_', names(freqtable)[end])) {
			start = which(gsub('^end_', '', names(freqtable)[end]) == names(freqtable) & names(freqtable) != names(freqtable)[end]);
			
		}
	}
	
	
	
	
	
	
	
	# freqtable <- table(df$behavior);
	# for (end in 1:length(names(freqtable))) {
		# if (grepl('^end_', names(freqtable)[end])) {
			# start = which(gsub('^end_', '', names(freqtable)[end]) == names(freqtable) & names(freqtable) != names(freqtable)[end]);
			# if (freqtable[start] != freqtable[end]) {
				# # differenceIsEnd = FALSE;
				# # if (freqtable[start] == freqtable[end] + 1) {
					# # behvec = df$behavior
					# # iLastStart = length(behvec);
					# # while (behvec[iLastStart] != start && iLastStart > 0) {iLastStart = iLastStart - 1}
					# # for (i in iLastStart:length(behvec)) {
						# # if (behvec[i] == end) break;
						# # if (i == length(behvec)) differenceIsEnd = TRUE;
					# # }
				# # }
				
				# # if (differenceIsEnd) {
					# warning(paste('THERE ARE ', freqtable[start], ' OCCURANCES OF ', names(freqtable[start]), ' AND ',
								# freqtable[end], ' OCCURANCES OF ', names(freqtable[end]), '!', sep = ''), immediate. = TRUE)
				# # }
			# }
		# }
	# }
}



# NOT_REIMPLEMENTED
# Reads in a data list (given by .getData) and returns a data list with all characters uppercase.
# The purpose of this is to combine identical behaviors performed to different fish. To keep a letter
#  from being combined, put the uppercase version in the vector "charsToKeepSeparate".
.mergeUppercaseLowercase = function (dataList, charsToKeepSeparate = "i") {
	for (i in 1:dim(dataList$data)[1]) {
		if (!(dataList$data$behavior[i] %in% charsToKeepSeparate)) {
			dataList$data$behavior[i] <- toupper(dataList$data$behavior[i]);
		}
	}
	for (i in 1:dim(dataList$codes)[1]) {
		if (!(dataList$codes[i,1] %in% charsToKeepSeparate)) {
			dataList$codes[i,1] <- toupper(dataList$codes[i,1]);
		}
	}
	
	behaviorVec <- str_extract(dataList$codes[,2], "[^(]*");
	subjectVec <- str_extract(dataList$codes[,2], "at [a-zA-Z]*");
	for (i in dim(dataList$codes)[1]:1) {
		if (i < dim(dataList$codes)[1] && dataList$codes[i,1] %in% dataList$codes[(i+1):dim(dataList$codes)[1],1]) {
			dataList$codes <- dataList$codes[-i,];
		} else {
			dataList$codes[i,2] = if (is.na(subjectVec[i])) {behaviorVec[i]}
			else {paste(behaviorVec[i], subjectVec[i], sep = "")}
		}
	}
	return(dataList);
}


.groupLevelPs = function(data, byTotal = FALSE, leader, follower) {
	probMats <- lapply(data, function(d) {.getProbabilityMatrix(d$behavior, byTotal=byTotal)});
	counts <- lapply(data, function(d) {table(d$behavior)})
	behnames = names(table(unlist(lapply(probMats, function(d){dimnames(d)[[1]]}))));
	groupNames = names(table(gsub("/.+", "", names(data))));
	probMatsByGroup = list();
	for (i in 1:length(groupNames)) {
		groupData = list(probMats = probMats[grepl(paste("^", groupNames[i], sep = ""), names(probMats))],
		                 counts = counts[grepl(paste("^", groupNames[i], sep = ""), names(counts))]);
		probMatsByGroup[[i]] = .TP_VEC(groupData, behnames, byTotal, leader, follower);
		names(probMatsByGroup)[i] <- groupNames[i];
	}
	return(probMatsByGroup);
}

.TP_VEC = function(data, behnames = NULL, byTotal, leader, follower) {
	if (is.null(behnames)) {behnames = names(table(unlist(lapply(probmas, function(d){dimnames(d)[[1]]}))));}
	probmas = data$probMats;
	counts = data$counts;
	probMat = matrix(data = numeric(length(behnames) * length(behnames)),
					 nrow = length(behnames), ncol = length(behnames), dimnames = list(behnames, behnames));
	probVec = numeric();
	countVec = numeric(length(behnames));
	names(countVec) <- behnames;
	nfish <- countVec;
	for (row in rownames(probMat)) {
		for (subject in 1:length(probmas)) {
			if (row %in% rownames(probmas[[subject]]) && row == leader) {
				countVec[row] <- countVec[row] + counts[[subject]][row];
				nfish[row] <- nfish[row] + 1;
				for (col in colnames(probmas[[subject]])) {
					if (col == follower) {probVec = c(probVec, probmas[[subject]][row,col]);}
					else {probVec = c(probVec, 0);}
				}
			}
		}
	}
	
	if (byTotal) {
		probMat = probMat / length(probmas);
	} else {
		for (col in colnames(probMat)) {
			probMat[,col] <- probMat[,col] / nfish;
		}
	}
	return(probVec);
}