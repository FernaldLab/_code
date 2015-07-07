library(Rgraphviz);# library(WGCNA);
options(stringsAsFactors=F);

.getDescriptionTableFromRawData = function (data0) {
	desc_row = which(grepl('^start.*description$', data0[,1])) + 2;
	desc_table = matrix(nrow=1, ncol=2);
	while (TRUE) {
		if (grepl('^-.*-$', data0[desc_row, ])) {
			break;
		} else {
			tmp = unlist(strsplit(data0[desc_row, ], '0'));
			code = gsub(' ','',tmp[1]);
			description = gsub('  ','',tmp[length(tmp)]);
			desc_table = rbind(desc_table, c(code, description));
			desc_row = desc_row + 1;
		}
	}
	desc_table = desc_table[-1, ];
	return(desc_table);
}

.getData = function (filename) {
	data0 = read.table(filename, fill=T, colClasses='character', sep='\t', header=F, quote='', blank.lines.skip=T, strip.white=T);
	desc_table = .getDescriptionTableFromRawData(data0);
	
	start_raw = which(data0=='RAW LOG') + 4;
	start_full = which(data0=='FULL LOG') + 4;
	end_raw = which(data0=='FULL LOG') - 2;
	end_full = which(data0=='NOTES') - 2;
	
	raw = strsplit(data0[start_raw:end_raw, ], ' ');
	raw = lapply(raw, function(f) f[c(1, length(f))]);
	df = data.frame(matrix(nrow=1, ncol=4));
	for (entry in raw) {
		if (nchar(entry)[2]==1) {
			df = rbind(df, c(as.character(entry[2]), as.numeric(entry[1]), NA, NA));
		}
	}
	df = df[-1, ];
	names(df) = c('behavior', 'start_frame', 'end_frame', 'duration');

	full = strsplit(data0[start_full:end_full, ], ' ');
	full = lapply(full, function(f) f[c(1, length(f))]);
	for (i in 1:length(full)) {
		entry = full[[i]]
		if (entry[2]=='stop') {
			this_frame = as.numeric(entry[1]);
			last_frame = as.numeric(full[[(i-1)]][1]);
			rrow = which(as.numeric(df[,2])==last_frame);
			df[rrow, 3:4] = c(as.numeric(this_frame), as.numeric(this_frame)-as.numeric(last_frame));
		}
	}
	if (nrow(df)<1) {
		warning(paste('No data in ', filename, sep=''));
		return(list(data=df, codes=NULL));
	}
	check1 = which(as.numeric(df$end_frame)<as.numeric(df$start_frame));
	check2 = which(as.numeric(df$end_frame)==as.numeric(df$duration));
	df[unique(c(check1, check2)), 3:4] = NA;
	desc_table = desc_table[as.character(desc_table[, 1]) %in% as.character(df$behavior), ];
	return(list(data=df, codes=desc_table));
}


#

.getDataBatch = function (filenames) {
	data = list();
	for (f in 1:length(filenames)) {
		print(filenames[f])
		data[[f]] = .getData(filenames[f]);
		names(data)[f] = filenames[f];
	}
	return(data);
}

#

.getProjectIdFromFilenames = function (filenames) {
	rawids = suppressWarnings(as.numeric(unlist(strsplit(filenames, '_'))));
	return(rawids[!is.na(rawids)]);
}

#

.countTransitionsInVec = function (dataVec, leader, follower) {
	total = 0;
	for (i in 1:length(dataVec)) {
		if (dataVec[i] == leader) {
			if (i == length(dataVec)) {
				break;
			} else if (dataVec[i+1] == follower) {
				total = total + 1;
			}
		}
	}
	return(total);
}

.checkForEndSingletonAndRemove = function (dataVec) {
	counts = table(dataVec);
	for (i in 1:length(counts)) {
		if (counts[i] == 1) {
			if (names(counts)[i] == dataVec[length(dataVec)]) {
				dataVec = dataVec[-length(dataVec)];
			}
		}
	}
	return(dataVec);
}

.checkInputDataVecOK = function (dataVec, thresh1=3, thresh2=10) {
	if (is.null(dataVec)) {
		warning('Input vector was empty, returning NULL')
		return(NULL);
	}
	if (length(dataVec) < thresh1) {
		warning('Input vector contains fewer than 3 entries, returning NULL');
		return(NULL);
	}
	if (length(dataVec) <= thresh2) {
		warning('Input vector contains 10 or fewer entries, results may not be meaningful');
		return(TRUE);
	}
	return(TRUE);
}

.getTransitionCountMatrixFromVec = function (dataVec, removeEndSingletons=T) {
	if (is.null(.checkInputDataVecOK(dataVec))) {
		return(NULL);
	} 
	behs = table(dataVec);
	if (removeEndSingletons) {
		dataVec = .checkForEndSingletonAndRemove(dataVec);
		behs = table(dataVec);
	}
	cMat = matrix(nrow=length(behs), ncol=length(behs), dimnames=list(names(behs), names(behs)));
	for (leader in 1:nrow(cMat)) {
		for (follower in 1:ncol(cMat)) {
			cMat[leader, follower] = .countTransitionsInVec(dataVec, rownames(cMat)[leader], colnames(cMat)[follower]);
		}
	}
	return(cMat);
}

.getTransitionCountMatrixFromList = function (getDataBatchOutputList, removeEndSingletons=T) {
	dat = getDataBatchOutputList;
	cList = list();
	for (f in 1:length(dat)) {
		cList[[f]] = .getTransitionCountMatrixFromVec(dat[[f]]$data$behavior, removeEndSingletons=removeEndSingletons);
	}
	names(cList) = names(dat);
	return(cList);
}

#

.getBehaviorDurationStats = function (dataMatFrom.getData) {
	data = dataMatFrom.getData;
	behaviors = names(table(data$behavior));
	dMat = matrix(nrow=length(behaviors), ncol=5, dimnames=list(behaviors, c('min','max','mean','median','total')));
	for (b in behaviors) {
		d = as.numeric(data$duration[data$behavior==b]);
		stats = summary(d[!is.na(d)]);
		dMat[match(b, rownames(dMat)), ] = c(stats[1], stats[6], stats[4], stats[3], sum(d[!is.na(d)]));
	}
	return(dMat);
}

.getBehaviorDurationMatrixFromFishList = function (getDataBatchOutput, stat='total', NAto0=T) {
	if (!all(stat %in% c('min','max','mean','median','total'))) {
		stop('Invalid stat argument, must be one of: min, max, mean, median, total');
	}
	dat = getDataBatchOutput;
	codes = .combineBehaviorCodesFromFishList(dat);
	dMat = matrix(nrow=nrow(codes), ncol=length(dat), dimnames=list(codes[, 1], names(dat)));
	for (i in 1:length(dat)) {
		f = dat[[i]];
		if (nrow(f$data)<1) {
			next;
		} else {
			stats = .getBehaviorDurationStats(f$data);
			for (beh in 1:nrow(stats)) {
				row = match(rownames(stats)[beh], rownames(dMat));
				col = match(names(dat)[i], colnames(dMat));
				dMat[row, col] = stats[beh, match(stat, colnames(stats))];
			}
		}
	}
	if (NAto0) {
		dMat[is.na(dMat)] = 0;
	}
	return(dMat);
}

#

.dataSummary = function (getDataOutput, plots=T, ...) {
	data = getDataOutput$data$behavior;
	if (is.null(.checkInputDataVecOK(data))) {
		return(NULL);
	} 
	codes = getDataOutput$codes;
	frames = as.numeric(getDataOutput$data$start_frame);
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
			 axes=F, xlab='time (frame #)', ylab='', 
			 col='blue', cex=3, pch=3, 
			 ...);
		axis(2, at=1:num_beh, labels=codes[,2], tick=F, las=2);
		axis(1, yaxp=c(0, max(frames), 10), col='white', col.ticks='black');
		for (i in 1:num_beh) { 
			abline(h=i, col='darkgrey');
		}
		
		boxplot(frames ~ as.factor(data), frame.plot=F,
				col='grey', notch=F, width=table(data)/length(data), 
				horizontal=T, names=codes[,2], las=1
				);
	}
	probMat = .getProbabilityMatrix(data);
	return(list(beh_counts=counts, 
				frame_diffs=diffs, 
				frame_diffsAvg=avg_diffs, 
				frames_per_beh=max(frames)/sum(counts),
				ethogram=.computeEntropyProbMatrix(probMat),
				data=getDataOutput
				)
			);
}

#

.combineBehaviorCodesFromFishList = function (getDataBatchOutput) {
	dat = getDataBatchOutput;
	for (f in dat) {
		if (!is.null(f$codes)) {
			codeMat = f$codes;
			break;
		}
	}; 
	for (f in dat) {
		if (is.null(f$codes)) {next}
		if (is.vector(f$codes)) {
			if (f$codes[1] %in% codeMat[, 1]) {
				next;
			} else {
				codeMat = rbind(codeMat, f$codes);
			}
		} else if (is.matrix(f$codes)) {
			if (all(f$codes[, 1] %in% codeMat[, 1])) {
				next;
			} else {
				codeMat = rbind(codeMat, f$codes[!(f$codes[, 1] %in% codeMat[, 1]), ]);
			}
		} else {
			stop('error');
		}
	}
	return(codeMat[order(codeMat[, 1]), ]);
}

#

.getBehaviorCountMatrixFromFishList = function (getDataBatchOutput, percent=T, clean=T) {
	dat = getDataBatchOutput;
	codes = .combineBehaviorCodesFromFishList(dat);
	cMat = matrix(nrow=nrow(codes), ncol=length(dat), dimnames=list(codes[, 1], names(dat)));
	for (i in 1:length(dat)) {
		f = dat[[i]];
		if (nrow(f$data)<1) {
			next;
		} else {
			counts = table(f$data$behavior);
			for (beh in 1:length(names(counts))) {
				row = match(names(counts)[beh], rownames(cMat));
				col = match(names(dat)[i], colnames(cMat));
				cMat[row, col] = counts[beh];
			}
		}
	}
	if (clean) {
		cat('removing fish with no behaviors and replacing NAs with 0s\n');
		checkNA = apply(cMat, 2, function(ff) sum(is.na(ff))) == nrow(cMat);
		cMat = cMat[, !checkNA];
		cMat[is.na(cMat)] = 0;
	}
	if (percent) {
		pMat = cMat;
		colTotals = apply(pMat, 2, sum, na.rm=T);
		for (f in 1:length(colTotals)) {
			pMat[, f] = pMat[, f] / colTotals[f];
		}
		out = list(cMat=cMat, pMat=pMat);
	} else {
		out = list(cMat=cMat, pMat=NULL);
	}
	return(out);
}

#

.computeTransitionProbability = function (data, leader, follower) {
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
	total_leader = sum(data == leader);
	prob = count / total_leader;
	return(list(probability=prob, termination=termination, count_transitions=count, count_leader=total_leader));
}

#

.getProbabilityMatrix = function (data, removeZeroCol=F) {
	beh = names(table(data));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeTransitionProbability(data=data, leader=leader, follower=follower);
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

#

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

#

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

#

.analyzeBehavior = function (filename, probMat=NULL) {
	if (!is.null(filename) & is.null(probMat)) {
		data = .getData(filename);
		probMat = .getProbabilityMatrix(data$data$behavior);
	} else if (!is.null(probMat) & is.null(filename)) {
		probMat = probMat;
	} 
	h = .computeEntropyProbMatrix(probMat);
	return(list(analysis=h, data=data));
}















# node_labels = data$codes[, 1];
# g = graphNEL(nodes=node_labels, edgeL=list(), edgemode='directed');
# for (row in 1:nrow(probMat)) {
	# for (col in 1:ncol(probMat)) {
		# g = addEdge(node_labels[row],node_labels[col],g,probMat[row,col])
	# }
# }
# g=layoutGraph(g,recipEdges='distinct')
# w=as.vector(t(probMat))
# names(w) = names(edgeRenderInfo(g)$direction)
# edgeRenderInfo(g) = list(lwd=w*5)
# nodesize = nodeRenderInfo(g)$height*2;
# nodeRenderInfo(g) = list(fixedsize=F, rWidth=nodesize, lWidth=nodesize, fontsize=5)
# #g=layoutGraph(g,recipEdges='distinct')
# renderGraph(g);



.corHeatmapWithPvalLabels = function(data, mar=c(8,15,3,3), col=blueWhiteRed(50), main='', cex.text=.5, thresh=.05/(ncol(data)^2/2-(ncol(data)/2)), ...) {
	tmp = corAndPvalue(data, ...);
	cors = tmp$cor; 
	pvals = tmp$p;
	textMat = paste(signif(cors, 2), '\n(', signif(pvals, 1), ')', sep='');
	dim(textMat) = dim(cors);
	textMat[pvals>thresh] = '';
	par(mar=mar);
	labeledHeatmap(Matrix=cors, xLabels=names(data), yLabels=names(data), ySymbols=names(data), colorLabels=F, colors=col, textMatrix=textMat, setStdMargins=F, cex.text=cex.text, zlim=c(-1, 1), main=main, ...)
}

.combineMats = function (mat1, mat2) {
	all_beh = sort(unique(c(rownames(mat1), rownames(mat2))));
	if (nrow(mat1) > 1) {
		mat1 = mat1[order(rownames(mat1)), ];
	}
	if (nrow(mat2) > 1) {
		mat2 = mat2[order(rownames(mat2)), ];
	}
	
	if (suppressWarnings(all(rownames(mat1) == all_beh))) {
		newMat = mat1;
		smallMat = mat2;
	} else if (suppressWarnings(all(rownames(mat2) == all_beh))) {
		newMat = mat2;
		smallMat = mat1;
	} else {
		newMat = matrix(rep(0, length(all_beh)^2),
						ncol=length(all_beh),
						dimnames=list(all_beh, all_beh)
						);
		smallMat = NULL;
	}
	
	if (!is.null(smallMat)) {
		cols = match(colnames(smallMat), colnames(newMat));
		for (row in 1:nrow(smallMat)) {
			rN = match(rownames(smallMat)[row], rownames(newMat));
			newMat[rN, cols] = newMat[rN, cols] + smallMat[row, ];
		}
	} else {
		rc1 = match(rownames(mat1), rownames(newMat));
		rc2 = match(rownames(mat2), rownames(newMat));
		newMat[rc1, rc1] = newMat[rc1, rc1] + mat1;
		newMat[rc2, rc2] = newMat[rc2, rc2] + mat2;
	}

	return(newMat);
}

# depends on combineMats
.collapseCountMatList = function (countMatList) {
	# initialize with largest count matrix
	dimNums = sapply(countMatList, nrow);
	checkNULL = lapply(dimNums, is.null);
	dimNums[unlist(lapply(checkNULL, which))] = 1;
	bigInd = which(unlist(dimNums) == max(unlist(dimNums)));
	if (length(bigInd) > 1) {
		bigInd = bigInd[1];
	}
	comboMat = countMatList[[bigInd]];
	countMatList = countMatList[-bigInd];
	# loop through reminaing matrices
	for (m in 1:length(countMatList)) {
		this_mat = countMatList[[m]];
		if (is.null(this_mat)) { next };
		comboMat = .combineMats(comboMat, this_mat);
	}
	return(comboMat);
}

# countsToProbs = function (countMat) {
	# probMat = countMat / apply(countMat, 1, sum);
	# if (any(is.nan(probMat))) {
		# probMat[is.nan(probMat)] = 0;
		# warning('Changed NaNs to 0s, check probMatrix')
	# }
	# return(probMat);
# }

# buildDotFile = function(probMatrix, countMat,
						# #originalDataVec, 
						# file='', title='untitled', fontsize=24) {
	# # write top line to file
	# cat('digraph', title, '\n', '	{\n', file=file);
		
	# # # get behavior frequencies
	# freqs = apply(countMat, 1, sum);
	# #freqs=rep(0,nrow(probMatrix));
	# #names(freqs) = rownames(probMatrix); print(freqs)
	# # # check that behaviors are in same order in freqs and probMatrix
	# #if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
	# # # compute proportions of behaviors for relative node size
	# for (beh in 1:length(freqs))
	# {
		# prop = freqs[beh] / sum(freqs) * 10; #print(prop)
		# cat('		', names(freqs)[beh], ' [width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
	# }
	
	# # loop through probMatrix to get probabilities
	# probMatrix=probMatrix*10;
	# for (row in 1:nrow(probMatrix))
	# {
		# for (col in 1:ncol(probMatrix))
		# {
			# val = probMatrix[row,col];
			# if (val > 0)
			# {
				# leader = rownames(probMatrix)[row];#print(leader);
				# follower = colnames(probMatrix)[col];#print(follower);
		
				# cat('		', leader, ' -> ', follower, ' [label="", style="setlinewidth(', val, ')", arrowsize=1];','\n' ,sep='', file=file, append=T);	
		 	# }
		# }
	# }
	
	# # write last line of file
	# cat('	}', file=file, append=T);
	
	# return(NULL)
# }

# #####################

# buildComboGraphFromFiles = function(files, filebase='x') {
	# data = list();
	# for (f in 1:length(files)) {
		# data[[f]] = getDataRA(files[f]);
		# names(data)[f] = files[f];
	# }; rm(f);

	# counts = list();
	# for (s in 1:length(data)) {
		# this_dat = data[[s]]$rawData;
		# counts[[s]] = getCountMat(this_dat);
		# names(counts)[s] = names(data)[s];
	# }; rm(s,this_dat);

	# countsCombined = collapseCountMatList(counts);#print(countsCombined);print(countsToProbs(countsCombined))
	# buildDotFile(countsToProbs(countsCombined),countsCombined,paste(filebase, '.dot', sep=''));
# }


#########################
# currently requires WGCNA::verboseBoxplot

.verboseBoxplotFromMatrixRow = function (mat, row, factor) {
	verboseBoxplot(mat[row, ], as.factor(factor),
				   xlab='',
				   )
}

# .boxplotRowsCompareColumnSubsets = function (mat, factor, ) {
	
	
# }

# .boxplotColSumsWithFactors = function (mat, factors) {
	# toPlot = apply(mat, 2, sum);
	
# }