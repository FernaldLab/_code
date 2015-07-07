# count number of times a given sequence is observed in a vector of behaviors
.countBehaviorSequences = function (behaviorVec, seq) {
	seq = strsplit(seq, '')[[1]];
	beh1 = seq[1];
	toTest = paste('behaviorVec[b+1]==\"', seq[2], '\"', sep='');
	if (length(seq) > 2) {
		for (b in 3:length(seq)) {
			toTest = paste(toTest, ' && behaviorVec[b+', b-1, ']==\"', seq[b], '\"', sep='');
		}
	}
	count = 0;
	firstpos = c(0);
	for (b in 1:(length(behaviorVec)-length(seq)+1)) {
		if (behaviorVec[b] == beh1) {
			if (b < (firstpos[length(firstpos)] + length(seq)) ) {
				next;
			} else {
				if (eval(parse(text=toTest))) {
					count = count+1;
					firstpos = c(firstpos, b);
				}
			}
		}
	}
	firstpos = firstpos[-1];
	return(list(count=count,pos=firstpos));
}

# shuffle behavior vector and count occurences of a given sequence to generate a null distribution
.generateNull = function (behaviorVec, seq, runs=1000) {
	nullDist = c();
	for (r in 1:runs) {
		if (r %% 1000 == 0) {
			cat('run ', r, '\n', sep='');
		}
		pVec = sample(behaviorVec);
		nullDist = c(nullDist, .countBehaviorSequences(behaviorVec=pVec, seq=seq)$count);
	}
	return(nullDist);
}

# compute a pvalue for how often a given sequence is observed
.computeOverrepPval = function (behaviorVec, seq, runs=1000, plot=T) {
	actual = .countBehaviorSequences(behaviorVec=behaviorVec, seq=seq);
	nullDist = .generateNull(behaviorVec=behaviorVec, seq=seq, runs=runs);
	pval = sum(nullDist >= actual$count) / runs;
	if (plot) {
		hist(nullDist, col='grey', border='darkgrey', main=paste('p=',pval,sep=''), xlab='');
		abline(v=actual$count, lty='dashed');
	}
	return(list(pval=pval, nullDist=nullDist, actual=actual));
}

# generate a vector of all possible behavior sequences of a given length

### missing some, e.g. no bb ###
.getPossibleCombos = function (behaviorVec, len) {
	l = unique(combn(behaviorVec, len, simplify=F));
	l = unlist(lapply(l, function(f) paste(f, collapse='')));
	return(sort(l));
}


.getOverrepCombos = function (behaviorVec, len, runs) {
	possibleSeqs = .getPossibleCombos(behaviorVec=behaviorVec, len=len);
	l = list();
	for (i in 1:length(possibleSeqs)) {
		cat(possibleSeqs[i], ' ');
		l[[i]] = .computeOverrepPval(behaviorVec=behaviorVec, seq=possibleSeqs[i], runs=runs, plot=F);
	}
	names(l) = possibleSeqs;
	m = cbind(count=sapply(l, function(f) f$actual$count), 
			  mean_exp=sapply(l, function(f) mean(f$nullDist)),
			  pval=sapply(l, function(f) f$pval)
			  );
	return(list(m[order(m[,3], -m[,1]), ], l));
}


################################################################################

.getOverrepCombos2 = function (behaviorVec, len, runs) {
	possibleSeqs = .getPossibleCombos(behaviorVec=behaviorVec, len=len);
	pCountMat = as.data.frame(matrix(ncol=length(possibleSeqs), nrow=runs));
	names(pCountMat) = possibleSeqs;
	for (r in 1:runs) {
		pCounts = c();
		if (r %% 10 == 0){cat('Run ', r, '\n', sep='')}
		pVec = sample(behaviorVec);
		for (b in 1:length(possibleSeqs)) {
			pCounts = c(pCounts, .countBehaviorSequences(behaviorVec=pVec, seq=possibleSeqs[b])$count);
			names(pCounts)[b] = possibleSeqs[b];
		}
		pCountMat[r, ] = pCounts;
	}
	resMat = as.data.frame(matrix(ncol=3, nrow=length(possibleSeqs), 
								  dimnames=list(possibleSeqs, c('count','mean_exp','pval'))
								  )
						   );
	for (b in 1:ncol(pCountMat)) {
		#cat(names(pCountMat)[b], ' ');
		actual = .countBehaviorSequences(behaviorVec=behaviorVec, seq=names(pCountMat)[b])$count;
		mean_exp = mean(pCountMat[,b]);
		pval = sum(pCountMat[,b] >= actual) / runs;
		resMat[b, ] = c(actual, mean_exp, pval)
	}
	return(list(pvals=resMat[order(resMat[,3],-resMat[,1]), ], nullDists=pCountMat));
}