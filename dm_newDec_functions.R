# helper functions

.getRowsAboveCountThresh = function (dat, thresh=20, cols=c(6,7,9,10)) {
	goodRows = vector(length=nrow(dat));
	for (row in 1:nrow(dat)) {
		c1 = sum(as.numeric(dat[row, cols[1]:cols[2]]));
		c2 = sum(as.numeric(dat[row, cols[3]:cols[4]]));
		if (c1>thresh && c2>thresh) {
			goodRows[row] = TRUE;
		} else {
			goodRows[row] = FALSE;
		}
	}
	return(goodRows)
}

.getRowsAboveCountThreshSingle = function (dat, thresh=20, cols=c(6,7)) {
	goodRows = vector(length=nrow(dat));
	for (row in 1:nrow(dat)) {
		c1 = sum(as.numeric(dat[row, cols[1]:cols[2]]));
		if (c1>thresh) {
			goodRows[row] = TRUE;
		} else {
			goodRows[row] = FALSE;
		}
	}
	return(goodRows)
}

.addFisherPvals = function (dat, cols=c(6,7,9,10)) {
	p = apply(dat, 1, function(f) fisher.test(matrix(as.numeric(c(f[cols[1]:cols[2]], f[cols[3]:cols[4]])), ncol=2))$p.value);
	return(as.data.frame(cbind(dat, p=p)));
}

.getComplementNuc = function (nuc) {
	nucs = c('A','C','G','T');
	mat = matrix(c(nucs, 'N', nucs[length(nucs):1], 'N'), ncol=2);
	if (!(nuc %in% mat[,1])) {
		stop('invalid nucleotide')
	} else {
		return(mat[mat[,1]==nuc,2]);
	}
}

.getReverseComplementSeq = function (seq) {
	seqvec = strsplit(seq,'')[[1]];
	if (all(seqvec %in% c('A','C','G','T','N'))) {
		for (n in 1:length(seqvec)) {
			seqvec[n] = .getComplementNuc(seqvec[n]);
		}
		return(paste(seqvec[length(seqvec):1],collapse=''));
	} else {
		stop('invalid nucleotide')
	}
}

.addReverseComplementColumn = function (dat, nucsCol=4, strandCol=3) {
	dat = as.data.frame(cbind(dat, nucsRC=rep(NA,nrow(dat))));#print(head(dat))
	for (r in 1:nrow(dat)) {
		if (dat[r, strandCol]=='-') {#print(r)
			dat$nucsRC[r] = .getReverseComplementSeq(dat[r,nucsCol]);
		}
	}
	return(dat);
}

.addMerCol = function (dat, strandCol=3, nucsCol=4, rcCol=12, l=2, colname=NULL) {
	dat = as.data.frame(cbind(dat, rep(NA,nrow(dat))));
	if (is.null(colname)) {
		names(dat)[ncol(dat)] = paste(l, 'mer', sep='');
	} else {
		names(dat)[ncol(dat)] = colname;
	}
	for (row in 1:nrow(dat)) {
		if (dat[row,strandCol]=="+") {
			dat[row,ncol(dat)]=substr(dat[row,nucsCol],3,3+l-1);
		}
		if (dat[row,strandCol]=="-") {
			dat[row,ncol(dat)]=substr(dat[row,rcCol],3,3+l-1);
		}
	}
	return(dat);
}

.flipNonGtoH = function (seq, checkPos) {
	seqvec = strsplit(seq,'')[[1]];
	if (seqvec[checkPos] != 'G') {
		seqvec[checkPos] = 'H';
		return(as.character(paste(seqvec, collapse='')));
	} else {
		return(seq);
	}
}

# # # .getHigherValue = function (dat, cols=c(11,12)) {
	# # # vals = c();
	# # # for (row in 1:nrow(dat)) {
		# # # v1 = as.numeric(dat[row, cols[1]]);
		# # # v2 = as.numeric(dat[row, cols[2]]);
		# # # if (v1 > v2) {
			# # # vals = c(vals, v1);
			# # # names(vals)[row] = colnames(dat)[cols[1]];
		# # # } else {
			# # # vals = c(vals, v2);
			# # # names(vals)[row] = colnames(dat)[cols[2]];
		# # # } 
	# # # }
	# # # return(vals);
# # # }

.makeRaggedMatrix = function (listOfVecs) {
	numvecs = length(listOfVecs);
	maxveclength = max(sapply(listOfVecs,length));
	mat = matrix(rep(NA, (maxveclength*numvecs)), ncol=numvecs);
	for (v in 1:numvecs) {
		thisvec = listOfVecs[[v]];
		mat[1:length(thisvec), v] = thisvec;
	}
	return(mat);
}

.barplotStackedFromVecs = function (listOfVecs, cex.text=1, labOffset=.2, ...) {
	mat = .makeRaggedMatrix(listOfVecs);
	mids = barplot(mat, beside=F, ...);
	hvec = listOfVecs;
	for (v in 1:length(hvec)) {
		vec = hvec[[v]];
		text(mids[v]+labOffset, vec[1] / 2, 
			 paste(names(vec)[1],' - ',round(vec[1],2)*100,'%',sep=''), 
			 cex=cex.text);
		for (i in 2:length(vec)) {
			h = sum(vec[1:i]);
			d = h - sum(vec[1:(i-1)]);
			#text(mids[v]+labOffset, h-(d/2), names(vec)[i], cex=cex.text);
			text(mids[v]+labOffset, h-(d/2), 
				 paste(names(vec)[i],' - ',round(vec[i],2)*100,'%',sep=''),
				 cex=cex.text);
		}
	}
}

# .checkForDuplicatedRowsAcrossColumns = function (mat, cols) {
	
# }

.addStreamColumn = function (dat, strandCol, distCol) {
	dat = cbind(dat, stream=rep(NA,nrow(dat)));
	for (r in 1:nrow(dat)) {
		dist = as.numeric(dat[r,distCol]);
		if (dist == 0) {
			dat$stream[r] = 'hit';
		} else if (dist > 0) {
			if (dat[r,strandCol]=='+') {
				dat$stream[r] = 'up';
			} else {
				dat$stream[r] = 'down';
			}
		} else if (dist < 0) {
			if (dat[r,strandCol]=='-') {
				dat$stream[r] = 'up';
			} else {
				dat$stream[r] = 'down';
			}
		}
	}
	return(dat);
}


# .annotateDMwithFeatures = function (datC, ov) {
	# d = datC;
	# for (o in 1:length(ov)) {
		# d = as.data.frame(cbind(d, rep(NA, nrow(d))));
		# names(d)[ncol(d)] = names(ov)[o];
	# }
	# return(d)
# }

.seqLogoFromSeqVec = function (seq, ...) {
	if (length(unique(nchar(seq))) > 1) {
		stop('all seqs must be same length');
	}
	m = matrix(rep(0,4*nchar(seq[1])), ncol=nchar(seq[1]));
	nt = c('A','C','G','T');
	for (j in 1:ncol(m)) {
		for (i in 1:4) {
			m[i,j] = sum(substr(seq,j,j)==nt[i]) / length(seq);
		}
	}
	seqLogo(makePWM(m), ...);
	return(m);
}




.addClosestHitColumn = function (dat, distCols) {
	#dat = as.data.frame(cbind(dat, closest=rep(NA,nrow(dat)), closestdist=rep(0,nrow(dat))));
	dt = dat[, distCols];
	maxhits = max(apply(dt,1,function(f) sum(!is.na(f))));
	rmat = as.data.frame(matrix(nrow=nrow(dat),ncol=(maxhits*2)));
	for (r in 1:nrow(dt)) {
		good = !is.na(dt[r,]);
		tmp = dt[r, good];#print(tmp)
		if (length(tmp)==1) {
			rmat[r,1:2] = c(gsub('dist','',names(dt)[good]), abs(as.numeric(tmp)));
		} else {
			inds = order(abs(as.numeric(tmp)));
			rmat[r, seq(1,length(tmp)*2,2)] = gsub('dist','',names(tmp)[inds]);
			rmat[r, seq(2,length(tmp)*2,2)] = sort(abs(as.numeric(tmp)));
		}		
		#mindist = min(abs(as.numeric(dt[r, ])), na.rm=T); print(mindist)
		#dat$closest[r] = paste(names(dt)[which(abs(as.numeric(dt[r, ]))==mindist)],collapse=':');
		#print(names(dt)[which(abs(as.numeric(dt[r, ]))==mindist)])
		#dat$closestdist[r] = mindist;
		#print(rmat)
	}
	names(rmat)[1:2] = c('closest','closestdist');
	return(as.data.frame(cbind(dat, rmat)));
}         

.getDistanceToUpstreamAndHits = function (dat, streamCols, distCols) {
	stream = dat[, streamCols];
	dists = dat[, distCols];
	names(stream) = gsub('stream','',names(stream));
	names(dists) = gsub('dist','',names(dists));
	if (any(names(stream) != names(dists))) {
		stop('check streamCols and distCols');
	}
	res = list();
	for (i in 1:ncol(stream)) {
		ind = !is.na(stream[, i]);
		res[[names(stream)[i]]] = data.frame(type=stream[ind, i], dist=abs(as.numeric(dists[ind, i])));
	}
	return(res);
}


.percentHits = function (dat, distCols, interval) {
	tdat = dat[, distCols];
	res = vector(length=length(ncol(tdat)));
	for (i in 1:ncol(tdat)) {
		res[i] = sum(abs(as.numeric(tdat[,i])) %in% interval[1]:interval[2]) / nrow(tdat);
	}
	return(res);
}

.generateIntervals = function (start=0, stop=100, intervalSize=10) {
	
}









.slideWindowSum = function (dat, interval=c(0,5000), window=10, step=1, ...) {
	total = interval[2];
	spots = seq(0, total-window, step);
	result = vector(length=length(spots));
	for (i in 1:length(spots)) {
		result[i] = sum(dat %in% spots[i]:(spots[i]+window), ...);
	}
	return(result);
}

.plotWindowSums = function (dlist, ylim=NULL, normBy=137, ...) {
	linecols = colors(T); 
	linecols = linecols[!grepl('black|white|light',linecols)];
	if (is.null(ylim)) {
		ylim = c(0,1);
	} else {
		ylim = ylim;
	}
	
	toplot = .slideWindowSum(dlist[[1]], ...) / normBy;
	plot(0:(length(toplot)-1), toplot, type='l', ylim=ylim);
	for (i in dlist[-1]) {
		par(new=T);
		toplot = .slideWindowSum(i, ...) / normBy;
		thiscol = linecols[sample(1:length(linecols), 1, F)];
		plot(0:(length(toplot)-1), toplot, type='l', ylim=ylim, col=thiscol, axes=F);
	}
}
