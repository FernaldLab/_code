source('~/Documents/_Fernald_lab/_code/behavior_syntax.R');
setwd('~/Documents/_Fernald_lab/matLabBoldnessRA/');
files = list.files()[grepl('clean$', list.files())];

for (f in 1:length(files)) {
	print(files[f]);
	buildGraphFromFile(filename=files[f],dotFileBase=files[f],datVersion='RA')
}; rm(f);

#####################

getDataRA = function (filename) {
	dat = as.character(read.table(filename, sep=' ', colClasses='character'));
	freqs = table(dat);
	dat_num = match(dat, names(freqs));
	return(list(rawData=dat, numData=dat_num, freqs=freqs));
}

getCounts = function (dataVec, leader, follower) {
	counts = 0;
	for (b in 1:length(dataVec)) {
		if (dataVec[b] == leader) {
			if (b == length(dataVec)) {
				break;
			} 
			if (dataVec[b+1] == follower) {
				counts = counts + 1;
			}
		}
	}
	return(counts);
}

# depends on getCounts
getCountMat = function (dataVec) {
	freqs = table(dataVec);
	counts = matrix(nrow=length(freqs), ncol=length(freqs), 
					dimnames=list(names(freqs), names(freqs))
					);
	for (leader in 1:nrow(counts)) {
		this_l = rownames(counts)[leader];
		for (follower in 1:ncol(counts)) {
			this_f = colnames(counts)[follower];
			counts[leader, follower] = getCounts(dataVec, this_l, this_f);
		}
	}
	return(counts);
}

combineMats = function (mat1, mat2) {
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
collapseCountMatList = function (countMatList) {
	# initialize with largest count matrix
	dimNums = sapply(countMatList, nrow);
	bigInd = which(dimNums == max(dimNums));
	if (length(bigInd) > 1) {
		bigInd = bigInd[1];
	}
	comboMat = countMatList[[bigInd]];
	countMatList = countMatList[-bigInd];
	# loop through reminaing matrices
	for (m in 1:length(countMatList)) {
		this_mat = countMatList[[m]];
		comboMat = combineMats(comboMat, this_mat);
	}
	return(comboMat);
}

countsToProbs = function (countMat) {
	probMat = countMat / apply(countMat, 1, sum);
	return(probMat);
}

buildDotFile = function(probMatrix, 
						#originalDataVec, 
						file='', title='untitled', fontsize=24) {
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# # get behavior frequencies
	#freqs = table(originalDataVec);
	freqs=rep(0,nrow(probMatrix));
	names(freqs) = rownames(probMatrix)
	# # check that behaviors are in same order in freqs and probMatrix
	#if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
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

#####################

nov6files = files[grepl('nov6',files)];
nov6data = list();
for (f in 1:length(nov6files)) {
	nov6data[[f]] = getDataRA(nov6files[f]);
	names(nov6data)[f] = nov6files[f];
}; rm(f);

nov6counts = list();
for (s in 1:length(nov6data)) {
	this_dat = nov6data[[s]]$rawData;
	nov6counts[[s]] = getCountMat(this_dat);
	names(nov6counts)[s] = names(nov6data)[s];
}; rm(s,this_dat);

nov6countsCombined = collapseCountMatList(nov6counts);
buildDotFile(countsToProbs(nov6countsCombined),'nov6.dot');

