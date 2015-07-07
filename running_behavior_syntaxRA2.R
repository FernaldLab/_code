## first run python script in terminal to prep raw files:
# ls | while read stdin; do echo $stdin; /Volumes/fishstudies/_scripts/cleanBehaviorLogs.py $stdin; done

#source('~/Documents/_Fernald_lab/_code/behavior_syntax.R');
source('/Volumes/fishstudies/_code/behavior_syntax.R');
#setwd('~/Documents/_Fernald_lab/matLabBoldnessRA2/');
setwd('/Volumes/fishstudies/_behaviorRA/Nov2013Behavior');
files = list.files()[grepl('clean$', list.files())];

for (f in 1:length(files)) {
	print(files[f]);
	buildGraphFromFile(filename=files[f],dotFileBase=files[f],datVersion='RA')
}; rm(f);

apr23files = files[grepl('Apr23',files)];
apr23filesD = files[grepl('Apr23.*_D',files)];
apr23filesS = files[grepl('Apr23.*_S',files)];
buildComboGraphFromFiles(apr23files, 'apr23')
buildComboGraphFromFiles(apr23filesD, 'apr23_D')
buildComboGraphFromFiles(apr23filesS, 'apr23_S')

apr24files = files[grepl('Apr24',files)];
apr24filesD = files[grepl('Apr24.*_D',files)];
apr24filesS = files[grepl('Apr24.*_S',files)];
buildComboGraphFromFiles(apr24files, 'apr24')
buildComboGraphFromFiles(apr24filesD, 'apr24_D')
buildComboGraphFromFiles(apr24filesS, 'apr24_S')


filesD=files[grepl('_D_',files)]
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
	if (any(is.nan(probMat))) {
		probMat[is.nan(probMat)] = 0;
		warning('Changed NaNs to 0s, check probMatrix')
	}
	return(probMat);
}

buildDotFile = function(probMatrix, countMat,
						#originalDataVec, 
						file='', title='untitled', fontsize=24) {
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# # get behavior frequencies
	freqs = apply(countMat, 1, sum);
	#freqs=rep(0,nrow(probMatrix));
	#names(freqs) = rownames(probMatrix); print(freqs)
	# # check that behaviors are in same order in freqs and probMatrix
	#if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
	# # compute proportions of behaviors for relative node size
	for (beh in 1:length(freqs))
	{
		prop = freqs[beh] / sum(freqs) * 10; #print(prop)
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

buildComboGraphFromFiles = function(files, filebase='x') {
	data = list();
	for (f in 1:length(files)) {
		data[[f]] = getDataRA(files[f]);
		names(data)[f] = files[f];
	}; rm(f);

	counts = list();
	for (s in 1:length(data)) {
		this_dat = data[[s]]$rawData;
		counts[[s]] = getCountMat(this_dat);
		names(counts)[s] = names(data)[s];
	}; rm(s,this_dat);

	countsCombined = collapseCountMatList(counts);#print(countsCombined);print(countsToProbs(countsCombined))
	buildDotFile(countsToProbs(countsCombined),countsCombined,paste(filebase, '.dot', sep=''));
}

nov6files = files[grepl('nov6',files,ignore.case=T)];
nov7files = files[grepl('nov7',files,ignore.case=T)];
nov25files = files[grepl('nov25',files,ignore.case=T)];
nov26files = files[grepl('nov26',files,ignore.case=T)];

buildComboGraphFromFiles(nov6files, 'nov6')


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
buildDotFile(countsToProbs(nov6countsCombined),nov6countsCombined,'nov6.dot');

########################

# T vs PKT group graphs for each day

source('~/Documents/_Fernald_lab/_code/behavior_syntax.R');
setwd('~/Documents/_Fernald_lab/matLabBoldnessRA2/');
files = list.files()[grepl('clean$', list.files())];

ids_t = '^3|^5|^7|^9';
ids_pkt = '^2|^12|^13|^15|^18|^21';
days = c('nov6', 'nov7', 'nov25', 'nov26');

for (d in days) {
	assign(paste('t_', d, sep=''), files[grepl(ids_t, files) & grepl(d, files)]);
}; rm(d);

for (d in days) {
	assign(paste('pkt_', d, sep=''), files[grepl(ids_pkt, files) & grepl(d, files)]);
}; rm(d);

t_days = ls()[grepl('^t_nov',ls())];
for (d in t_days) {
	buildComboGraphFromFiles(get(d), d);
}; rm(d);

pkt_days = ls()[grepl('^pkt_nov',ls())];
pkt_days = pkt_days[-which(pkt_days=='pkt_nov26')];
for (d in pkt_days) {
	buildComboGraphFromFiles(get(d), d);
}; rm(d);



################


pkt_nov7data = list();
for (f in 1:length(pkt_nov7)) {
	#print(pkt_nov7[f])
	pkt_nov7data[[f]] = getDataRA(pkt_nov7[f]);
	names(pkt_nov7data)[f] = pkt_nov7[f];
}; rm(f);

pkt_nov7counts = list();
for (s in 1:length(pkt_nov7data)) {
	this_dat = pkt_nov7data[[s]]$rawData;
	pkt_nov7counts[[s]] = getCountMat(this_dat);
	names(pkt_nov7counts)[s] = names(pkt_nov7data)[s];
}; rm(s,this_dat);

pkt_nov7countsCombined = collapseCountMatList(pkt_nov7counts);


###############


