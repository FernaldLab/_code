rm(list=ls());
setwd('/Volumes/fishstudies/_methylation/new3157v3165overlapsCH');
files = list.files();
l = list();
for (f in 1:length(files)) {
	l[[f]] = read.table(files[f], header=F, sep='\t', colClasses='character');
	names(l)[f] = files[f];
}; rm(f);
names(l) = c('TE','LNC','SNP_ass','UTR','CDS','intron','SNP_ND','CNV','SNP_D','microsat');

# create version of list with duplicate methylation hits removed from each matrix
lm = l;
lm = lapply(lm, function(f) f[!duplicated(paste(f$V1, f$V2, f$V3)), ]);

# function to add column with methylation level of animal that had higher level
# will replace last column
.addHighLevelCol = function(mat, col1=11, col2=12) {
	lastCol = ncol(mat);
	mat[, lastCol] = apply(mat, 1, function(f) if(as.numeric(f[col1])>as.numeric(f[col2])){as.numeric(f[col1])} else {as.numeric(f[col2])});
	return(mat);
}
lm = lapply(lm, .addHighLevelCol);

# function to add column denoting whether methylation hit is on sense or antisense strand to feature
.addSenseAntisenseCol = function(mat, col1, col2) {
	mat = cbind(mat, sense=mat[,col1]==mat[,col2]);
	mat$sense[mat$sense==T] = 'S';
	mat$sense[mat$sense==F] = 'AS';
	return(mat);
}



# investigate TEs
lm$TE = .addSenseAntisenseCol(lm$TE, 4, 23);
TEtype = levels(as.factor(lm$TE$V21));
TEtype = TEtype[-4];	# remove RC since it has only one hit

par(mfrow=c(2,3));
for (te in TEtype) {
	verboseBoxplot(lm$TE$V24[lm$TE$V21==te], 
				   as.factor(paste(lm$TE$V14[lm$TE$V21==te], lm$TE$V4[lm$TE$V21==te])),
				   xlab='', ylab='meth.lvl', main=te, col='grey'
				   );
}; rm(te);

par(mfrow=c(2,3));
for (te in TEtype) {
	verboseBoxplot(lm$TE$V24[lm$TE$V21==te], 
				   as.factor(paste(lm$TE$V14[lm$TE$V21==te], lm$TE$sense[lm$TE$V21==te])),
				   xlab='', ylab='meth.lvl', main=te, col='grey'
				   );
}; rm(te);

# within TE type, compare levels of dinucleotides on sense vs antisense
te = 'DNA';
te_ct = levels(as.factor(lm$TE$V14[lm$TE$V21==te]));
par(mfrow=c(1,4));
for (j in 1:length(te_ct)) {
	g = as.factor(lm$TE$sense[lm$TE$V21==te & lm$TE$V14==te_ct[j]]);
	num = lm$TE$V24[lm$TE$V21==te & lm$TE$V14==te_ct[j]];
	if (length(unique(g))<2) {
		num = c(num, 0);
		g = c(g, 'X');
	}
	verboseBoxplot(num, g, xlab='',ylab='meth.lvl',main=paste(te, te_ct[j]),col='grey');
}; rm(j,te_ct);

# plot levels of trinucleotides centered on methylated position
centered_tri = apply(lm$TE, 1, function(f) if(f[4]=='+'){substr(f[5],2,4)}else{substr(f[10],2,4)});
te = 'DNA';
te_ct = levels(as.factor(centered_tri[lm$TE$V21==te & lm$TE$V16=='CH']));
par(mfrow=c(2,6));
for (j in 1:length(te_ct)) {
	g = as.factor(lm$TE$V4[lm$TE$V21==te & centered_tri==te_ct[j] & lm$TE$V16=='CH']);
	num = lm$TE$V24[lm$TE$V21==te & centered_tri==te_ct[j] & lm$TE$V16=='CH'];
	if (length(unique(g))<2) {
		num = c(num, 0);
		g = c(g, 'X');
	}
	verboseBoxplot(num, g, xlab='',ylab='meth.lvl',main=paste(te, te_ct[j]),col='grey');
}; rm(te_ct,j,num,g);


# function to plot frequencies of trinucleotides centered on methylated position
.plotTriCenFreqs = function(te, ...) {
	tmpp = sort(table(centered_tri[lm$TE$V21==te & lm$TE$V4=='+' & lm$TE$V16=='CH']));
	tmpp = tmpp / sum(table(centered_tri[lm$TE$V21==te & lm$TE$V4=='+']));
	tmpm = sort(table(centered_tri[lm$TE$V21==te & lm$TE$V4=='-' & lm$TE$V16=='CH'])); 
	tmpm = tmpm / sum(table(centered_tri[lm$TE$V21==te & lm$TE$V4=='-']));
	if (length(tmpp)==length(tmpm)) {
		tmpm = tmpm[match(names(tmpp), names(tmpm))];
	} else if (length(tmpp)<length(tmpm)) {
		tmpp = tmpp[match(names(tmpm), names(tmpp))];
		tmpp[is.na(tmpp)] = 0;
		names(tmpp)[is.na(names(tmpp))] = names(tmpm)[is.na(names(tmpp))];
	} else {
		tmpm = tmpm[match(names(tmpp), names(tmpm))];
		tmpm[is.na(tmpm)] = 0;
		names(tmpm)[is.na(names(tmpm))] = names(tmpp)[is.na(names(tmpm))];
	}
	names(tmpp) = paste(names(tmpp), '+', sep='');
	names(tmpm) = paste(names(tmpm), '-', sep='');
	tmp = c();
	for (i in 1:length(tmpp)) {
		tmp = c(tmp, tmpp[i], tmpm[i]);
	}
	barplot(tmp,las=2,horiz=T,col=rep(c('black','grey'),length(tmp)/2),main=te, ...);
}


.plotTriCenFreqsSAS = function(te, ...) {
	tmpp = sort(table(centered_tri[lm$TE$V21==te & lm$TE$sense=='S' & lm$TE$V16=='CH']));
	#tmpp = tmpp / sum(table(centered_tri[lm$TE$V21==te & lm$TE$sense=='S']));
	tmpp = tmpp / sum(tmpp);
	tmpm = sort(table(centered_tri[lm$TE$V21==te & lm$TE$sense=='AS' & lm$TE$V16=='CH'])); 
	#tmpm = tmpm / sum(table(centered_tri[lm$TE$V21==te & lm$TE$sense=='AS']));
	tmpm = tmpm / sum(tmpm);
	if (length(tmpp)==length(tmpm)) {
		tmpm = tmpm[match(names(tmpp), names(tmpm))];
	} else if (length(tmpp)<length(tmpm)) {
		tmpp = tmpp[match(names(tmpm), names(tmpp))];
		tmpp[is.na(tmpp)] = 0;
		names(tmpp)[is.na(names(tmpp))] = names(tmpm)[is.na(names(tmpp))];
	} else {
		tmpm = tmpm[match(names(tmpp), names(tmpm))];
		tmpm[is.na(tmpm)] = 0;
		names(tmpm)[is.na(names(tmpm))] = names(tmpp)[is.na(names(tmpm))];
	}
	names(tmpp) = paste(names(tmpp), ':A', sep='');
	names(tmpm) = paste(names(tmpm), ':AS', sep='');
	tmp = c();
	for (i in 1:length(tmpp)) {
		tmp = c(tmp, tmpp[i], tmpm[i]);
	}
	barplot(tmp,las=2,horiz=T,col=rep(c('black','grey'),length(tmp)/2),main=te, ...);
}


par(mfrow=c(2,3), oma=c(0,0,2,0));
for (te in TEtype) {
	.plotTriCenFreqs(te);
}; rm(te);
title('DM, NCH trinucleotides - % of all trinucleotides in each TE type', outer=T);











verboseBoxplot(lm$TE$V24[lm$TE$V21==te & lm$TE$V16=='CH'], 
				   as.factor(paste(centered_tri[lm$TE$V21==te & lm$TE$V16=='CH'], lm$TE$V4[lm$TE$V21==te & lm$TE$V16=='CH'])),
				   xlab='', ylab='meth.lvl', main=te, col='grey'
				   );






li = l;
li = lapply(li, function(f) f[!duplicated(paste(f$V1, f$V2, f$V3)), ]);



###
# 

.getDiTriByMeth = function(mat, cols=14:17) {
	return(mat[!duplicated(paste(mat$V1, mat$V2, mat$V3)), cols]);
}


