rm(list=ls());
setwd('/Volumes/fishstudies/_methylation/');

# function to add column of dinucleotides 
.addDiCol = function (mat) {
	mat = cbind(mat, di=rep('N', nrow(x2)));
	mat$di = as.character(mat$di);
	for (row in 1:nrow(mat)) {
		if (mat[row,4]=="+") {
			mat$di[row] = substr(mat[row,5],3,4);
		}
		if (mat[row,4]=="-") {
			mat$di[row] = substr(mat[row,10],3,4);
		}
	}
	return(mat);
}

# function to add column of trinucleotides
.addTriCol = function (mat) {
	mat = cbind(mat, tri=rep("N",nrow(mat)));
	mat$tri = as.character(mat$tri);
	for (row in 1:nrow(mat)) {
		if (mat[row,4]=="+") {
			mat$tri[row]=substr(mat[row,5],3,5);
		}
		if (mat[row,4]=="-") {
			mat$tri[row]=substr(mat[row,10],3,5);
		}
	}
	return(mat);
}

########
# read in data
x = read.table('3157_vs_3165.filtered.50percent.significant.bed_Gflip2', colClasses='character');

# add Fisher's exact test p-vals
x2 = cbind(x, p=apply(x, 1, function(f) fisher.test(matrix(as.numeric(c(f[6:7],f[8:9])),ncol=2))$p.value));

# add dinucleotide column and clean
x2 = .addDiCol(x2);
x2 = x2[!(x2$di %in% c('CN', '00')), ];

# add trinucleotide column
x2 = .addTriCol(x2);

# add dinucleotide CH column
x2 = cbind(x2, diH=x2$di);
x2$diH[x2$diH %in% c('CA','CC','CT')] = 'CH';

# add trinucleotide CH column
x2 = cbind(x2, triH=x2$tri);
x2$triH[!(x2$triH %in% c('CGA','CGC','CGG','CGT'))] = paste(substr(x2$triH[!(x2$triH %in% c('CGA','CGC','CGG','CGT'))],1,1), 
															'H', 
															substr(x2$triH[!(x2$triH %in% c('CGA','CGC','CGG','CGT'))],3,3), 
															sep='');


# check how often D or ND had higher methylation level
sum(x2$V12>x2$V11);
sum(x2$V12<x2$V11);

### analyze ###
library(WGCNA);
source('/Volumes/fishstudies/_code/bootstrapFunctions_6-16-13.R');

# check if D or ND had higher average level
verboseBoxplot(c(as.numeric(x2$V12), as.numeric(x2$V11)), 
			   c(rep('D',nrow(x2)), rep('ND',nrow(x2))), 
			   xlab='', ylab='methylation level',
			   main=paste('diff.meth positions(n=',nrow(x2),')\n',sep=''),
			   col='grey'
			   );
a=bootstrap.2independent(as.numeric(x2$V12), as.numeric(x2$V11),
								  Func='median', trials=10000,
								  dataDescriptor='methylation level', 
								  groupNames=c('D','ND')
								  );
par(mfrow=c(1,2), oma=c(0,0,2,0));
hist(as.numeric(x2$V11), breaks=nrow(x2), xlab='meth.lvl', main='ND');
hist(as.numeric(x2$V12), breaks=nrow(x2), xlab='meth.lvl', main='D');
title('diff.meth.positions', outer=T);
## Kruskal test in verboseBoxplot says D has significantly higher level, 
## but bootstrap says no difference, and distributions are clearly at least bimodal, thus...

# compare levels of each animal only from hits where that animal had higher level
x2D = x2[x2$V12>x2$V11, ];
x2ND = x2[x2$V11>x2$V12, ];


write.table(x2,file="3157_vs_3165.filtered.50percent.significant.bed_Gflip2_R",sep="\t",quote=F,row.names=F,col.names=F);
write.table(x2ND,file="3157_vs_3165.filtered.50percent.significant.bed_Gflip2_R_ND",sep="\t",quote=F,row.names=F,col.names=F);
write.table(x2D,file="3157_vs_3165.filtered.50percent.significant.bed_Gflip2_R_D",sep="\t",quote=F,row.names=F,col.names=F);

# also check each strand
# Kruskal
par(mfrow=c(1,3), oma=c(0,0,2,0));
verboseBoxplot(c(as.numeric(x2D$V12), as.numeric(x2ND$V11)), 
			   c(rep(paste('D (n=',nrow(x2D),')',sep=''),nrow(x2D)), 
			   	 rep(paste('ND (n=',nrow(x2ND),')',sep=''),nrow(x2ND))),
			   xlab='', ylab='methylation level', main='both strands\n', col='grey');
			   
verboseBoxplot(c(as.numeric(x2D$V12[x2D$V4=='+']),as.numeric(x2ND$V11[x2ND$V4=='+'])),
			   c(rep(paste('D (n=',sum(x2D$V4=='+'),')',sep=''), sum(x2D$V4=='+')),
			     rep(paste('ND (n=',sum(x2ND$V4=='+'),')',sep=''), sum(x2ND$V4=='+'))),
			   xlab='', ylab='methylation level', main='fwd strand\n', col='grey');

verboseBoxplot(c(as.numeric(x2D$V12[x2D$V4=='-']),as.numeric(x2ND$V11[x2ND$V4=='-'])),
			   c(rep(paste('D (n=',sum(x2D$V4=='-'),')',sep=''), sum(x2D$V4=='-')),
			     rep(paste('ND (n=',sum(x2ND$V4=='-'),')',sep=''), sum(x2ND$V4=='-'))),
			   xlab='', ylab='methylation level', main='rev strand\n', col='grey');
title('only levels sig. higher than in other', outer=T);

# bootstrap
a=bootstrap.2independent(as.numeric(x2D$V12), as.numeric(x2ND$V11),
								   Func='median', trials=10000,
								   dataDescriptor='methylation level', groupNames=c('D','ND')
								   );
									   
a=bootstrap.2independent(as.numeric(x2D$V12[x2D$V4=='+']), as.numeric(x2ND$V11[x2ND$V4=='+']),
									  Func='median', trials=10000,
									  dataDescriptor='methylation level', groupNames=c('D-fwd','ND-fwd')
									  );
									  
a=bootstrap.2independent(as.numeric(x2D$V12[x2D$V4=='-']), as.numeric(x2ND$V11[x2ND$V4=='-']),
									  Func='median', trials=10000,
									  dataDescriptor='methylation level', groupNames=c('D-rev','ND-rev')
									  );
									  
# check levels for different dinucleotides
par(mfrow=c(2,2));
verboseBoxplot(as.numeric(c(x2$V11,x2$V12)), as.factor(c(x2$di,x2$di)),xlab='',ylab='meth.lvl',main=paste('DM,ND&D all (n=',nrow(x2)*2,')\n',sep=''),col='grey');
verboseBoxplot(as.numeric(c(x2$V11,x2$V12)), as.factor(c(x2$diH,x2$diH)),xlab='',ylab='meth.lvl',main=paste('DM,ND&D all (n=',nrow(x2)*2,')\n',sep=''),col='grey');
verboseBoxplot(as.numeric(c(x2ND$V11,x2D$V12)), as.factor(c(x2ND$di,x2D$di)),xlab='',ylab='meth.lvl',main=paste('DM,ND&D high (n=',nrow(x2)*2,')\n',sep=''),col='grey');
verboseBoxplot(as.numeric(c(x2ND$V11,x2D$V12)), as.factor(c(x2ND$diH,x2D$diH)),xlab='',ylab='meth.lvl',main=paste('DM,ND&D high (n=',nrow(x2)*2,')\n',sep=''),col='grey');


CG_hi = as.numeric(c(x2ND$V11[x2ND$diH=='CG'], x2D$V12[x2D$diH=='CG']));
CH_hi = as.numeric(c(x2ND$V11[x2ND$diH=='CH'], x2D$V12[x2D$diH=='CH']));
di_hiMat = matrix(nrow=max(length(CG_hi),length(CH_hi)), ncol=2, dimnames=list(NULL, c('CG','CH')));
di_hiMat[1:length(CG_hi), 1] = CG_hi;
di_hiMat[1:length(CH_hi), 2] = CH_hi;
a=bootstrap.2independent(CG_hi, CH_hi, 
					     Func='median', trials=10000, 
					     dataDescriptor='meth.lvl',groupNames=c('CG','CH')
					     );
					     
# check levels for different trinucleotides	
par(mfrow=c(1,3));			     
verboseBoxplot(as.numeric(c(x2D$V12[x2D$diH=='CH'], x2ND$V11[x2ND$diH=='CH'])), 
			   as.factor(c(x2D$triH[x2D$diH=='CH'], x2ND$triH[x2ND$diH=='CH'])),
			   xlab='', ylab='meth.lvl', main='DM,ND & D high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2D$V12[x2D$diH=='CH']), 
			   as.factor(x2D$triH[x2D$diH=='CH']),
			   xlab='', ylab='meth.lvl', main='DM, D high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2ND$V11[x2ND$diH=='CH']), 
			   as.factor(x2ND$triH[x2ND$diH=='CH']),
			   xlab='', ylab='meth.lvl', main='DM, ND high\n', col='grey'
			   );
			   
			   
			   
par(mfrow=c(2,3));		   
verboseBoxplot(as.numeric(x2D$V12[x2D$diH=='CH' & x2D$V4=='+']), 
			   as.factor(x2D$triH[x2D$diH=='CH' & x2D$V4=='+']),
			   xlab='', ylab='meth.lvl', main='DM, D fwd high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2D$V12[x2D$diH=='CH' & x2D$V4=='-']), 
			   as.factor(x2D$triH[x2D$diH=='CH' & x2D$V4=='-']),
			   xlab='', ylab='meth.lvl', main='DM, D rev high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2D$V12[x2D$diH=='CH']), 
			   as.factor(x2D$triH[x2D$diH=='CH']),
			   xlab='', ylab='meth.lvl', main='DM, D high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2ND$V11[x2ND$diH=='CH' & x2ND$V4=='+']), 
			   as.factor(x2ND$triH[x2ND$diH=='CH' & x2ND$V4=='+']),
			   xlab='', ylab='meth.lvl', main='DM, ND fwd high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2ND$V11[x2ND$diH=='CH' & x2ND$V4=='-']), 
			   as.factor(x2ND$triH[x2ND$diH=='CH' & x2ND$V4=='-']),
			   xlab='', ylab='meth.lvl', main='DM, ND rev high\n', col='grey'
			   );
verboseBoxplot(as.numeric(x2ND$V11[x2ND$diH=='CH']), 
			   as.factor(x2ND$triH[x2ND$diH=='CH']),
			   xlab='', ylab='meth.lvl', main='DM, ND high\n', col='grey'
			   );
			   
			   
			   
par(mfrow=c(2,4));
for (tri in c('CHA','CHC','CHG','CHT')) {
	verboseBoxplot(as.numeric(x2D$V12[x2D$triH==tri]), 
				   as.factor(x2D$V4[x2D$triH==tri]), 
				   col='grey', xlab=tri, ylab='meth.lvl', main='D, DM high\n'
				   );
}; rm(tri);
for (tri in c('CHA','CHC','CHG','CHT')) {
	verboseBoxplot(as.numeric(x2ND$V11[x2ND$triH==tri]), 
				   as.factor(x2ND$V4[x2ND$triH==tri]), 
				   col='grey', xlab=tri, ylab='meth.lvl', main='ND, DM high\n'
				   );
}; rm(tri);