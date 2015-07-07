rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('~/Downloads/RemoveOutliers/RemoveOutliers_0.9');

#library(affy)library(cluster)library(impute)library(preprocessCore)

# read in raw data
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];
rownames(dat) = ids$hsa;
dat=cbind(rownames(dat),dat);
names(dat)[1]='id';
dat0 = dat;

# try removing heart
remove = 'ts';
dat = dat[, !grepl(remove,names(dat))];

#numzeros = apply(dat[, grep('hs',names(dat))], 1, function(f) sum(f==0));
numzeros = apply(dat[, 2:61], 1, function(f) sum(f==0));
dat = dat[numzeros==0, ];

sampleInfo = as.data.frame(matrix(nrow=ncol(dat)-1, ncol=4));
names(sampleInfo) = c('sampleID','species','tissue','sex');
sampleInfo$sampleID = names(dat)[-1];
sampleInfo$species = substr(names(dat)[-1],1,3);
sampleInfo$tissue = substr(names(dat)[-1],5,6);
sampleInfo$sex = substr(names(dat)[-1],8,8);
sampleInfo = cbind(sampleInfo, fake_gp=c(rep('ok',60), rep('sucks',11)))
#tmp = substr(names(dat)[-1],1,3);
tmp = substr(names(dat)[-1],5,6);

indices = list();
sp = names(table(tmp));
for (s in 1:length(sp)) {
	indices[[s]] = grep(sp[s], names(dat));
}; rm(s);
#indices=indices[c(2,6,4,1,5,3)];
indices=indices[c(2,1)]

RemoveOutliers(datexpr1=dat,
			   subset1=NULL,
			   impute1=FALSE,
			   skip1=1,
			   indices1=indices,
			   sampleinfo1=sampleInfo,
			   samplelabels1=1,
			   grouplabels1=5,
			   trait1=c(2,3,4),
			   btrait1=c(2,3,4),
			   asfactors1=c(2,3,4),
			   projectname1='combatTest',
			   cexlabels=.9,
			   normalize1=TRUE,
			   replacenegs1=TRUE,
			   fitmodels1=TRUE,
			   whichfit1='mean',
			   exportfigures1=TRUE
			   )
			   
			   
#############################################
			   

			   