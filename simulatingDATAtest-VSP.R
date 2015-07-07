# clear workspace and set working directory
rm(list = ls());
setwd(dir = '~/Documents/_tutorial-Neuron_paper');

# load libraries and set options
library(WGCNA);
options(stringsAsFactors = F);
allowWGCNAThreads();

setwd('~/Documents/_analysis_compare/_workspaces');

load('WORKSPACE-VSP_NET27-areaXcompareDAVIDnewestNEWpgNOV-PAINT');

# 

ls(); dim(DATA.VSP);

nSingers = sum(traitsVSP$motifs2>10);
nNonSingers = sum(traitsVSP$motifs2<10);

highSinger = rownames(traitsVSP)[traitsVSP$motifs2 == max(traitsVSP$motifs2)];
nonSinger = rownames(traitsVSP)[traitsVSP$motifs2 == 0][1];

DATA.VSP.highSinger = as.vector(DATA.VSP[rownames(DATA.VSP)==highSinger, ], mode='numeric');
names(DATA.VSP.highSinger) = names(DATA.VSP);
DATA.VSPnorm.highSinger = as.vector(DATA.VSPnorm[rownames(DATA.VSPnorm)==highSinger, ], mode='numeric');
names(DATA.VSPnorm.highSinger) = names(DATA.VSPnorm);

DATA.VSP.nonSinger = as.vector(DATA.VSP[rownames(DATA.VSP)==nonSinger, ], mode='numeric');
names(DATA.VSP.nonSinger) = names(DATA.VSP);
DATA.VSPnorm.nonSinger = as.vector(DATA.VSPnorm[rownames(DATA.VSPnorm)==nonSinger, ], mode='numeric');
names(DATA.VSPnorm.nonSinger) = names(DATA.VSPnorm);

norm.test = apply(DATA.VSP, 2, shapiro.test);
norm.test.pvals = 1:ncol(DATA.VSP);
for (g in 1:length(norm.test))
{
	norm.test.pvals[g] = norm.test[[g]][[2]]
}

cvs = 1:ncol(DATA.VSP);
for (g in 1:ncol(DATA.VSP))
{
	cvs[g] = sd(DATA.VSP[,g]) / mean(DATA.VSP[,g])
}

sds = 1:ncol(DATA.VSP);
for (g in 1:ncol(DATA.VSP))
{
	sds[g] = sd(DATA.VSP[,g]);
}

######

DATA.both = data.frame(DATA.VSP.highSinger, DATA.VSP.nonSinger);
meansAll = apply(DATA.VSP, 2, mean);
meansBoth = apply(DATA.both, 1, mean);


















############


DATA.VSP.highSinger.norm = DATA.VSP.highSinger / max(DATA.VSP.highSinger,na.rm=T);
DATA.VSP.nonSinger.norm = DATA.VSP.nonSinger / max(DATA.VSP.nonSinger,na.rm=T);

DATA.VSP.boot = DATA.VSP;

for (s in 1:nSingers)
{
	
}