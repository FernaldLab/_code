rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');

# read in raw data
datGTEX = read.table('Extened_nonHuman_GTEx_Data/Raw Data/GTEx_RPKM_Network_Matrix.txt', header=T, sep='\t');


dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat0 = dat0[, 6:23];
rownames(dat0) = ids$hsa; rm(ids);

# strip transcript numbers from human Ensembl ids
x = unlist(strsplit(datGTEX$Name,'\\.'));
datGTEX$Name = x[grep('^E', x)]; rm(x);

# restrict to common genes
common = intersect(datGTEX$Name, rownames(dat0));
datGTEX = datGTEX[match(common, datGTEX$Name), ];
dat0 = dat0[match(common, rownames(dat0)), ];

dat = cbind(datGTEX, dat0);

### will change some after getting GTEX M/F info

names(dat) = gsub('Frontal_Cortex|br', 'frontalCtx', names(dat));
names(dat) = gsub('Cerebellum|cb', 'cerebellum', names(dat));
names(dat) = gsub('Heart|ht', 'heart', names(dat));
names(dat) = gsub('Kidney|kd', 'kidney', names(dat));
names(dat) = gsub('Liver|lv', 'liver', names(dat));

names(dat) = gsub('hsa.', '', names(dat));

# remove genes with zeros
numzeros = apply(dat[, 2:ncol(dat)], 1, function(f) sum(f==0));
dat = dat[numzeros==0, ];

# log???
#



### build sampleInfo table for RemoveOutliers
sampleInfo = as.data.frame(matrix(nrow=ncol(dat)-1, ncol=3));
names(sampleInfo) = c('sampleID','tissue','study');

# add sample IDs
sampleInfo$sampleID = names(dat)[-1];

# add study info
gtexSamples = grepl('GTEX', sampleInfo$sampleID);
sampleInfo$study[gtexSamples] = 'GTEX';
sampleInfo$study[!gtexSamples] = 'HSA';

# add tissue info
tmp = unlist(strsplit(sampleInfo$sampleID[gtexSamples],'_'));
sampleInfo$tissue[gtexSamples] = tmp[!grepl('^GTEX',tmp)];

tmp = unlist(strsplit(sampleInfo$sampleID[!gtexSamples],'\\.'));
sampleInfo$tissue[!gtexSamples] = tmp[seq(1,length(tmp),3)];

# limit to tissues in both datasets???
common_tissue = intersect(names(table(sampleInfo$tissue[gtexSamples])), names(table(sampleInfo$tissue[!gtexSamples])));
common_tissueInd = (1:nrow(sampleInfo))[sampleInfo$tissue %in% common_tissue];
sampleInfo = sampleInfo[common_tissueInd, ];
dat = dat[, c(1, (common_tissueInd + 1))];

# build group indices
grp = names(table(sampleInfo$tissue));

indices = list();
for (s in grp) {
	indices[[s]] = grep(s, names(dat));
}; rm(s);
indices = indices[c(2,1,4,5,3)];


source('~/Downloads/RemoveOutliers/RemoveOutliers_0.9');
#library(affy)
library(cluster)
library(impute)
library(preprocessCore)


RemoveOutliers(datexpr1=dat,
			   subset1=NULL,
			   impute1=FALSE,
			   skip1=1,
			   indices1=indices,
			   sampleinfo1=sampleInfo,
			   samplelabels1=1,
			   grouplabels1=2,
			   trait1=3,
			   btrait1=3,
			   asfactors1=3,
			   projectname1='combatTest',
			   cexlabels=.9,
			   normalize1=TRUE,
			   replacenegs1=TRUE,
			   fitmodels1=TRUE,
			   whichfit1='mean',
			   exportfigures1=TRUE
			   )