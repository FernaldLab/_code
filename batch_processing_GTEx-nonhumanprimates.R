rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');

# read in raw data
dat.hs0 = read.table('Extened_nonHuman_GTEx_Data/Raw Data/GTEx_RPKM_Network_Matrix.txt', header=T, sep='\t');
dat.nhp0 = read.table('Extened_nonHuman_GTEx_Data/Raw Data/NonHuman_FPKM_Network_Matrx_ENSEMBL.txt', header=T, sep='\t');

# strip transcript numbers from human Ensembl ids
x = unlist(strsplit(dat.hs0$Name,'\\.'));
dat.hs0$Name = x[grep('^E', x)]; rm(x);

# restrict to common genes
common = intersect(dat.hs0$Name, rownames(dat.nhp0));
dat.hs = dat.hs0[match(common, dat.hs0$Name), ];
dat.nhp = dat.nhp0[match(common, rownames(dat.nhp0)), ];

dat = cbind(dat.hs, dat.nhp);
save(dat, file='Extened_nonHuman_GTEx_Data/Raw Data/GTEx-NonHuman_combined_commonGenes.RData');
rm(dat.hs0, dat.nhp0, dat.hs, dat.nhp);gc();

###########################################################################################################################

rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('~/Downloads/RemoveOutliers/RemoveOutliers_0.9');
#library(affy)
library(cluster)
library(impute)
library(preprocessCore)

load(file='Extened_nonHuman_GTEx_Data/Raw Data/GTEx-NonHuman_combined_commonGenes.RData');

# remove genes with zeros
numzeros = apply(dat[, 2:ncol(dat)], 1, function(f) sum(f==0));
dat = dat[numzeros==0, ];

# log???
#

### build sampleInfo table for RemoveOutliers
sampleInfo = as.data.frame(matrix(nrow=ncol(dat)-1, ncol=4));
names(sampleInfo) = c('sampleID','species','tissue','study');

# add sample IDs
sampleInfo$sampleID = names(dat)[-1];

# add species info
hsSamples = grepl('GTEX', sampleInfo$sampleID);
sampleInfo$species[hsSamples] = 'GTEX';

tmp = sampleInfo$sampleID[!hsSamples];
tmp = unlist(strsplit(tmp, '_'));
sampleInfo$species[!hsSamples] = tmp[seq(2, length(tmp), 2)];

# add tissue info
sampleInfo$tissue[hsSamples] = gsub('_.*$', '', sampleInfo$sampleID[hsSamples]);
sampleInfo$tissue[!hsSamples] = tmp[seq(1, length(tmp), 2)];

# make sure everything matches up
sampleInfo$tissue = gsub('Cerebelum', 'Cerebellum', sampleInfo$tissue);
sampleInfo$sampleID = gsub('Cerebelum', 'Cerebellum', sampleInfo$sampleID);
names(dat) = gsub('Cerebelum', 'Cerebellum', names(dat));

sampleInfo$tissue = gsub('Frontal|FrontalCortex', 'FrontalCtx', sampleInfo$tissue);
sampleInfo$sampleID = gsub('Frontal_Cortex|FrontalCortex', 'FrontalCtx', sampleInfo$sampleID);
names(dat) = gsub('Frontal_Cortex|FrontalCortex', 'FrontalCtx', names(dat));

# add study info
sampleInfo$study[hsSamples] = 'GTEX';
sampleInfo$study[!hsSamples] = 'NHP';

# limit to tissues in both datasets???
common_tissue = intersect(names(table(sampleInfo$tissue[hsSamples])), names(table(sampleInfo$tissue[!hsSamples])));
common_tissueInd = (1:nrow(sampleInfo))[sampleInfo$tissue %in% common_tissue];
sampleInfo = sampleInfo[common_tissueInd, ];
dat = dat[, c(1, (common_tissueInd + 1))];

### build list of group indices
grpInd = 3;
grp = names(table(sampleInfo[,grpInd]));

indices = list();
for (s in grp) {
	indices[[s]] = grep(s, names(dat));
}; rm(s);

# # if only common tissues
# indices = indices[c(3,1,6,5,2,7,4)];

# # if all tissues included   
# indices = indices[c(4,2,10,6,7,8,3,11,5,12,1,9)];

# # if using species
# indices = indices[c(4,2,3,1,6,7,5)];

tr = (1:ncol(sampleInfo))[!(1:ncol(sampleInfo) %in% grpInd)];
samplelabels1 = 1;
tr = tr[!(tr %in% samplelabels1)];

RemoveOutliers(datexpr1=dat,
			   subset1=NULL,
			   impute1=FALSE,
			   skip1=1,
			   indices1=indices,
			   sampleinfo1=sampleInfo,
			   samplelabels1=samplelabels1,
			   grouplabels1=grpInd,
			   trait1=tr,
			   btrait1=tr,
			   asfactors1=tr,
			   projectname1='combatTest',
			   cexlabels=.9,
			   normalize1=TRUE,
			   replacenegs1=TRUE,
			   fitmodels1=TRUE,
			   whichfit1='mean',
			   exportfigures1=TRUE
			   )









