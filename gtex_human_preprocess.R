rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');

# read in raw data
datGTEX = read.table('Extened_nonHuman_GTEx_Data/Raw Data/GTEx_RPKM_Network_Matrix.txt', header=T, sep='\t');

# strip transcript numbers from human Ensembl ids
x = unlist(strsplit(datGTEX$Name,'\\.'));
datGTEX$Name = x[grep('^E', x)]; rm(x);

rownames(datGTEX) = datGTEX$Name;
datGTEX = datGTEX[, -1];

# remove genes with too many zeros
#zthresh = ceiling(ncol(datGTEX)/3);
zthresh = 1;
datGTEX = datGTEX[apply(datGTEX, 1, function(f) sum(f==0)) < zthresh, ];

# remove noise 
datGTEX[datGTEX < .5] = NA;

countNAs = apply(datGTEX, 1, function(f) sum(is.na(f)));
datGTEX = datGTEX[countNAs < zthresh, ];

# log2 transform, don't know why but this makes sample clustering work
datGTEXlog2 = log2(datGTEX);


#####################################

source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
# remove outliers and quantile normalize
#keepMe = names(datGTEX)!='hsa.br.M.4' #& names(dat.hs)!='hsa.kd.F.1';

#X = dat.hs[,keepMe];
X = datGTEXlog2
out.hs = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
###########################################################################################

source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');
library(WGCNA); allowWGCNAThreads();

DATAhs = as.data.frame(t(out.hs$data_Qnorm));

sft = exn.plotPowers(DATAhs, networkType='signed', blockSize=ncol(DATAhs)+1, verbose=3);
k16 = exn.computeConnectivityAndPlotScaleFreeness(DATAhs, networkType='signed', power=16);

stem = 'GTEX_no0_remove.5_noNA_log2_default_preProc';
dir.create(stem);
setwd(stem);
assign(stem, DATAhs);

blockwiseModulesEnrichedIterate(DATA=get(stem), power=16, densityPermTest=T, skipThresh=200, deepSplitVec=c(2, 4), permTestPvalThresh=0.005);