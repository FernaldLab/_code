rm(list=ls());options(stringsAsFactors=F);setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
library(WGCNA); allowWGCNAThreads();
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');



load(file='dat.hs_dat.nhp_noZeros_removeNoise_higher-hs-brainVar.Rdata')

dat.hs = dat.hs[rownames(dat.hs) %in% rownames(dat.nhp), ];

keepMe = names(dat.hs)!='hsa.br.M.4';

##### Allow up to 1 missing value per gene
zthresh = 1;
X = dat.hs[,keepMe];
out.hs = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);


DATAhs = as.data.frame(t(out.hs$data_Qnorm));

sft = exn.plotPowers(DATAhs, networkType='signed', blockSize=ncol(DATAhs)+1, verbose=3);

k = exn.computeConnectivityAndPlotScaleFreeness(DATAhs, networkType='signed', power=18);

blockwiseModulesEnrichedIterate(DATA=DATAhs, power=18, 								                                densityPermTest=T, skipThresh=300,
							    deepSplitVec=c(2,4),
							    mergeCutHeightVec=c(.15,.2),
							    minModuleSizeVec=c(10,20,40,80,100));
