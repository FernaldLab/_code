rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');

# read in raw data
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];

# separate human samples
dat.hs0 = dat[, grep('hs', names(dat))];
rownames(dat.hs0) = ids$hsa;

# try removing ts samples
#dat.hs0 = dat.hs0[, 1:16];

# separate non-human primate samples
keep = 'ptr|ppa'
dat.nhp0 = dat[, grep(keep, names(dat))];
rownames(dat.nhp0) = ids$hsa;

# clean up
rm(dat0, dat);

#########################################################################################

# process human data

dat.hs = dat.hs0;

# remove genes with too many zeros
zthresh = ceiling(ncol(dat.hs)/3);
dat.hs = dat.hs[apply(dat.hs, 1, function(f) sum(f==0)) < zthresh, ];

# remove noise 
dat.hs[dat.hs < .5] = NA;

# log2 transform, don't know why but this makes sample clustering work
dat.hs = log2(dat.hs);

# ttypes = unique(substr(gsub('hsa.','',names(dat.hs)),1,2));
# for (ttype in ttypes) {
	# assign(paste(ttype, 'Var', sep=''), rank(apply(dat.hs[,grepl(ttype,names(dat.hs))], 1, var)));
# }; rm(ttype);

# avgVarBrain = c();
# for (g in 1:length(brVar)) {
	# avgVarBrain = c(avgVarBrain, mean(c(brVar[g], cbVar[g])));
# }; rm(g);
# names(avgVarBrain) = names(brVar);

# avgVarNonBrain = c();
# for (g in 1:length(brVar)) {
	# avgVarNonBrain = c(avgVarNonBrain, mean(c(htVar[g], kdVar[g], lvVar[g], tsVar[g])));
# }; rm(g);
# names(avgVarNonBrain) = names(brVar);

# higherVar = avgVarNonBrain - avgVarBrain;
# dat.hs = dat.hs[higherVar < 0, ];



source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
# remove outliers and quantile normalize
keepMe = names(dat.hs)!='hsa.br.M.4' #& names(dat.hs)!='hsa.kd.F.1';

X = dat.hs[,keepMe];
#X = dat.hs
out.hs = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
				 
				 
# process nonhuman data

dat.nhp = dat.nhp0;

# remove genes with too many zeros
zthresh = ceiling(ncol(dat.nhp)/3);
dat.nhp = dat.nhp[apply(dat.nhp, 1, function(f) sum(f==0)) < zthresh, ];

# remove noise 
dat.nhp[dat.nhp < .5] = NA;

# log2 transform, don't know why but this makes sample clustering work
dat.nhp = log2(dat.nhp);

source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
# remove outliers and quantile normalize
keepMe = names(dat.nhp)!='ppa.ht.M.1' #& names(dat.nhp)!='hsa.kd.F.1';

X = dat.nhp[,keepMe];
#X = dat.nhp
out.nhp = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=NULL, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
save(out.hs, out.nhp, file='hsa_nhp_.33zeros_.5rpkm_log2_nohs.br.M.4_defaultPreProc.Rdata');				 
				 
###########################################################################################

source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');
library(WGCNA); allowWGCNAThreads();

DATAhs = as.data.frame(t(out.hs$data_Qnorm));
DATAnhp = as.data.frame(t(out.nhp$data_Qnorm));

DATAhs = DATAhs[, apply(DATAhs, 2, function(f) sum(is.na(f))) == 0];
common_genes = intersect(names(DATAhs), names(DATAnhp));
DATAhs = DATAhs[, names(DATAhs) %in% common_genes];
DATAnhp = DATAnhp[, names(DATAnhp) %in% common_genes];

#cv = function (x) { return(sd(x, na.rm=T) / mean(x, na.rm=T)) };
#DATAhs = DATAhs[, apply(DATAhs,2,cv) > 0];

sft = exn.plotPowers(DATAhs, networkType='signed', blockSize=ncol(DATAhs)+1, verbose=3);
k18 = exn.computeConnectivityAndPlotScaleFreeness(DATAhs, networkType='signed', power=18);

stem = 'hsa_1zeros_.5rpkm_log2_higherBrainVar_nohs.br.M.4_defaultPreProc_noNAs_common';
dir.create(stem);
setwd(stem);
assign(stem, DATAhs);
# blockwiseModulesEnrichedIterate(DATA=get(stem), power=18, densityPermTest=F, deepSplitVec=c(2, 3, 4), minModuleSizeVec=c(10, 25, 50, 75, 100), mergeCutHeightVec=c(0.05, 0.10, 0.15, 0.20));

blockwiseModulesEnrichedIterate(DATA=get(stem), power=18, densityPermTest=T, skipThresh=200, deepSplitVec=c(2, 3, 4), permTestPvalThresh=0.005);

########################

rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');
library(WGCNA); allowWGCNAThreads();


stem='hsa_.33zeros_.5rpkm_log2_nohs.br.M.4_defaultPreProc/get(stem)_signed_p18_ds4_mm40_mch0.15run11'
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));
DATAhs=DATA;NEThs=net;rm(DATA,net)
load(file='hsa_nhp_.33zeros_.5rpkm_log2_nohs.br.M.4_defaultPreProc.Rdata');
DATAnhp = as.data.frame(t(out.nhp$data_Qnorm));

commonGenes = intersect(names(DATAhs),names(DATAnhp));
DATAhs = DATAhs[, names(DATAhs) %in% commonGenes];
DATAnhp = DATAnhp[, names(DATAnhp) %in% commonGenes];

NEThs = blockwiseModulesEnriched(DATA=DATAhs, maxBlockSize=(ncol(DATAhs)+1), 
						 power=18, networkType='signed', 
						 deepSplit=4, minModuleSize=40, mergeCutHeight=.15,
						 verbose=3, densityPermTest=F
						 );
load('run17DATA.RData')
load('run17NET.RData')		
DATAhs0=DATAhs
DATAhs=DATA		
NEThs0=NEThs
NEThs=net

DATAnhp = DATAnhp[, names(DATAnhp) %in% names(DATAhs)];
		 
NETnhp = blockwiseModules(datExpr=DATAnhp, maxBlockSize=(ncol(DATAnhp)+1), 
						 power=18, networkType='signed', 
						 deepSplit=4, minModuleSize=40, mergeCutHeight=.15,
						 verbose=3
						 );

NETnhp$colors = matchLabels(NETnhp$colors, NEThs$colors);
NETnhp$MEs = moduleEigengenes(DATAnhp, NETnhp$colors)$eigengenes;
#################################################################################

exn.getNetworkBasics(NEThs, DATAhs, '.hs');
exn.getNetworkBasics(NETnhp, DATAnhp, '.nhp');


presRanks = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=T);
medRanks = presRanks$medRanks;
