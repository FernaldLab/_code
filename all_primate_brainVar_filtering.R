rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');

# read in raw data
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];
rownames(dat) = ids$hsa;

# remove unwanted samples
remove = 'ts|mml|ppy';
dat = dat[, !grepl(remove,names(dat))];

# try removing hsa.br.M.4 and ppa.lv.F.1 because of sample dendrogram
dat = dat[, -c(4,40,41)];

# remove genes with zeros across samples
zthresh = 1;
dat = dat[apply(dat, 1, function(f) sum(f==0)) < zthresh, ];

# remove singular tissue samples
nums = table(substr(names(dat),1,6));
remove = names(nums)[nums==1];
dat = dat[, !(substr(names(dat),1,6) %in% remove)];

# remove ts and ppy and mml samples
# # # keep = !grepl('ts|mml|ppy', names(dat));
# # # dat = dat[, keep];

# center genes to all have equivalent means
x = as.data.frame(t(apply(dat, 1, function(f) 200 + as.numeric(f) - mean(as.numeric(f)))));
names(x)=names(dat);
dat = x; rm(x); 

save(dat, file='allData_NOts_NOmmlppy.RData');
save(dat, file='allData_NOts_NOmmlppy_removed_earlier_USETHIS.RData');
save(dat, file='allData_NOts_NOmmlppy_removed_earlier_USETHIS_removed_hsabrM4_ppalvF1.RData');
###############################################
###############################################

load('allData_NOts_NOmmlppy_removed_earlier_USETHIS_removed_hsabrM4_ppalvF1.RData');
# get variance ranks in each tissue type
species = unique(substr(names(dat),1,3));

# compute species/tissue specific variance ranks
for (s in species) {
	print(s)
	datTMP = dat[, grep(s, names(dat))];
	ttypes = unique(substr(gsub(paste(s,'.',sep=''),'',names(datTMP)),1,2)); print(ttypes)
	for (ttype in ttypes) {
		assign(paste(s, '.', ttype, 'Var', sep=''), rank(apply(datTMP[,grepl(ttype,names(datTMP))], 1, var)));
	}
}

DATA = dat;
rm(list=ls()[!grepl('Var|species|DATA',ls())]);


#nonhsa.avgBrainVar = (ggo.brVar + ppa.brVar + ptr.brVar) / 3;
#hsaVSnonhsa.higherBrainVar = nonhsa.avgBrainVar - hsa.brVar;

# compute avg brain variance ranks
for (s in species) {
	brVar = paste(s,'.brVar',sep='');
	cbVar = paste(s,'.cbVar',sep='');
	if (exists(brVar) & exists(cbVar)) {
		assign(paste(s,'.avgBrainVar',sep=''), (get(brVar)+get(cbVar))/2);
	} else if (exists(brVar)) {
		assign(paste(s,'.avgBrainVar',sep=''), get(brVar));
	} else {
		assign(paste(s,'.avgBrainVar',sep=''), get(cbVar));
	}
}; rm(s,brVar,cbVar);

nonhsa.avgBrainVar = (ggo.avgBrainVar + ppa.avgBrainVar + ptr.avgBrainVar) / 3;
hsaVSnonhsa.higherBrainVar = nonhsa.avgBrainVar - hsa.avgBrainVar;

DATA.hsaVSnonhsa.higherBrainVar = DATA[hsaVSnonhsa.higherBrainVar < 0, ];
save(DATA.hsaVSnonhsa.higherBrainVar, file='DATA.hsaVSnonhsa.higherBrainVar_removed_hsabrM4_ppalvF1.RData');

# # # compute avg non-brain variance ranks
# # for (s in species) {
	# # assign(paste(s, '.avgNonBrainVar', sep=''), 
		  # # ( get(paste(s,'.htVar',sep=''))+get(paste(s,'.kdVar',sep=''))+get(paste(s,'.lvVar',sep='')) ) / 3
		  # # );
# # }; rm(s);

# # # compute variance differences for each species
# # for (s in species) {
	# # nonbr = get(paste(s, '.avgNonBrainVar', sep=''));
	# # br = get(paste(s, '.avgBrainVar', sep=''));
	# # assign(paste(s, '.higherBrainVar', sep=''),  nonbr - br);
# # }; rm(s,nonbr,br);

# # filter data for each species
# for (s in species) {
	# spDATA = DATA[, grepl(s, names(DATA))];
	# spVAR = get(paste(s, '.higherBrainVar', sep=''));
	# assign(paste('DATA.', s, '.higherBrainVar', sep=''), spDATA[spVAR < 0, ]);
# }; rm(s,spDATA,spVAR);

# # get genes common to all species
# tmp = ls()[grepl('DATA.',ls())];
# commonGenes = rownames(get(tmp[1]));
# for (s in tmp[2:length(tmp)]) {
	# tmpGenes = rownames(get(s));
	# commonGenes = commonGenes[commonGenes %in% tmpGenes];
# }; rm(s,tmpGenes);
# print(length(commonGenes))

# # filter species data to common genes
# for (s in tmp) {
	# assign(paste(s, '.Common', sep=''), get(s)[rownames(get(s)) %in% commonGenes, ]);
# }; rm(s);

# tmp = ls()[grepl('.Common',ls())];
# DATA.Common = get(tmp[1]);
# for (s in tmp[2:length(tmp)]) {
	# DATA.Common = cbind(DATA.Common, get(s));
# }; rm(s);

# tmp = ls()[grepl('DATA.',ls())];
# save(list=tmp, file='all_primate_brainVar_filtering_speciesDATA.RData');

######

# hsa.higherBrainVar = hsa.avgNonBrainVar - hsa.avgBrainVar;
# DATA.hsa.higherBrainVar = DATA[hsa.higherBrainVar < 0, ];
# save(DATA.hsa.higherBrainVar, file='all_primate_brainVar_filtering_DATA.hsa.higherBrainVar.RData');

######

# compute overall average brain variance
# tmp = ls()[grepl('avgBrainVar',ls())];
# total = get(tmp[1]);
# for (s in tmp[2:length(tmp)]) {
	# total = total + get(s);
# }; rm(s);
# all.avgBrainVar = total / length(tmp);
# rm(tmp);

# # compute overall average non-brain variance
# tmp = ls()[grepl('avgNonBrainVar',ls())];
# total = get(tmp[1]);
# for (s in tmp[2:length(tmp)]) {
	# total = total + get(s);
# }; rm(s);
# all.avgNonBrainVar = total / length(tmp);
# rm(tmp);

# all.higherBrainVar = all.avgNonBrainVar - all.avgBrainVar;
# DATA.all.higherBrainVar = DATA[all.higherBrainVar < 0, ];

# save.image(file='all_primate_brainVar_filtering_WORKSPACE.RData');
# save(DATA.all.higherBrainVar, file='all_primate_brainVar_filtering_DATA.all.higherBrainVar.RData');



###################
#####################
rm(list=ls());options(stringsAsFactors=F); 
library(WGCNA);
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
load(file='DATA.hsaVSnonhsa.higherBrainVar_removed_hsabrM4_ppalvF1.RData');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

DATA = DATA.hsaVSnonhsa.higherBrainVar;
out = preProcess(datIN=DATA)

DATA = as.data.frame(t(out$data_Qnorm));
sft = pickSoftThreshold(DATA, networkType='signed', verbose=3);
k = softConnectivity(DATA, type='signed', power=14, blockSize=ncol(DATA));

par(mfrow=c(1,2));
scaleFreePlot(k); hist(k, breaks=length(k),col='grey',border='darkgrey');

# # # > sft = pickSoftThreshold(DATA, networkType='signed', verbose=3);
# # # pickSoftThreshold: will use block size 4102.
 # # # pickSoftThreshold: calculating connectivity for given powers...
   # # # ..working on genes 1 through 4102 of  4102
   # # # Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# # # 1      1   0.0260  2.64          0.932 2100.00   2100.00 2340.0
# # # 2      2   0.0358 -1.97          0.867 1170.00   1180.00 1460.0
# # # 3      3   0.0816 -1.90          0.888  696.00    703.00  965.0
# # # 4      4   0.2000 -2.38          0.916  438.00    440.00  694.0
# # # 5      5   0.2700 -2.13          0.921  288.00    285.00  519.0
# # # 6      6   0.3520 -2.01          0.932  197.00    192.00  402.0
# # # 7      7   0.4400 -1.96          0.942  139.00    133.00  319.0
# # # 8      8   0.5220 -1.86          0.957  100.00     93.40  259.0
# # # 9      9   0.5920 -2.01          0.959   74.40     67.30  219.0
# # # 10    10   0.6300 -2.11          0.960   56.20     49.00  187.0
# # # 11    12   0.7040 -2.20          0.973   33.70     27.20  142.0
# # # 12    14   0.7620 -2.23          0.980   21.30     15.80  110.0
# # # 13    16   0.7980 -2.24          0.982   14.10      9.55   88.1
# # # 14    18   0.7930 -2.32          0.974    9.66      5.93   71.6
# # # 15    20   0.8070 -2.29          0.972    6.83      3.79   59.1

########################

source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');
# source('/Volumes/fishstudies/_code/exploreNetwork.R');
# source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');




net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=(ncol(DATA)+1), power=14, deepSplit=2, minModuleSize=10, skipThresh=200, saveFileBase='hsaVSnonhsa.higherBrainVar_removed_hsabrM4_ppalvF1_b14_mm10_mch.2_ds2', mergeCutHeight=.2);


rm(list=ls());options(stringsAsFactors=F); 
library(WGCNA);
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
stem = 'hsaVSnonhsa.higherBrainVar_removed_hsabrM4_ppalvF1_b14_mm10_mch.2_ds2run4';
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));

NET=net; rm(net)

dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

traits = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA))));
rownames(traits) = rownames(DATA);
#ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
#traits = cbind(traits, age.index=ages[match(rownames(DATA), gsub(' ','.',ages$Sample)), 5:7][,2]);
traits = cbind(traits, ptr_rest=as.numeric(grepl('ptr', rownames(DATA))));
traits = cbind(traits, ppa_rest=as.numeric(grepl('ppa', rownames(DATA))));
traits = cbind(traits, ggo_rest=as.numeric(grepl('ggo', rownames(DATA))));

traits = cbind(traits, br=as.numeric(grepl('br', rownames(DATA))));
traits = cbind(traits, cb=as.numeric(grepl('cb', rownames(DATA))));
traits = cbind(traits, ht=as.numeric(grepl('ht', rownames(DATA))));
traits = cbind(traits, kd=as.numeric(grepl('kd', rownames(DATA))));
traits = cbind(traits, lv=as.numeric(grepl('lv', rownames(DATA))));


#traits = cbind(traits, ppy_rest=as.numeric(grepl('ppy', rownames(DATA))));
#traits = cbind(traits, mml_rest=as.numeric(grepl('mml', rownames(DATA))));
#traits = traits[, c(1,3:7,2)];

traitCors = exn.computeAndPlotMETraitCors(traits, MEs, main=stem);
GS = exn.computeGS(traits, DATA);
kME = exn.computekME(DATA, MEs)$all;

TRAIT = 'human_rest';
exn.plotAllModsGSkME(colors,TRAIT,GS,kME,mfrow=c(4,7));

tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
primates_hiv1_interactions = primates_hiv1_interactions[primates_hiv1_interactions==1]
mod.genes=exn.getModuleGenes(DATA, colors);
checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA));
hiv = checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA))$pvals

tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);
mammals_pubmed_interactions = mammals_pubmed_interactions[mammals_pubmed_interactions==1]
pubmed = checkGeneListEnrichmentList(names(mammals_pubmed_interactions),mod.genes,names(DATA))$pvals

load(file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
disease = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease[[s]] = checkGeneListEnrichmentList(setsEns[[s]], mod.genes, names(DATA))$pvals
	names(disease)[s] = names(setsEns)[s];
}; rm(s);lapply(disease, head)


load(file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData')
celltype = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype[[s]] = checkGeneListEnrichmentList(ctEns[[s]], mod.genes, names(DATA))$pvals
	names(celltype)[s] = names(ctEns)[s];
}; rm(s); lapply(celltype, head);


# ps = read.table('journal.pgen.1000840.s009.txt');
# ps.sig = ps[ps$V5<.05, ];
# selection = checkGeneListEnrichmentList(ps.sig[,1], mod.genes, names(DATA))$pvals;selection

# kscores0 = read.table('~/Downloads/primates_positive_selection',sep='\t',header=T,row.names=1);
# kscores = kscores0[,16];
# names(kscores) = rownames(kscores0);
# kscores = kscores[names(kscores) %in% names(DATA)];
# kscores = kscores[match(names(DATA),names(kscores))];
# names(kscores) = names(DATA);

# .plotModkMEAndkscore = function(kME, colors, kscores, module, thresh=0, ...) {
	# modkME = exn.getModulekME(module,colors,kME);
	# modkscore = kscores[names(kscores) %in% rownames(modkME)];	
	# modkscore = modkscore[match(rownames(modkME),names(modkscore))];#print(head(modkscore))
	# modkscore[is.na(modkscore)] = 0;
	# modkscore = modkscore[modkscore>thresh];#print(head(modkscore))
	# verboseScatterplot(modkME[rownames(modkME) %in% names(modkscore), 1], modkscore,frame.plot=F,abline=T,xlab=names(modkME)[1], ylab='kscores', col=module,pch=20, ...);
	# return(modkscore);
# }

# par(mfrow=MFROW)
# for ( m in 1:length(names(table(colors)))) {
	# tmp = .plotModkMEAndkscore(kME,colors,kscores,names(table(colors))[m])
# }

# ks.sig = kscores[kscores>quantile(kscores,.95,na.rm=T)]
# ks.sig = ks.sig[!is.na(ks.sig)]
