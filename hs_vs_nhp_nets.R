rm(list=ls()); options(stringsAsFactors=F); 
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
library(RDAVIDWebService); library(biomaRt);
library(WGCNA); allowWGCNAThreads();

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
collectGarbage();

# remove genes with zeros across samples
zthresh = 1;
dat.hs = dat.hs0[apply(dat.hs0, 1, function(f) sum(f==0)) < zthresh, ];
dat.nhp = dat.nhp0[apply(dat.nhp0, 1, function(f) sum(f==0)) < zthresh, ];

dat.hs = dat.hs[rownames(dat.hs) %in% rownames(dat.nhp), ];
dat.nhp = dat.nhp[rownames(dat.nhp) %in% rownames(dat.hs), ];


# center genes to all have equivalent means
x = as.data.frame(t(apply(dat.hs, 1, function(f) 200 + as.numeric(f) - mean(as.numeric(f)))));
names(x)=names(dat.hs);
dat.hs = x; rm(x); 
dat.hs = as.data.frame(t(scale(t(dat.hs)))) + 10;

# # x = as.data.frame(t(apply(dat.nhp, 1, function(f) 200 + as.numeric(f) - mean(as.numeric(f)))));
# # names(x)=names(dat.nhp);
# # dat.nhp = x; rm(x);
dat.nhp = as.data.frame(t(scale(t(dat.nhp)))) + 10;

######## get variance ranks in each tissue type
ttypes = unique(substr(gsub('hsa.','',names(dat.hs)),1,2));
for (ttype in ttypes) {
	assign(paste(ttype, 'Var', sep=''), rank(apply(dat.hs[,grepl(ttype,names(dat.hs))], 1, var)));
}; rm(ttype);

avgVarBrain = c();
for (g in 1:length(brVar)) {
	avgVarBrain = c(avgVarBrain, mean(c(brVar[g], cbVar[g])));
}; rm(g);
names(avgVarBrain) = names(brVar);

avgVarNonBrain = c();
for (g in 1:length(brVar)) {
	avgVarNonBrain = c(avgVarNonBrain, mean(c(htVar[g], kdVar[g], lvVar[g], tsVar[g])));
}; rm(g);
names(avgVarNonBrain) = names(brVar);

higherVar = avgVarNonBrain - avgVarBrain;
dat.hs = dat.hs[higherVar < 0, ];

dat.nhp = dat.nhp[rownames(dat.nhp) %in% rownames(dat.hs), ];

# # filter out genes with higher variance across non-brain samples
# brainVar = apply(dat.hs[,grepl('br|cb',names(dat.hs))], 1, var);
# nonBrainVar = apply(dat.hs[,!grepl('br|cb',names(dat.hs))], 1, var);
# hsHigherVar = rank(nonBrainVar) - rank(brainVar);
# dat.hs = dat.hs[hsHigherVar < 0, ];

# remove outliers and quantile normalize
keepMe = names(dat.hs)!='hsa.br.M.4' & names(dat.hs)!='hsa.kd.F.1';

X = dat.hs[,keepMe];
#X = dat.hs
out = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);

DATA = as.data.frame(t(out$data_Qnorm));
#BLOCKSIZE = 1000;
#TYPE = 'signed';
#sft = pickSoftThreshold(DATA, networkType=TYPE, verbose=3, blockSize=BLOCKSIZE);
sft = exn.plotPowers(DATA, networkType='signed', blockSize=ncol(DATA)+1, verbose=3);

POWER = 18;
# DS = 2;
# MM = 10;
# MCH = 0.15;
# NPERM = 10000
# ST = 100;
# k = softConnectivity(DATA, type=TYPE, power=POWER, blockSize=BLOCKSIZE);
# par(mfrow=c(1,2));
# scaleFreePlot(k); hist(k,col='grey',border='darkgrey');

k = exn.computeConnectivityAndPlotScaleFreeness(DATA, networkType='signed', power=18);
	
# net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=(ncol(DATA)+1), power=POWER, 
							   # deepSplit=DS, minModuleSize=MM, mergeCutHeight=MCH,
							   # densityPermTest=T, skipThresh=ST, nPerm=NPERM,
							   # saveFileBase=paste('hs_all-tissue_',TYPE,'_p',POWER,'_ds',DS,'_mm',MM,'_mch',MCH,'_st',ST,'_perm',NPERM,sep=''),
							   # verbose=3
							   # );
							   
POWER = 18;
TYPE = 'signed';
for (DS in c(2,4)) {
	for (MM in c(20, 40, 80, 100)) {
		for (MCH in c(.15, .2)) {
			net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=(ncol(DATA)+1), power=POWER, networkType=TYPE,
			deepSplit=DS, minModuleSize=MM, mergeCutHeight=MCH,
			densityPermTest=F, verbose=3,
			saveFileBase=paste('hs_all-tissueSCALED_',TYPE,'_p',POWER,'_ds',DS,'_mm',MM,'_mch',MCH,sep='')
			);
		}
	}
}							   
							   
							   

							   
							   
####################					   
							   
							   
							   
							   
							   



# unsigned
# # # sft = pickSoftThreshold(DATA, networkType='unsigned', verbose=3, blockSize=BLOCKSIZE);
   # # # # Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# # # # 1      1    0.329  2.830          0.988  2580.0   2530.00   3700
# # # # 2      2    0.026  0.368          0.988  1160.0   1130.00   2150
# # # # 3      3    0.076 -0.484          0.986   626.0    600.00   1410
# # # # 4      4    0.298 -0.896          0.983   381.0    355.00    995
# # # # 5      5    0.470 -1.120          0.966   251.0    226.00    734
# # # # 6      6    0.577 -1.210          0.948   176.0    151.00    563
# # # # 7      7    0.648 -1.210          0.919   129.0    105.00    445
# # # # 8      8    0.835 -1.110          0.992    98.0     75.10    360
# # # # 9      9    0.869 -1.280          0.995    76.7     55.20    328
# # # # 10    10    0.891 -1.390          0.990    61.5     41.40    303
# # # # 11    12    0.922 -1.520          0.981    41.9     24.30    263
# # # # 12    14    0.938 -1.560          0.966    30.2     15.20    234
# # # # 13    16    0.948 -1.580          0.962    22.8      9.95    212
# # # # 14    18    0.955 -1.580          0.959    17.8      6.73    193
# # # # 15    20    0.958 -1.570          0.955    14.3      4.69    180

# # # POWER = 9
# # # k = softConnectivity(DATA, type='unsigned', power=POWER, blockSize=BLOCKSIZE);
# # # par(mfrow=c(1,2));
# # # scaleFreePlot(k); hist(k, breaks=length(k),col='grey',border='darkgrey');
# # # net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=(ncol(DATA)+1), power=POWER, networkType='unsigned',
							   # # # deepSplit=2, minModuleSize=40, mergeCutHeight=.2,
							   # # # densityPermTest=F, 
							   # # # saveFileBase=paste('hs_all-tissue_p',POWER,'_ds2_mm40_mch.2_UN',sep='')
							   # # # );
# # # rm(list=ls());
# # # load('hs_all-tissue_p9_ds2_mm40_mch.2_UNrun10NET.RData')
# # # load('hs_all-tissue_p9_ds2_mm40_mch.2_UNrun10DATA.RData')
#################
################

setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls()); options(stringsAsFactors=F);
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
library(RDAVIDWebService); library(biomaRt);
library(WGCNA); allowWGCNAThreads();
#stem = 'hs_all-tissue_signed_p18_ds2_mm40_mch0.15run6';
stem = 'hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10'
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));



#net= blockwiseModules(datExpr=DATA, maxBlockSize=(ncol(DATA)+1), power=18, deepSplit=2, minModuleSize=20, mergeCutHeight=.15,verbose=3,networkType='signed')

NET=net;rm(net)
DATA=DATA;
dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

exn.plotDendroAndColors(dendro, colors, block=block, blockGenes=blockGenes);
exn.plotEigengeneNetworks2(MEs)

kME = exn.computekME(DATA, MEs)$all;
mod.genes=exn.getModuleGenes(DATA, colors);


traits = data.frame(br=as.numeric(grepl('br', rownames(DATA))));
rownames(traits) = rownames(DATA);

traits = cbind(traits, cb=as.numeric(grepl('cb', rownames(DATA))));
traits = cbind(traits, br_cb=as.numeric(grepl('br|cb', rownames(DATA))));
traits = cbind(traits, ht=as.numeric(grepl('ht', rownames(DATA))));
traits = cbind(traits, kd=as.numeric(grepl('kd', rownames(DATA))));
traits = cbind(traits, lv=as.numeric(grepl('lv', rownames(DATA))));
traits = cbind(traits, ts=as.numeric(grepl('ts', rownames(DATA))));


ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits = cbind(traits, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATA)), gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);

traitCors = exn.computeAndPlotMETraitCors(traits, MEs, main=stem);

GS = exn.computeGS(traits, DATA);

.getModGenesRankedBykME = function(module_names,colors,kME) {
	outList = list();
	for (m in 1:length(module_names)) {
		outList[[m]] = exn.getModulekME(module_names[m],colors,kME)
	}
	names(outList) = module_names;
	return(outList);
}
modkMEs = .getModGenesRankedBykME(names(table(colors)),colors,kME);


length(table(colors))
MFROW = c(4,9);
factors = unlist(strsplit(rownames(DATA),'\\.'))[seq(2,length(unlist(strsplit(rownames(DATA),'\\.'))),4)];
par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(factors), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=.8,cex.main=.8);
}

factors2 = factors; factors2[6:7] = 'br';
par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(factors2), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=1);
}

par(mfrow=MFROW)
for(i in 1:ncol(MEs)) {
	verboseBoxplot(MEs[,i], as.factor(c(rep('br/cb',7), rep('rest',9))), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i],cex.axis=1);
}


###############################

tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
primates_hiv1_interactions = primates_hiv1_interactions[primates_hiv1_interactions==1]
checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA));
hiv = checkGeneListEnrichmentList(names(primates_hiv1_interactions),mod.genes,names(DATA))$pvals;hiv

tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);
mammals_pubmed_interactions = mammals_pubmed_interactions[mammals_pubmed_interactions==1]
pubmed = checkGeneListEnrichmentList(names(mammals_pubmed_interactions),mod.genes,names(DATA))$pvals;pubmed



# dir = '/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/';
# files = list.files(dir)[grepl('noHead$', list.files(dir))];
# sets = list();
# for (f in 1:length(files)) {
	# print(files[f]);
	# sets[[f]] = read.table(paste(dir, files[f], sep=''));
	# names(sets)[f] = files[f];
# }; rm(f);

# # convert gene symbols to ensembl ids
# mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
# setsEns = list();
# for (s in 1:length(sets)) {
	# print(names(sets)[s]);
	# setsEns[[s]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=sets[[s]][,1], mart=mart)[,1];
	# names(setsEns)[s] = names(sets)[s];
# }; rm(s);
#save(setsEns,file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
load(file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
disease = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease[[s]] = checkGeneListEnrichmentList(setsEns[[s]], mod.genes, names(DATA))$pvals
	names(disease)[s] = names(setsEns)[s];
}; rm(s);lapply(disease, head)


common_disease_genes = intersect(intersect(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", setsEns$Autistic.noHead), setsEns$Bipolar.noHead) ;

checkGeneListEnrichmentList(common_disease_genes, mod.genes, names(DATA))$pvals

# dir = '_cahoyMaterials';
# files = paste(dir, list.files(dir)[grepl('csv$', list.files(dir))], sep='/');
# ct = list();
# for (f in 1:length(files)) {
	# ct[[f]] = read.csv(files[f]);
	# names(ct)[f] = gsub(dir, '', gsub('.csv', '', files[f]));
# }; rm(f);
# ct = ct[unlist(lapply(ct,nrow))!=0];
# ct0 = ct;

# ct = ct0;
# THRESH = 5;
# ct = lapply(ct, function(f) f[f[,2]>THRESH, ]);

# ctEns = list();
# for (ctl in 1:length(ct)) {
	# print(names(ct)[ctl]);
	# ctEns[[ctl]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=ct[[ctl]][,1], mart=mart)[,1];
	# names(ctEns)[ctl] = names(ct)[ctl];
# }; rm(ctl);
# save(ctEns,file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData');

load(file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData')
celltype = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype[[s]] = checkGeneListEnrichmentList(ctEns[[s]], mod.genes, names(DATA))$pvals
	names(celltype)[s] = names(ctEns)[s];
}; rm(s); lapply(celltype, head);


ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
selection = checkGeneListEnrichmentList(ps.sig[,1], mod.genes, names(DATA))$pvals;selection

kscores0 = read.table('~/Downloads/primates_positive_selection',sep='\t',header=T,row.names=1);
kscores = kscores0[,16];
names(kscores) = rownames(kscores0);
kscores = kscores[names(kscores) %in% names(DATA)];
kscores = kscores[match(names(DATA),names(kscores))];
names(kscores) = names(DATA);

.plotModkMEAndkscore = function(kME, colors, kscores, module, thresh=0, ...) {
	modkME = exn.getModulekME(module,colors,kME);
	modkscore = kscores[names(kscores) %in% rownames(modkME)];	
	modkscore = modkscore[match(rownames(modkME),names(modkscore))];#print(head(modkscore))
	modkscore[is.na(modkscore)] = 0;
	modkscore = modkscore[modkscore>thresh];#print(head(modkscore))
	verboseScatterplot(modkME[rownames(modkME) %in% names(modkscore), 1], modkscore,frame.plot=F,abline=T,xlab=names(modkME)[1], ylab='kscores', col=module,pch=20, ...);
	return(modkscore);
}

par(mfrow=MFROW)
for ( m in 1:length(names(table(colors)))) {
	tmp = .plotModkMEAndkscore(kME,colors,kscores,names(table(colors))[m])
}

ks.sig = kscores[kscores>quantile(kscores,.99,na.rm=T)];
ks.sig = ks.sig[!is.na(ks.sig)];
selection_k = checkGeneListEnrichmentList(names(ks.sig), mod.genes, names(DATA))$pvals;selection_k

################################

MODS = names(table(colors));
d = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d, inputIds=names(DATA), idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
for (m in 1:length(MODS)) {
	print(MODS[m]);
	modGenes = names(DATA)[colors==MODS[m]];
	addList(d, inputIds=modGenes, idType='ENSEMBL_GENE_ID', listName=MODS[m], listType='Gene');
}; rm(m, modGenes);
setCurrentBackgroundPosition(d, which(getBackgroundListNames(d) == 'BG'));
setCurrentGeneListPosition(d, 1);
modDAVID = exn.getModChartsFromDAVID(d);

.getGenesAvgkMEInMod = function(genes, modkME) {
	return(mean(modkME[rownames(modkME) %in% genes, 1]));
}
.addkMEColToDAVID = function(modDAVIDmod, modkME) {
	modDAVIDmod = cbind(modDAVIDmod, avgkME=1:nrow(modDAVIDmod));
	for (r in 1:nrow(modDAVIDmod)) {
		modDAVIDmod[r, ]$avgkME = .getGenesAvgkMEInMod(.getTermGenes(modDAVIDmod,r), modkME);
	}
	return(modDAVIDmod);
}
.addkMEColsToDAVIDList = function(modDAVID, modkME, order=T) {
	for (m in 1:length(modDAVID)) {
		if (nrow(modDAVID[[m]])==0) {
			modDAVID[[m]] = modDAVID[[m]];
		} else {
			modDAVID[[m]] = .addkMEColToDAVID(modDAVID[[m]], modkME[[m]]);
			if (order) {
				modDAVID[[m]] = modDAVID[[m]][order(modDAVID[[m]]$avgkME, decreasing=T), ];
			}
		}
	}
	return(modDAVID);
}
modDAVID = .addkMEColsToDAVIDList(modDAVID,modkMEs);

.getTermGenes = function(modDAVIDmod,term) {
	if (is.numeric(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[term], ', ')));
	} else if (is.character(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[modDAVIDmod$Term==term], ', ')));
	} else {
		stop()
	}
}
# assumes nums has gene names
.getAvgNumForGenes = function(genes, nums) {
	return(mean(nums[names(nums) %in% genes], na.rm=T));
}

.addGeneAvgColToDAVID = function(modDAVIDmod, nums, colname='colname') {
	modDAVIDmod = cbind(modDAVIDmod, colname=1:nrow(modDAVIDmod));
	for (row in 1:nrow(modDAVIDmod)) {
		g = .getTermGenes(modDAVIDmod, row);
		modDAVIDmod[row, ncol(modDAVIDmod)] = .getAvgNumForGenes(g, nums);
	}
	names(modDAVIDmod)[ncol(modDAVIDmod)] = colname;
	return(modDAVIDmod)
}

for(m in 1:length(modDAVID)) {
	modDAVID[[m]] = .addGeneAvgColToDAVID(modDAVID[[m]],kscores,'kscores');
}; rm(m)






#lapply(modDAVID,function(f) head(f[, c(1,2,3,4,5,11:13)]));
#lapply(modDAVID,function(f) head(f[grepl('GOTERM',f$Category), c(1,2,3,4,5,11:13)]));

for (m in 1:length(modDAVID)) {
	modDAVID[[m]][,13] = as.numeric(modDAVID[[m]][,13]);
}; rm(m);
modDAVIDf = exn.filterModDAVIDList(modDAVID,13,10,T);
lapply(modDAVIDf,function(f) (f[1:30,c(1,2,3,4,5,11:ncol(f))]));

save(modDAVID,modDAVIDf,file='/Volumes/fishstudies/_mammalian_RNAseq/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_DAVID.RData');

terms = exn.modDAVIDByTerm(modDAVID);
termsU = unlist(terms$uniqueToMod);

sort(table(termsU));
sort(table(termsU) / table(unlist(terms[[1]])));

PAT = 'synap|nerve|neuro|neura'
sort(table(unlist(terms[[1]][grepl(PAT,names(terms[[1]]))])));
sort(table(termsU[grepl(PAT,names(termsU))]));


#############
MOD = 'turquoise';
TERM = 'GO:0016570~histone modification'

MODDAVID = modDAVIDf[names(modDAVIDf)==MOD][[1]]
MODKME = modkMEs[names(modkMEs)==MOD][[1]];
xg = .getTermGenes(MODDAVID, TERM);
xk = MODKME[rownames(MODKME) %in% xg, ]
xkk = kscores[names(kscores) %in% xg];
tmp = cbind(xk, kscore=xkk[match(rownames(xk),names(xkk))]);

par(mfrow=c(1,3));
verboseScatterplot(tmp[,1], tmp[,3], abline=T,xlab='kME',ylab='kscore')
verboseBoxplot(tmp[,1], as.factor(rownames(tmp) %in% names(mammals_pubmed_interactions)),xlab='',ylab='kME')
verboseBoxplot(tmp[,3], as.factor(rownames(tmp) %in% names(mammals_pubmed_interactions)),xlab='',ylab='kscore')

###################################

modPvals = hiv;
modPvals = cbind(modPvals, pubmed=pubmed[match(rownames(modPvals), rownames(pubmed)), 2]);
rownames(modPvals) = modPvals[,1];
modPvals = modPvals[, -1];
names(modPvals)[1] = 'hiv';

for(s in 1:length(disease)) {
	tmp = disease[[s]];
	tmp = tmp[match(rownames(modPvals), tmp[,1]),];
	modPvals = cbind(modPvals, tmp[,2]);
}; rm(s);
names(modPvals)[3:ncol(modPvals)] = names(disease);

for(s in 1:length(celltype)) {
	tmp = celltype[[s]];
	tmp = tmp[match(rownames(modPvals), tmp[,1]),];
	modPvals = cbind(modPvals, tmp[,2]);
}; rm(s);
names(modPvals)[14:ncol(modPvals)] = names(celltype);

modPvals = cbind(modPvals, ps=selection[match(rownames(modPvals),selection[,1]), 2]);

modPvals = matrix(as.numeric(as.matrix(modPvals)), ncol=ncol(modPvals),dimnames=dimnames(modPvals));



tmp = matrix(as.numeric(modPvals<(.05/(sum(lower.tri(cor(modPvals)))))),ncol=ncol(modPvals))
#tmp = matrix(as.numeric(modPvals<.01),ncol=ncol(modPvals)))
dimnames(tmp)=dimnames(modPvals)
heatmap(tmp,Rowv=NA,Colv=NA,scale='none',main='');

heatmap(cor(modPvals),main='',symm=T);

########################################




.modExprHeatmap = function(module, DATA, colors, orderByExpr=T, ...) {
	modExpr = as.matrix(DATA[, colors==module]);
	if (orderByExpr) {
		modExpr = modExpr[, order(apply(modExpr, 2, mean, na.rm=T),decreasing=T)];
	}
	heatmap(t(as.matrix(modExpr)), Rowv=NA, Colv=NA, ...)
}






.checkTermGeneEnrichment = function(modDAVIDmod, row, genes) {
	termGenes = unlist(strsplit(modDAVIDmod[row,]$Genes, ', '));
	return(checkGeneListEnrichment(genes, termGenes, unique(unlist(strsplit(modDAVIDmod$Genes, ', ')))));
}

.printTermGeneEnrichment = function(modDAVIDmod, genes, thresh=.05) {
	for(row in 1:nrow(modDAVIDmod)){  
		tmp = .checkTermGeneEnrichment(modDAVIDmod,row,genes); 
		if(tmp[[2]]$p.value<thresh){
			print(modDAVIDmod[row,c(1,2,3,4,5,11:ncol(modDAVIDmod))]);
			print(tmp)
		}   
	}
}

####################################


.plotTermCors = function(modDAVIDmod,cols=c(14,15), main='',bg='',col='black') {
	verboseScatterplot(modDAVIDmod[,cols[1]], modDAVIDmod[,cols[2]],abline=T,xlab=colnames(modDAVIDmod)[cols[1]],ylab=colnames(modDAVIDmod)[cols[2]],frame.plot=F,abline.col='red',pch=20,main=main,col=col)
}

par(mfrow=c(3,5));
for(m in 1:length(modDAVIDf)){
	if(nrow(modDAVIDf[[m]])<10) {next}
	.plotTermCors(modDAVIDf[[m]],main=names(modDAVIDf)[m],col=names(modDAVIDf)[m])
}



####################


hivgenes = names(primates_hiv1_interactions);


par(mfrow=MFROW)
for(m in 1:length(modkMEs)) {
	verboseBoxplot(modkMEs[[m]][,1], as.factor(rownames(modkMEs[[m]]) %in% hivgenes), 
				   xlab='', ylab=names(modkMEs[[m]])[1], col=names(modkMEs)[m], main=names(modkMEs)[m]
				   );
}; rm(m)


par(mfrow=MFROW)
for(m in 1:length(mod.genes)) {
	verboseBoxplot(kscores[names(kscores) %in% mod.genes[[m]]], as.factor(names(kscores[names(kscores) %in% mod.genes[[m]]]) %in% hivgenes), 
				   xlab='', ylab='kscores', col=names(mod.genes)[m], main=names(mod.genes)[m]
				   );
}; rm(m)

par(mfrow=MFROW)
for (mod in names(mod.genes)) {
	MOD = mod;
	modkME0 = modkMEs[names(modkMEs)==MOD][[1]];
	modkME = rank(modkME0[, 1]);
	modkscores = kscores[match(rownames(modkME0), names(kscores))];
	modkME = modkME[!is.na(modkscores)];
	modkscores = modkscores[!is.na(modkscores)];
	modkscores = rank(modkscores);
	verboseBoxplot(apply(cbind(modkME,modkscores), 1, mean), as.factor(names(apply(cbind(modkME,modkscores), 1, mean)) %in% hivgenes), 
					   xlab='', ylab='', col=MOD, main=''
					   )
}


source('/Volumes/fishstudies/_code/circlePlot.R');
MOD = 'blue';
TERM = 'GO:0045202~synapse';
MODDAVID = modDAVID;

termgenes = .getTermGenes(MODDAVID[names(MODDAVID)==MOD][[1]], TERM);
termgenesDATA = DATA[, names(DATA) %in% termgenes];

termgenesADJ = adjacency(termgenesDATA, type='signed', power=18);



#MOD = '';
#ADJ = adjacency(DATA[,colors==MOD],type='signed',power=18);
ADJ = termgenesADJ
COL = rep(MOD, nrow(ADJ));
COL[rownames(ADJ) %in% hivgenes] = 'black'
circlePlot(ADJ,
		   colnames(ADJ),
		   order(-apply(ADJ,1,sum)),
		   lineColors=grey2red(50,.9,1),
		   pointBg=COL,
		   max.cex.labels=1.5,
		   min.cex.labels=0.75
		   )


# min.line.width=1;
# max.line.width=10;
# xLabelOffset=0.01;
# yLabelOffset=0;
# min.cex.labels=0.75;  #use min.cex.labels=1 for smaller plots
# max.cex.labels=1.5;
# min.cex.points=1;
# max.cex.points=3;
# circlePlot(bluePDadj[[1]],
			# bluePDadj[[2]][,2],
			# order(-bluePDk),
			# min.line.width=min.line.width,
			# max.line.width=max.line.width, 
			# lineColors=grey2red(50,.9,1),
			# min.cex.labels=min.cex.labels,
			# max.cex.labels=max.cex.labels,
			# xLabelOffset=xLabelOffset,
			# yLabelOffset=yLabelOffset,
			# pointBg=pointBg,
			# pointColors=pointColors,
			# min.cex.points=min.cex.points,
			# max.cex.points=max.cex.points
			# );

### yellow module: common_disease_genes enrichment
# "ENSG00000152413" "ENSG00000091129" "ENSG00000148053" "ENSG00000125675" : HOMER1, NCAM, NTRK2, GRIA3



save.image(file=paste('WORKSPACE_', stem, '_AUG13.RData', sep=''));
#######################
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];
dat.nhp0 = dat[, grep('ptr|ppa|ggo', names(dat))];
rownames(dat.nhp0) = ids$hsa;

dat.nhp = dat.nhp0[rownames(dat.nhp0) %in% names(DATA), ];

zthresh=1
dat.nhp = dat.nhp[apply(dat.nhp, 1, function(f) sum(f==0)) < zthresh, ];


x = as.data.frame(t(apply(dat.nhp, 1, function(f) 200 + as.numeric(f) - mean(as.numeric(f)))));
names(x)=names(dat.nhp);
dat.nhp = x; rm(x);

source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
zthresh = 1;
out = preProcess(datIN=dat.nhp, 
				 removeOutlierProbes=F, deviate=3, 
				 removeTooManyNAs=F, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=F, IACthresh=2, 
				 Qnorm=T);
				 
				 
				 
DATAp = as.data.frame(t(out$data_Qnorm));			 
				 

							   
net= blockwiseModules(datExpr=DATAp, maxBlockSize=(ncol(DATAp)+1), power=18, deepSplit=2, minModuleSize=20, mergeCutHeight=.15,verbose=3,networkType='signed')

NETp=net;rm(net)
dendrop=NETp$dendrograms;
blockp = 1;
blockGenesp = NETp$blockGenes;
MEsp = NETp$MEs;

NETp$colors = matchLabels(NETp$colors, colors);
colorsp = NETp$colors;

MEsp=moduleEigengenes(DATAp,colorsp)$eigengenes




exn.plotEigengeneNetworks2(MEsp)

kMEp = exn.computekME(DATAp, MEsp)$all;
mod.genesp=exn.getModuleGenes(DATAp, colorsp);






exn.plotDendroAndColors(dendrop, colorsp, block=blockp, blockGenes=blockGenesp);

plotDendroAndColors(dendro[[block]], data.frame(colors[blockGenes[[1]]], colorsp[blockGenesp[[1]]]), dendroLabels=F, groupLabels=c('hs','rest'));
exn.plotModuleOverlaps(colors,colorsp);


traitsp = data.frame(br=as.numeric(grepl('br', rownames(DATAp))));
rownames(traitsp) = rownames(DATAp);

traitsp = cbind(traitsp, cb=as.numeric(grepl('cb', rownames(DATAp))));
traitsp = cbind(traitsp, br_cb=as.numeric(grepl('br|cb', rownames(DATAp))));
traitsp = cbind(traitsp, ht=as.numeric(grepl('ht', rownames(DATAp))));
traitsp = cbind(traitsp, kd=as.numeric(grepl('kd', rownames(DATAp))));
traitsp = cbind(traitsp, lv=as.numeric(grepl('lv', rownames(DATAp))));
traitsp = cbind(traitsp, ts=as.numeric(grepl('ts', rownames(DATAp))));
traitsp = cbind(traitsp, ptr=as.numeric(grepl('ptr', rownames(DATAp))));
traitsp = cbind(traitsp, ppa=as.numeric(grepl('ppa', rownames(DATAp))));
traitsp = cbind(traitsp, ggo=as.numeric(grepl('ggo', rownames(DATAp))));
# ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
# traits = cbind(traits, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATA)), gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);

traitCorsp = exn.computeAndPlotMETraitCors(traitsp, MEsp, main='');