setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
library(WGCNA); library(RDAVIDWebService); library(biomaRt);

rm(list=ls());
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
#stem = 'NET.br.s.ds4_mm80_mch.15run12';
stem = 'NET/NET.br.s.ds2_mm40_mch.2run12';
load(paste(stem,'NET.RData', sep=''));
load(paste(stem,'DATA.RData', sep=''));

NET=net; DATA=DATA; rm(net);

dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

exn.plotDendroAndColors(dendro, colors, block=block, blockGenes=blockGenes);

traits = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA))), 
					primate_nonprimate=as.numeric(grepl('hsa|ptr|ppa|ggo|ppy|mml', rownames(DATA))),
					human_nonhumanprimate=c(rep(1,5), rep(0,15), rep(NA,9)),
					male_female=as.numeric(grepl('M', rownames(DATA)))
					);
rownames(traits) = rownames(DATA);

ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits = cbind(traits, age.index=ages[match(rownames(DATA), gsub(' ','.',ages$Sample)), 5:7][,2])

traitCors = exn.computeAndPlotMETraitCors(traits, MEs);
MEnet = exn.plotEigengeneNetworks2(MEs, returnCors=T);
kME = exn.computekME(DATA, MEs)$all;
GS = exn.computeGS(traits, DATA);

TRAIT = 'primate_nonprimate';
exn.plotAllModsGSkME(colors,TRAIT,GS,kME,mfrow=c(3,7));
MOD = 'blue';
exn.plotModuleAllGSkME(MOD,colors,GS,kME);

GENES = rownames(GS[GS$q.GS.human_rest<.05, ]);
d = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d, inputIds=names(DATA), idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
addList(d, inputIds=GENES, idType='ENSEMBL_GENE_ID', listName='GENES', listType='Gene');
setCurrentBackgroundPosition(d, which(getBackgroundListNames(d) == 'BG'));
ch = getFunctionalAnnotationChart(d);
tmp = c();
for (row in 1:nrow(ch)) {
	row_genes = unlist(strsplit(ch$Genes[row], ', '));
	row_mods = colors[match(row_genes, names(DATA))];
	tmp = c(tmp, paste(row_mods, collapse=', '));
}; rm(row, row_genes, row_mods);
ch = data.frame(ch, modules=tmp);
rm(tmp);

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
lapply(modDAVID,function(f) head(f[, c(1,2,3,4,5,11:13)]));
lapply(modDAVID,function(f) head(f[grepl('GOTERM',f$Category), c(1,2,3,4,5,11:13)]));
modDAVIDf = exn.filterModDAVIDList(modDAVID,13,10,T);
lapply(modDAVIDf,function(f) f[,c(1,2,3,4,5,11:13)]);

terms = exn.modDAVIDByTerm(modDAVID);
termsU = unlist(terms$uniqueToMod);



dir = '/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Disease Sets/';
files = list.files(dir)[grepl('noHead$', list.files(dir))];
sets = list();
for (f in 1:length(files)) {
	print(files[f]);
	sets[[f]] = read.table(paste(dir, files[f], sep=''));
	names(sets)[f] = files[f];
}; rm(f);

# convert gene symbols to ensembl ids
mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
setsEns = list();
for (s in 1:length(sets)) {
	print(names(sets)[s]);
	setsEns[[s]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=sets[[s]][,1], mart=mart)[,1];
	names(setsEns)[s] = names(sets)[s];
}; rm(s);


modGenes = list();
for (m in 1:length(table(colors))) {
	modGenes[[m]] = names(DATA)[colors==names(table(colors))[m]];
	names(modGenes)[m] = names(table(colors))[m];
}; rm(m);

for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	print(checkGeneListEnrichmentList(setsEns[[s]], modGenes, names(DATA))$pvals);
}; rm(s)


# get genes common to lists
thresh = 50;
tmpsets = sets[which(sapply(sets,nrow) > thresh)];
cgenes = tmpsets[[1]][, 1];
for (s in 1:length(tmpsets)) {
	cgenes = cgenes[cgenes %in% tmpsets[[s]][,1]];
}; rm(s);
mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
cgenesEns = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=cgenes, mart=mart)[,1]
checkGeneListEnrichmentList(cgenesEns, modGenes, names(DATA))$pvals;




dir = '/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_cahoyMaterials';
files = paste(dir, list.files(dir)[grepl('csv$', list.files(dir))], sep='/');
ct = list();
for (f in 1:length(files)) {
	ct[[f]] = read.csv(files[f]);
	names(ct)[f] = gsub(dir, '', gsub('.csv', '', files[f]));
}; rm(f);

ct0 = ct;

ct = ct0;
THRESH = 5;
ct = lapply(ct, function(f) f[f[,2]>THRESH, ]);

ctEns = list();
for (ctl in 1:length(ct)) {
	print(names(ct)[ctl]);
	ctEns[[ctl]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=ct[[ctl]][,1], mart=mart)[,1];
	names(ctEns)[ctl] = names(ct)[ctl];
}; rm(ctl);

for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	print(checkGeneListEnrichmentList(ctEns[[s]], modGenes, names(DATA))$pvals);
}; rm(s)









####################
####################
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt', header=T, sep='\t');

ids = dat0[, 1:9];
dat = dat0[, 10:141];

dat.cb = dat[, grep('cb',names(dat))];
rownames(dat.cb) = ids$hsa;
rm(dat0,dat,ids);


out = preProcess(datIN=dat.cb, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.cb)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);

DATA.cb = as.data.frame(t(out$data_Qnorm));
DATA.cb = DATA.cb[, match(names(DATA), names(DATA.cb))];
#.ds2_mm40_mch.2run12
NET.cb = blockwiseModules(DATA.cb, power=14, networkType='signed',deepSplit=2,minModuleSize=40,mergeCutHeight=.2);
NET.cb$colors = matchLabels(NET.cb$colors, colors);
dendro.cb=NET.cb$dendrograms;
block.cb = 1;
blockGenes.cb = NET.cb$blockGenes;
colors.cb = NET.cb$colors;

# need to recompute with new colors
MEs.cb=moduleEigengenes(DATA.cb,colors.cb)$eigengenes


exn.plotDendroAndColors(dendro.cb, colors.cb, block=block.cb, blockGenes=blockGenes.cb);

plotDendroAndColors(dendro[[block]], data.frame(colors[blockGenes[[1]]], colors.cb[blockGenes.cb[[1]]]), dendroLabels=F, groupLabels=c('br','cb'));
exn.plotModuleOverlaps(colors,colors.cb);


traits.cb = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA.cb))), 
					primate_nonprimate=as.numeric(grepl('hsa|ptr|ppa|ggo|ppy|mml', rownames(DATA.cb))),
					human_nonhumanprimate=c(rep(1,2), rep(0,9), rep(NA,9)),
					male_female=as.numeric(grepl('M', rownames(DATA.cb)))
					);
rownames(traits.cb) = rownames(DATA.cb);
traitCors.cb = exn.computeAndPlotMETraitCors(traits.cb, MEs.cb);

traits.cb = cbind(traits.cb, age.index=ages[match(rownames(DATA.cb), gsub(' ','.',ages$Sample)), 5:7][,2])


setLabels=c('br','cb');
multiExpr=list();
multiExpr[[1]]=list(data=DATA);
multiExpr[[2]]=list(data=DATA.cb);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=colors;
colorList[[2]]=colors.cb;
names(colorList)=setLabels;

system.time( {				# should only take a minute or two
 	
 	mp = modulePreservation(multiExpr, 
 				   colorList,
 				   networkType = 'signed', 
 				   nPermutations = 0,			# no stats, just get ranks
 				   maxModuleSize = max( table(colors) ),
 				   verbose = 3,
 				   savePermutedStatistics = F
 				   );
 	}
 );


# compute ranks
tmp = cbind(mp$observed$br$intra$br_vs_cb, mp$observed$br$inter$br_vs_cb);
tmp = tmp[!(rownames(tmp)=='gold'), ];
tmp2 = apply(-tmp, 2, rank);
 
# con - cor.kIM, cor.kME, cor.cor
medRankCon = apply( tmp2[ , c(8, 9, 11) ], 1, median);
# density - propVarExplained, meanSignAwareKME, meanSignAwareCorDat, meanAdj
medRankDen = apply(tmp2[ , c(1, 2, 4, 5) ], 1, median);
 
medRankPres = (medRankCon + medRankDen) / 2;

medRanks = cbind(medRankPres,medRankCon,medRankDen);
medRanks = medRanks[order(medRanks[,1]),];


cor(traitCors[[1]], medRanks[match(gsub('ME','',rownames(traitCors[[1]])), rownames(medRanks)),]);
corPvalueFisher(cor(traitCors[[1]], medRanks[match(gsub('ME','',rownames(traitCors[[1]])), rownames(medRanks)),]),nrow(traitCors[[1]]));
verboseScatterplot(traitCors[[1]][,3], medRanks[match(gsub('ME','',rownames(traitCors[[1]])), rownames(medRanks)),1], abline=T, col=gsub('ME','',rownames(traitCors[[1]])),pch=19,cex=3,xlab='MEcor:human_nonhumanprimate',ylab='preservation rank');


##########################

# human specificity
hs = apply(DATA[1:6,],2,mean,na.rm=T)
rest = apply(DATA[7:29,],2,mean,na.rm=T)
fc=cbind(hs,rest)

sum(apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2)
fc[apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2, ]
hgenes = rownames(fc[apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2, ])
###########################

## species permutations

sp = substr(rownames(traits),1,3);
sp.names = unique(sp);

mat = matrix(ncol=nrow(DATA))
for (s in 1:length(sp.names)) {
	testmat = combn(sp.names, s);
	for (col in 1:ncol(testmat)) {
		mat = rbind(mat, as.numeric(sp %in% testmat[,col]));
	}
}; rm(s, col);
#######################



### other tissue

source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt', header=T, sep='\t');

ids = dat0[, 1:9];
dat = dat0[, 10:141];
## kidney network
dat.kd = dat[, grep('kd',names(dat))];
rownames(dat.kd) = ids$hsa;
rm(dat0,dat,ids);


out = preProcess(datIN=dat.kd, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.kd)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);


DATA.kd = as.data.frame(t(out$data_Qnorm));
DATA.kd = DATA.kd[, match(names(DATA), names(DATA.kd))];
#.ds2_mm40_mch.2run12
NET.kd = blockwiseModules(DATA.kd, power=14, networkType='signed',deepSplit=2,minModuleSize=40,mergeCutHeight=.2);
NET.kd$colors = matchLabels(NET.kd$colors, colors);
dendro.kd=NET.kd$dendrograms;
block.kd = 1;
blockGenes.kd = NET.kd$blockGenes;
colors.kd = NET.kd$colors;

# need to recompute with new colors
MEs.kd=moduleEigengenes(DATA.kd,colors.kd)$eigengenes


exn.plotDendroAndColors(dendro.kd, colors.kd, block=block.kd, blockGenes=blockGenes.kd);