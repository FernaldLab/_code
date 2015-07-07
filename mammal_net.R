rm(list=ls()); options(stringsAsFactors=F); library(WGCNA);
#setwd('~/Desktop/_mammalian_RNAseq/');
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
#source('~/Desktop/Austin/__fromMBP/_codeOct4/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
#source('~/Desktop/Austin/_code/blockwiseModulesEnriched-Feb2013.R');

dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt', header=T, sep='\t');

ids = dat0[, 1:9];
dat = dat0[, 10:141];

dat.br = dat[, grep('br',names(dat))];

# if remove gallus
dat.br = dat.br[, -c(30,31)]
dat.br = dat.br[, !(names(dat.br) %in% c('hsa.br.M.4' ,'ppa.br.F.2'))];
out1 = preProcess(datIN=dat.br, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.br)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);
# removed hsa.br.M.4 ppa.br.F.2



dat.br0=dat.br;
dat.br = as.data.frame(t(out1$data_Qnorm));

######




###################
source('/Volumes/handsfiler$/FishStudies/_code/exploreNetwork.R');
load('blockwiseModulesEnriched_DS4.permTest.run15NET.RData');
load('blockwiseModulesEnriched_DS4.permTest.run15DATA.RData');

dendro=net$dendrograms;
block = 1;
blockGenes = net$blockGenes;
colors = net$colors;
MEs = net$MEs;

exn.plotDendroAndColors(dendro, colors, block=block, blockGenes=blockGenes);

# get gene ids
genes = ids$hsa[as.numeric(names(DATA))];
names(DATA) = genes; rm(genes);

# 
traits = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA))), 
					primate_nonprimate=as.numeric(grepl('hsa|ptr|ppa|ggo|ppy|mml', rownames(DATA))),
					human_nonhumanprimate=c(rep(1,5), rep(0,15), rep(NA,9)),
					male_female=as.numeric(grepl('M', rownames(DATA)))
					);
rownames(traits) = rownames(DATA);

traitCors = exn.computeAndPlotMETraitCors(traits, MEs);
MEnet = exn.plotEigengeneNetworks2(MEs, returnCors=T);
kME = exn.computekME(DATA, MEs)$all;
GS = exn.computeGS(traits, DATA);

exn.plotAllModsGSkME(colors,'primate_nonprimate',GS,kME,mfrow=c(2,6));

library(RDAVIDWebService);
david = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(david, inputIds=names(DATA), idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
mods = names(table(colors));
for (m in 1:length(mods)) {
	print(mods[m]);
	modGenes = names(DATA)[colors==mods[m]];
	addList(david, inputIds=modGenes, idType='ENSEMBL_GENE_ID', listName=mods[m], listType='Gene');
}; rm(m, modGenes);
setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
setCurrentGeneListPosition(david, 1);

modDAVID = exn.getModChartsFromDAVID(david);
modDAVIDf = exn.filterModDAVIDList(modDAVID,13,10,T);

terms = exn.modDAVIDByTerm(modDAVID);
termsU = unlist(terms$uniqueToMod);

CHART = modDAVID$turquoise;
MOD = 'turquoise';
TRAIT = 'primate_nonprimate';
avgkMEvec = c();
avgGSvec = c();
for (row in 1:nrow(CHART)) {
	genes = strsplit(CHART[row, 6], ', ')[[1]];
	modcol = match(paste('kME', MOD, sep=''), names(kME));
	avgkME = mean(kME[rownames(kME) %in% genes, modcol]);
	avgkMEvec = c(avgkMEvec, avgkME);
	tcol = match(paste('GS.', TRAIT, sep=''), names(GS));
	avgGS = mean(GS[rownames(GS) %in% genes, tcol]);
	avgGSvec = c(avgGSvec, avgGS);
}; rm(row, genes);
tDAVID = cbind(CHART, kME=avgkMEvec, GS=avgGSvec);

###############

# check gene set enrichment
library(biomaRt);
#source('/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_analysis_compare/_code/checkGeneListEnrichment.R');
source('~/Documents/_Fernald_lab/_code/checkGeneListEnrichment2014.R');
#alz=read.table('/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Disease Sets/Alzheimer',skip=5);
#mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
#hits = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='hgnc_symbol', values=alz[,1], mart=mart);

# remove headers with:
# ls | awk '{print "awk '\''BEGIN{FS=\"\\t\"}$2 ~/[0-9]+/'\'' "$1" > "$1".noHead "}' | bash

# read in genes
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
for (m in 1:length(table(net$colors))) {
	modGenes[[m]] = names(DATA)[net$colors==names(table(net$colors))[m]];
	names(modGenes)[m] = names(table(net$colors))[m];
}; rm(m);

for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	print(checkGeneListEnrichmentList(setsEns[[s]], modGenes, names(DATA))$pvals);
}; rm(s)

# res = list();
# MOD = 'blue';
# for (s in 1:length(setsEns)) {
	# print(names(setsEns)[s]);
	# res[[s]] = checkGeneListEnrichment(setsEns[[s]], names(DATA)[net$colors==MOD], names(DATA));
	# names(res)[s] = names(setsEns)[s];
# }; rm(s);

metamodGenes = list(mm1=c(modGenes$turquoise, modGenes$blue, modGenes$greenyellow),
					mm2=c(modGenes$magenta, modGenes$pink, modGenes$brown, modGenes$green),
					mm3=c(modGenes$yellow, modGenes$black, modGenes$red, modGenes$purple, modGenes$tan)
					);

qts = c(.5, .8, .9, .95);
setsFilt = vector(mode='list', length=length(sets));
for (sF in 1:length(setsFilt)) {
	setList = list();
	this_set = sets[[sF]];
	for (qt in 1:length(qts)) {
		setList[[qt]] = this_set[this_set[,2] > quantile(this_set[,2], qts[qt]), ];
		names(setList)[qt] = qts[qt];
	}
	setsFilt[[sF]] = setList;
	names(setsFilt)[sF] = names(sets)[sF];
}
rm(sF,qt,this_set);

mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
setsFiltEns = list();
for (s in 1:length(setsFilt)) {
	print(names(setsFilt)[s]);
	filtList = list();
	for (qt in 1:length(setsFilt[[1]])) {
		print(names(setsFilt[[s]])[qt]);
		filtList[[qt]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=setsFilt[[s]][[qt]][,1], mart=mart)[,1];
		names(filtList)[qt] = names(setsFilt[[s]])[qt];
	}
	setsFiltEns[[s]] = filtList;
	names(setsFiltEns)[s] = names(setsFilt)[s];
}; rm(s);

for (s in 1:length(setsFiltEns)) {
	print(names(setsFiltEns)[s]);
	for (qt in 1:length(setsFiltEns[[1]])) {
		print(names(setsFiltEns[[1]])[qt]);
		print(checkGeneListEnrichmentList(setsFiltEns[[s]][[qt]], modGenes, names(DATA))$pvals);
	}
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



# cell type lists
# astro = read.csv('/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_cahoyMaterials/s4astro.csv');
# astroEns = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=astro[,1], mart=mart)[,1];
# thresh = 3
# astroFilt = astro[astro[,2]>=thresh, ];
# astroFiltEns = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=astroFilt[,1], mart=mart)[,1];

dir = '/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_cahoyMaterials';
files = paste(dir, list.files(dir)[grepl('csv$', list.files(dir))], sep='/');
ct = list();
for (f in 1:length(files)) {
	ct[[f]] = read.csv(files[f]);
	names(ct)[f] = gsub(dir, '', gsub('.csv', '', files[f]));
}; rm(f);

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

################


pickSoft = pickSoftThreshold(dat.br, networkType='signed', blockSize=5500);
pickSoftUn = pickSoftThreshold(dat.br, networkType='unsigned', blockSize=5500);

set = pickSoft$fitIndices;par(mfrow=c(1,3));
	x = set$Power;
	y = set$SFT.R.sq;
	plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '');
	text(x, y, labels = x);
	abline( h = 0.8, col = 'red', lty = 'dashed' );
	y = set$mean.k.;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '');
	text(x, y, labels = x);
	abline( h = 50, col = 'red', lty = 'dashed' );
	y = set$slope;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = '');
	text(x, y, labels = x);
	abline( h = -2, col = 'red', lty = 'dashed' );
	abline( h = -1, col = 'red', lty = 'dashed' )
	
k12 = softConnectivity(dat.br, type='signed', power=12, blockSize=5500);
par(mfrow=c(1,2))
hist(k12); scaleFreePlot(k12);


# try unsigned beta=5
kUn5 = softConnectivity(dat.br, type='unsigned', power=5, blockSize=5500);
par(mfrow=c(1,2))
hist(kUn5); scaleFreePlot(kUn5);
source('/Volumes/handsfiler$/FishStudies/_code/exploreNetwork.R');

						   
						   
#################
dat.cb = dat[, grep('cb',names(dat))];
rownames(dat.cb) = ids$hsa

zeros=apply(dat.cb,1,function(f) sum(f==0,na.rm=T));
dat.cb = dat.cb[zeros<=15, ];


#out = preProcess(datIN=dat.cb, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.cb)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);

out = preProcess(datIN=dat.cb, removeOutlierProbes=T, deviate=2.5, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.cb)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);
# removed mmu.cb.M.1 ggo.cb.M.1

# mmu.cb.M.1 ggo.cb.M.1 removed beforehand
#out2 = preProcess(datIN=dat.cb, removeOutlierProbes=T, deviate=2.5, removeTooManyNAs=T, probe_thresh=floor(ncol(dat.cb)/3), sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);

DATA.cb = as.data.frame(t(out$data_Qnorm));


DATA.cb = DATA.cb[, names(DATA.cb) %in% names(DATA)];



NET.cb = blockwiseModules(DATA.cb,deepSplit=4);

#pickSoftUn = pickSoftThreshold(DATA.cb, networkType='unsigned', blockSize=5500);
# try unsigned beta=6
#kUn6 = softConnectivity(DATA.cb, type='unsigned', power=6, blockSize=5500);
#par(mfrow=c(1,2))
#hist(kUn6); scaleFreePlot(kUn6);


#NET.cb = blockwiseModulesEnriched(DATA=DATA.cb, networkType='unsigned', power=6, skipThresh=200,saveFileBase='cerebellum_un6_matchCtx');

NET.cb$colors = matchLabels(NET.cb$colors, NET$colors[names(DATA) %in% names(DATA)]);

setLabels=c('ctx','cb');
multiExpr=list();
multiExpr[[1]]=list(data=DATA);
multiExpr[[2]]=list(data=DATA.cb2);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=NET$colors;
colorList[[2]]=NET.cb$colors;
names(colorList)=setLabels;

system.time( {				# should only take a minute or two
 	
 	mp = modulePreservation(multiExpr, 
 				   colorList,
 				   networkType = 'unsigned', 
 				   nPermutations = 0,			# no stats, just get ranks
 				   maxModuleSize = max( table(NET$colors) ),
 				   verbose = 3,
 				   savePermutedStatistics = F
 				   );
 	}
 )


















#######
###
# get genes only expressed in humans

