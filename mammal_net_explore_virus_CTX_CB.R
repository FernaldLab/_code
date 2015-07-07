########################################################################
##### ANALYZE CTX AND CEREBELLUM NETS, INCLUDING MAMMAL VIRUS DATA #####
########################################################################

### LOAD NETS AND LIBS ###
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls());
library(WGCNA); library(RDAVIDWebService); library(biomaRt);
allowWGCNAThreads();
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');


stem = 'NET/NET.br.s.ds2_mm40_mch.2run12';
load(paste(stem,'NET.RData', sep=''));
load(paste(stem,'DATA.RData', sep=''));
# cerebellum network constructed in mammal_net_explore.R
load(paste('NET/NET.cb.matchTo_NET.br.s.ds2_mm40_mch.2run12','.RData',sep=''));
load(paste('DATA.cb.matchTo_NET.br.s.ds2_mm40_mch.2run12','.RData',sep=''));

### SET VARIABLES ###
NET=net; DATA=DATA; rm(net);
dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

dendro.cb=NET.cb$dendrograms;
blockGenes.cb=NET.cb$blockGenes;
colors.cb=NET.cb$colors;
MEs.cb=NET.cb$MEs;

### GET TRAITS ###
# ctx
traits = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA))), 
					primate_nonprimate=as.numeric(grepl('hsa|ptr|ppa|ggo|ppy|mml', rownames(DATA))),
					human_nonhumanprimate=c(rep(1,5), rep(0,15), rep(NA,9)),
					male_female=as.numeric(grepl('M', rownames(DATA)))
					);
rownames(traits) = rownames(DATA);

ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits = cbind(traits, age.index=ages[match(rownames(DATA), gsub(' ','.',ages$Sample)), 5:7][,2])

# cerebellum
traits.cb = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA.cb))), 
					primate_nonprimate=as.numeric(grepl('hsa|ptr|ppa|ggo|ppy|mml', rownames(DATA.cb))),
					human_nonhumanprimate=c(rep(1,2), rep(0,9), rep(NA,9)),
					male_female=as.numeric(grepl('M', rownames(DATA.cb)))
					);
rownames(traits.cb) = rownames(DATA.cb);
traits.cb = cbind(traits.cb, age.index=ages[match(rownames(DATA.cb), gsub(' ','.',ages$Sample)), 5:7][,2])

### TRAIT CORRELATIONS ###
# modules
# capture images and save in /Volumes/fishstudies/_mammalian_RNAseq/__analysis
traitCors = exn.computeAndPlotMETraitCors(traits, MEs, main=stem);
traitCors.cb = exn.computeAndPlotMETraitCors(traits.cb, MEs.cb, main=paste('cerebellum: ', stem, sep=''));

# genes
GS = exn.computeGS(traits, DATA);
GS.cb = exn.computeGS(traits.cb, DATA.cb);

### kME ###
kME = exn.computekME(DATA, MEs)$all;
kME.cb = exn.computekME(DATA.cb, MEs.cb)$all;

### PLOTS ###
TRAIT = 'human_nonhumanprimate';
exn.plotAllModsGSkME(colors.cb,TRAIT,GS.cb,kME.cb,mfrow=c(3,7));		# NaNs in cerebellum GS
MOD = 'blue';
exn.plotModuleAllGSkME(MOD,colors,GS,kME);

plotDendroAndColors(dendro[[block]], data.frame(colors[blockGenes[[1]]], colors.cb[blockGenes.cb[[1]]]), dendroLabels=F, groupLabels=c('br','cb'));
exn.plotModuleOverlaps(colors,colors.cb);

### LOAD VIRUS INTERACTION DATA ###
tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
tmp = read.table('mammals_hiv1_interactions',header=F,sep='\t',row.names=1);
mammals_hiv1_interactions = tmp[,1]; names(mammals_hiv1_interactions) = rownames(tmp); rm(tmp);
tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);

### FILTER TO GENES IN NETWORK ###
primates.HIV.network.INTERSECT = names(primates_hiv1_interactions) %in% colnames(DATA);
mammal.HIV.network.INTERSECT = names(mammals_hiv1_interactions) %in% colnames(DATA);
mammal.PUBMED.network.INTERSECT = names(mammals_pubmed_interactions) %in% colnames(DATA);

# Intersection of interaction data with mammals
temp1 = mammals_hiv1_interactions[mammal.HIV.network.INTERSECT];
temp2 = mammals_pubmed_interactions[mammal.PUBMED.network.INTERSECT];
temp3 = primates_hiv1_interactions[primates.HIV.network.INTERSECT];

#Subset intersection for only genes with real interactions
mammal.HIV.interaction.GENES = temp1[temp1==1]; rm(temp1);
mammal.PUBMED.interaction.GENES = temp2[temp2==1]; rm(temp2);
primate.HIV.interaction.GENES = temp3[temp3==1]; rm(temp3);

#All mammal/primate interaction ids in network
mammal.HIV.module.genes = table(colors[colnames(DATA) %in% names(mammal.HIV.interaction.GENES)]);
mammal.HIV.module.genes.cb = table(colors.cb[colnames(DATA.cb) %in% names(mammal.HIV.interaction.GENES)]);

#Get gene names in module
mod.genes=exn.getModuleGenes(DATA, colors)
mod.genes.cb=exn.getModuleGenes(DATA.cb, colors.cb)

#Check for enrichments of interaction genes
mammal.HIV.enrichments.1=checkGeneListEnrichmentList(names(mammal.HIV.interaction.GENES),mod.genes,names(DATA))
mammal.PUBMED.enrichments.1 =checkGeneListEnrichmentList(names(mammal.PUBMED.interaction.GENES),mod.genes,names(DATA)) 
#primate.HIV.enrichments.1 =checkGeneListEnrichmentList(names(primate.HIV.interaction.GENES),mod.genes,names(DATA)) 

mammal.HIV.enrichments.1.cb=checkGeneListEnrichmentList(names(mammal.HIV.interaction.GENES),mod.genes.cb,names(DATA.cb));
mammal.PUBMED.enrichments.1.cb =checkGeneListEnrichmentList(names(mammal.PUBMED.interaction.GENES),mod.genes.cb,names(DATA.cb));

### MAKE VIRUS INTERACTION TABLE ###
mammal.HIV.interaction.GENES = names(mammal.HIV.interaction.GENES);
mammal.PUBMED.interaction.GENES = names(mammal.PUBMED.interaction.GENES);

m.hiv.net = mammals_hiv1_interactions[mammal.HIV.network.INTERSECT];
m.pm.net = mammals_pubmed_interactions[mammal.PUBMED.network.INTERSECT];

cols= c('module','kME',names(GS)[seq(1,length(names(GS)),3)]);

m.hiv = as.data.frame(matrix(nrow=length(m.hiv.net), ncol=7, dimnames=list(names(m.hiv.net), cols)));
for (r in 1:nrow(m.hiv)) {
	g = rownames(m.hiv)[r];
	m.hiv[r, 1] = colors[names(DATA)==g];
	m.hiv[r, 2] = kME[rownames(kME)==g, match(m.hiv[r, 1], gsub('kME','',names(kME)))];
	m.hiv[r, 3:7] = GS[rownames(GS)==g, match(names(m.hiv)[3:7], names(GS))];
}; rm(r);
m.hiv = as.data.frame(cbind(m.hiv.net, m.hiv));
names(m.hiv)[1] = 'int';


m.pm = as.data.frame(matrix(nrow=length(m.pm.net), ncol=7, dimnames=list(names(m.pm.net), cols)));
for (r in 1:nrow(m.pm)) {
	g = rownames(m.pm)[r];
	m.pm[r, 1] = colors[names(DATA)==g];
	m.pm[r, 2] = kME[rownames(kME)==g, match(m.pm[r, 1], gsub('kME','',names(kME)))];
	m.pm[r, 3:7] = GS[rownames(GS)==g, match(names(m.pm)[3:7], names(GS))];
}; rm(r);
m.pm = as.data.frame(cbind(m.pm.net, m.pm));
names(m.pm)[1] = 'int';

MODS = names(table(colors));
hits = mammal.HIV.interaction.GENES;
par(mfrow=c(4,7))
for (m in 1:length(MODS)) {
	mod = MODS[m];
	tmp = exn.getModulekME(mod, colors, kME);
	verboseBoxplot(tmp$kME, as.factor(rownames(tmp) %in% hits), col='grey',xlab='',ylab=names(tmp)[1]);
}; rm(m, mod, tmp);



### GENE ONTOLOGY ###
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

for (m in 1:length(modDAVID)) {
	modDAVID[[m]][,13] = as.numeric(modDAVID[[m]][,13]);
}; rm(m);
modDAVIDf = exn.filterModDAVIDList(modDAVID,13,10,T);
lapply(modDAVIDf,function(f) f[,c(1,2,3,4,5,11:13)]);

terms = exn.modDAVIDByTerm(modDAVID);
termsU = unlist(terms$uniqueToMod);

dir = '/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/';
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


#############

### POSITIVE SELECTION ###

ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];

checkGeneListEnrichmentList(ps.sig[,1], modGenes, names(DATA))$pvals;
red=exn.getModulekME('red', colors, kME);


.kscorePlots = function(module, modGenes, ks, kME) {
	
	ks = ks[ks[,1] %in% modGenes[[which(names(modGenes)==module)]], ];
	ks = ks[match(rownames(kME.mod), ks[,1]),];
	plot(1:nrow(ks), ks[,5], type='l',xlab='genes ranked by kME',ylab='kscore',main=module)
}

.quantileKscores = function(module, modGenes, colors, ks, kME) {print(module)
	ks = ks[ks[,1] %in% modGenes[[which(names(modGenes)==module)]], ];
	kME.mod = exn.getModulekME(module,colors,kME);
	ks = ks[match(rownames(kME.mod), ks[,1]),];
	
	qM = quantile(kME.mod[,1])
	qma = c()
	for (r in 1:nrow(kME.mod)) {
		x = kME.mod[r,1];
		if (x>=qM[1] & x<qM[2]) {
			qma=c(qma,1)
		} else if (x>=qM[2] & x<qM[3]) {
			qma=c(qma,2)
		} else if (x>=qM[3] & x<qM[4]) {
			qma=c(qma,3)
		} else if (x>qM[4]){
			qma=c(qma,4)
		} else {
			qma=c(qma,0)
		}
	}
	print(length(qma)); print(dim(ks))
	verboseBoxplot(ks[,5],as.factor(qma),main=module,col=module,xlab='',ylab='kscore')
}
module='green'
kME.mod = exn.getModulekME(module,colors,kME);
qM = quantile(kME.mod[,1])
qma = c()
for (r in 1:nrow(kME.mod)) {
	x = kME.mod[r,1];
	if (x>=qM[1] & x<qM[2]) {
		qma=c(qma,1)
	} else if (x>=qM[2] & x<qM[3]) {
		qma=c(qma,2)
	} else if (x>=qM[3] & x<qM[4]) {
		qma=c(qma,3)
	} else if (x>qM[4]){
		qma=c(qma,4)
	}
}; rm(r,x)




modGenes.cb = list();
for (m in 1:length(table(colors.cb))) {
	modGenes.cb[[m]] = names(DATA.cb)[colors.cb==names(table(colors.cb))[m]];
	names(modGenes.cb)[m] = names(table(colors.cb))[m];
}; rm(m);


















MODS.cb = names(table(colors.cb));
d.cb = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d.cb, inputIds=names(DATA.cb), idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
for (m in 1:length(MODS.cb)) {
	print(MODS.cb[m]);
	modGenes.cb = names(DATA.cb)[colors.cb==MODS.cb[m]];
	addList(d, inputIds=modGenes, idType='ENSEMBL_GENE_ID', listName=MODS[m], listType='Gene');
}; rm(m, modGenes);
setCurrentBackgroundPosition(d, which(getBackgroundListNames(d) == 'BG'));
setCurrentGeneListPosition(d, 1);
modDAVID = exn.getModChartsFromDAVID(d);
lapply(modDAVID,function(f) head(f[, c(1,2,3,4,5,11:13)]));
lapply(modDAVID,function(f) head(f[grepl('GOTERM',f$Category), c(1,2,3,4,5,11:13)]));


