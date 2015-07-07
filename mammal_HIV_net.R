DATA.hiv = DATA[, names(DATA) %in% mammal.HIV.interaction.GENES];
NET.hiv = blockwiseModules(DATA.hiv, power=14, networkType='signed',deepSplit=2,minModuleSize=10,mergeCutHeight=.2);

DATA.cb.hiv = DATA.cb[, names(DATA.cb) %in% names(mammal.HIV.interaction.GENES)];
NET.cb.hiv = blockwiseModules(DATA.cb.hiv, power=14, networkType='signed',deepSplit=2,minModuleSize=10,mergeCutHeight=.2);

dendro.hiv=NET.hiv$dendrograms;
block.hiv = 1;
blockGenes.hiv = NET.hiv$blockGenes;
colors.hiv = NET.hiv$colors;
MEs.hiv = NET.hiv$MEs;

NET.cb.hiv$colors = matchLabels(NET.cb.hiv$colors, colors.hiv);
dendro.cb.hiv=NET.cb.hiv$dendrograms;
block.cb.hiv = 1;
blockGenes.cb.hiv = NET.cb.hiv$blockGenes;
colors.cb.hiv = NET.cb.hiv$colors;
# need to recompute with new colors
MEs.cb.hiv=moduleEigengenes(DATA.cb.hiv,colors.cb.hiv)$eigengenes

exn.plotDendroAndColors(dendro.hiv, colors.hiv, block=block.hiv, blockGenes=blockGenes.hiv);

plotDendroAndColors(dendro.hiv[[block.hiv]], data.frame(colors.hiv[blockGenes.hiv[[1]]], colors.cb.hiv[blockGenes.cb.hiv[[1]]]), dendroLabels=F, groupLabels=c('br','cb'));
exn.plotModuleOverlaps(colors.hiv,colors.cb.hiv);


traitCors.hiv = exn.computeAndPlotMETraitCors(traits, MEs.hiv);
traitCors.cb.hiv = exn.computeAndPlotMETraitCors(traits.cb, MEs.cb.hiv, main=paste('cerebellum: ', sep=''));

MEnet.hiv = exn.plotEigengeneNetworks2(MEs.hiv, returnCors=T);
kME.hiv = exn.computekME(DATA.hiv, MEs.hiv)$all;
GS.hiv = exn.computeGS(traits, DATA.hiv);

MEnet.cb.hiv = exn.plotEigengeneNetworks2(MEs.cb.hiv, returnCors=T);
kME.cb.hiv = exn.computekME(DATA.cb.hiv, MEs.cb.hiv)$all;
GS.cb.hiv = exn.computeGS(traits.cb, DATA.cb.hiv);

TRAIT = 'primate_nonprimate';
exn.plotAllModsGSkME(colors.hiv,TRAIT,GS.hiv,kME.hiv,mfrow=c(3,7));

MOD = 'turquoise';
exn.plotModuleAllGSkME(MOD,colors.hiv,GS.hiv,kME.hiv);


MODS.hiv = names(table(colors.hiv));
d.hiv = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d.hiv, inputIds=names(DATA.hiv), idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
for (m in 1:length(MODS.hiv)) {
	print(MODS.hiv[m]);
	modGenes.hiv = names(DATA.hiv)[colors.hiv==MODS.hiv[m]];
	addList(d.hiv, inputIds=modGenes.hiv, idType='ENSEMBL_GENE_ID', listName=MODS.hiv[m], listType='Gene');
}; rm(m, modGenes.hiv);
setCurrentBackgroundPosition(d.hiv, which(getBackgroundListNames(d.hiv) == 'BG'));
setCurrentGeneListPosition(d.hiv, 1);
modDAVID.hiv = exn.getModChartsFromDAVID(d.hiv);
lapply(modDAVID.hiv,function(f) head(f[, c(1,2,3,4,5,11:13)]));
lapply(modDAVID.hiv,function(f) head(f[grepl('GOTERM',f$Category), c(1,2,3,4,5,11:13)]));

for (m in 1:length(modDAVID.hiv)) {
	modDAVID.hiv[[m]][,13] = as.numeric(modDAVID.hiv[[m]][,13]);
}; rm(m);

modDAVIDf.hiv = exn.filterModDAVIDList(modDAVID.hiv,13,10,T);
lapply(modDAVIDf.hiv,function(f) f[,c(1,2,3,4,5,11:13)]);

terms.hiv = exn.modDAVIDByTerm(modDAVID.hiv);
termsU.hiv = unlist(terms.hiv$uniqueToMod);



############


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

modGenes.hiv = list();
for (m in 1:length(table(colors.hiv))) {
	modGenes.hiv[[m]] = names(DATA.hiv)[colors.hiv==names(table(colors.hiv))[m]];
	names(modGenes.hiv)[m] = names(table(colors.hiv))[m];
}; rm(m);

for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	print(checkGeneListEnrichmentList(setsEns[[s]], modGenes.hiv, names(DATA.hiv))$pvals);
}; rm(s)



##############


dir = '/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_cahoyMaterials';
files = paste(dir, list.files(dir)[grepl('csv$', list.files(dir))], sep='/');
ct = list();
for (f in 1:length(files)) {
	ct[[f]] = read.csv(files[f]);
	names(ct)[f] = gsub(dir, '', gsub('.csv', '', files[f]));
}; rm(f);

ct0 = ct;

ct = ct0;
THRESH = 2.5;
ct = lapply(ct, function(f) f[f[,2]>THRESH, ]);

ctEns = list();
for (ctl in 1:length(ct)) {
	print(names(ct)[ctl]);
	ctEns[[ctl]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=ct[[ctl]][,1], mart=mart)[,1];
	names(ctEns)[ctl] = names(ct)[ctl];
}; rm(ctl);

for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	print(checkGeneListEnrichmentList(ctEns[[s]], modGenes.hiv, names(DATA.hiv))$pvals);
}; rm(s)

##################
ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
checkGeneListEnrichmentList(ps.sig[,1], modGenes.hiv, names(DATA.hiv))$pvals;

ps = read.table('journal.pgen.1000840.s009.txt');
psi = ps[, 5]; names(psi)=ps[,1]
psi=psi[names(psi) %in% names(DATA.hiv)]
tmpGS=exn.makeTmpGS(psi,'psi',GS.hiv)

ps.br = ps[ps[,1]%in%names(DATA.hiv), ]
rownames(ps.br)=ps.br[,1];
ps.br=ps.br[match(rownames(GS.hiv), rownames(ps.br)),]
tmpGS=cbind(GS.hiv,ps.br[,5])
names(tmpGS)[16]='GS.psi';
exn.plotAllModsGSkME(colors.hiv,'psi',tmpGS,kME.hiv,mfrow=c(3,7));

checkGeneListEnrichment(ps[ps[,5]<.05,1], names(DATA.hiv), ps[ps[,1]%in%names(DATA),1]);

########################

# find human specific genes

####################
# module preservation

setLabels=c('br','cb');
multiExpr=list();
multiExpr[[1]]=list(data=DATA.hiv);
multiExpr[[2]]=list(data=DATA.cb.hiv);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=colors.hiv;
colorList[[2]]=colors.cb.hiv;
names(colorList)=setLabels;

system.time( {				# should only take a minute or two
 	
 	mp = modulePreservation(multiExpr, 
 				   colorList,
 				   networkType = 'signed', 
 				   nPermutations = 0,			# no stats, just get ranks
 				   maxModuleSize = max( table(colors.hiv) ),
 				   verbose = 3,
 				   savePermutedStatistics = F
 				   );
 	}
 )
 
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