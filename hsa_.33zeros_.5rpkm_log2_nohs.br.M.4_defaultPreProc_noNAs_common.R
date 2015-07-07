setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls()); options(stringsAsFactors=F);
library(WGCNA); allowWGCNAThreads();
library(RDAVIDWebService);library(biomaRt);
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/exploreNetwork_testing.R');

# load hs network and data
stem = 'hsa_.33zeros_.5rpkm_log2_higherBrainVar_nohs.br.M.4_defaultPreProc_noNAs_common/get(stem)_signed_p18_ds2_mm100_mch0.15run9';
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));
NEThs = net; rm(net);
DATAhs = DATA; rm(DATA);

load('hsa_nhp_.33zeros_.5rpkm_log2_nohs.br.M.4_defaultPreProc.Rdata');
DATAnhp = as.data.frame(t(out.nhp$data_Qnorm)); rm(out.nhp);
DATAnhp = DATAnhp[, names(DATAnhp) %in% names(DATAhs)];
collectGarbage();
#################################################################################

stemsplit = unlist(strsplit(unlist(strsplit(stem,'/',fixed=T))[2], '_'));
TYPE = stemsplit[grep('signed', stemsplit)];
POWER = stemsplit[grep('p[0-9]*', stemsplit)];
POWER = as.numeric(substr(POWER,2,nchar(POWER)));
DEEPSPLIT = as.numeric(substr(stemsplit[grep('^ds', stemsplit)], 3, 3));
MINMODULESIZE = stemsplit[grep('^mm', stemsplit)];
MINMODULESIZE = as.numeric(substr(MINMODULESIZE, 3, nchar(MINMODULESIZE)));
MERGECUTHEIGHT = stemsplit[grep('^mch', stemsplit)];
MERGECUTHEIGHT = as.numeric(substr(gsub('run[0-9]*$','',MERGECUTHEIGHT),4,nchar(gsub('run[0-9]*$','',MERGECUTHEIGHT))));
#################################################################################

# build nhp network with same settings
NETnhp = blockwiseModules(datExpr=DATAnhp, maxBlockSize=(ncol(DATAnhp)+1), 
						 power=POWER, networkType=TYPE, 
						 deepSplit=DEEPSPLIT, minModuleSize=MINMODULESIZE, mergeCutHeight=MERGECUTHEIGHT,
						 verbose=3
						 );
NETnhp$colors = matchLabels(NETnhp$colors, NEThs$colors);
NETnhp$MEs = moduleEigengenes(DATAnhp, NETnhp$colors)$eigengenes;

#################################################################################

exn.getNetworkBasics(NEThs, DATAhs, '.hs');
exn.getNetworkBasics(NETnhp, DATAnhp, '.nhp');

plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, colors.nhp), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs', 'nhp'), ylab="Distance (1-TO)", main = '',
					hang=.05, colorHeight=.5, cex.axis=.8
					);

#################################################################################					
# test module preservation and robustness
# to just get ranks
presRanks = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=T);
medRanks = presRanks$medRanks;

presStats = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=F);




#################################################################################

overlap = exn.plotModuleOverlaps(colors.hs, colors.nhp);
exn.plotEigengeneNetworks2(MEs.hs);
exn.plotEigengeneNetworks2(MEs.nhp);

# build traits table and compute ME correlations
traits.hs = data.frame(br=as.numeric(grepl('br', rownames(DATAhs))));
rownames(traits.hs) = rownames(DATAhs);
traits.hs = cbind(traits.hs, cb=as.numeric(grepl('cb', rownames(DATAhs))));
traits.hs = cbind(traits.hs, br_cb=as.numeric(grepl('br|cb', rownames(DATAhs))));
traits.hs = cbind(traits.hs, ht=as.numeric(grepl('ht', rownames(DATAhs))));
traits.hs = cbind(traits.hs, kd=as.numeric(grepl('kd', rownames(DATAhs))));
traits.hs = cbind(traits.hs, lv=as.numeric(grepl('lv', rownames(DATAhs))));
traits.hs = cbind(traits.hs, ts=as.numeric(grepl('ts', rownames(DATAhs))));
ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits.hs = cbind(traits.hs, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATAhs)), 
												  gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);

traits.nhp = data.frame(br=as.numeric(grepl('br', rownames(DATAnhp))));
rownames(traits.nhp) = rownames(DATAnhp);
traits.nhp = cbind(traits.nhp, cb=as.numeric(grepl('cb', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, br_cb=as.numeric(grepl('br|cb', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, ht=as.numeric(grepl('ht', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, kd=as.numeric(grepl('kd', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, lv=as.numeric(grepl('lv', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, ts=as.numeric(grepl('ts', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATAnhp)), 
													gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);
traits.nhp = cbind(traits.nhp, ptr=as.numeric(grepl('ptr', rownames(DATAnhp))));
traits.nhp = cbind(traits.nhp, ppa=as.numeric(grepl('ppa', rownames(DATAnhp))));
#traits.nhp = cbind(traits.nhp, ggo=as.numeric(grepl('ggo', rownames(DATAnhp))));

traitCors.hs = exn.computeAndPlotMETraitCors(traits.hs, MEs.hs, main=stem);
traitCors.nhp = exn.computeAndPlotMETraitCors(traits.nhp, MEs.nhp, main=paste(stem,'_NHP',sep=''));

# par(mfrow=c(2,4))
# for (tr in 1:ncol(traits.hs)) {
	# verboseScatterplot(medRanks[,1], 
					   # abs(traitCors.hs$cor[match(rownames(medRanks), gsub('ME','',rownames(traitCors.hs$cor))), tr]),
					   # abline=T,
					   # pch=19,
					   # col=rownames(medRanks),
					   # frame.plot=F,
					   # xlab=colnames(medRanks)[1],
					   # ylab=colnames(traitCors.hs$cor)[tr],
					   # cex=2,
					   # ylim=c(0,1), 
					   # corOptions="method='s'"
					   # );
# }; rm(tr)

exn.plotMETraitCorsAgainstNums(traitCors=abs(traitCors.hs$cor), modNums=medRanks[,1], 
							   mfrow=c(2,4), ylim=c(0,1), xlab='medRank');
							   
# plot kruskal.test of ME expression as function of traits
factors = unlist(strsplit(rownames(DATAhs),'\\.'))[seq(2,length(unlist(strsplit(rownames(DATAhs),'\\.'))),4)];
factors2 = factors; factors2[6:7] = 'br';
factors3 = c(rep('br/cb',7), rep('rest',9));
exn.plotMEsExprWithTraits(MEs.hs, factors, c(4,6), cex.axis=.8, cex.main=1.1);

GS.hs = exn.computeGS(traits.hs, DATAhs);
GS.nhp = exn.computeGS(traits.nhp, DATAnhp);

colorFrame = data.frame(colors.hs, 
							   colors.nhp, 
							   numbers2colors(GS.hs$GS.br), 
							   numbers2colors(GS.hs$GS.cb), 
							   numbers2colors(GS.hs$GS.br_cb), 
							   numbers2colors(GS.hs$GS.ht),
							   numbers2colors(GS.hs$GS.kd),
							   numbers2colors(GS.hs$GS.lv),
							   numbers2colors(GS.hs$GS.ts));
plotDendroAndColors(dendro.hs[[1]], colorFrame,
					dendroLabels=F, addGuide=F,
					groupLabels=gsub('numbers2colors.GS.hs.','',names(colorFrame),fixed=T), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);
# colorFrame = data.frame(colors.hs, 
							   # colors.nhp, 
							   # numbers2colors(GS.hs$GS.br), 
							   # numbers2colors(GS.hs$GS.cb), 
							   # numbers2colors(GS.hs$GS.br_cb), 
							   # numbers2colors(GS.hs$GS.ht),
							   # numbers2colors(GS.hs$GS.kd),
							   # numbers2colors(GS.hs$GS.lv),
							   # numbers2colors(GS.hs$GS.ts), labels2colors(as.numeric(names(DATAhs) %in% names(GS.hsth$numsig)[GS.hsth$numsig==2]),colorSeq='black')
							   # );	
# plotDendroAndColors(dendro.hs[[1]], colorFrame,
					# dendroLabels=F, addGuide=F,
					# groupLabels=gsub('numbers2colors.GS.hs.','',names(colorFrame),fixed=T), ylab="Distance (1-TO)", main = '',
					# hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					# );				
					

#################################################################################

tf = read.table('Homo_sapiens_TF_EnsemblID.txt')[,1];
numTFs = sapply(modGenes.hs, function(f) sum(f %in% tf));
pcTFs = numTFs/table(colors.hs);
par(mfrow=c(2,4))
for (tr in 1:ncol(traits.hs)) {
	verboseScatterplot(pcTFs, 
					   traitCors.hs$cor[match(names(pcTFs), gsub('ME','',rownames(traitCors.hs$cor))), tr],
					   abline=T,
					   pch=19,
					   col=names(pcTFs),
					   frame.plot=F,
					   xlab='pcTFs',
					   ylab=colnames(traitCors.hs$cor)[tr],
					   cex=2,
					   ylim=c(-1,1),
					   corOptions="method='s'"
					   );
}; rm(tr)

par(mfrow=c(2,4))
for (tr in 1:ncol(traits.hs)) {
	verboseScatterplot(pcTFs, 
					   abs(traitCors.hs$cor[match(names(pcTFs), gsub('ME','',rownames(traitCors.hs$cor))), tr]),
					   abline=T,
					   pch=19,
					   col=names(pcTFs),
					   frame.plot=F,
					   xlab='pcTFs',
					   ylab=colnames(traitCors.hs$cor)[tr],
					   cex=2,
					   ylim=c(0,1), corOptions="method='s'"
					   );
}; rm(tr)

# tf.hs = .checkGeneListEnrichmentList(tf, modGenes.hs, names(DATAhs));
# x = tf.hs$pvals$ratio;
# names(x) = rownames(tf.hs$pvals)
# .plotMETraitCorsAgainstNums(traitCors.hs$cor,x,c(2,4),c(-1,1),corOptions="method='s'", xlab='TF odds ratio');
# .plotMETraitCorsAgainstNums(abs(traitCors.hs$cor),x,c(2,4),c(0,1),corOptions="method='s'", xlab='TF odds ratio');
# rm(x);

tf.hs = .plotMETraitCorsAgainstGeneEnrichmentOddsRatio(tf, modGenes.hs, names(DATAhs), traitCors.hs$cor, ylim=c(-1,1));
.plotMETraitCorsAgainstGeneEnrichmentOddsRatio(tf, modGenes.hs, names(DATAhs), abs(traitCors.hs$cor), ylim=c(0,1), return=F);

#################################################################################

load(file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
disease.hs = list();
disease.nhp = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease.hs[[s]] = .checkGeneListEnrichmentList(setsEns[[s]], modGenes.hs, names(DATAhs))$pvals;
	disease.nhp[[s]] = .checkGeneListEnrichmentList(setsEns[[s]], modGenes.nhp, names(DATAnhp))$pvals
	names(disease.hs)[s] = names(setsEns)[s];
	names(disease.nhp)[s] = names(setsEns)[s];
}; rm(s);lapply(disease.hs, head); lapply(disease.nhp,head)
common_disease_genes = intersect(intersect(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", setsEns$Autistic.noHead), setsEns$Bipolar.noHead);
common_disease_genes2 = intersect(setsEns$"Alzheimer.noHead", setsEns$"Frontotemporal.noHead");

.plotMETraitCorsAgainstGeneEnrichmentOddsRatio(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", 
											   modGenes.hs, names(DATAhs), traitCors.hs$cor, 
											   ylim=c(-1,1), 
											   return=F
											   );
.checkGeneListEnrichmentList(common_disease_genes, modGenes.hs, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(common_disease_genes2, modGenes.hs, names(DATAhs))$pvals;

#################################################################################

load(file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData')
celltype.hs = list();
celltype.nhp = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype.hs[[s]] = .checkGeneListEnrichmentList(ctEns[[s]], modGenes.hs, names(DATAhs))$pvals;
	celltype.nhp[[s]] = .checkGeneListEnrichmentList(ctEns[[s]], modGenes.nhp, names(DATAnhp))$pvals;
	names(celltype.hs)[s] = names(ctEns)[s];
	names(celltype.nhp)[s] = names(ctEns)[s];
}; rm(s); lapply(celltype.hs, head); lapply(celltype.nhp, head);

.plotMETraitCorsAgainstGeneEnrichmentOddsRatio(ctEns$"/s6neuron", 
											   modGenes.hs, names(DATAhs), traitCors.hs$cor, 
											   ylim=c(-1,1), 
											   return=F
											   );
save.image(file=paste('WORKSPACE',gsub('/get(stem)', '_', stem,fixed=T),'.Rdata',sep=''));
#################################################################################

kscores0 = read.table('~/Downloads/primates_ryan_matrix.txt',sep='\t',header=F,row.names=1);
kscoresTotal = kscores0[,16];
kscoresHuman = kscores0[,10];
names(kscoresTotal) = rownames(kscores0);
names(kscoresHuman) = rownames(kscores0);

kscoresTotal = kscoresTotal[names(kscoresTotal) %in% names(DATAhs)];
kscoresTotal = kscoresTotal[match(names(DATAhs),names(kscoresTotal))];
names(kscoresTotal) = names(DATAhs);

kscoresHuman = kscoresHuman[names(kscoresHuman) %in% names(DATAhs)];
kscoresHuman = kscoresHuman[match(names(DATAhs),names(kscoresHuman))];
names(kscoresHuman) = names(DATAhs);

MFROW = c(4,6);
kToPlot = kscoresHuman;
par(mfrow=MFROW)
for ( m in 1:length(names(table(colors.hs)))) {
	tmp = exn.plotModkMEAndkscore(kME.hs,colors.hs,kToPlot,names(table(colors.hs))[m], abline.col='red',corOptions="method='s',use='p'");
	abline(h=quantile(kToPlot[names(kToPlot) %in% names(DATAhs)],.9,na.rm=T), lty='dashed')
}

# par(mfrow=MFROW)
# for ( m in 1:length(names(table(colors.hs)))) {
	# tmp = exn.plotModkMEAndkscoreNoThresh(kME.hs,colors.hs,kToPlot,names(table(colors.hs))[m], abline.col='red',corOptions="method='s',use='p'");
	# abline(h=quantile(kToPlot[names(kToPlot) %in% names(DATAhs)],.9,na.rm=T), lty='dashed')
# }

tmp = kscoresTotal
mod.k.hs = .getAvgNumForModules(modGenes.hs, tmp);
mod.k.nhp = .getAvgNumForModules(modGenes.nhp, tmp);

par(mfrow=c(1,3));
for (i in 1:3) {
	verboseScatterplot(mod.k.hs[match(rownames(medRanks),names(mod.k.hs))],
						medRanks[,i], 
					   abline=T,
					   pch=19,
					   col=rownames(medRanks),
					   frame.plot=F,
					   ylab=colnames(medRanks)[i],
					   xlab='avg.kscore',
					   cex=2, corOptions="method='s'"
					   );
}


tmp = kscoresTotal
ks.sig = tmp[tmp>quantile(tmp,.95,na.rm=T)];
ks.sig = ks.sig[!is.na(ks.sig)];
selection_k.hs = .checkGeneListEnrichmentList(names(ks.sig), modGenes.hs, names(DATAhs))$pvals;selection_k.hs
selection_k.nhp = .checkGeneListEnrichmentList(names(ks.sig), modGenes.nhp, names(DATAnhp))$pvals;selection_k.nhp

stem











#################################################################################
# DAVID enrichments
d.hs = exn.startAndUploadBgAndModListsToDAVID2(colors.hs, 'ahilliar@stanford.edu', 'ENSEMBL_GENE_ID', names(DATAhs));
modDAVID.hs0 = exn.getModChartsFromDAVID(d.hs);
#modDAVID.hs = modDAVID.hs0[!names(modDAVID.hs0)=='grey'];
modDAVID.hs = .addkMEColsToDAVIDList(modDAVID.hs0,modkMEs.hs);
for(m in 1:length(modDAVID.hs)) {
	if (names(modDAVID.hs)[m]=='grey') {modDAVID.hs[[m]]=cbind(modDAVID.hs[[m]], kscores=rep(NA,1)); next}
	modDAVID.hs[[m]] = .addGeneAvgColToDAVID(modDAVID.hs[[m]],kscores,'kscores');
}; rm(m)
for (m in 1:length(modDAVID.hs)) {
	modDAVID.hs[[m]][,13] = as.numeric(modDAVID.hs[[m]][,13]);
}; rm(m);
modDAVID.hsf = exn.filterModDAVIDList(modDAVID.hs,13,10,T);
lapply(modDAVID.hsf,function(f) (f[,c(1,2,3,4,5,11:ncol(f))]));


#################################################################################

### hCondels

hcondels0 = read.table('hCONDELS.txt', header=T, sep='\t');

mart = useMart('ensembl', 'hsapiens_gene_ensembl');
attributes = c('hgnc_symbol','ensembl_gene_id','entrezgene','description');

hcondels0In = unique(hcondels0$Within);
hcondels0In = unlist(strsplit(hcondels0In, ', '));
hcondels0InBM = getBM(attributes=attributes, filters='hgnc_symbol', values=hcondels0In, mart=mart);
hcondels0InBM_ENS = unique(hcondels0InBM$ensembl_gene_id);
hcondels0InBM_ENSinNet = hcondels0InBM_ENS[hcondels0InBM_ENS %in% names(DATAhs)];
hcondels0InBM_ENSinNetInfo = .exploreGenes(hcondels0InBM_ENSinNet, modGenes.hs, modkMEs.hs, modDAVID.hs);

.checkGeneListEnrichmentList(hcondels0InBM_ENSinNet, ctEns, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(hcondels0InBM_ENSinNet, setsEns, names(DATAhs))$pvals;
.verboseBoxplotColumns(GS.hs[,grep('^GS', names(GS.hs))], rownames(GS.hs) %in% hcondels0InBM_ENSinNet, c(2,4));
verboseBoxplot(kscoresTotal, names(kscoresTotal) %in% hcondels0InBM_ENSinNet);
verboseBoxplot(kscoresHuman, names(kscoresHuman) %in% hcondels0InBM_ENSinNet);

hcondels0Down = unique(hcondels0$Downstream.gene);
hcondels0Down = unlist(strsplit(hcondels0Down, ', '));
hcondels0DownBM = getBM(attributes=attributes, filters='hgnc_symbol', values=hcondels0Down, mart=mart);
hcondels0DownBM_ENS = unique(hcondels0DownBM$ensembl_gene_id);
hcondels0DownBM_ENSinNet = hcondels0DownBM_ENS[hcondels0DownBM_ENS %in% names(DATAhs)];
hcondels0DownBM_ENSinNetInfo = .exploreGenes(hcondels0DownBM_ENSinNet, modGenes.hs, modkMEs.hs, modDAVID.hs);

.checkGeneListEnrichmentList(hcondels0DownBM_ENSinNet, ctEns, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(hcondels0DownBM_ENSinNet, setsEns, names(DATAhs))$pvals;
.verboseBoxplotColumns(GS.hs[,grep('^GS', names(GS.hs))], rownames(GS.hs) %in% hcondels0DownBM_ENSinNet, c(2,4));
verboseBoxplot(kscoresTotal, names(kscoresTotal) %in% hcondels0DownBM_ENSinNet);
verboseBoxplot(kscoresHuman, names(kscoresHuman) %in% hcondels0DownBM_ENSinNet);

hcondels0Up = unique(hcondels0$Upstream.gene);
hcondels0Up = unlist(strsplit(hcondels0Up, ', '));
hcondels0UpBM = getBM(attributes=attributes, filters='hgnc_symbol', values=hcondels0Up, mart=mart);
hcondels0UpBM_ENS = unique(hcondels0UpBM$ensembl_gene_id);
hcondels0UpBM_ENSinNet = hcondels0UpBM_ENS[hcondels0UpBM_ENS %in% names(DATAhs)];
hcondels0UpBM_ENSinNetInfo = .exploreGenes(hcondels0UpBM_ENSinNet, modGenes.hs, modkMEs.hs, modDAVID.hs);

.checkGeneListEnrichmentList(hcondels0UpBM_ENSinNet, ctEns, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(hcondels0UpBM_ENSinNet, setsEns, names(DATAhs))$pvals;
.verboseBoxplotColumns(GS.hs[,grep('^GS', names(GS.hs))], rownames(GS.hs) %in% hcondels0UpBM_ENSinNet, c(2,4));
verboseBoxplot(kscoresTotal, names(kscoresTotal) %in% hcondels0UpBM_ENSinNet);
verboseBoxplot(kscoresHuman, names(kscoresHuman) %in% hcondels0UpBM_ENSinNet);

hcondelsENSinNet = unique(c(hcondels0InBM_ENSinNet, hcondels0DownBM_ENSinNet, hcondels0UpBM_ENSinNet));
.checkGeneListEnrichmentList(hcondelsENSinNet, ctEns, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(hcondelsENSinNet, setsEns, names(DATAhs))$pvals;
.verboseBoxplotColumns(GS.hs[,grep('^GS', names(GS.hs))], rownames(GS.hs) %in% hcondelsENSinNet, c(2,4));
verboseBoxplot(kscoresTotal, names(kscoresTotal) %in% hcondelsENSinNet);

hcondelsENSinNetInfo = .exploreGenes(hcondelsENSinNet, modGenes.hs, modkMEs.hs, modDAVID.hs);



tmp = GS.hs$GS.br;
names(tmp) = rownames(GS.hs);
.verboseBoxplotAcrossModules(modGenes.hs, tmp, hcondelsENSinNet, c(3,8));

.verboseBoxplotAcrossModules(modGenes.hs, kscoresTotal, hcondelsENSinNet, c(3,8));
.verboseBoxplotAcrossModules(modGenes.hs, kscoresHuman, hcondelsENSinNet, c(3,8));
