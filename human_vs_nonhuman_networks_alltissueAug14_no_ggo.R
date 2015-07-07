setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls()); options(stringsAsFactors=F);
library(WGCNA); allowWGCNAThreads();
library(RDAVIDWebService);library(biomaRt);
#source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');
#################################################################################
# load hs network and data
# hs_all_tissue_sd3_densityPermTest_signed_p20_ds2_mm20_mch0.15run5
stem = 'hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10';
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));
NET = net; rm(net);

#################################################################################
load('hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_nonhumanDATAnhp_no_ggo.RData');
common_genes = intersect(names(DATA), names(DATAnhp));
DATAhs = DATA[, names(DATA) %in% common_genes];
DATAnhp = DATAnhp[, names(DATAnhp) %in% common_genes];
load('hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltNEThs_no_ggo.RData');
load('hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs_no_ggo.RData');
load('hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtNETnhp_no_ggo.RData');
load('hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtDATAnhp_no_ggo.RData');
NETnhp$colors = matchLabels(NETnhp$colors, NEThs$colors);
NETnhp$MEs = moduleEigengenes(DATAnhp, NETnhp$colors)$eigengenes;
#################################################################################

exn.getNetworkBasics(NEThs, DATAhs, '.hs');
exn.getNetworkBasics(NETnhp, DATAnhp, '.nhp');

################ merge modules ########################
MEs.hs_merge = mergeCloseModules(DATAhs, colors.hs, MEs=MEs.hs, cutHeight=.2);
MEs.hs = MEs.hs_merge$newMEs;
colors.hs = MEs.hs_merge$colors;
#dendro.hs = MEs.hs_merge$dendro;

MEs.nhp_merge = mergeCloseModules(DATAnhp, colors.nhp, MEs=MEs.nhp, cutHeight=.2);
MEs.nhp = MEs.nhp_merge$newMEs;
dendro.nhp = MEs.nhp_merge$dendro;
colors.nhp = matchLabels(MEs.nhp_merge$colors, colors.hs);
MEs.nhp = moduleEigengenes(DATAnhp, colors.nhp)$eigengenes;
#################################################################


kME.hs = exn.computekME(DATAhs, MEs.hs)$all;
modGenes.hs = exn.getModuleGenes(DATAhs, colors.hs);
modkMEs.hs = .getModGenesRankedBykME(names(table(colors.hs)), colors.hs, exn.computekME(DATAhs, MEs.hs)$all);
homekMEs.hs = .getkMEsFromAssignedModules(modkMEs.hs, data=DATAhs);

kME.nhp = exn.computekME(DATAnhp, MEs.nhp)$all;
modGenes.nhp = exn.getModuleGenes(DATAnhp, colors.nhp);
modkMEs.nhp = .getModGenesRankedBykME(names(table(colors.nhp)), colors.nhp, exn.computekME(DATAnhp, MEs.nhp)$all);

# test module preservation and robustness
# to just get ranks
presRanks = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=T);
medRanks = presRanks$medRanks;
#save(presRanks, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs_modulePreservationQuality_postMerge.RData');

#presStats = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=F);
#save(presStats, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs_modulePreservationStats.RData');

#presStatsObs = presStats$preservation$observed$ref.hs$inColumnsAlsoPresentIn.nhp;
#presStatsObs = presStatsObs[!(rownames(presStatsObs)=='gold'), ];

#presStats$referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.nhp
# # # $referenceSeparability$Z
# # # $referenceSeparability$Z$ref.hs
# # # $referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.hs
# # # [1] NA

# # # $referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.nhp



################################################################################

exn.plotDendroAndColors(dendro.hs, colors.hs, block=1, blockGenes=blockGenes.hs, hang=.05);
plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, colors.nhp), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs', 'nhp'), ylab="Distance (1-TO)", main = '',
					hang=.05, colorHeight=.5, cex.axis=.8
					);
					
					

					
# # # plotDendroAndColors(dendro.hs[[1]], 
					# # # data.frame(colors.hs, labels2colors(as.numeric(names(DATAhs)%in%tf),colorSeq='black'), numbers2colors(kscores[match(names(DATAhs),names(kscores))],colors=blueWhiteRed(1000)[501:1000]), labels2colors(as.numeric(names(DATAhs)%in%setsEns$Alzheimer.noHead),colorSeq='black'), labels2colors(as.numeric(names(DATAhs)%in%setsEns$'Schizophrenia_phenopedia_03-04-2014.txt.noHead'),colorSeq='black')), 
					# # # dendroLabels=F, addGuide=F,
					# # # groupLabels=c('hs', 'tf','kscore','alz','sz'), ylab="Distance (1-TO)", main = '',
					# # # hang=.05, colorHeight=.5, cex.axis=.8
					# # # );
					
plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, colors.nhp, numbers2colors(GS.hs$GS.br), numbers2colors(GS.hs$GS.cb), numbers2colors(GS.hs$GS.br_cb), labels2colors(as.numeric(names(DATAhs)%in%tf),colorSeq='black'), labels2colors(as.numeric(names(DATAhs)%in%names(ks.sig)),colorSeq='black')), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs','nhp','GS.br', 'GS.cb','GS.br_cb','tf','ks.sig'), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.5, cex.axis=.8
					);
					
plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, 
							   colors.nhp, 
							   numbers2colors(GS.hs$GS.br), 
							   numbers2colors(GS.hs$GS.cb), 
							   numbers2colors(GS.hs$GS.br_cb), 
							   numbers2colors(GS.hs$GS.ht),
							   numbers2colors(GS.hs$GS.kd),
							   numbers2colors(GS.hs$GS.lv),
							   numbers2colors(GS.hs$GS.ts),
							   labels2colors(as.numeric(names(DATAhs)%in%tf),colorSeq='black'), 
							   labels2colors(as.numeric(names(DATAhs)%in%names(ks.sig)),colorSeq='black'),
							   labels2colors(as.numeric(names(DATAhs)%in%setsEns$'Schizophrenia_phenopedia_03-04-2014.txt.noHead'),colorSeq='black')), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs','nhp','GS.br', 'GS.cb','GS.br_cb','GS.ht','GS.kd','GS.lv','GS.ts','tf','ks.sig','sz'), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);
					
overlap = exn.plotModuleOverlaps(colors.hs, colors.nhp);
MEnet.hs = exn.plotEigengeneNetworks2(MEs.hs, returnCors=T);
MEnet.nhp = exn.plotEigengeneNetworks2(MEs.nhp, returnCors=T);


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

GS.hs = exn.computeGS(traits.hs, DATAhs);
GS.nhp = exn.computeGS(traits.nhp, DATAnhp);

plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, 
							   colors.nhp, 
							   numbers2colors(GS.hs$GS.br), 
							   numbers2colors(GS.hs$GS.cb), 
							   numbers2colors(GS.hs$GS.br_cb), 
							   numbers2colors(GS.hs$GS.ht),
							   numbers2colors(GS.hs$GS.kd),
							   numbers2colors(GS.hs$GS.lv),
							   numbers2colors(GS.hs$GS.ts)), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs','nhp','GS.br', 'GS.cb','GS.br_cb','GS.ht','GS.kd','GS.lv','GS.ts'), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);
					
					
					


# # # par(mfrow=c(2,4))
# # # for (tr in 1:ncol(traits.hs)) {
	# # # verboseScatterplot(medRanks[,1], 
					   # # # traitCors.hs$cor[match(rownames(medRanks), gsub('ME','',rownames(traitCors.hs$cor))), tr],
					   # # # abline=T,
					   # # # pch=19,
					   # # # col=rownames(medRanks),
					   # # # frame.plot=F,
					   # # # xlab=colnames(medRanks)[1],
					   # # # ylab=colnames(traitCors.hs$cor)[tr],
					   # # # cex=2,
					   # # # ylim=c(-1,1)
					   # # # );
# # # }; rm(tr)

exn.plotMETraitCorsAgainstNums(traitCors=traitCors.hs$cor, modNums=medRanks[,1], mfrow=c(2,4), ylim=c(-1,1), xlab='medRank');




tf = read.table('Homo_sapiens_TF_EnsemblID.txt')[, 1];
numTFs = sapply(modGenes.hs, function(f) sum(f %in% tf));
pcTFs = numTFs/table(colors.hs);
par(mfrow=c(2,4))
for (tr in 1:ncol(traits.hs)) {
	verboseScatterplot(pcTFs, 
					   traitCors.hs$cor[match(names(numTFs), gsub('ME','',rownames(traitCors.hs$cor))), tr],
					   abline=T,
					   pch=19,
					   col=names(pcTFs),
					   frame.plot=F,
					   xlab='pcTFs',
					   ylab=colnames(traitCors.hs$cor)[tr],
					   cex=2,
					   ylim=c(-1,1)
					   );
}; rm(tr)

x=.checkGeneListEnrichmentList(tf,modGenes.hs,names(DATAhs))$pvals$ratio
names(x)=rownames(.checkGeneListEnrichmentList(tf,modGenes.hs,names(DATAhs))$pvals)
.plotMETraitCorsAgainstNums(traitCors.hs$cor,x,c(2,4),c(-1,1),corOptions="method='s'", xlab='TF odds ratio');

# plot kruskal.test of ME expression as function of traits
factors = unlist(strsplit(rownames(DATAhs),'\\.'))[seq(2,length(unlist(strsplit(rownames(DATAhs),'\\.'))),4)];
factors2 = factors; factors2[6:7] = 'br';
factors3 = c(rep('br/cb',7), rep('rest',9));
exn.plotMEsExprWithTraits(MEs.hs, factors, c(3,9), cex.axis=.8, cex.main=.9);

#################################################################################
# check for gene list enrichments in modules
# virus interactions
tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
primates_hiv1_interactions = primates_hiv1_interactions[primates_hiv1_interactions==1];
tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);
mammals_pubmed_interactions = mammals_pubmed_interactions[mammals_pubmed_interactions==1];

phiv = names(primates_hiv1_interactions)[names(primates_hiv1_interactions) %in% DATAhs]

hiv.hs = .checkGeneListEnrichmentList(names(primates_hiv1_interactions),modGenes.hs,names(DATAhs))$pvals;hiv.hs
hiv.nhp = .checkGeneListEnrichmentList(names(primates_hiv1_interactions),modGenes.nhp,names(DATAnhp))$pvals;hiv.nhp

pubmed.hs = .checkGeneListEnrichmentList(names(mammals_pubmed_interactions),modGenes.hs,names(DATAhs))$pvals;pubmed.hs
pubmed.nhp = .checkGeneListEnrichmentList(names(mammals_pubmed_interactions),modGenes.nhp,names(DATAnhp))$pvals;pubmed.nhp

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


x=disease.hs$'Schizophrenia_phenopedia_03-04-2014.txt.noHead'$ratio
names(x)=rownames(disease.hs$'Schizophrenia_phenopedia_03-04-2014.txt.noHead')
.plotMETraitCorsAgainstNums(traitCors.hs$cor,x,c(2,4),c(-1,1),corOptions="method='s'")


.checkGeneListEnrichmentList(common_disease_genes, modGenes.hs, names(DATAhs))$pvals;
.checkGeneListEnrichmentList(common_disease_genes, modGenes.nhp, names(DATAnhp))$pvals;

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

ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
selection.hs = .checkGeneListEnrichmentList(ps.sig[,1], modGenes.hs, names(DATAhs))$pvals;selection.hs
selection.nhp = .checkGeneListEnrichmentList(ps.sig[,1], modGenes.nhp, names(DATAnhp))$pvals;selection.nhp

kscores0 = read.table('~/Downloads/primates_positive_selection',sep='\t',header=T,row.names=1);
kscores = kscores0[,16];
names(kscores) = rownames(kscores0);
kscores = kscores[names(kscores) %in% names(DATAhs)];
kscores = kscores[match(names(DATAhs),names(kscores))];
names(kscores) = names(DATAhs);

ks.sig = kscores[kscores>quantile(kscores,.95,na.rm=T)];
ks.sig = ks.sig[!is.na(ks.sig)];
selection_k.hs = .checkGeneListEnrichmentList(names(ks.sig), modGenes.hs, names(DATAhs))$pvals;selection_k.hs
selection_k.nhp = .checkGeneListEnrichmentList(names(ks.sig), modGenes.nhp, names(DATAnhp))$pvals;selection_k.nhp

MFROW = c(3,9);
par(mfrow=MFROW)
for ( m in 1:length(names(table(colors.hs)))) {
	tmp = exn.plotModkMEAndkscore(kME.hs,colors.hs,kscores,names(table(colors.hs))[m], abline.col='red');
	abline(h=quantile(kscores[names(kscores) %in% names(DATAhs)],.9,na.rm=T), lty='dashed')
}

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
save(modDAVID.hs0, modDAVID.hs, modDAVID.hsf, file=paste(stem, '_DAVID.hs.RData_no_ggo_postMerge'));
load(paste(stem, '_DAVID.hs.RData_no_ggo_postMerge'));

d.nhp = exn.startAndUploadBgAndModListsToDAVID2(colors.nhp, 'ahilliar@stanford.edu', 'ENSEMBL_GENE_ID', names(DATAnhp));
modDAVID.nhp0 = exn.getModChartsFromDAVID(d.nhp);
#modDAVID.nhp = modDAVID.nhp0[!names(modDAVID.nhp0)=='grey'];
modDAVID.nhp = .addkMEColsToDAVIDList(modDAVID.nhp0,modkMEs.nhp);
for(m in 1:length(modDAVID.nhp)) {
	modDAVID.nhp[[m]] = .addGeneAvgColToDAVID(modDAVID.nhp[[m]],kscores,'kscores');
}; rm(m)
for (m in 1:length(modDAVID.nhp)) {
	modDAVID.nhp[[m]][,13] = as.numeric(modDAVID.nhp[[m]][,13]);
}; rm(m);
modDAVID.nhpf = exn.filterModDAVIDList(modDAVID.nhp,13,10,T);
lapply(modDAVID.nhpf,function(f) (f[,c(1,2,3,4,5,11:ncol(f))]));
# # # save(modDAVID.nhp0, modDAVID.nhp, modDAVID.nhpf, file=paste(stem, '_DAVID.nhp.RData_no_ggo_postMerge'));
load(paste(stem, '_DAVID.nhp.RData_no_ggo_postMerge'));

term_k.quant = quantile(unlist(sapply(modDAVID.hs, function(f) f$kscores)), .9, na.rm=T);
termCorskME_k.hs = .plotTermCorsAllModulesThresh(modDAVID.hs, mfrow=c(4,7), verbose=F, cex=1.3, cex.main=1, cex.lab=1, cex.axis=1, thresh=20, horizontal=term_k.quant);
cor(termCorskME_k.hs[match(gsub('ME','',rownames(traitCors.hs$cor)), names(termCorskME_k.hs))], traitCors.hs$cor, use='p');
cor(termCorskME_k.hs[match(rownames(medRanks), names(termCorskME_k.hs))], medRanks, use='p');
verboseScatterplot(termCorskME_k.hs[match(rownames(medRanks), names(termCorskME_k.hs))], 
				   medRanks[,1], corOptions="use='p'", 
				   xlab='kME-k_cor',ylab='preservation', 
				   abline=T, frame.plot=F, 
				   bg=gsub('ME','',rownames(traitCors.hs$cor)), pch=21, col='black', cex=1.5
				   );

termCorskME_k.hsf = .plotTermCorsAllModulesThresh(modDAVID.hsf, mfrow=c(3,5), verbose=F, cex=2, thresh=10);
cor(termCorskME_k.hsf[match(gsub('ME','',rownames(traitCors.hs$cor)), names(termCorskME_k.hsf))], traitCors.hs$cor, use='p');
cor(termCorskME_k.hsf[match(rownames(medRanks), names(termCorskME_k.hsf))], medRanks, use='p');
verboseScatterplot(termCorskME_k.hsf[match(gsub('ME','',rownames(traitCors.hs$cor)), names(termCorskME_k.hsf))], 
				   traitCors.hs$cor[,1], corOptions="use='p'", 
				   xlab='kME-k_cor',ylab='trait_cor_br', 
				   abline=T, frame.plot=F, 
				   bg=gsub('ME','',rownames(traitCors.hs$cor)), pch=21, col='black', cex=1.5
				   );
#################

terms.hs = exn.modDAVIDByTerm(modDAVID.hs);
terms.hsU = unlist(terms.hs$uniqueToMod);
sort(table(unlist(terms.hs$uniqueToMod)));
sort(table(unlist(terms.hs$uniqueToMod)) / table(unlist(terms.hs$modsByTerm)));

.countTermsWithPatternInModules('synap|nerve|neuro|neura', modDAVID.hs);
.printTermGeneEnrichment(modDAVID.hsf$yellow, setsEns$'Schizophrenia_phenopedia_03-04-2014.txt.noHead', .05);
.comparekMEWithNumsForFilteredGenesInModTerm(modDAVID.hsf$pink,'GO:0050877~neurological system process',kscores,modkMEs.hs$pink,setsEns$'Schizophrenia_phenopedia_03-04-2014.txt.noHead',notch=F,cex=2);

#################

# enrichments for hs module subsets in nhp modules
g.hsDAVID0 = exn.getDAVIDchartsForModuleSubsetsInOtherNetwork('green', c('green','grey','saddlebrown'), modGenes.hs, modGenes.nhp);
g.hsDAVID = .addkMEColsToDAVIDList2(g.hsDAVID0, modkMEs.hs$green);
for (m in 1:length(g.hsDAVID)) {
	g.hsDAVID[[m]] = .addGeneAvgColToDAVID(g.hsDAVID[[m]],kscores,'kscores');
}; rm(m);
for (m in 1:length(g.hsDAVID)) {
	g.hsDAVID[[m]][,13] = as.numeric(g.hsDAVID[[m]][,13]);
}; rm(m);
g.hsDAVIDf = exn.filterModDAVIDList(g.hsDAVID,13,10,T);
lapply(g.hsDAVIDf,function(f) (f[,c(1,2,3,4,5,11:ncol(f))]));





#################




tf3 = toupper(read.table('TF_Binding_Enrichments_BAYES3.txt')[,1]);
tf6 = toupper(read.table('TF_Binding_Enrichments_BAYES6.txt')[,1]);

mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
tf3b = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters='hgnc_symbol', values=tf3, mart=mart);
tf3b = tf3b[grep('^ENS',tf3b[,1]), ];
tf3b = cbind(tf3b,kscore=kscores[match(tf3b[,1],names(kscores))]);
tf3bnet = tf3b[tf3b[,1]%in%names(DATAhs),];
tf3bnet = cbind(tf3bnet, module=colors.hs[match(tf3bnet[,1], names(DATAhs))]);

exn.getAllTermsWithGene('ENSG00000196628', modDAVID.hs$royalblue);
exn.getAllTermsWithGene('ENSG00000172845', modDAVID.hs$brown);


tf6b = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters='hgnc_symbol', values=tf6, mart=mart);
tf6b = tf6b[grep('^ENS',tf6b[,1]), ];
tf6b = cbind(tf6b,kscore=kscores[match(tf6b[,1],names(kscores))]);				   
tf6bnet = tf6b[tf6b[,1]%in%names(DATAhs),];			   




x='ABHD6 ALDH5A1 ANXA6 AP2S1 AP3D1 ARNTL C6orf55 CBR1 CLN8 COMMD7 COMT CPOX DHX33 DOLPP1 EBSP ENTH EXOSC2 EXOSC5 FADS1 FADS2 FLJ10521 FLJ12270 FLJ12788 GPI GTPBP3 H17 HIP2 HMGCR HSPC121 HTF9C ICK IMP4 JMJD2B KIAA0182 LOC133308 LOC90637 MGC11335 MGLL MMAB MRPS34 NMT1 PCANAP6 PIGL POP7 PQLC2 PXMP2 RAD23B RBM14 RNF170 SCO2 SFXN1 SRD5A1 SS18L1 SSB3 STAF42 STIM1 STRA13 SYT7 THOP1 TIMM8A TMEM23 TNFAIP1 TOLLIP TOR1A TRAF2 VHL XPMC2H YIF1 ZDHHC18 ZHX1';
e2f1bpink = .convertIDsWithBiomaRt(fromID='hgnc_symbol',toID='ensembl_gene_id',ids=strsplit(x,' ',fixed=T)[[1]])


plotDendroAndColors(dendro.hs[[1]], 
					data.frame(colors.hs, colors.nhp, labels2colors(as.numeric(names(DATAhs)%in%e2f1q3),colorSeq='black')), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('hs', 'nhp', 'e2f1_q3'), ylab="Distance (1-TO)", main = '',
					hang=.05, colorHeight=.5, cex.axis=.8
					);


#############################

files = paste('MSigDB_files/', list.files('MSigDB_files'), sep='');
msigdb0 = .parseMSigDBgenesFromMultipleFiles(files);
msigdb = .convertIDsWithBiomaRtList(fromID='hgnc_symbol', toID='ensembl_gene_id', modGenes=msigdb0);
names(msigdb) = gsub('.txt$', '', gsub('MSigDB_files/', '', names(msigdb)));
				   
#############################################



			   