#################################################################################
#
# data for hs network was preprocessed in "hs_vs_nhp_nets.R"
# briefly:
# 	1) removed genes with any 0s
#	2) centered genes to have equivalent means
# 	3) computed variance for each gene in each tissue type and ranked them
#		- averaged br+cb ranks and ht+kd+lv+ts ranks
# 		- subtracted average brain variance from average non-brain variance
#		- kept genes with higher variance in brain 
#	4) removed outliers, then genes with >1 NA, then q.normed using preProcess 
#	5) tested networks with different parameters using blockwiseModulesEnriched
#		- went forward with hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10
#
#################################################################################

setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls()); options(stringsAsFactors=F);
library(WGCNA); allowWGCNAThreads();
library(RDAVIDWebService);
#source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');

#################################################################################
# load hs network and data
stem = 'hs_all-tissue_signed_p18/hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10';
#hs_all-tissueSCALED_signed_p18_ds4_mm80_mch0.15run8dendro-block1
#hs_all-tissueSCALED_signed_p18_ds2_mm80_mch0.15run10dendro-block1
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));
NET = net; rm(net);

#################################################################################
# load non-human data to process for comparison to human
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];
dat.nhp0 = dat[, grep('ptr|ppa', names(dat))];
rownames(dat.nhp0) = ids$hsa;
rm(dat0);

# remove genes with any 0s
zthresh = 1
dat.nhp = dat.nhp0[apply(dat.nhp0, 1, function(f) sum(f==0)) < zthresh, ];

# center genes to all have equivalent means
x = as.data.frame(t(apply(dat.nhp, 1, function(f) 200 + as.numeric(f) - mean(as.numeric(f)))));
names(x)=names(dat.nhp);
dat.nhp = x; rm(x);

# # # dat.nhp = as.data.frame(t(scale(t(dat.nhp)))) + 10

# # # # preprocess to remove outliers and normalize
keepMe = !(names(dat.nhp) %in% c('ptr.kd.M.1', 'ggo.cb.M.1'));
X = dat.nhp[,keepMe];
out = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
# # # DATAnhp = as.data.frame(t(out$data_Qnorm));
# # # save(DATAnhp, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_nonhumanDATAnhp.RData');

load('hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_nonhumanDATAnhp.RData');

#################################################################################
# limit each dataset to genes in both
common_genes = intersect(names(DATA), names(DATAnhp));
DATAhs = DATA[, names(DATA) %in% common_genes];
DATAnhp = DATAnhp[, names(DATAnhp) %in% common_genes];

# rebuild human network with same settings as before				### SHOULD I REBUILD OR JUST FILTER EXISTING NETWORK TO BE ABLE TO COMPARE TO NHP?
NEThs = blockwiseModules(datExpr=DATAhs, maxBlockSize=(ncol(DATAhs)+1), 
						 power=18, networkType='signed', 
						 deepSplit=2, minModuleSize=20, mergeCutHeight=.15,
						 verbose=3
						 );
# # # save(NEThs, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltNEThs.RData');					 
# # # save(DATAhs, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs.RData');	
load('hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltNEThs.RData');
load('hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs.RData');

# build nhp network with same settings
NETnhp = blockwiseModules(datExpr=DATAnhp, maxBlockSize=(ncol(DATAnhp)+1), 
						 power=18, networkType='signed', 
						 deepSplit=2, minModuleSize=20, mergeCutHeight=.15,
						 verbose=3
						 );
# # # save(NETnhp, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtNETnhp.RData');					 
# # # save(DATAnhp, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtDATAnhp.RData');
load('hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtNETnhp.RData');
load('hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_builtDATAnhp.RData');
NETnhp$colors = matchLabels(NETnhp$colors, NEThs$colors);
NETnhp$MEs = moduleEigengenes(DATAnhp, NETnhp$colors)$eigengenes;
#################################################################################


#################################################################################

exn.getNetworkBasics(NEThs, DATAhs, '.hs');
exn.getNetworkBasics(NETnhp, DATAnhp, '.nhp');

# test module preservation and robustness
# to just get ranks
presRanks = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=T);
medRanks = presRanks$medRanks;
#save(mp, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs_modulePreservationQuality.RData');

presStats = exn.computeModulePreservation(DATAhs, DATAnhp, colors.hs, colors.nhp, c('hs','nhp'), ranksOnly=F);
save(presStats, file='hs_all-tissue_signed_p18_ds2_mm20_mch0.15run10_rebuiltDATAhs_modulePreservationStats.RData');

presStatsObs = presStats$preservation$observed$ref.hs$inColumnsAlsoPresentIn.nhp;
presStatsObs = presStatsObs[!(rownames(presStatsObs)=='gold'), ];

presStats$referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.nhp
# # # $referenceSeparability$Z
# # # $referenceSeparability$Z$ref.hs
# # # $referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.hs
# # # [1] NA

# # # $referenceSeparability$Z$ref.hs$inColumnsAlsoPresentIn.nhp


modQual.hs = presStats$quality$Z$ref.hs$inColumnsAlsoPresentIn.nhp

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
					
# # # plotDendroAndColors(dendro.hs[[1]], 
					# # # data.frame(colors.hs, colors.nhp, numbers2colors(GS.hs$GS.br), numbers2colors(GS.hs$GS.cb), numbers2colors(GS.hs$GS.br_cb), labels2colors(as.numeric(names(DATAhs)%in%tf),colorSeq='black'), labels2colors(as.numeric(names(DATAhs)%in%names(ks.sig)),colorSeq='black')), 
					# # # dendroLabels=F, addGuide=F,
					# # # groupLabels=c('hs','nhp','GS.br', 'GS.cb','GS.br_cb','tf','ks.sig'), ylab="Distance (1-TO)", main = '',
					# # # hang=.05, autoColorHeight=F,colorHeight=.5, cex.axis=.8
					# # # );
					
exn.plotModuleOverlaps(colors.hs, colors.nhp);
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

par(mfrow=c(2,4))
for (tr in 1:ncol(traits.hs)) {
	verboseScatterplot(medRanks[,1], 
					   traitCors.hs$cor[match(rownames(medRanks), gsub('ME','',rownames(traitCors.hs$cor))), tr],
					   abline=T,
					   pch=19,
					   col=rownames(medRanks),
					   frame.plot=F,
					   xlab=colnames(medRanks)[1],
					   ylab=colnames(traitCors.hs$cor)[tr],
					   cex=2,
					   ylim=c(-1,1)
					   );
}; rm(tr)


tf = read.table('Homo_sapiens_TF_EnsemblID.txt')[,1];
numTFs = sapply(modGenes.hs, function(f) sum(f %in% tf));
pcTFs = numTFs/table(colors.hs);
par(mfrow=c(2,4))
for (tr in 1:ncol(traits.hs)) {
	verboseScatterplot(numTFs, 
					   traitCors.hs$cor[match(names(numTFs), gsub('ME','',rownames(traitCors.hs$cor))), tr],
					   abline=T,
					   pch=19,
					   col=names(numTFs),
					   frame.plot=F,
					   xlab='numTFs',
					   ylab=colnames(traitCors.hs$cor)[tr],
					   cex=2,
					   ylim=c(-1,1)
					   );
}; rm(tr)
checkGeneListEnrichmentList(tf, modGenes.hs, names(DATAhs))

# plot kruskal.test of ME expression as function of traits
factors = unlist(strsplit(rownames(DATAhs),'\\.'))[seq(2,length(unlist(strsplit(rownames(DATAhs),'\\.'))),4)];
factors2 = factors; factors2[6:7] = 'br';
factors3 = c(rep('br/cb',7), rep('rest',9));
exn.plotMEsExprWithTraits(MEs.hs, factors, c(4,9), cex.axis=.8, cex.main=1.1);

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

hiv.hs = checkGeneListEnrichmentList(names(primates_hiv1_interactions),modGenes.hs,names(DATAhs))$pvals;hiv.hs
hiv.nhp = checkGeneListEnrichmentList(names(primates_hiv1_interactions),modGenes.nhp,names(DATAnhp))$pvals;hiv.nhp

pubmed.hs = checkGeneListEnrichmentList(names(mammals_pubmed_interactions),modGenes.hs,names(DATAhs))$pvals;pubmed.hs
pubmed.nhp = checkGeneListEnrichmentList(names(mammals_pubmed_interactions),modGenes.nhp,names(DATAnhp))$pvals;pubmed.nhp

load(file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
disease.hs = list();
disease.nhp = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease.hs[[s]] = checkGeneListEnrichmentList(setsEns[[s]], modGenes.hs, names(DATAhs))$pvals;
	disease.nhp[[s]] = checkGeneListEnrichmentList(setsEns[[s]], modGenes.nhp, names(DATAnhp))$pvals
	names(disease.hs)[s] = names(setsEns)[s];
	names(disease.nhp)[s] = names(setsEns)[s];
}; rm(s);lapply(disease.hs, head); lapply(disease.nhp,head)
common_disease_genes = intersect(intersect(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", setsEns$Autistic.noHead), setsEns$Bipolar.noHead) ;

checkGeneListEnrichmentList(common_disease_genes, modGenes.hs, names(DATAhs))$pvals;
checkGeneListEnrichmentList(common_disease_genes, modGenes.nhp, names(DATAnhp))$pvals;

load(file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData')
celltype.hs = list();
celltype.nhp = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype.hs[[s]] = checkGeneListEnrichmentList(ctEns[[s]], modGenes.hs, names(DATAhs))$pvals;
	celltype.nhp[[s]] = checkGeneListEnrichmentList(ctEns[[s]], modGenes.nhp, names(DATAnhp))$pvals;
	names(celltype.hs)[s] = names(ctEns)[s];
	names(celltype.nhp)[s] = names(ctEns)[s];
}; rm(s); lapply(celltype.hs, head); lapply(celltype.nhp, head);

ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
selection.hs = checkGeneListEnrichmentList(ps.sig[,1], modGenes.hs, names(DATAhs))$pvals;selection.hs
selection.nhp = checkGeneListEnrichmentList(ps.sig[,1], modGenes.nhp, names(DATAnhp))$pvals;selection.nhp

kscores0 = read.table('~/Downloads/primates_positive_selection',sep='\t',header=T,row.names=1);
kscores = kscores0[,16];
names(kscores) = rownames(kscores0);
kscores = kscores[names(kscores) %in% names(DATAhs)];
kscores = kscores[match(names(DATAhs),names(kscores))];
names(kscores) = names(DATAhs);

ks.sig = kscores[kscores>quantile(kscores,.95,na.rm=T)];
ks.sig = ks.sig[!is.na(ks.sig)];
selection_k.hs = checkGeneListEnrichmentList(names(ks.sig), modGenes.hs, names(DATAhs))$pvals;selection_k.hs
selection_k.nhp = checkGeneListEnrichmentList(names(ks.sig), modGenes.nhp, names(DATAnhp))$pvals;selection_k.nhp

MFROW = c(4,9);
par(mfrow=MFROW)
for ( m in 1:length(names(table(colors.hs)))) {
	tmp = exn.plotModkMEAndkscore(kME.hs,colors.hs,kscores,names(table(colors.hs))[m], abline.col='red');
	abline(h=quantile(kscores[names(kscores) %in% names(DATAhs)],.9,na.rm=T), lty='dashed')
}


### TO DO ###

# diff.k against kscore
# avg kscore in modules
# kscore / kME ratio
# kME of highest kscore

tf3 = toupper(read.table('TF_Binding_Enrichments_BAYES3.txt')[,1]);
tf6 = toupper(read.table('TF_Binding_Enrichments_BAYES6.txt')[,1]);

mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
tf3b = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters='hgnc_symbol', values=tf3, mart=mart);
tf3b = tf3b[grep('^ENS',tf3b[,1]), ];
tf3b = cbind(tf3b,kscore=kscores[match(tf3b[,1],names(kscores))]);
tf3bnet = tf3b[tf3b[,1]%in%names(DATAhs),]

tf6b = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters='hgnc_symbol', values=tf6, mart=mart);
tf6b = tf6b[grep('^ENS',tf6b[,1]), ];
tf6b = cbind(tf6b,kscore=kscores[match(tf6b[,1],names(kscores))]);

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
# # # save(modDAVID.hs0, modDAVID.hs, modDAVID.hsf, file=paste(stem, '_DAVID.hs.RData'));
load(paste(stem, '_DAVID.hs.RData'));

# # # d.nhp = exn.startAndUploadBgAndModListsToDAVID2(colors.nhp, 'ahilliar@stanford.edu', 'ENSEMBL_GENE_ID', names(DATAnhp));
# # # modDAVID.nhp0 = exn.getModChartsFromDAVID(d.nhp);
# # # #modDAVID.nhp = modDAVID.nhp0[!names(modDAVID.nhp0)=='grey'];
# # # modDAVID.nhp = .addkMEColsToDAVIDList(modDAVID.nhp0,modkMEs.nhp);
# # # for(m in 1:length(modDAVID.nhp)) {
	# # # modDAVID.nhp[[m]] = .addGeneAvgColToDAVID(modDAVID.nhp[[m]],kscores,'kscores');
# # # }; rm(m)
# # # for (m in 1:length(modDAVID.nhp)) {
	# # # modDAVID.nhp[[m]][,13] = as.numeric(modDAVID.nhp[[m]][,13]);
# # # }; rm(m);
# # # modDAVID.nhpf = exn.filterModDAVIDList(modDAVID.nhp,13,10,T);
# # # lapply(modDAVID.nhpf,function(f) (f[1:30,c(1,2,3,4,5,11:ncol(f))]));
# # # save(modDAVID.nhp0, modDAVID.nhp, modDAVID.nhpf, file=paste(stem, '_DAVID.nhp.RData'));
load(paste(stem, '_DAVID.nhp.RData'));

term_k.quant = quantile(unlist(sapply(modDAVID.hs, function(f) f$kscores)), .9);
termCorskME_k.hs = .plotTermCorsAllModulesThresh(modDAVID.hs, mfrow=c(4,9), verbose=F, cex=1.5, cex.main=1, cex.lab=1, cex.axis=1, thresh=20, horizontal=term_k.quant);
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

.plotTermCorsAllModulesThresh(modDAVID.nhpf, mfrow=c(2,3), verbose=F, cex=2, thresh=10);

terms.hs = exn.modDAVIDByTerm(modDAVID.hs);
terms.hsU = unlist(terms.hs$uniqueToMod);
sort(table(unlist(terms.hs$uniqueToMod)));
sort(table(unlist(terms.hs$uniqueToMod)) / table(unlist(terms.hs$modsByTerm)));

.countTermsWithPatternInModules('synap|nerve|neuro|neura', modDAVID.hs);
.printTermGeneEnrichment(modDAVID.hsf$red, names(primates_hiv1_interactions), .1);
.comparekMEWithNumsForFilteredGenesInModTerm(modDAVID.hsf$red,'membrane',kscores,modkMEs.hs$red,names(primates_hiv1_interactions),notch=F,cex=2);

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

exn.getTermInfoAcrossMods('GO:0045202~synapse', modDAVID.hs, terms.hs);
rr = exn.getAllTermsWithGenesFromTerm('GO:0045202~synapse', modDAVID.hs$red);
pt = exn.getAllTermsWithGenesFromTerm('GO:0045202~synapse', modDAVID.hs$paleturquoise);
exn.getDAVIDchartOverlapByTerm(pt[-6], rr[-6]);


rs = modGenes.hs$red[modGenes.hs$red %in% setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead"];
exn.circlePlotCompareGenesAcrossNetworks(rs, DATAhs, DATAnhp, colors.hs, colors.nhp,linecol.gamma=2,cex.labels=c(.4,1),linecol.base=.99,max.cex.points=4);

tmpcol = rep('red', length(modGenes.hs$red));
tmpcol[match(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", modGenes.hs$red)] = 'black';

#########################################################
library(GOFunction);library(biomaRt);
allgenesEntrez = .convertIDsWithBiomaRt(ids=names(DATAhs));
modGenesEntrez.hs = .convertIDsWithBiomaRtList(modGenes=modGenes.hs);


load('human_GO2ALLEGSlist.RData'); 
z=list();
zc=1; 
for (i in 1:length(GO2ALLEGSlist)) { 
	if(!is.na(match('9310',GO2ALLEGSlist[[i]]))) {
		print(names(GO2ALLEGSlist)[i])
		z[[zc]]=GO2ALLEGSlist[[i]]; 
		names(z)[zc]=names(GO2ALLEGSlist)[i]
		zc=zc+1
	} 
}; rm(i,zc)

comm = z[[1]];
for (tt in 2:length(z)) {
	comm = intersect(comm, z[[tt]])
}; rm(tt);

comm1 = comm[comm %in% allgenesEntrez[,1]];
comm1ens = .convertIDsWithBiomaRt(fromID='entrezgene',toID='ensembl_gene_id',ids=comm1);
comm1ens = comm1ens[comm1ens[,1]%in%names(DATAhs), 1];

comm2 = .getkMERankAndQuantileForGenes(comm1ens,modkMEs.hs);
comm2 = comm2[order(comm2$module),]
comm2[comm2$quantile<=.2,]
#########################################################
mod.k.hs = .getAvgNumForModules(modGenes.hs, kscores);
mod.k.nhp = .getAvgNumForModules(modGenes.nhp, kscores);

par(mfrow=c(1,3));
for (i in 1:3) {
	verboseScatterplot(medRanks[,i], 
					   mod.k.hs[match(rownames(medRanks),names(mod.k.hs))],
					   abline=T,
					   pch=19,
					   col=rownames(medRanks),
					   frame.plot=F,
					   xlab=colnames(medRanks)[i],
					   ylab='avg.kscore',
					   cex=2
					   );
}


high.kscoreRank.hs = c();
for (m in 1:length(modGenes.hs)) {
	ks = kscores[names(kscores) %in% modGenes.hs[[m]]];
	tmp=cbind(modkMEs.hs[[m]], kscore=ks[match(rownames(modkMEs.hs[[m]]), names(ks))]);print(head(tmp))
	high.kscoreRank.hs = c(high.kscoreRank.hs, which(tmp$kscore==max(as.numeric(tmp$kscore),na.rm=T)))
}
names(high.kscoreRank.hs) = names(modGenes.hs);


topTermAvgkscore.hs = c();
for (m in 1:length(modDAVID.hs)) {
	tmp=modDAVID.hs[[m]][1:5,];
	tmp = tmp[match('kscores', names(tmp))][[1]];
	topTermAvgkscore.hs = c(topTermAvgkscore.hs, mean(tmp))
}
names(topTermAvgkscore.hs) = names(modDAVID.hs)

topTermAvgkscore.hs = c();
for (m in 1:length(modDAVID.hs)) {
	tmp=modDAVID.hs[[m]][1:50,];
	tmp = tmp[match('kscores', names(tmp))][[1]];
	topTermAvgkscore.hs = c(topTermAvgkscore.hs, mean(unique(tmp)[1:5]))
}
names(topTermAvgkscore.hs) = names(modDAVID.hs)

verboseScatterplot(medRanks[,1], 
					   topTermAvgkscore.hs[match(rownames(medRanks),names(topTermAvgkscore.hs))],
					   abline=T,
					   pch=19,
					   col=rownames(medRanks),
					   frame.plot=F,
					   xlab=colnames(medRanks)[1],
					   ylab='avg.kscore.topTerm',
					   cex=2
					   );
					   
					   
					   
###############################################					   
# get hs module subsets across nhp network
ot = exn.plotModuleOverlaps(colors.hs, colors.nhp, plot=F)$pTable;		
modSubsets = .buildModuleSubsetsFromOverlap(ot, modGenes.hs, modGenes.nhp);		
modSubsets.kscores = .getNumsForModuleOverlapSubsets(modSubsets, kscores); 
modSubsets.kME = .getOriginalkMEForModuleOverlapSubsets(modSubsets, modkMEs.hs);  

par(mfrow=c(3,7))
for (mod in 1:length(modSubsets.kscores)) {
	if (!('grey' %in% names(modSubsets.kscores[[mod]])) | length(modSubsets.kscores[[mod]])==1) {next}; 
	.verboseBoxplotCompareModuleSubsets(modSubsets.kscores[[mod]],'grey',main=names(modSubsets.kscores)[mod],cex.main=1, col=names(modSubsets.kscores)[mod], cex.axis=.8, ylab='')
}; rm(mod);

par(mfrow=c(4,4))
for (mod in c(1:22,24:36)) {
	if (length(modSubsets.kscores[[mod]][[1]])==0 | length(modSubsets.kscores[[mod]])==1) {next}
	.verboseBoxplotCompareModuleSubsets(modSubsets.kscores[[mod]],names(modSubsets.kscores[[mod]]),main=names(modSubsets.kscores)[mod],cex.main=1, col=names(modSubsets.kscores)[mod], cex.axis=.8, ylab='')
}; rm(mod);


biglist = list();
for (row in 1:nrow(ot)) {
	mod = rownames(ot)[row];print(mod)
	omods = colnames(ot)[ot[row, ] < (.05/(nrow(ot)*ncol(ot)))];#print(omods)
	if (length(omods)==0){next}
	modgenes = modGenes.hs[match(mod, names(modGenes.hs))][[1]];#print(length(modgenes))
	modlist = list();
	for (m in 1:length(omods)) {
		print(omods[m])
		omodgenes = modGenes.nhp[match(omods[m], names(modGenes.nhp))][[1]];
		testgenes = intersect(modgenes, omodgenes);
		testchart = exn.getDAVIDchartForGeneList(testgenes,unlist(modGenes.hs))
		modlist[[m]] = testchart;
		names(modlist)[m]=omods[m];
	}
	biglist[[row]] = modlist;
}
names(biglist) = rownames(ot);

###################################

homekMEs.hs = .getkMEsFromAssignedModules(modkMEs.hs, data=DATAhs);
homekMEs.nhp = .getkMEsFromAssignedModules(modkMEs.nhp, data=DATAnhp);

verboseScatterplot(homekMEs.hs[,1], kscores[match(rownames(homekMEs.hs),names(kscores))],
				   col=colors.hs, pch=19, frame.plot=F, 
				   xlab='kME.hs', ylab='kscore', 
				   abline=T);
verboseScatterplot(homekMEs.nhp[,1], kscores[match(rownames(homekMEs.nhp),names(kscores))],
				   col=colors.hs, pch=19, frame.plot=F, 
				   xlab='kME.nhp', ylab='kscore', 
				   abline=T);
				   
	
.plotComparekMEAndExpressionAcrossNetworks(homekMEs.hs[,1],homekMEs.nhp[,1],DATAhs,DATAnhp,colors1=colors.hs,subsetModules=names(modGenes.hs)[!(names(modGenes.hs) %in% c('red','yellow','cyan','turquoise'))])	
				   
.verboseBoxplotComparekMEOfGeneSetsWithinAllModules(modkMEs.hs, tf, c(4,9),names=c('!TF','TF'))

#################################


k.hs = softConnectivity(DATAhs, type='signed', power=18, blockSize=(ncol(DATAhs)+1));
k.hs = k.hs/max(k.hs);
names(k.hs) = names(DATAhs);
k.nhp = softConnectivity(DATAnhp, type='signed', power=18, blockSize=(ncol(DATAnhp)+1));
k.nhp = k.nhp/max(k.nhp);
names(k.nhp) = names(DATAnhp);
k.diff = k.nhp - k.hs;
verboseScatterplot(k.hs, k.nhp, col=colors.hs, pch=19, abline=T, frame.plot=F);
.verboseScatterplotAcrossModules(k.diff,ks,colors.hs,c(4,9),cex.main=.9)

par(mfrow=c(4,9));
for (m in 1:length(modGenes.hs)) {
	mod = names(modGenes.hs)[m];
	g = colors.hs==mod
	verboseScatterplot(k.hs[g], k.nhp[g], col=colors.hs[g], pch=19, abline=T, frame.plot=F, main=mod,xlab='k.hs',ylab='k.nhp',cex.main=1);
}

kIN.hs = .computekINAllModules(DATAhs, colors.hs);
kIN.nhp = .computekINAllModules(DATAnhp, colors.nhp);
kIN.hsVec = c();
for (m in kIN.hs) {
	kIN.hsVec = c(kIN.hsVec, m);
}; rm(m)
kIN.hsVec = kIN.hsVec[match(names(DATAhs), names(kIN.hsVec))];
kIN.nhpVec = c();
for (m in kIN.nhp) {
	kIN.nhpVec = c(kIN.nhpVec, m);
}; rm(m)
kIN.nhpVec = kIN.nhpVec[match(names(DATAhs), names(kIN.nhpVec))];
kIN.diff = kIN.nhpVec - kIN.hsVec;
.verboseScatterplotAcrossModules(kIN.diff,ks,colors.hs,c(4,9),cex.main=.9)

par(mfrow=c(2,2))
verboseScatterplot(sapply(kIN.hs, function(f) median(f/max(f))),
				   sapply(modGenes.hs, function(f) sum(f %in% tf)),
				   col=names(modGenes.hs), 
				   cex=1.8, frame.plot=F, pch=19, abline=T, corOptions="method='p'",
				   ylab='# of TFs',
				   xlab='median kIN');
verboseScatterplot(table(colors.hs),
				   sapply(modGenes.hs, function(f) sum(f %in% tf)),
				   col=names(modGenes.hs), 
				   cex=1.8, frame.plot=F, pch=19, abline=T, corOptions="method='p'",
				   ylab='# of TFs',
				   xlab='# of genes');
verboseScatterplot(traitCors.hs$cor[,1],
				   sapply(modGenes.hs, function(f) sum(f %in% tf)),
				   col=names(modGenes.hs), 
				   cex=1.8, frame.plot=F, pch=19, abline=T, corOptions="method='p'",
				   ylab='# of TFs',
				   xlab=colnames(traitCors.hs$cor)[1]); 
verboseScatterplot(abs(traitCors.hs$cor[,1]),
				   sapply(modGenes.hs, function(f) sum(f %in% tf)),
				   col=names(modGenes.hs), 
				   cex=1.8, frame.plot=F, pch=19, abline=T, corOptions="method='p'",
				   ylab='# of TFs',
				   xlab=colnames(traitCors.hs$cor)[1]);



par(mfrow=c(4,9));
for (m in 1:length(modGenes.hs)) {
	mod = names(modGenes.hs)[m];
	g = colors.hs==mod
	verboseScatterplot(homekMEs.hs[g,1], homekMEs.nhp[g,1], col=colors.hs[g], pch=19, abline=T, frame.plot=F, main=mod,xlab='kME.hs',ylab='kME.nhp',cex.main=1,xlim=c(.2,1),ylim=c(-1,1));
}

par(mfrow=c(2,2));verboseScatterplot(ks[names(ks) %in% tf],k.hs[names(k.hs) %in% tf]); verboseScatterplot(ks[names(ks) %in% tf],kIN.hsVec[names(kIN.hsVec) %in% tf]);verboseScatterplot(ks[!(names(ks) %in% tf)],k.hs[!(names(k.hs) %in% tf)]); verboseScatterplot(ks[!(names(ks) %in% tf)],kIN.hsVec[!(names(kIN.hsVec) %in% tf)])


med.hs = apply(DATAhs/max(DATAhs, na.rm=T), 2, median, na.rm=T);
med.nhp = apply(DATAnhp/max(DATAnhp, na.rm=T), 2, median, na.rm=T);


ks=kscores[names(kscores) %in% names(DATAhs)]
ks=ks[match(names(DATAhs),names(ks))]
names(ks) = names(DATAhs)
modkscores.hs = .getNumsForModulesList(ks, colors.hs, normalize=T)
modkscores.nhp = .getNumsForModulesList(ks, colors.nhp, normalize=T);







keep = colors.hs!='grey';
par(mfrow=c(2,2))
#verboseBoxplot(homekMEs.hs[colors.hs!='grey',1],rownames(homekMEs.hs)[colors.hs!='grey'] %in% tf)
verboseBoxplot(kIN.hsVec[keep],names(kIN.hsVec[keep]) %in% tf, frame.plot=F, xlab='', ylab='intramodular connectivity (kIN.hs)', names=c('!TF','TF'), col='grey')
verboseBoxplot(k.hs[keep],names(k.hs[keep]) %in% tf, frame.plot=F, xlab='', ylab='network connectivity (k.hs)', names=c('!TF','TF'), col='grey')
verboseBoxplot(kIN.diff[keep],names(kIN.diff[keep]) %in% tf, frame.plot=F, xlab='', ylab='differential.kIN.hs', names=c('!TF','TF'),ylim=c(-1,1), col='grey')
verboseBoxplot(k.diff[keep],names(k.diff[keep]) %in% tf, frame.plot=F, xlab='', ylab='differential.k.hs', names=c('!TF','TF'),ylim=c(-1,1), col='grey')



par(mfrow=c(4,9))
for (m in 1:length(modkscores.hs)){
	if(sum(names(modkscores.hs[[m]]) %in% tf)==0){next}
	verboseBoxplot(modkscores.hs[[m]], names(modkscores.hs[[m]]) %in% tf, col=names(modkscores.hs)[m])
}


#########################

gnums = c();
gfact = c();
mcols = c();
for (m in 1:length(modkMEs.hs)) {
	this = modkMEs.hs[[m]];
	if (sum(rownames(this) %in% tf)==0) {next}
	gmean = mean(this[rownames(this) %in% tf, 1]);
	gnums = c(gnums, this[rownames(this) %in% tf, 1]);
	gfact = c(gfact, rep(names(modkMEs.hs)[m], sum(rownames(this) %in% tf)));
	mcols = c(mcols, names(modkMEs.hs)[m])
}
verboseBoxplot(gnums, gfact, col=mcols, cex.axis=.6)