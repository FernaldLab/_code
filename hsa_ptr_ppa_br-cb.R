rm(list=ls());options(stringsAsFactors=F);setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/exploreNetwork_testing.R');

##### read in raw data
dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');
ids = dat0[, 1:5];
dat = dat0[, 6:83];

##### separate human/chimp/bonobo samples
dat = dat[, grepl('hsa|ptr|ppa', names(dat))];
rownames(dat) = ids$hsa;

##### remove genes with any zeros
zthresh = 1;
dat = dat[apply(dat, 1, function(f) sum(f==0)) < zthresh, ];

##### remove values less than 0.5
dat[dat < 0.5] = NA;

######## get variance ranks in each tissue type
ttypes = unique(substr(gsub('hsa.|ppa.|ptr.','',names(dat)),1,2));
for (ttype in ttypes) {
	assign(paste(ttype, 'Var', sep=''), rank(apply(dat[,grepl(ttype,names(dat))], 1, var)));
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
dat = dat[higherVar < 0, ];

##### separate br/cb samples
dat = dat[, grepl('br|cb', names(dat))];

# remove outliers and quantile normalize
keepMe = names(dat)!='hsa.br.M.4'# & names(dat)!='hsa.cb.F.1';

X = dat[,keepMe];
#X = dat
out = preProcess(datIN=X, 
				 removeOutlierProbes=T, deviate=3, 
				 removeTooManyNAs=T, probe_thresh=zthresh, 
				 sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
				 Qnorm=T);
				 
				 
DATA = as.data.frame(t(out$data_Qnorm));
save(DATA, file='dat_hsa-ptr-ppa_br-cb_readyForNet.Rdata')
##### check scale-freeness and connectivity of human networks from different power values
sft = exn.plotPowers(DATAhs, networkType='signed', blockSize=ncol(DATAhs)+1, verbose=3);

##### check degree distribution for power=18
k = exn.computeConnectivityAndPlotScaleFreeness(DATA, networkType='signed', power=18);


# load hs network and data

#stem = 'DATA_hsa_ptr_ppa_br_cb_signed_p18_ds2_mm40_mch0.15run10';
stem='DATA_hsa_ptr_ppa_br_cb_signed_p18_ds4_mm80_mch0.15run12';
load(paste(stem, 'NET.RData', sep=''));
load(paste(stem, 'DATA.RData', sep=''));
NET = net; rm(net);

exn.getNetworkBasics(NET, DATA);


# build traits table and compute ME correlations
traits.1 = data.frame(br=as.numeric(grepl('br', rownames(DATA))));
rownames(traits.1) = rownames(DATA);
traits.1 = cbind(traits.1, cb=as.numeric(grepl('cb', rownames(DATA))));
traits.1 = cbind(traits.1, hsa=as.numeric(grepl('hsa', rownames(DATA))));
traits.1 = cbind(traits.1, ptr=as.numeric(grepl('ptr', rownames(DATA))));
traits.1 = cbind(traits.1, ppa=as.numeric(grepl('ppa', rownames(DATA))));
ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits.1 = cbind(traits.1, age.index=ages[match(gsub('br|cb|ht|lv|kd|ts','',rownames(DATA)), 
												  gsub('br|cb|ht|lv|kd|ts','',gsub(' ','.',ages$Sample))), 5:7][,2]);
												  
traitCors.1 = exn.computeAndPlotMETraitCors(traits.1, MEs.1, main=stem);

GS.1 = exn.computeGS(traits.1, DATA);

plotDendroAndColors(dendro.1[[1]], 
					data.frame(colors.1, numbers2colors(GS.1$GS.br), numbers2colors(GS.1$GS.cb), numbers2colors(GS.1$GS.hsa), numbers2colors(GS.1$GS.ptr),numbers2colors(GS.1$GS.ppa)
							  ), 
					dendroLabels=F, addGuide=F,
					groupLabels=c('module','GS.br', 'GS.cb','GS.hsa','GS.ptr','GS.ppa'), ylab="Distance (1-TO)", main = '',
					hang=.05, autoColorHeight=F,colorHeight=.6, cex.axis=.8
					);
					
					
tf = read.table('Homo_sapiens_TF_EnsemblID.txt')[, 1];
.checkGeneListEnrichmentList(tf,modGenes.1,names(DATA))$pvals

tmp = read.table('primates_hiv1_interactions',header=F,sep='\t',row.names=1); 
primates_hiv1_interactions = tmp[,1]; names(primates_hiv1_interactions) = rownames(tmp); rm(tmp);
primates_hiv1_interactions = primates_hiv1_interactions[primates_hiv1_interactions==1];
tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);
mammals_pubmed_interactions = mammals_pubmed_interactions[mammals_pubmed_interactions==1];

phiv = names(primates_hiv1_interactions)[names(primates_hiv1_interactions) %in% DATA]

hiv.1 = .checkGeneListEnrichmentList(names(primates_hiv1_interactions),modGenes.1,names(DATA))$pvals;hiv.1

pubmed.1 = .checkGeneListEnrichmentList(names(mammals_pubmed_interactions),modGenes.1,names(DATA))$pvals;pubmed.1

load(file='/Volumes/fishstudies/_mammalian_RNAseq/Disease Sets/ensemblIDs.RData');
disease.1 = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease.1[[s]] = .checkGeneListEnrichmentList(setsEns[[s]], modGenes.1, names(DATA))$pvals;
	names(disease.1)[s] = names(setsEns)[s];
}; rm(s);lapply(disease.1, head); 
common_disease_genes = intersect(intersect(setsEns$"Schizophrenia_phenopedia_03-04-2014.txt.noHead", setsEns$Autistic.noHead), setsEns$Bipolar.noHead);
common_disease_genes2 = intersect(setsEns$"Alzheimer.noHead", setsEns$"Frontotemporal.noHead");

.checkGeneListEnrichmentList(common_disease_genes, modGenes.1, names(DATA))$pvals;
.checkGeneListEnrichmentList(common_disease_genes2, modGenes.1, names(DATA))$pvals;

load(file='/Volumes/fishstudies/_mammalian_RNAseq/_cahoyMaterials/ensemblIDs_celltype.RData')
celltype.1 = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype.1[[s]] = .checkGeneListEnrichmentList(ctEns[[s]], modGenes.1, names(DATA))$pvals;
	names(celltype.1)[s] = names(ctEns)[s];
}; rm(s); lapply(celltype.1, head);

ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
selection.1 = .checkGeneListEnrichmentList(ps.sig[,1], modGenes.1, names(DATA))$pvals;selection.1

kscores0 = read.table('~/Downloads/primates_positive_selection',sep='\t',header=T,row.names=1);
kscores = kscores0[,16];
names(kscores) = rownames(kscores0);
kscores = kscores[names(kscores) %in% names(DATA)];
kscores = kscores[match(names(DATA),names(kscores))];
names(kscores) = names(DATA);

ks.sig = kscores[kscores>quantile(kscores,.95,na.rm=T)];
ks.sig = ks.sig[!is.na(ks.sig)];
selection_k.1 = .checkGeneListEnrichmentList(names(ks.sig), modGenes.1, names(DATA))$pvals;selection_k.1

MFROW = c(4,7);
par(mfrow=MFROW)
for ( m in 1:length(names(table(colors.1)))) {
	tmp = exn.plotModkMEAndkscore(kME.1,colors.1,kscores,names(table(colors.1))[m], abline.col='red');
	abline(h=quantile(kscores[names(kscores) %in% names(DATA)],.9,na.rm=T), lty='dashed')
}

# DAVID enrichments
d.1 = exn.startAndUploadBgAndModListsToDAVID2(colors.1, 'ahilliar@stanford.edu', 'ENSEMBL_GENE_ID', names(DATA));
modDAVID0 = exn.getModChartsFromDAVID(d.1);
#modDAVID = modDAVID0[!names(modDAVID0)=='grey'];
modDAVID = .addkMEColsToDAVIDList(modDAVID0,modkMEs.1);
for(m in 1:length(modDAVID)) {
	if (names(modDAVID)[m]=='grey') {modDAVID[[m]]=cbind(modDAVID[[m]], kscores=rep(NA,1)); next}
	modDAVID[[m]] = .addGeneAvgColToDAVID(modDAVID[[m]],kscores,'kscores');
}; rm(m)
for (m in 1:length(modDAVID)) {
	modDAVID[[m]][,13] = as.numeric(modDAVID[[m]][,13]);
}; rm(m);
modDAVIDf = exn.filterModDAVIDList(modDAVID,13,10,T);
lapply(modDAVIDf,function(f) (f[,c(1,2,3,4,5,11:ncol(f))]));




term_k.quant = quantile(unlist(sapply(modDAVID, function(f) f$kscores)), .9, na.rm=T);
termCorskME_k = .plotTermCorsAllModulesThresh(modDAVID, mfrow=c(4,7), verbose=F, cex=1.3, cex.main=1, cex.lab=1, cex.axis=1, thresh=20, horizontal=term_k.quant);
cor(termCorskME_k.hs[match(gsub('ME','',rownames(traitCors.hs$cor)), names(termCorskME_k.hs))], traitCors.hs$cor, use='p');
cor(termCorskME_k.hs[match(rownames(medRanks), names(termCorskME_k.hs))], medRanks, use='p');
verboseScatterplot(termCorskME_k.hs[match(rownames(medRanks), names(termCorskME_k.hs))], 
				   medRanks[,1], corOptions="use='p'", 
				   xlab='kME-k_cor',ylab='preservation', 
				   abline=T, frame.plot=F, 
				   bg=gsub('ME','',rownames(traitCors.hs$cor)), pch=21, col='black', cex=1.5
				   );