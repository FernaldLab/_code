rm(list=ls()); options(stringsAsFactors=F); library(WGCNA); allowWGCNAThreads();
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');
source('/Volumes/fishstudies/_code/blockwiseModulesEnriched-Feb2013.R');

dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');

ids = dat0[, 1:5];
dat = dat0[, 6:83];
dat.br0 = dat[, grep('br',names(dat))];
rownames(dat.br0) = ids$hsa;
rm(dat0, dat);
collectGarbage();

# remove F samples
dat.br0 = dat.br0[, grepl('M', names(dat.br0))]


# remove genes with too many 0's
zthresh = floor(ncol(dat.hs0)*2/3);
dat.hs = dat.hs0[apply(dat.hs0, 1, function(f) sum(f==0))  < zthresh, ];

# remove genes in bottom 5% expression level
#dat.br = dat.br[apply(dat.br, 1, sum) > quantile(apply(dat.br, 1, sum), .05), ];

out = preProcess(datIN=dat.br, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=zthresh, sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);


DATA = as.data.frame(t(out$data_Qnorm));
# try without removing any samples
sft = pickSoftThreshold(DATA, networkType='signed', verbose=3);
k = softConnectivity(DATA, type='signed', power=16, blockSize=ncol(DATA));

par(mfrow=c(1,2));
scaleFreePlot(k); hist(k, breaks=length(k),col='grey',border='darkgrey');

# pickSoftThreshold: will use block size 11808.
 # pickSoftThreshold: calculating connectivity for given powers...
   # ..working on genes 1 through 11808 of  11808
   # Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1 0.098300  8.230          0.912  5980.0    5950.0   6510
# 2      2 0.000155 -0.157          0.905  3430.0    3420.0   4180
# 3      3 0.009430 -0.646          0.936  2160.0    2160.0   2980
# 4      4 0.049100 -0.989          0.943  1450.0    1450.0   2280
# 5      5 0.112000 -1.160          0.945  1030.0    1020.0   1820
# 6      6 0.196000 -1.310          0.948   756.0     742.0   1510
# 7      7 0.294000 -1.450          0.955   573.0     556.0   1270
# 8      8 0.366000 -1.480          0.965   445.0     427.0   1100
# 9      9 0.434000 -1.540          0.971   353.0     334.0    955
# 10    10 0.498000 -1.610          0.976   285.0     265.0    841
# 11    12 0.600000 -1.710          0.983   194.0     173.0    669
# 12    14 0.671000 -1.790          0.987   138.0     119.0    545
# 13    16 0.725000 -1.840          0.992   101.0      83.6    453
# 14    18 0.769000 -1.880          0.995    76.7      60.7    381
# 15    20 0.798000 -1.900          0.996    59.3      45.0    325


net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=(ncol(DATA)+1), power=16, deepSplit=4, minModuleSize=100, skipThresh=300, saveFileBase='_primate/primate_M_b16_mm100_mch.2_ds4', mergeCutHeight=.2);


#################################
setwd('/Volumes/fishstudies/_mammalian_RNAseq/');
rm(list=ls());
library(WGCNA); library(RDAVIDWebService); library(biomaRt);
allowWGCNAThreads();
source('/Volumes/fishstudies/_code/exploreNetwork.R');
source('/Volumes/fishstudies/_code/checkGeneListEnrichment.R');


stem = '_primate/primate_M_b16_mm100_mch.2_ds4run16';
load('_primate/primate_M_b16_mm100_mch.2_ds4run16DATA.RData');
load('_primate/primate_M_b16_mm100_mch.2_ds4run16NET.RData');

NET=net;rm(net)


dendro=NET$dendrograms;
block = 1;
blockGenes = NET$blockGenes;
colors = NET$colors;
MEs = NET$MEs;

traits = data.frame(human_rest=as.numeric(grepl('hsa', rownames(DATA))));
rownames(traits) = rownames(DATA);
ages = read.table('Mammalian Samples.txt',header=T,sep='\t');
traits = cbind(traits, age.index=ages[match(rownames(DATA), gsub(' ','.',ages$Sample)), 5:7][,2]);
traits = cbind(traits, ptr_rest=as.numeric(grepl('ptr', rownames(DATA))));
traits = cbind(traits, ppa_rest=as.numeric(grepl('ppa', rownames(DATA))));
traits = cbind(traits, ggo_rest=as.numeric(grepl('ggo', rownames(DATA))));
traits = cbind(traits, ppy_rest=as.numeric(grepl('ppy', rownames(DATA))));
traits = cbind(traits, mml_rest=as.numeric(grepl('mml', rownames(DATA))));
traits = traits[, c(1,3:7,2)];

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
          # module      pval
# 13    lightgreen 4.633e-07
# 12     lightcyan 1.841e-06
# 21     royalblue 8.499e-05
# 20           red 0.0001567
# 1          black 0.0001867
# 8  darkturquoise 0.0004537
# 5      darkgreen  0.001442
# 4           cyan  0.002359
# 11        grey60  0.002727
# 22        salmon  0.006935
# 14   lightyellow  0.007516
# 17        orange  0.008628
# 25        yellow   0.01649
# 2           blue   0.02551
# 19        purple   0.03658
# 24     turquoise   0.05533
# 7        darkred   0.06417
# 6       darkgrey    0.1014
# 18          pink    0.1195
# 16  midnightblue    0.2254
# 10   greenyellow     0.283
# 9          green    0.5751
# 3          brown    0.6824
# 23           tan    0.7196
# 15       magenta         1

tmp = read.table('mammals_pubmed_interactions',header=F,sep='\t',row.names=1);
mammals_pubmed_interactions = tmp[,1]; names(mammals_pubmed_interactions) = rownames(tmp); rm(tmp);
mammals_pubmed_interactions = mammals_pubmed_interactions[mammals_pubmed_interactions==1]
pubmed = checkGeneListEnrichmentList(names(mammals_pubmed_interactions),mod.genes,names(DATA))$pvals
          # module     pval
# 2           blue 7.55e-05
# 1          black  0.00264
# 22        salmon 0.005953
# 13    lightgreen 0.006388
# 4           cyan 0.007041
# 10   greenyellow  0.02853
# 16  midnightblue  0.05401
# 20           red  0.07722
# 12     lightcyan  0.08018
# 19        purple   0.2021
# 24     turquoise   0.2034
# 11        grey60   0.2095
# 14   lightyellow     0.28
# 17        orange   0.3234
# 8  darkturquoise   0.3332
# 3          brown    0.477
# 9          green   0.4944
# 6       darkgrey   0.5853
# 5      darkgreen   0.6417
# 15       magenta    0.673
# 25        yellow   0.6848
# 21     royalblue   0.7816
# 23           tan   0.8268
# 18          pink   0.8418
# 7        darkred        1
#################

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
disease = list();
for (s in 1:length(setsEns)) {
	print(names(setsEns)[s]);
	disease[[s]] = checkGeneListEnrichmentList(setsEns[[s]], mod.genes, names(DATA))$pvals
	names(disease)[s] = names(setsEns)[s];
}; rm(s)
# [1] "Alzheimer.noHead"
        # module      pval
# 9        green 0.0001504
# 14 lightyellow  0.009973
# 18        pink   0.01204
# 23         tan    0.0173
# 15     magenta   0.06772
# 7      darkred    0.2411
# [1] "Autistic.noHead"
          # module      pval
# 1          black 7.826e-05
# 10   greenyellow  0.001173
# 5      darkgreen  0.002146
# 2           blue  0.004182
# 25        yellow   0.01631
# 8  darkturquoise    0.0232
# [1] "Bipolar.noHead"
        # module      pval
# 1        black 2.953e-06
# 10 greenyellow 1.993e-05
# 2         blue 3.978e-05
# 5    darkgreen 0.0001698
# 11      grey60 0.0005314
# 20         red 0.0008719
# [1] "Creutzfeldt-Jakob.noHead"
        # module    pval
# 3        brown 0.05369
# 25      yellow   0.095
# 24   turquoise 0.09873
# 17      orange   0.161
# 14 lightyellow  0.3225
# 12   lightcyan  0.3684
# [1] "Depressive.noHead"
        # module      pval
# 1        black 3.341e-09
# 10 greenyellow  0.003297
# 2         blue  0.004452
# 15     magenta  0.007473
# 14 lightyellow   0.02439
# 23         tan    0.0496
# [1] "Frontotemporal.noHead"
        # module      pval
# 9        green 0.0002032
# 14 lightyellow   0.01919
# 15     magenta   0.04423
# 23         tan   0.04629
# 18        pink   0.05698
# 10 greenyellow   0.07235
# [1] "Huntington.noHead"
      # module    pval
# 25    yellow 0.03407
# 3      brown 0.04103
# 5  darkgreen   0.137
# 9      green  0.1702
# 11    grey60   0.263
# 20       red  0.2687
# [1] "Lewy.noHead"
        # module    pval
# 2         blue 0.04017
# 11      grey60 0.05247
# 23         tan 0.08545
# 14 lightyellow  0.2771
# 12   lightcyan  0.3181
# 22      salmon  0.3947
# [1] "Parkinson.noHead"
          # module     pval
# 8  darkturquoise 0.005054
# 9          green  0.01453
# 15       magenta  0.02708
# 1          black  0.02927
# 5      darkgreen  0.06502
# 24     turquoise    0.132
# [1] "Pick.noHead"
        # module    pval
# 14 lightyellow 0.06279
# 22      salmon 0.09547
# 3        brown  0.2255
# 1        black       1
# 2         blue       1
# 4         cyan       1
# [1] "Schizophrenia_phenopedia_03-04-2014.txt.noHead"
         # module      pval
# 1         black 1.529e-12
# 10  greenyellow 1.558e-06
# 11       grey60 0.0002675
# 2          blue 0.0003567
# 20          red  0.001564
# 16 midnightblue   0.00204

#################
dir = '_cahoyMaterials';
files = paste(dir, list.files(dir)[grepl('csv$', list.files(dir))], sep='/');
ct = list();
for (f in 1:length(files)) {
	ct[[f]] = read.csv(files[f]);
	names(ct)[f] = gsub(dir, '', gsub('.csv', '', files[f]));
}; rm(f);
ct = ct[unlist(lapply(ct,nrow))!=0];
ct0 = ct;

ct = ct0;
THRESH = 10;
ct = lapply(ct, function(f) f[f[,2]>THRESH, ]);

ctEns = list();
for (ctl in 1:length(ct)) {
	print(names(ct)[ctl]);
	ctEns[[ctl]] = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=ct[[ctl]][,1], mart=mart)[,1];
	names(ctEns)[ctl] = names(ct)[ctl];
}; rm(ctl);
celltype = list();
for (s in 1:length(ctEns)) {
	print(names(ctEns)[s]);
	celltype[[s]] = checkGeneListEnrichmentList(ctEns[[s]], mod.genes, names(DATA))$pvals
	names(celltype)[s] = names(ctEns)[s];
}; rm(s)

# [1] "/s15astro_dev"
      # module    pval
# 1      black 0.01675
# 2       blue 0.06767
# 25    yellow  0.2721
# 3      brown       1
# 4       cyan       1
# 5  darkgreen       1
# [1] "/s16astro_mat"
      # module   pval
# 9      green 0.1075
# 24 turquoise 0.2961
# 1      black      1
# 2       blue      1
# 3      brown      1
# 4       cyan      1
# [1] "/s17oligo_pc"
         # module      pval
# 2          blue 3.043e-07
# 5     darkgreen 1.272e-05
# 14  lightyellow 0.0002563
# 1         black 0.0005161
# 24    turquoise 0.0007726
# 16 midnightblue  0.002148
# [1] "/s18oligo_pcVSmyelin"
         # module      pval
# 7       darkred 1.899e-19
# 9         green 4.435e-05
# 2          blue   0.01634
# 23          tan   0.07732
# 16 midnightblue    0.1198
# 4          cyan    0.1234
# [1] "/s4astro"
          # module      pval
# 14   lightyellow 2.734e-07
# 18          pink   0.03642
# 8  darkturquoise   0.06485
# 9          green   0.07779
# 12     lightcyan    0.1897
# 22        salmon     0.266
# [1] "/s5oligo"
    # module      pval
# 7  darkred 7.556e-19
# 9    green 1.816e-06
# 2     blue   0.02818
# 25  yellow   0.03332
# 20     red   0.08044
# 19  purple    0.1748
# [1] "/s6neuron"
      # module      pval
# 2       blue 3.002e-11
# 1      black  1.83e-10
# 5  darkgreen  3.16e-10
# 24 turquoise 6.433e-09
# 18      pink   0.00137
# 15   magenta  0.003088

#################

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

##################

ps = read.table('journal.pgen.1000840.s009.txt');
ps.sig = ps[ps$V5<.05, ];
selection = checkGeneListEnrichmentList(ps.sig[,1], mod.genes, names(DATA))$pvals;


# .kscorePlots = function(module, modGenes, ks, kME, colors) {
	# kME.mod = exn.getModulekME(module, colors, kME)
	# ks = ks[ks[,1] %in% modGenes[[which(names(modGenes)==module)]], ];
	# ks = ks[match(rownames(kME.mod), ks[,1]),];
	# plot(1:nrow(ks), ks[,5], type='l',xlab='genes ranked by kME',ylab='kscore',main=module)
# }


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


#############
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

##################
modFNC = exn.fundamentalModuleConcepts(colors,DATA,power=16);
modDCH = matrix(nrow=length(modFNC),ncol=3, dimnames=list(names(modFNC), c('density','centralization','heterogeneity')));
for (m in 1:length(modFNC)) {
	modDCH[m, ] = c(modFNC[[m]]$Density, modFNC[[m]]$Centralization, modFNC[[m]]$Heterogeneity);
}; rm(m);

modInfo = cbind(modDCH, modPvals[match(rownames(modDCH), rownames(modPvals)), ]);
modInfo = cbind(modInfo, traitCors$cor[match(rownames(modInfo), gsub('ME','',rownames(traitCors$cor))), ]);

par(mfrow=c(3,7));for ( x in 4:24){verboseBoxplot(modInfo[,x],as.factor(modInfo[,25]<0),xlab='',ylab=colnames(modInfo)[x])}
###################


tmp = matrix(as.numeric(modPvals<(.05/450)),ncol=ncol(modPvals))
dimnames(tmp)=dimnames(modPvals)
heatmap(tmp,Rowv=NA,Colv=NA,scale='none',main=stem)


source('/Volumes/fishstudies/_code/preProcATH-for_web_noVSN.R');

dat0 = read.table('/Volumes/fishstudies/Mammalian RNA-seq Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Primate1to1Orthologues.txt', header=T, sep='\t');

ids = dat0[, 1:5];
dat = dat0[, 6:83];
## kidney network
dat.kd = dat[, grep('kd',names(dat))];
rownames(dat.kd) = ids$hsa;
rm(dat0,dat,ids);
dat.kd = dat.kd[rownames(dat.kd) %in% names(DATA), ]
zthresh = floor(ncol(dat.kd)*2/3);
out.kd = preProcess(datIN=dat.kd, removeOutlierProbes=T, deviate=3, removeTooManyNAs=T, probe_thresh=zthresh, sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, Qnorm=T);


DATA.kd = as.data.frame(t(out.kd$data_Qnorm));


#######################

.getModGenesRankedBykME = function(module_names,colors,kME) {
	outList = list();
	for (m in 1:length(module_names)) {
		outList[[m]] = exn.getModulekME(module_names[m],colors,kME)
	}
	names(outList) = module_names;
	return(outList);
}

modkMEs = .getModGenesRankedBykME(names(table(colors)),colors,kME);

lookup = cbind(names(table(colors)),sample(1:length(names(table(colors))),length(names(table(colors)))));
.writekMErownamesRandomNames = function(modkMEs,lookup) {
	for (m in 1:length(modkMEs)) {
		tmp = rownames(modkMEs[[m]])
		write.table(tmp,file=paste('mod_', lookup[match(names(modkMEs)[m], lookup[,1]),2],'_genes_ENS.txt',sep=''),quote=F,col.names=F,row.names=F);
	}
}


#######################
# human specificity
# hs = apply(DATA[1:5,],2,mean,na.rm=T)
# rest = apply(DATA[6:ncol(DATA),],2,mean,na.rm=T)
# fc=cbind(hs,rest)
fold=2
hgenesUp = fc[apply(fc,1,function(f) f[1]>(fold*f[2])), ]
hgenesDown = fc[apply(fc,1,function(f) f[2]>(fold*f[1])), ]
h = rbind(hgenesUp,hgenesDown)



hgenes=h
sort(table(colors[names(DATA) %in% rownames(hgenes)]));
checkGeneListEnrichmentList(rownames(hgenes), mod.genes, names(DATA))$pvals
MOD = 'turquoise'
g=names(DATA)[colors==MOD & names(DATA) %in% rownames(hgenes)]
gd=strsplit(modDAVIDf$turquoise$Genes, ', '); names(gd) = modDAVIDf$turquoise$Term
checkGeneListEnrichmentList(g,gd,names(DATA)[colors==MOD])$pvals





hs5fold=checkGeneListEnrichmentList(rownames(hgenes), mod.genes, names(DATA))$pvals

# sum(apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2)
# fc[apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2, ]
# hgenes = rownames(fc[apply(fc,1,function(f) if(f[1]>f[2]){f[1]/f[2]}else{f[2]/f[1]}) > 2, ])



#####################################

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


avgkscores = c();
par(mfrow=c(5,5));
for(mod in 1:length(mod.genes)) {
	tmp=.plotModkMEAndkscore(kME,colors,kscores,names(mod.genes)[mod],abline.lty='dashed',abline.color='grey');
	avgkscores = c(avgkscores, mean(tmp,na.rm=T));
};
names(avgkscores) = names(mod.genes);




#####################################

.addkscoreToTerms = function(kscores, modDAVID) {
	modDAVID = cbind(modDAVID, avgkscore=1:nrow(modDAVID));
	for (term in 1:nrow(modDAVID)) {
		termIDs = unlist(strsplit(modDAVID[term, ]$Genes,', '));
		avgk = mean(kscores[names(kscores) %in% termIDs], na.rm=T);
		modDAVID[term, ncol(modDAVID)] = avgk;
	}
	return(modDAVID);
}