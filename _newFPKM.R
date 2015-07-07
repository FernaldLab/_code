# isoform network, based on data in _newFPKM
rm(list=ls());
setwd('~/Documents/_Fernald_lab/_hybridRNAseq/_newFPKM');
source('../../_code/preProcATH-for_web_noVSN.R');
options('stringsAsFactors'=F);

##############################################################
###### read in data and get genes common to all samples ######
##############################################################

# read in fpkm
dirs = list.files();
dirs=dirs[grepl('^CV',dirs) | grepl('^MC',dirs)];
files = paste(dirs,'/isoforms.fpkm_tracking.noZeros',sep='');

DAT0 = list();
for (file in 1:length(files)) {
	DAT0[[file]] = read.table(files[file], header=T, row.names=1);
	cat(files[file], '\n');
	print(head(DAT0[[file]]));
	cat('\n');
}; rm(file);
names(DAT0) = gsub('Cufflinkswithannotation', '', dirs);
lapply(DAT0,dim);

# FPKM_status filter
for (set in 1:length(DAT0)) {
	cat(names(DAT0)[set], ':\n', sep='');
	print(dim(DAT0[[set]]));
	print(sum(DAT0[[set]]$FPKM_status!='OK'));
	DAT0[[set]] = DAT0[[set]][DAT0[[set]]$FPKM_status=='OK',]
	print(dim(DAT0[[set]]));
}; rm(set);

# get genes common to all datasets
genes = rownames(DAT0[[1]]);
for (set in 2:length(DAT0)) {
	genes = genes[genes %in% rownames(DAT0[[set]])];
	#genes = intersect(genes, rownames(DAT0[[set]]));
}; rm(set);

# build dataframe of common genes
DAT = data.frame(matrix(nrow=length(genes), ncol=length(DAT0)));
rownames(DAT) = genes;
for (set in 1:length(DAT0)) {
	thisSet = DAT0[[set]];
	thisSet = thisSet[match(genes, rownames(thisSet)), ];
	DAT[, set] = thisSet$FPKM;
	names(DAT)[set] = names(DAT0)[set];
}; rm(set, thisSet);

collectGarbage();

##############################################################
###### remove outliers and normalize #########################
##############################################################

out = preProcess(DAT, removeOutlierProbes=T, deviate=2, removeTooManyNAs=T, removeOutlierSamples=T, Qnorm=T);
DATA=as.data.frame(t(out$data_Qnorm));

save(DAT0, DAT, dirs, files, out, file='_newFPKM_isoforms_rawAndPreprocData.RData');
save(DATA, file='_newFPKM_isoforms_DATA-for_net.RData');

##############################################################
###### test parameters for network construction ##############
##############################################################

rm(list=ls());
library(WGCNA); allowWGCNAThreads(); options('stringsAsFactors'=F);
load('_newFPKM_isoforms_DATA-for_net.RData');
load('_newFPKM_isoforms_rawAndPreprocData.RData')

###################################################################
##### optional: collapse to one transcript per annotated gene #####

rowGroup = gsub('.[0-9]+$', '', gsub('mz.mrna.', '', names(DATA)));
DATAcollapsed = collapseRows(t(DATA), rowGroup=rowGroup, rowID=names(DATA));
save(DATAcollapsed, file='_newFPKM_isoforms_DATAcollapsed-for_net.RData');

DATA2 = DATAcollapsed[[1]];
DATA2 = as.data.frame(t(DATA2));

annos=read.table('../../_broadftp/cichlid_geneNames_MZnoNONE.txt',header=F,sep='\t',stringsAsFactors=F,quote="");
DATA2 = DATA2[, names(DATA2) %in% gsub('mz.gene.','',annos[,1])];

powers = seq(10,18,2);
klist=list();
for (k in 1:length(powers)) {
	klist[[k]] = softConnectivity(DATA2,type='signed',power=powers[k],blockSize=6000,verbose=3);
	names(klist[[k]]) = names(DATA2);
	names(klist)[k] = paste('k',powers[k],sep='');
}; rm(powers, k);

par(mfrow=c(2,5));
for (k in 1:length(klist)) {
	hist(klist[[k]],main=names(klist)[k]);
	};rm(k);
for (k in 1:length(klist)) {
	scaleFreePlot(klist[[k]],main=names(klist)[k]);
	};rm(k);
	
source('~/Documents/_Fernald_lab/_code/blockwiseModulesEnriched-Feb2013.R');
out=blockwiseModulesEnriched(DATA2, power=10, mergeCutHeight=.2, nPerm=5000, saveFileBase='iso_collapsed_s10_mch.2');
###################################################################


##############

# impute data, DATA should have genes in rows and samples in columns
DATAim = impute.knn(as.matrix(t(DATA)), k = min(10, ncol(DATA) - 1));
DATAim = as.data.frame(t(DATAim$data));

sft = pickSoftThreshold(DATAim,networkType='signed',blockSize=6000,verbose=3);

   # Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
# 1      1 0.599000  4.5500          0.682    3110      3130   3630
# 2      2 0.000267 -0.0234         -0.204    1960      1930   2630
# 3      3 0.401000 -0.8730          0.740    1380      1300   2090
# 4      4 0.588000 -1.1000          0.792    1050       937   1760
# 5      5 0.642000 -1.1400          0.767     836       704   1540
# 6      6 0.657000 -1.0800          0.679     690       556   1370
# 7      7 0.680000 -1.0500          0.661     585       461   1250
# 8      8 0.707000 -1.0500          0.677     506       388   1150
# 9      9 0.735000 -1.0600          0.699     444       343   1070
# 10    10 0.765000 -1.0600          0.728     395       298   1000
# 11    12 0.836000 -1.0700          0.808     323       234    896
# 12    14 0.867000 -1.0900          0.837     272       199    811
# 13    16 0.881000 -1.0900          0.855     234       165    745
# 14    18 0.899000 -1.1000          0.876     206       138    691
# 15    20 0.903000 -1.1000          0.879     183       118    647

par(mfrow=c(1,3));
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '');
text(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, labels = sft$fitIndices$Power);
abline( h = 0.8, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$mean.k, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '');
text(sft$fitIndices$Power, sft$fitIndices$mean.k, labels = sft$fitIndices$Power);
abline( h = 50, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$slope, type = 'n', xlab = 'Power', ylab = 'Slope', main = '');
text(sft$fitIndices$Power, sft$fitIndices$slope, labels = sft$fitIndices$Power);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

# compare potential power values
k12 = softConnectivity(DATAim,type='signed',power=12,blockSize=6000,verbose=3);
names(k12) = names(DATAim);
k16 = softConnectivity(DATAim,type='signed',power=16,blockSize=6000,verbose=3);
names(k16) = names(DATAim);
k17 = softConnectivity(DATAim,type='signed',power=17,blockSize=6000,verbose=3);
names(k17) = names(DATAim);
k18 = softConnectivity(DATAim,type='signed',power=18,blockSize=6000,verbose=3);
names(k18) = names(DATAim);

par(mfrow=c(2,4));
hist(k12); hist(k16); hist(k17); hist(k18);
scaleFreePlot(k12); scaleFreePlot(k16); scaleFreePlot(k17); scaleFreePlot(k18);

# try power=18
source('../../../_code/blockwiseModulesEnriched-Feb2013.R');
blockOut=blockwiseModulesEnriched(DATAim, power=18, mergeCutHeight=.2, nPerm=2000, skipThresh=300, saveFileBase='_newFPKM_isoforms_impute_p18mch.2');

######################
# or try removing any transcripts with NAs
numNA = apply(DATA, 2, function(x){sum(is.na(x))});
numNA = numNA > 0;
DATAnoNA = DATA[, !numNA];
sft = pickSoftThreshold(DATAnoNA, networkType='signed', blockSize=5000, verbose=3);

 # pickSoftThreshold: calculating connectivity for given powers...
   # ..working on genes 1 through 4640 of  4640
   # Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.6000  4.210          0.680    2440    2440.0   2840
# 2      2   0.0192 -0.160         -0.146    1540    1500.0   2060
# 3      3   0.4870 -0.953          0.614    1080    1010.0   1640
# 4      4   0.6060 -1.130          0.700     814     722.0   1380
# 5      5   0.6680 -1.230          0.754     643     543.0   1200
# 6      6   0.7120 -1.260          0.792     525     425.0   1070
# 7      7   0.7320 -1.280          0.809     441     345.0    977
# 8      8   0.7470 -1.310          0.814     377     286.0    900
# 9      9   0.7600 -1.310          0.831     328     241.0    837
# 10    10   0.7660 -1.310          0.838     290     206.0    784
# 11    12   0.7850 -1.310          0.867     232     156.0    704
# 12    14   0.7820 -1.300          0.885     193     122.0    643
# 13    16   0.7790 -1.300          0.891     164      98.1    595
# 14    18   0.7820 -1.300          0.899     142      80.8    556
# 15    20   0.7800 -1.290          0.909     125      67.6    524

par(mfrow=c(1,3));
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '');
text(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, labels = sft$fitIndices$Power);
abline( h = 0.8, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$mean.k, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '');
text(sft$fitIndices$Power, sft$fitIndices$mean.k, labels = sft$fitIndices$Power);
abline( h = 50, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$slope, type = 'n', xlab = 'Power', ylab = 'Slope', main = '');
text(sft$fitIndices$Power, sft$fitIndices$slope, labels = sft$fitIndices$Power);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

# compare potential power values
k15 = softConnectivity(DATA,type='signed',power=15,blockSize=6000,verbose=3);
names(k15) = names(DATA);
k16 = softConnectivity(DATA,type='signed',power=16,blockSize=6000,verbose=3);
names(k16) = names(DATA);
k17 = softConnectivity(DATA,type='signed',power=17,blockSize=6000,verbose=3);
names(k17) = names(DATA);
k18 = softConnectivity(DATA,type='signed',power=18,blockSize=6000,verbose=3);
names(k18) = names(DATA);

par(mfrow=c(2,4));
hist(k15); hist(k16); hist(k17); hist(k18);
scaleFreePlot(k15); scaleFreePlot(k16); scaleFreePlot(k17); scaleFreePlot(k18);

# try power=18
source('../../../_code/blockwiseModulesEnriched-Feb2013.R');
blockOut=blockwiseModulesEnriched(DATAnoNA, power=18, mergeCutHeight=.2, nPerm=5000, skipThresh=300, saveFileBase='_newFPKM_isoforms_noNA_p18mch.2');

load('_newFPKM_isoforms_noNA_p18mch.2run2DATA.RData');
DATAnoNA = DATA;

# to get quality stats
setLabels=c('hybrid1', 'hybrid2');
multiExpr=list();
multiExpr[[1]]=list(data=DATAnoNA);
multiExpr[[2]]=list(data=DATAnoNA);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=blockOut$net$colors;
colorList[[2]]=blockOut$net$colors;
names(colorList)=setLabels;

mp = modulePreservation(multiExpr, colorList, networkType='signed',nPermutations=100,verbose=3);
save(mp, file='_newFPKM_isoforms_noNA_p18mch.2_NET2_modulePreservationQuality.RData');

##############

# test powers for network construction
sft = pickSoftThreshold(DATA,networkType='signed',blockSize=6000,verbose=3);

   # Power   SFT.R.sq      slope truncated.R.sq   mean.k. median.k.   max.k.
# 1      1 0.59690454  3.9327881      0.6922037 3123.5833 3160.1076 3691.093
# 2      2 0.09627724 -0.5163478     -0.0744732 2007.3272 1980.4442 2787.693
# 3      3 0.46461833 -1.3068575      0.3660897 1447.5283 1368.1308 2318.152
# 4      4 0.62205840 -1.5433551      0.5763017 1122.0500 1007.7588 2043.932
# 5      5 0.66690427 -1.7221500      0.6637221  912.4955  795.3685 1855.030
# 6      6 0.69264587 -1.7300419      0.7328929  767.4465  645.2627 1718.890
# 7      7 0.71385545 -1.6924196      0.7911635  661.5934  541.5945 1611.573
# 8      8 0.73126833 -1.6543501      0.8254915  581.2014  473.4698 1523.587
# 9      9 0.74035921 -1.6387975      0.8441033  518.2232  417.0179 1449.425
# 10    10 0.75676432 -1.6093068      0.8665518  467.6494  373.3128 1385.623
# 11    12 0.77314254 -1.5960949      0.8865635  391.6658  298.5990 1280.517
# 12    14 0.76889451 -1.6354075      0.8736211  337.4771  244.3332 1196.616
# 13    16 0.77455700 -1.6174240      0.8814763  296.9961  208.0185 1127.457
# 14    18 0.77651719 -1.6007120      0.8912035  265.6683  180.2475 1069.092
# 15    20 0.77413128 -1.6086133      0.8840928  240.7420  157.3970 1018.937

par(mfrow=c(1,3));
plot(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '');
text(sft$fitIndices$Power, sft$fitIndices$SFT.R.sq, labels = sft$fitIndices$Power);
abline( h = 0.8, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$mean.k, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '');
text(sft$fitIndices$Power, sft$fitIndices$mean.k, labels = sft$fitIndices$Power);
abline( h = 50, col = 'red', lty = 'dashed' );

plot(sft$fitIndices$Power, sft$fitIndices$slope, type = 'n', xlab = 'Power', ylab = 'Slope', main = '');
text(sft$fitIndices$Power, sft$fitIndices$slope, labels = sft$fitIndices$Power);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

# compare potential power values
k15 = softConnectivity(DATA,type='signed',power=15,blockSize=6000,verbose=3);
names(k15) = names(DATA);
k16 = softConnectivity(DATA,type='signed',power=16,blockSize=6000,verbose=3);
names(k16) = names(DATA);
k17 = softConnectivity(DATA,type='signed',power=17,blockSize=6000,verbose=3);
names(k17) = names(DATA);
k18 = softConnectivity(DATA,type='signed',power=18,blockSize=6000,verbose=3);
names(k18) = names(DATA);

par(mfrow=c(2,4));
hist(k15); hist(k16); hist(k17); hist(k18);
scaleFreePlot(k15); scaleFreePlot(k16); scaleFreePlot(k17); scaleFreePlot(k18);

save(k15, k16, k17, k18, file='_newFPKM_isoforms_k15-18.RData');
rm(k15, k16, k17, k18);
collectGarbage();
##############################################################
###### build network #########################################
##############################################################

##### just use blockwiseModulesEnriched #####


# # try beta=18
# beta=18;
# type='signed';
# net0 = blockwiseModules(DATA,
						   # maxBlockSize = 6000,
						   # power = beta,
						   # networkType = type,
						   # deepSplit = 2,
						   # minModuleSize = 10,
						   # verbose = 3,
						   # saveTOMs = T,
						   # saveTOMFileBase = '_newFPKM_isoforms_run0TOM'
						   # );
# collectGarbage();
# save(net0, file='_newFPKM_isoforms_net0_s18.RData');

# # plot dendrograms
block = 1;
plotDendroAndColors(net0$dendrograms[[block]],
					net0$colors[ net0$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = net0$colors[ net0$blockGenes[[block]] ],
					main = paste('block', block),
					dendroLabels = F,
					#hang = 0.03,
					addGuide = T,
					guideHang = 0.05)
rm(block); collectGarbage();


# source('~/Documents/_analysis_compare/_code/getModDensitiesFromTOM.R');

# TOM0 = TOMsimilarityFromExpr(DATA,
							   # networkType = type,
						       # power = beta,
				   		       # verbose = 3
							   # );
# collectGarbage();
# save(TOM0, file='_newFPKM_isoforms_run0TOMmatrix-block.1.RData');

# modDensities = getModDensitiesFromTOM(TOM0, net0$colors);
# system.time( {		# probably take ~5min with nPerm=1000
	# permTest = modDensityPerm( TOM0, net0$colors, modDensities );
	# } );
# save(permTest, file='_newFPKM_isoforms_run0_permTest.RData');

# weakTOM = names( table(net0$colors)[ permTest$pvals > 0.01 ] );
# bgGenes = net0$colors %in% weakTOM;
# DATA1 = DATA[ , !bgGenes];
# save(DATA1, file='_newFPKM_isoforms_DATA1.RData');

##############################################################
###### examine modules #########################################
##############################################################
rm(list=ls());
setwd('~/Documents/_Fernald_lab/_hybridRNAseq/_newFPKM');
library(WGCNA); allowWGCNAThreads();
library(org.Hs.eg.db);
options('stringsAsFactors'=F);
source('~/Documents/_Fernald_lab/_code/exploreNetwork.R');

load('_newFPKM_isoforms_blockwiseModulesEnriched.run4DATA.RData');
load('_newFPKM_isoforms_blockwiseModulesEnriched.run4NET.RData');
annos=read.table('../../_broadftp/cichlid_geneNames_MZnoNONE.txt',header=F,sep='\t',stringsAsFactors=F,quote="");
annos[annos[,1]=='mz.gene.s6.92', 3] = 'SPINC';
IDs=as.list(org.Hs.egSYMBOL2EG);

mod = exn.buildModuleGeneInfo(DATA, net$colors, annos, IDs);
tr = exn.transcriptsAcrossModules(names(DATA), net$colors);
kMEs = exn.computekME(DATA, net$MEs);

###########################################
# to get quality stats
setLabels=c('hybrid1', 'hybrid2');
multiExpr=list();
multiExpr[[1]]=list(data=DATA);
multiExpr[[2]]=list(data=DATA);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=net$colors;
colorList[[2]]=net$colors;
names(colorList)=setLabels;

mp = modulePreservation(multiExpr, colorList, networkType='signed',nPermutations=100,verbose=3)
save(mp, file='_newFPKM_isoforms_modulePreservationQuality.RData');
####################################################

# getting GO terms using biomaRt

library(biomaRt);
library(GO.db);

# fish=c('gaculeatus_gene_ensembl', 'drerio_gene_ensembl', 'olatipes_gene_ensembl', 'tnigroviridis_gene_ensembl');
# marts = vector(mode='list', length=length(fish));
# for (f in 1:length(fish)) {
	# marts[[f]] = useMart('ensembl', dataset=fish[f])
# }; rm(f);
# names(marts) = c('gac', 'dar', 'orl', 'tni');
# hsa = useMart('ensembl','hsapiens_gene_ensembl');
# save(fish, marts, hsa, file='_newFPKM_isoforms_marts.RData');
load('_newFPKM_isoforms_marts.RData');


# allAnnosIDs = exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork(DATA,annos,IDs,F);
# allGO = getBM(attributes=c('go_id'), filters='entrezgene', values=unlist(allAnnosIDs$allGeneIDs), hsa);
# modGO = exn.getModGOfromEntrez(mod$geneIDs);
# save(allAnnosIDs, allGO, file='_newFPKM_isoforms_allAnnosAndGOforNetwork.RData');
# save(modGO, file='_newFPKM_isoforms_GOtermsByModule.RData');
load('_newFPKM_isoforms_allAnnosAndGOforNetwork.RData');
load('_newFPKM_isoforms_GOtermsByModule.RData');

# get GO terms by gene
# GOterms = as.list(GOTERM);
# save(GOterms, file='GOterms.RData');
load('GOterms.RData');

#allGenesGO = exn.getGOfromEntrezByGene(unlist(allAnnosIDs$allGeneIDs), GOterms, write=T);

####################################################
##### SNPs #####
####################################################

snps0 = exn.getIDsFromGTF('../CichlidBAMs/Metriaclima_zebra.BROADMZ2cp.V2.gtfSNPsASE.V3');


x=exn.simulateAlleleFreqs(DATA[,1:1000],10);
# net0=blockwiseModules(x$DATA,networkType='signed',power=18);
# net0new=blockwiseModules(x$DATAnew,networkType='signed',power=18);

# block = 1;
# plotDendroAndColors(net0$dendrograms[[block]],
					# net0$colors[ net0$blockGenes[[block]] ],
					# groupLabels = 'module',
					# rowText = net0$colors[ net0$blockGenes[[block]] ],
					# main = paste('block', block),
					# dendroLabels = F,
					# #hang = 0.03,
					# addGuide = T,
					# guideHang = 0.05)
# plotDendroAndColors(net0new$dendrograms[[block]],
					# net0new$colors[ net0new$blockGenes[[block]] ],
					# groupLabels = 'module',
					# rowText = net0new$colors[ net0new$blockGenes[[block]] ],
					# main = paste('block', block),
					# dendroLabels = F,
					# #hang = 0.03,
					# addGuide = T,
					# guideHang = 0.05)
# rm(block);

# names(x$DATA)[names(x$DATA) %in% x$tr];
# names(x$DATAnew)[grep('[AB]$', names(x$DATAnew))];
# net0$colors[names(x$DATA) %in% x$tr];
# net0new$colors[grepl('[AB]$', names(x$DATAnew))];

# as.data.frame(cbind(names(x$DATA)[names(x$DATA) %in% x$tr], net0$colors[names(x$DATA) %in% x$tr]));
# as.data.frame(cbind(names(x$DATAnew)[grep('[AB]$', names(x$DATAnew))], net0new$colors[grepl('[AB]$', names(x$DATAnew))]));


overlap = overlapTable(net1$colors, net2$colors);
 
numMat = -log10(overlap$pTable);
numMat[numMat > 50] = 50;
 
textMat = paste(overlap$countTable, '\n', signif(overlap$pTable, 2));
dim(textMat) = dim(numMat);
 
xlabels = paste('M', sort(unique(net2$colors)));
ylabels = paste('M', sort(unique(net1$colors)));
xSymbols = paste(sort(unique(net2$colors)), ': ', table(net2$colors), sep = '');
ySymbols = paste(sort(unique(net1$colors)), ': ', table(net1$colors), sep = '');

sizeGrWindow(7, 7); fp = FALSE
fcex = 1.00;
pcex = .7;
fcex1 = .7;
pcex1 = 1.00;
par(mar = c(6, 7, 2, 1.0));
labeledHeatmap(Matrix = numMat,
			   xLabels = xlabels, xSymbols = xSymbols,
			   yLabels = ylabels, ySymbols = ySymbols,
			   colorLabels = T, 
			   colors = greenWhiteRed(100)[50:100],
			   textMatrix = textMat, cex.text = pcex, setStdMargins = F,
			   cex.lab = fcex1,
			   xColorWidth = 0.08,
			   main = ''
			   );

####################################################
##### traits #####
####################################################
avgH0 = read.table('Total Entropy.txt', header=F, row.names=1);
avgH = avgH0[, 1];
names(avgH) = rownames(avgH0);
avgH = c(avgH[1:3], CV2=NA, avgH[4:5]);
avgH = avgH[order(names(avgH))];

bH0 = read.table('Behavioral Entropy.txt', header=T, row.names=1);
bH = rbind(bH0[1:3, ], rep(NA, 7), bH0[4:5, ]);
rownames(bH)[4] = 'CV2';
bH = bH[order(rownames(bH)), ];

files = list.files()[grepl('tMat.txt', list.files())];
tMats = list();
for (f in 1:length(files)) {
	tMats[[f]] = read.table(files[f], header=T, row.names=1);
	names(tMats)[f] = gsub('tMat.txt', '', files[f]);
}; rm(f);
tMats = c(tMats[1], tMats);
tMats[[2]] = NA;
names(tMats)[2] = 'CV2';

####################################################
##### DAVID #####
####################################################

all = exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork(DATA,annos,IDs,F);

library(RDAVIDWebService);

#david = exn.startAndUploadBgAndModListsToDAVID();
#modDAVID0 = exn.getModChartsFromDAVID(david);
#save(modDAVID0, file='_newFPKM_isoforms_rawDAVIDcharts.RData');
load('_newFPKM_isoforms_rawDAVIDcharts.RData');
modDAVID = exn.addGeneInfoToDAVIDChartList(modDAVID0);
termDAVID = exn.modDAVIDByTerm(modDAVID);

###################################################
##### simulated allele networks #####
###################################################
source('~/Documents/_Fernald_lab/_code/exploreNetwork.R');
sim = exn.simulateAlleleFreqs(DATA=DATA[,1:500],nTranscripts=10,rSamples=c(5,6),jitter=.2);
ans = exn.getAlleleNetworks(sim$DATA.padded,sim$DATA.alleles,annos,IDs);
overlap = exn.plotModuleOverlaps(ans$padded$net$colors,ans$alleles$net$colors);



