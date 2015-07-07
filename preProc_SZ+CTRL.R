rm(list=ls());
setwd('/Volumes/OSX/_macbookPro_athBACKUP_05-16-13/Documents/_Fernald_lab');
source('_code/preProcATH-for_web_noVSN.R');

load('SZ_DATA_LOOKUPs.RData');

lapply(DATAlist, dim); lapply(LOOKUPlist, dim);
summary(DATAlist); summary(LOOKUPlist);

dat = as.data.frame(cbind(DATAlist$dat.SZ, DATAlist$dat.CTRL));

outNAs=removeTooManyNAs(dat, probe_thresh=ncol(dat)/2);

outTrans=removeOutlierProbesIterate(outNAs$dataClean, deviate=3);

outNAs2=removeTooManyNAs(outTrans$dataClean, probe_thresh=ncol(dat)/2);

outSamps=outlierSamplesIterate(outNAs2$dataClean);

outNAs3=removeTooManyNAs(outSamps$dataClean, probe_thresh=ncol(dat)/2);

datQnorm = as.data.frame(normalize.quantiles(as.matrix(outNAs3$dataClean)));
dimnames(datQnorm) = dimnames(outNAs3$dataClean);

pickSoft = pickSoftThreshold(as.data.frame(t(datQnorm)), networkType='signed', blockSize=6000, verbose=2);

 # pickSoftThreshold: calculating connectivity for given powers...
   # ..working on genes 1 through 6000 of  27018
   # ..working on genes 6001 through 12000 of  27018
   # ..working on genes 12001 through 18000 of  27018
   # ..working on genes 18001 through 24000 of  27018
   # ..working on genes 24001 through 27018 of  27018
   # Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1  0.04030  5.170          0.901 13800.0   13800.0  15100
# 2      2  0.10600  4.350          0.930  7610.0    7530.0   9320
# 3      3  0.04330  1.520          0.902  4460.0    4410.0   6300
# 4      4  0.00226 -0.212          0.895  2760.0    2710.0   4550
# 5      5  0.11100 -1.050          0.922  1780.0    1740.0   3420
# 6      6  0.30500 -1.450          0.934  1200.0    1160.0   2660
# 7      7  0.45100 -1.620          0.923   837.0     790.0   2120
# 8      8  0.53100 -1.670          0.905   602.0     553.0   1720
# 9      9  0.57700 -1.650          0.887   445.0     394.0   1420
# 10    10  0.62500 -1.580          0.889   336.0     285.0   1190
# 11    12  0.74100 -1.440          0.942   205.0     156.0    867
# 12    14  0.74800 -1.570          0.942   133.0      89.5    701
# 13    16  0.77400 -1.620          0.959    90.8      53.1    581
# 14    18  0.79900 -1.660          0.971    64.5      32.4    490
# 15    20  0.81300 -1.690          0.973    47.3      20.3    418

k14 = softConnectivity(as.data.frame(t(datQnorm)), type='signed', blockSize=6000, verbose=3);

save.image(file='preProc_SZ+CTRL_WORKSPACE.RData');
##

DATA = as.data.frame(t(datQnorm));
net = blockwiseModules(DATA,
						maxBlockSize=6000,
						power=14,
						networkType='signed',
						deepSplit=2,
						minModuleSize=100,
						verbose=3
						);

collectGarbage(); 						
sort(table(net$colors));
						   
block = 1;
plotDendroAndColors(net$dendrograms[[block]],
					net$colors[ net$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = net$colors[ net$blockGenes[[block]] ],
					main = paste('block', block),
					dendroLabels = F,
					#hang = 0.03,
					addGuide = T,
					guideHang = 0.05);
rm(block); 


MEs = moduleEigengenes(DATA, net$colors);
MEs = orderMEs(MEs$eigengenes)
kME = as.data.frame(cor(DATA, MEs, use='p'));
GS = as.data.frame(cor(DATA, vec, use='p'));

lookup=LOOKUPlist$lookup0;
lookup=lookup[match(rownames(DATA), rownames(lookup)), ];
vec=rep(0,nrow(lookup));
vec[lookup[,7]=='SZ']=1;

MEcor = cor(MEs, vec, use='p');
MEcorP = corPvalueFisher(MEcor, nrow(DATA));


salmonGenes = names(DATA)[net$colors=='salmon'];
salmonkME = kME[rownames(kME) %in% salmongenes,];
salmonkME=salmonkME[order(salmonkME$MEsalmon, decreasing=T),];

salmonGenesTrim = salmonGenes;
tmp = strsplit(salmonGenes,'|',fixed=T);
for (gene in 1:length(tmp))
{
	salmonGenesTrim[gene] = tmp[[gene]][1];
} 
rm(gene);
write.table(salmonGenesTrim,'salmonGenesTEST.txt',quote=F,row.names=F,col.names=F);