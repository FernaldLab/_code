# clear workspace and set working directory, yours is probably different
rm( list = ls() );
setwd('~/Documents/_analysis_compare');

# load WGCNA functions and set options
library('WGCNA');
options(stringsAsFactors = F);

load(file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed');	# 
load(file = '_Rdata/VSPprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_datVSPfilt11-usedINnetVSPfinal10');

datVSP = datVSPfilt11;
tmp = datVSPcollapsedList$deviate3$group2row;
genes = tmp[ match( names(datVSP), tmp[, 1] ), 2];
rm(datVSPfilt11, datVSPcollapsedList, tmp);
sum( names(datVSP) == names(genes) ) == ncol(datVSP); # should be TRUE


load(file = '_Rdata/XnewMBPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm');
datX = outQnorm$deviate3;
rm(outQnorm);
datX = datX[ match( genes, rownames(datX) ) , ];
rownames(datX) = names(genes);
datX = as.data.frame( t(datX) );


######
#####
######

load(file = '_Rdata/VSPprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_NETfinal10_perm5000p0.01mod10');

netX = blockwiseModules(datX,
						maxBlockSize = 6000,
						power = 14,
						networkType = 'signed',
						deepSplit = 2,
						minModuleSize = 10,
						verbose = 3);
collectGarbage();
#save(netX, file = '_Rdata/Xprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_NETX-withVSPgenes');

netX$colors = matchLabels(netX$colors, netVSPfinal10$colors);

# plot dendrograms
block = 1;
plotDendroAndColors(netX$dendrograms[[block]],
					netX$colors[ netX$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = netX$colors[ netX$blockGenes[[block]] ],
					main = paste('block', block),
					dendroLabels = F,
					#hang = 0.03,
					addGuide = T,
					guideHang = 0.05)
rm(block);



#######
#####
#######

setLabels=c('VSP','X');
multiExpr=list();
multiExpr[[1]]=list(data=datVSP);
multiExpr[[2]]=list(data=datX);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=netVSPfinal10$colors;
colorList[[2]]=netX$colors;
names(colorList)=setLabels;

system.time( {				# should only take a minute or two
 	
 	mp = modulePreservation(multiExpr, 
 				   colorList,
 				   networkType = 'signed', 
 				   nPermutations = 0,			# no stats, just get ranks
 				   maxModuleSize = max( table(netVSPfinal10$colors) ),
 				   verbose = 3,
 				   savePermutedStatistics = F
 				   );
 	}
 )

save(mp, file='_Rdata/modulePreservation_netVSPfinal10');
# compute ranks
tmp = cbind(mp$observed$VSP$intra$VSP_vs_X, mp$observed$VSP$inter$VSP_vs_X);
tmp = tmp[!(rownames(tmp)=='gold'), ];
tmp2 = apply(-tmp, 2, rank);
 
# con - cor.kIM, cor.kME, cor.cor
medRankCon = apply( tmp2[ , c(8, 9, 11) ], 1, median);
# density - propVarExplained, meanSignAwareKME, meanSignAwareCorDat, meanAdj
medRankDen = apply(tmp2[ , c(1, 2, 4, 5) ], 1, median);
 
medRankPres = (medRankCon + medRankDen) / 2;

medRanks = cbind(medRankPres,medRankCon,medRankDen);
medRanks = medRanks[order(medRanks[,1]),];
 
####
 
overlap = overlapTable(netVSPfinal10$colors, netX$colors);
 
numMat = -log10(overlap$pTable);
numMat[numMat > 50] = 50;
 
textMat = paste(overlap$countTable, '\n', signif(overlap$pTable, 2));
dim(textMat) = dim(numMat);
 
xlabels = paste('M', sort(unique(netX$colors)));
ylabels = paste('M', sort(unique(netVSPfinal10$colors)));
xSymbols = paste(sort(unique(netX$colors)), ': ', table(netX$colors), sep = '');
ySymbols = paste(sort(unique(netVSPfinal10$colors)), ': ', table(netVSPfinal10$colors), sep = '');

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