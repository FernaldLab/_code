rm(list = ls());
setwd('~/Documents/_analysis_compare');
library('WGCNA');
options(stringsAsFactors=F);

# load VSP network and data
load('VSP_3sd_1S3S_Qnorm_collapsed_NET27final.RData');
load('VSP_3sd_1S3S_Qnorm_collapsed_DATA27.RData');

# load area X data and network built using same settings as in VSP
load('VSP_3sd_1S3S_Qnorm_collapsed_NET27-areaXgenes.RData');
load('VSP_3sd_1S3S_Qnorm_collapsed_DATA27-areaX.RData');

#####

# re-name for consistency below
DATA.VSP = DATA27;
NET.VSP = NET27;
rm(DATA27, NET27);

# match area X module colors to those in VSP as closely as possible
NET.X$colors = matchLabels(NET.X$colors, NET.VSP$colors);

# re-name important aspects of network and define other variable for below
dendrogramVSP = NET.VSP$dendrograms[[1]];
MEsVSP = NET.VSP$MEs;
colorsVSP = NET.VSP$colors;
nGenesVSP = ncol(DATA.VSP); 
nSamplesVSP = nrow(DATA.VSP);

dendrogramX = NET.X$dendrograms[[1]];
MEsX = NET.X$MEs;
colorsX = NET.X$colors;
nGenesX = ncol(DATA.X); 
nSamplesX = nrow(DATA.X);

# compare module assignments in area X to VSP dendrogram
plotDendroAndColors(dendrogramVSP, 
					data.frame(colorsVSP, colorsX), 
					dendroLabels = F, 
					addGuide = T, 
					guideHang = 0.1, 
					groupLabels = c('VSP', 'X'), 
					rowText = data.frame(colorsVSP, colorsX), 
					main = paste('VSP: ', length(colorsVSP), ' genes', sep = '') 
					);
					
# re-compute eigengenes to make sure colors are consistent
MEs0VSP = moduleEigengenes( DATA.VSP, colorsVSP );
MEsVSP = orderMEs( MEs0VSP$eigengenes );
modNamesVSP = substring( names(MEsVSP), 3 ); 

# compute module membership (MM) in VSP
geneMM.VSP = as.data.frame( cor( DATA.VSP, MEsVSP, use = 'p' ) ); 
MMPvalueVSP = as.data.frame( corPvalueStudent( as.matrix(geneMM.VSP), nSamplesVSP ) ); 
names(geneMM.VSP) = paste('MM', modNamesVSP, sep = ''); 
names(MMPvalueVSP) = paste('p.MM', modNamesVSP, sep = '');

# compute and plot eigengene network
# function plotMEpairs() is also useful
MET.VSP = orderMEs(MEsVSP);
plotEigengeneNetworks(MET.VSP, 
					  setLabels = '',
					  marDendro = c(0, 4, 1, 2),
					  marHeatmap = c(3, 4, 1, 2)
					  );
					  
# re-compute eigengenes to make sure colors are consistent
MEs0X = moduleEigengenes( DATA.X, colorsX );
MEsX = orderMEs( MEs0X$eigengenes );
modNamesX = substring( names(MEsX), 3 ); 

# compute MM in area X
geneMM.X = as.data.frame( cor( DATA.X, MEsX, use = 'p' ) ); 
MMPvalueX = as.data.frame( corPvalueStudent( as.matrix(geneMM.X), nSamplesX ) ); 
names(geneMM.X) = paste('MM', modNamesX, sep = ''); 
names(MMPvalueX) = paste('p.MM', modNamesX, sep = '');

# compute and plot area X eigengene network
MET.X = orderMEs(MEsX);
plotEigengeneNetworks(MET.X, 
					  setLabels = '',
					  marDendro = c(0, 4, 1, 2),
					  marHeatmap = c(3, 4, 1, 2)
					  );
#####

# plot heatmap of module
mod = 'blue';
topOnly = T;	geneMM = geneMM.VSP; colors = colorsVSP; num = 100;	# if T, need to specify extras
DAT = DATA.VSP;

if( topOnly ) {
	modGenes = colors == mod;
	if( sum(modGenes) < num ) { num = sum(modGenes) }
	modGenes = rownames(geneMM)[modGenes];
	modMM = geneMM[modGenes, ];
	MMcol = match(mod, substring(names(geneMM), 3) );
	modMM = modMM[order( -modMM[, MMcol] ), ];
	modMM = modMM[1:num, ];
	modGenes = rownames(modMM);
	rm(geneMM, colors, num, modMM, MMcol);
} else {
	modGenes = colorsVSP == mod;
	modGenes = names(DAT)[modGenes];
}
# for subset of module genes, where subset is vector of gene symbols
# modGenes = modGenes[modGenes %in% subset];
modDAT = DAT[, match(modGenes, names(DAT))];
modTOM = TOMsimilarityFromExpr(modDAT, networkType = 'signed', power = 14);
modTOMdist = 1 - modTOM;
diag(modTOMdist) = NA;
if( ncol(modTOMdist) > 100 ) {labels = ''} else {labels = names(modDAT)}
heatmap(modTOMdist, 
		symm = T, 
		revC = T, 
		labRow = labels, 
		labCol = labels,
		cexRow = 0.5, cexCol = 0.5
		);	
rm(mod, modGenes, modDAT, modTOM, modTOMdist);

#####
# trait correlations
			  
# load singing and age data
traits = c('age', 
		   'SingingVsQuiet', 
		   'motifs2', 
		   'mean.amplitude', 
		   'mean.pitch', 
		   'mean.FM', 
		   'mean.AM.2', 
		   'mean.entropy', 
		   'mean.pitch.goodness',
		   'mean.mean.freq'
		   );
traits0 = read.csv('../_writeUp/traits_inclFeatures.csv', row.names = 1);
traitsVSP = traits0;
rownames(traitsVSP) = gsub('X', 'S', rownames(traitsVSP));
samplesVSP = rownames(DATA.VSP);
traitsVSP = traitsVSP[ match( samplesVSP, rownames(traitsVSP) ), ];
traitsVSP = traitsVSP[ , match( traits, names(traitsVSP) )];

traitsX = traits0;
samplesX = rownames(DATA.X);
traitsX = traitsX[ match( samplesX, rownames(traitsX) ), ];
traitsX = traitsX[ , match( traits, names(traitsX) )];

# compute ME-trait cors in VSP
moduleTraitCorVSP = cor( MEsVSP, traitsVSP, use = 'p' ); 
moduleTraitPvalueVSP = corPvalueStudent( moduleTraitCorVSP, nSamplesVSP );
textMatrixVSP = paste( signif(moduleTraitCorVSP, 2), '\n(', 
					   signif(moduleTraitPvalueVSP, 1), ')', sep = ''
					   ); 
dim(textMatrixVSP) = dim(moduleTraitCorVSP);

par( mar = c(8, 9, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCorVSP,
				xLabels = names(traitsVSP),
				yLabels = names(MEsVSP),
				ySymbols = names(MEsVSP),
				colorLabels = F,
				colors = greenWhiteRed(50),
				textMatrix = textMatrixVSP,
				setStdMargins = F,
				cex.text = 0.6,
				zlim = c(-1, 1),
				main = paste("VSP module-trait relationships")
				); 
sigVSP = moduleTraitPvalueVSP < (.05 / (ncol(moduleTraitPvalueVSP) * nrow(moduleTraitPvalueVSP)));

# compute ME-trait cors in X
moduleTraitCorX = cor( MEsX, traitsX, use = 'p' ); 
moduleTraitPvalueX = corPvalueStudent( moduleTraitCorX, nSamplesX );
textMatrixX = paste( signif(moduleTraitCorX, 2), '\n(', 
					   signif(moduleTraitPvalueX, 1), ')', sep = ''
					   ); 
dim(textMatrixX) = dim(moduleTraitCorX);

par( mar = c(8, 9, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCorX,
				xLabels = names(traitsX),
				yLabels = names(MEsX),
				ySymbols = names(MEsX),
				colorLabels = F,
				colors = greenWhiteRed(50),
				textMatrix = textMatrixX,
				setStdMargins = F,
				cex.text = 0.6,
				zlim = c(-1, 1),
				main = paste("X module-trait relationships")
				); 
sigX = moduleTraitPvalueX < (.05 / (ncol(moduleTraitPvalueX) * nrow(moduleTraitPvalueX)));

# compute GS in VSP, add qvals
for( t in 1:ncol(traitsVSP) ) {
	
	cat( 'working on', names(traitsVSP)[t], '\n' );
	temp = as.data.frame( cor( DATA.VSP, traitsVSP[ , t ], use = 'p' ) );
	tempPval = as.data.frame( corPvalueStudent( as.matrix(temp), nrow(DATA.VSP) ) );
	tempQval = as.data.frame( qvalue( tempPval[ , 1] )$qvalues );
	names(temp) = paste( 'GS.', names(traitsVSP)[t], '.VSP', sep='' );
	names(tempPval) = paste( 'p.GS.', names(traitsVSP)[t], '.VSP', sep='' );
	names(tempQval) = paste( 'q.GS.', names(traitsVSP)[t], '.VSP', sep='' );
	assign(names(temp), temp);
	assign(names(tempPval), tempPval);
	assign(names(tempQval), tempQval);
	
	cat('# of pvals < 0.05 ...', sum( tempPval[, 1] < 0.05 ), '\n');
	cat('# of qvals < 0.05 ...', sum( tempQval[, 1] < 0.05 ), '\n');
	
	}; 
	rm(t, temp, tempPval, tempQval);
	

# compute GS in X, add qvals
for( t in 1:ncol(traitsX) ) {
	
	cat( 'working on', names(traitsX)[t], '\n' );
	temp = as.data.frame( cor( DATA.X, traitsX[ , t ], use = 'p' ) );
	tempPval = as.data.frame( corPvalueStudent( as.matrix(temp), nrow(DATA.X) ) );
	tempQval = as.data.frame( qvalue( tempPval[ , 1] )$qvalues );
	names(temp) = paste( 'GS.', names(traitsX)[t], '.X', sep='' );
	names(tempPval) = paste( 'p.GS.', names(traitsX)[t], '.X', sep='' );
	names(tempQval) = paste( 'q.GS.', names(traitsX)[t], '.X', sep='' );
	assign(names(temp), temp);
	assign(names(tempPval), tempPval);
	assign(names(tempQval), tempQval);
	
	cat('# of pvals < 0.05 ...', sum( tempPval[, 1] < 0.05 ), '\n');
	cat('# of qvals < 0.05 ...', sum( tempQval[, 1] < 0.05 ), '\n');
	
	}; 
	rm(t, temp, tempPval, tempQval);	
	
#####
# plotting examples

#view relationship of MM to traits for single module
par( mfrow = c(2, 4) );
module = 'lightcyan';
column = match( module, modNamesVSP );
moduleGenes = colorsVSP == module;
for( t in 1:ncol(traitsVSP) ) {
	temp = get( paste('GS.', names(traitsVSP)[t], '.VSP', sep='') );
	verboseScatterplot(geneMM.VSP[moduleGenes, column],
					   temp[moduleGenes, 1],
					   abline = T,
					   abline.color = 'red',
					   abline.lty = 'dashed',
				   	   col = 'black',
				   	   bg = colorsVSP[moduleGenes],
				   	   pch = 21,
				   	   cex = 1.5,
					   xlab = paste('MM.', module, '.VSP', sep = ''),
					   ylab = names(temp),
					   xlim=c(0,1),
					   ylim=c(-1,1),
					   corOptions = "method = 'p'"
					   );
	abline(h=0, col='darkgrey');
	#abline(v=0, col='darkgrey');
	};
	rm(module, column, moduleGenes, t, temp);

# GS against MM for all modules	
par( mfrow = c(3, 6) );
temp0 = 'GS.mean.pitch.VSP';
temp = get(temp0);
for( m in 1:length(modNamesVSP) ) {
	module = modNamesVSP[m];
	moduleGenes = colorsVSP == module;
	MMcol = match( paste('MM', module, sep = ''), names(geneMM.VSP) );
	verboseScatterplot(geneMM.VSP[moduleGenes, MMcol],
						temp[moduleGenes, 1],
						abline = T,
						abline.color = 'red',
					    abline.lty = 'dashed',
						bg = colorsVSP[moduleGenes],
						col='black',
						pch = 21,
				   	    cex = 1,
						xlab = paste('MM.', module, sep = ''),
						ylab = temp0,
						xlim=c(0,1),
					    ylim=c(-1,1),
						);
	abline(h=0, col='darkgrey');
	
	};
	rm(temp, temp0, moduleGenes, module, MMcol);
	
#####
# module preservation

setLabels=c('VSP','X');
multiExpr=list();
multiExpr[[1]]=list(data=DATA.VSP);
multiExpr[[2]]=list(data=DATA.X);
names(multiExpr)=setLabels;

colorList=list();
colorList[[1]]=colorsVSP;
colorList[[2]]=colorsX;
names(colorList)=setLabels;

mp = modulePreservation(multiExpr, 
 				   colorList,
 				   networkType = 'signed', 
 				   nPermutations = 0,			# no stats, just get ranks
 				   maxModuleSize = max( table(colorsVSP) ),
 				   verbose = 3,
 				   savePermutedStatistics = F
 				   );
 				   
overlap = overlapTable(colorsVSP, colorsX);
 
numMat = -log10(overlap$pTable);
numMat[numMat > 50] = 50;
 
textMat = paste(overlap$countTable, '\n', signif(overlap$pTable, 2));
dim(textMat) = dim(numMat);
 
xlabels = paste('M', sort(unique(NET.X$colors)));
ylabels = paste('M', sort(unique(NET.VSP$colors)));
xSymbols = paste(sort(unique(NET.X$colors)), ': ', table(NET.X$colors), sep = '');
ySymbols = paste(sort(unique(NET.VSP$colors)), ': ', table(NET.VSP$colors), sep = '');

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
			   
# compute ranks
tmp = cbind(mp$observed$VSP$intra$VSP_vs_X, mp$observed$VSP$inter$VSP_vs_X);
tmp = tmp[!(rownames(tmp)=='gold'), ];
tmp2 = apply(-tmp, 2, rank);
tmp3 = tmp[,c(1,2,4,5,8,9,11)]
tmp3 = tmp3[ match( substring( rownames(moduleTraitCorVSP), 3) , rownames(tmp3)), ];
 
# con - cor.kIM, cor.kME, cor.cor
medRankCon = apply( tmp2[ , c(8, 9, 11) ], 1, median);
# density - propVarExplained, meanSignAwareKME, meanSignAwareCorDat, meanAdj
medRankDen = apply(tmp2[ , c(1, 2, 4, 5) ], 1, median);
 
medRankPres = (medRankCon + medRankDen) / 2;

medRanks = data.frame( cbind(medRankPres,medRankCon,medRankDen) );
medRanks = medRanks[order(medRanks[,1]),];

#####
# plots of preservation against ME-trait cors

# preservation ranks against ME-trait cors
medRanks = medRanks[ match( substring( rownames(moduleTraitCorVSP), 3) , rownames(medRanks)), ];
par(mfrow=c(3,7));
for( column in 1:3 ) {
	for( t in 1:ncol(traitsVSP) ) {
		verboseScatterplot(medRanks[,column],
				   (moduleTraitCorVSP[, t]),
				   xlab = colnames(medRanks)[column],
				   ylab = names(traitsVSP)[t],
				   col = substring( rownames(moduleTraitCorVSP), 3),
				   pch=20,
				   cex=3,
				   abline=T, abline.color='red', abline.lty='dashed'
				   )
	}
};
rm(column, t)

# density preservation stats against ME-trait cors
par(mfrow=c(4,7),oma=c(0,0,2,0));
for( column in 1:4 ) {
	for( t in 1:7 ){
		verboseScatterplot(tmp3[,column],
				   (moduleTraitCorVSP[, t]),
				   xlab = colnames(tmp3)[column],
				   ylab = names(traitsVSP)[t],
				   col = substring( rownames(moduleTraitCorVSP), 3),
				   pch=20,
				   cex=2, cex.lab=1.2,
				   abline=T, abline.color='red', abline.lty='dashed'
				   )
	}
};title('VSP: 5368 genes', outer=T)
rm(column, t)

# meta-module average ME-trait correlations
par(mfrow=c(2,4), oma=c(0,0,7,0))
for(t in 1:7){
	verboseBarplot(moduleTraitCorVSP[,t],
				   c(rep('A',4),rep('B',4),rep('C',4),rep('D',5)),
				   ylab=colnames(moduleTraitCorVSP)[t],
				   xlab='meta-modules in VSP',
				   numberStandardErrors=2,
				   ylim=c(-0.75,0.75)
				   )
};rm(t);
title('VSP: 5368 genes\nA = purple, pink, red, turquoise\nB = midnightblue, blue, green, greenyellow\nC = cyan, grey60, salmon, yellow\nD = black, brown, lightcyan, magenta, tan',outer=T)


####################################
# build module gene lists for enrichment tests

modGenesVSP=list();
for (mod in 1:length(modNamesVSP)) {
	modGenesVSP[[mod]]=names(DATA.VSP)[colorsVSP==modNamesVSP[mod]]; 
	names(modGenesVSP)[mod]=modNamesVSP[mod]
	}; rm(mod)
	
modGenesX=list();
for (mod in 1:length(modNamesX)) {
	modGenesX[[mod]]=names(DATA.X)[colorsX==modNamesX[mod]]
	names(modGenesX)[mod]=modNamesX[mod]
	}; rm(mod)

# get genes in same mod in both regions
modGenesBoth=list();
for (mod in 1:length(modNamesVSP)){
	if(modNamesVSP[mod] %in% modNamesX){
		modGenesBoth[[mod]] = names(DATA.VSP)[colorsVSP==modNamesVSP[mod] & colorsX==modNamesVSP[mod]];
		names(modGenesBoth)[mod]=modNamesVSP[mod];
	}
}

# use all genes on chip as reference
load(file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed');
geneRefSet = rownames(datVSPcollapsedList$deviate3[[2]]);
rm(datVSPcollapsedList);

# test for enrichment of PD associated genes
source('_code/checkGeneListEnrichment.R');
PD = read.csv('PDassociated.csv',header=F)[,1];
PD = toupper(PD);
PD = gsub(' ', '', PD);

checkGeneListEnrichment(PD, modGenesVSP$blue, geneRefSet);
checkGeneListEnrichmentList(PD, modGenesVSP, geneRefSet)$pvals;

# convert gene symbols to human Entrez IDs for DAVID/Ingenuity
# load library
library(org.Hs.eg.db);
IDs=as.list(org.Hs.egSYMBOL2EG);

# convert for each VSP module and store in list
modIDsVSP=list();
for (mod in 1:length(modGenesVSP)) {
	temp=IDs[names(IDs) %in% modGenesVSP[[mod]]];
	IDvec=c();
	for (gene in 1:length(temp)) {
		IDvec=c(IDvec,temp[[gene]]);
	}
	modIDsVSP[[mod]]=IDvec;
	names(modIDsVSP)[mod]=names(modGenesVSP)[mod];
}
rm(mod,temp,IDvec,gene);
# set filenames for outputing each module to own file
filenamesVSP=c();
for (mod in 1:length(modIDsVSP)){
	filenamesVSP=c(filenamesVSP,paste('VSP.', names(modIDsVSP)[mod], 'Entrez.txt', sep=''));
}
rm(mod);
# write files
for (mod in 1:length(modIDsVSP)) {
	write.table(modIDsVSP[[mod]], file=filenamesVSP[mod], quote=F, row.names=F, col.names=F)
}
rm(mod);

# convert for each X module and store in list	
modIDsX=list();
for (mod in 1:length(modGenesX)) {
	temp=IDs[names(IDs) %in% modGenesX[[mod]]];
	IDvec=c();
	for (gene in 1:length(temp)) {
		IDvec=c(IDvec,temp[[gene]]);
	}
	modIDsX[[mod]]=IDvec;
	names(modIDsX)[mod]=names(modGenesX)[mod];
} 
rm(mod,temp,IDvec,gene);
# set filenames for outputing each module to own file
filenamesX=c();
for (mod in 1:length(modIDsX)){
	filenamesX=c(filenamesX,paste('X.', names(modIDsX)[mod], 'Entrez.txt', sep=''));
}
rm(mod);
# write files
for (mod in 1:length(modIDsX)) {
	write.table(modIDsX[[mod]], file=filenamesX[mod], quote=F, row.names=F, col.names=F)
}
rm(mod);

# convert for each shared module and store in list	
modIDsBoth=list();
for (mod in 1:length(modGenesBoth)) {
	if(!(is.null(modGenesBoth[[mod]]))) {
		temp=IDs[names(IDs) %in% modGenesBoth[[mod]]];
		IDvec=c();
		for (gene in 1:length(temp)) {
			IDvec=c(IDvec,temp[[gene]]);
		}
		modIDsBoth[[mod]]=IDvec;
		names(modIDsBoth)[mod]=names(modGenesBoth)[mod];
	}	
} 
rm(mod,temp,IDvec,gene);
# set filenames for outputing each module to own file
filenamesBoth=c();
for (mod in 1:length(modIDsBoth)){
	filenamesBoth=c(filenamesBoth,paste('Both.', names(modIDsBoth)[mod], 'Entrez.txt', sep=''));
}
rm(mod);
# write files
for (mod in 1:length(modIDsBoth)) {
	write.table(modIDsBoth[[mod]], file=filenamesBoth[mod], quote=F, row.names=F, col.names=F)
}
rm(mod);

########################################

# build info table

load(file = '_RData/annotationInfo');
load(file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed');

geneToProbe = datVSPcollapsedList$deviate3$group2row;
geneToProbe = geneToProbe[ match( names(DATA.VSP), geneToProbe[, 1] ), 2];
cloneIDs = annos$cloneID[match(geneToProbe, annos$probeID)];
sequences = annos$sequence[match(geneToProbe, annos$probeID)];
geneNames = annos$geneName[match(geneToProbe, annos$probeID)];
rm(datVSPcollapsedList, annos);

info = data.frame(probeID = geneToProbe, 
				  cloneID = cloneIDs, 
				  sequence = sequences, 
				  geneName = geneNames,
				  moduleVSP = colorsVSP,
				  moduleX = colorsX
				  );
				  
# add GS with q-vals for X and S
for(t in 1:ncol(traitsVSP)){
	temp = get(paste('GS.',names(traitsVSP)[t],'.VSP',sep=''));
	temp2 = get(paste('q.GS.',names(traitsVSP)[t],'.VSP',sep=''));
	temp3 = get(paste('GS.',names(traitsVSP)[t],'.X',sep=''));
	temp4 = get(paste('q.GS.',names(traitsVSP)[t],'.X',sep=''));
	info = cbind(info,temp,temp2,temp3,temp4);
	};
	rm(t,temp,temp2,temp3,temp4);
	
# order modules by significance for pitch
modOrder = order( -abs(cor(MEsVSP, traitsVSP$mean.pitch, use="p")) );
# add module membership info in order
for (mod in 1:ncol(geneMM.VSP)) { 
	oldNames = names(info);
	info = data.frame(info, geneMM.VSP[, modOrder[mod]], MMPvalueVSP[, modOrder[mod]]); 
	names(info)=c(oldNames, paste('MM.', modNamesVSP[modOrder[mod]], sep = ''), paste('p.MM', modNamesVSP[modOrder[mod]], sep = '')); 
}
rm(mod, modOrder, oldNames);

