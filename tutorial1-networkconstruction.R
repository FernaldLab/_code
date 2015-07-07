################################
##### network construction #####
################################

# start up
rm(list = ls());
setwd('~/Documents/_analysis_compare');
library('WGCNA');
options(stringsAsFactors=F);
load(file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed');
load(file = '_Rdata/VSPprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_k14');

k14n = k14 / max(k14);
par(mfrow = c(1, 2));
hist(k14n);
hist(k14);

cutoff = 0.05;		# probably change for each dataset

DATA = DATAcollapsedList[[3]][[1]]; #DATA = datVSPcollapsedList[[3]][[1]];
rm(DATAcollapsedList);				#rm(datVSPcollapsedList);
DATA = as.data.frame( t(DATA) );
DATA = DATA[ , match( names(k14n)[k14n > cutoff] , names(DATA) )];

net0 = blockwiseModules(DATA,
						   maxBlockSize = 6000,
						   power = 14,
						   networkType = 'signed',
						   deepSplit = 2,
						   minModuleSize = 10,
						   verbose = 3,
						   saveTOMs = T,
						   saveTOMFileBase = 'VSPpostkfilt_run0'
						   );
collectGarbage();
#save(net0, file = '_Rdata/Xprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_NET0');
						   
# plot dendrograms
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

# compute average TO within each module
# #source('_code/getModuleDensities.R');
# #modDensities = getModuleDensities(datVSP, netVSP0$colors);
source('_code/getModDensitiesFromTOM.R');

TOM0 = TOMsimilarityFromExpr(DATA,
							   networkType = 'signed',
						       power = 14,
				   		       verbose = 3
							   );
collectGarbage();
#save(TOM0, file = '_Rdata/Xprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_TOM0');

modDensities = getModDensitiesFromTOM(TOM0, net0$colors);
nPerm = 5000;
# system.time( {	# skip grey for now to save time
	# permTest = modDensityPerm( TOM0, net0$colors[net0$colors!='grey'], modDensities[names(modDensities)!='grey'], nPerm = nPerm );
	# } )
# save(permTest, file = '_Rdata/Xprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_PERM5000');

system.time( {		# probably take 40-60min
	permTest = modDensityPerm( TOM0, net0$colors, modDensities, nPerm = nPerm );
	} )
#save(permTest, file = '_Rdata/Xprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_PERM5000');

weakTOM = names( table(net0$colors)[ permTest$pvals > 0.01 ] ); # includes grey
bgGenes = net0$colors %in% weakTOM;

DATAfilt = DATA[ , !bgGenes];

netfilt = blockwiseModules(DATAfilt,
						   maxBlockSize = 6000,
						   power = 14,
						   networkType = 'signed',
						   deepSplit = 2,
						   minModuleSize = 10,
						   verbose = 3);
collectGarbage();
#save(net0, file = '_Rdata/VSPprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_NETfiltperm5000p0.01mod10');