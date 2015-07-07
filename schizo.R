setwd('~/Documents/_Fernald_lab');
rm(list=ls());

# load data
# dat0 = as.matrix(read.csv('matrix.tab.csv', header=F, sep='\t'));
# rows = read.table('rowlabels.txt')[,1];
# cols = read.table('collabels.txt')[,1];
# rownames(dat0) = rows; colnames(dat0) = cols;
# rm(rows, cols);
# dat0[is.nan(dat0)] = NA;

# # load sample info
# lookup0 = as.matrix(read.csv('master_lookup.csv', header=T));
# lookup0 = gsub(' ', '', lookup0);
# rownames(lookup0) = lookup0[, 2];
# lookup0 = lookup0[, -2];
# lookup0[lookup0=='NaN'] = NA;

# save.image(file = 'dat0_lookup0.RData');

### filter data ###
load(file = 'dat0_lookup0.RData');

# remove data samples not in lookup
dat0 = dat0[, colnames(dat0) %in% rownames(lookup0)];
# remove lookup samples not in data
lookup0 = lookup0[rownames(lookup0) %in% colnames(dat0), ];

# re-order lookup rows to match data cols and check
lookup0 = lookup0[match(colnames(dat0), rownames(lookup0)), ];
nrow(lookup0) == ncol(dat0);
sum(rownames(lookup0)==colnames(dat0)) == nrow(lookup0);

# keep only schizophrenic and control samples
groupCol = match('My.Annotation', colnames(lookup0));
#toKeep = lookup0[, groupCol] %in% c('SZ', 'CTRL');
toKeepSZ = lookup0[, groupCol] == 'SZ';
toKeepCTRL = lookup0[, groupCol] == 'CTRL';
toKeepBPD = lookup0[, groupCol] == 'BPD';

#dat = dat0[, toKeep];
dat.SZ = dat0[, toKeepSZ];
dat.CTRL = dat0[, toKeepCTRL];
dat.BPD = dat0[, toKeepBPD];
#lookup = lookup0[toKeep, ];
lookup.SZ = lookup0[toKeepSZ, ];
lookup.CTRL = lookup0[toKeepCTRL, ];
lookup.BPD = lookup0[toKeepBPD, ];

DATAlist = list(dat0=dat0, dat.SZ=dat.SZ, dat.BPD=dat.BPD, dat.CTRL=dat.CTRL);
LOOKUPlist = list(lookup0=lookup0, lookup.SZ=lookup.SZ, lookup.BPD=lookup.BPD, lookup.CTRL=lookup.CTRL);

save(DATAlist, LOOKUPlist, file='SZ_DATA_LOOKUPs.RData');
#load('SZ_DATA_LOOKUPs.RData');
#################

# remove transcripts with too many missing values

DATAlist_filt = list();
for (set in 1:length(DATAlist))
{
	DATA = DATAlist[[set]];
	tooMany = floor(ncol(DATA) / 3);
	countNAs = c();
	
	for (transcript in 1:nrow(DATA))
	{
		countNAs[transcript] = sum(is.na(DATA[transcript, ]));
	}
	rm(transcript);
	
	DATAlist_filt[[set]] = DATA[countNAs < tooMany, ];
}
rm(set, tooMany, DATA, countNAs);
names(DATAlist_filt) = names(DATAlist);

# tooMany = floor(ncol(dat) / 3);
# countNAs = c();
# for (transcript in 1:nrow(dat))
# {
	# countNAs[transcript] = sum(is.na(dat[transcript, ]));
# }
# rm(transcript);

# dat = dat[countNAs < tooMany, ];

#########################################
### build preliminary net for testing ###
library(WGCNA);
options(stringsAsFactors = F);
allowWGCNAThreads();
collectGarbage();

DATA = as.data.frame(t(DATAlist_filt$dat.SZ));

system.time
({
	pickSoft = pickSoftThreshold(DATA,
							networkType = 'signed',
							blockSize = 6500,
							verbose = 2
							)	
})
collectGarbage();
							
k16 = softConnectivity(DATA, type='signed', power=16, blockSize=6500, verbose=3);
k18 = softConnectivity(DATA, type='signed', power=18, blockSize=6500, verbose=3);


net0 = blockwiseModules(dat,
 						   maxBlockSize = 6500,
 						   power = 14,
 						   networkType = 'signed',
 						   deepSplit = 2,
 						   minModuleSize = 10,
 						   verbose = 3,
 						   saveTOMs = T,
 						   saveTOMFileBase = 'SZ_CRTL_postkfilt_run0'
 						   );
 						   
block = 1;
plotDendroAndColors(net0$dendrograms[[block]],
					net0$colors[ net0$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = net0$colors[ net0$blockGenes[[block]] ],
					main = paste('block', block),
					dendroLabels = F,
					#hang = 0.03,
					addGuide = T,
					guideHang = 0.05);
rm(block); collectGarbage();

# compute average TO within each module
source('~/Documents/_analysis_compare/_code/getModDensitiesFromTOM.R');
TOM0 = TOMsimilarityFromExpr(dat,
							   networkType = 'signed',
						       power = 16,
				   		       verbose = 3
							   );
collectGarbage();