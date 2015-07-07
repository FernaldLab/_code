setwd('~/Documents/_Fernald_lab');
rm(list=ls());

load(file = 'dat0_lookup0.RData');

lookup = lookup0;
dat = dat0;

rm(dat0, lookup0);

dat = dat[, colnames(dat) %in% rownames(lookup)];

# remove lookup samples not in data
lookup = lookup[rownames(lookup) %in% colnames(dat), ];

# re-order lookup rows to match data cols and check
lookup = lookup[match(colnames(dat), rownames(lookup)), ];
nrow(lookup) == ncol(dat);
sum(rownames(lookup)==colnames(dat)) == nrow(lookup);

# keep only schizophrenic and contorl samples
groupCol = match('Liu.Group', colnames(lookup));
toKeep = lookup[, groupCol] %in% c('Schizophrenia', 'Control');

dat = dat[, toKeep];
lookup = lookup[toKeep, ];

# remove transcripts with too many missing values
tooMany = floor(ncol(dat) / 3);
countNAs = c();
for (transcript in 1:nrow(dat))
{
	countNAs[transcript] = sum(is.na(dat[transcript, ]));
}
rm(transcript);

dat = dat[countNAs < tooMany, ];
################


DATA=dat[1:300,];

library(WGCNA);
options(stringsAsFactors = F);
allowWGCNAThreads();

rm(countNAs, dat, groupCol, toKeep, tooMany);

DATA = as.data.frame(t(DATA));

net0 = blockwiseModules(DATA,
						   maxBlockSize = ncol(DATA)/2,
						   power = 16,
						   networkType = 'signed',
						   deepSplit = 2,
						   minModuleSize = 10,
						   verbose = 3,
						   saveTOMs = T,
						   saveTOMFileBase = 'testingTOMblocks'
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
					guideHang = 0.05)
rm(block);

##reload saved TOMs
# load(net0$TOMFiles[[1]]);
# TOM1=TOM;
# rm(TOM);

for (f in 1:3)
{
	load(net0$TOMFiles[[f]]);
	assign(paste('TOM', f, sep=''), TOM);
	rm(TOM);
}
rm(f);
