rm(list=ls());
setwd('~/Documents/_Fernald_lab');
library(WGCNA); allowWGCNAThreads();
options('stringsAsFactors'=F);
source('_code/preProcATH-for_web_noVSN.R');
source('_code/exploreNetwork.R');

dat0 = read.csv('GSE10624_series_matrix.csv', row.names=1, sep='\t');

files = paste('GSE10624_RAW/', list.files('GSE10624_RAW', '*clean'), sep='');
raw0 = list();
for (f in 1:length(files)) {
	this = gsub('GSE10624_RAW/', '', files[f]); print(this);
	raw0[[f]] = read.table(files[f], sep='\t', header=T, fill=T);
	names(raw0)[f] = this;
}; rm(f, this);
raw = list();
for (s in 1:length(raw0)) {
	tmp = raw0[[s]];
	tmp = tmp[tmp$ID %in% rownames(dat0), ];
	raw[[s]] = tmp[order(tmp$ID), ];
	names(raw)[s] = names(raw0)[s];
}; rm(s, tmp);

# both should return TRUE
all(unlist(lapply(raw,function(df){all(df$ID==rownames(dat0))})));
exn.checkNamesOfModListElements(raw);

for (s in 1:length(raw)) {
	raw[[s]] = raw[[s]][, -c(1:4, 6:8, 33:44, 54:56)]; 
}; rm(s);
names(raw) = gsub('.gpr.clean', '', names(raw));

raw_el = raw[[1]];
lim = 10;
for (row in 1:lim) {
	cat('B635: ', raw_el$B635[row], ', B635SD: ', raw_el$B635.SD[row], ', B635+2SD: ' , (2*raw_el$B635.SD[row])+raw_el$B635[row], ', F635med: ', raw_el$F635.Median[row], ', F635<B635+2SD: ', raw_el$F635.Median[row]<(2*raw_el$B635.SD[row])+raw_el$B635[row], '\n', sep='');
}; rm(raw_el, row, lim);


dat635 = as.data.frame(matrix(nrow=nrow(raw[[1]]), ncol=length(raw)));
dat532 = as.data.frame(matrix(nrow=nrow(raw[[1]]), ncol=length(raw)));
rownames(dat635) = raw[[1]]$ID;
rownames(dat532) = raw[[1]]$ID;
for (s in 1:length(raw)) {
	
	f635 = raw[[s]]$F635.Mean;
	b635 = raw[[s]]$B635.Mean;
	f635bgsub = f635 - b635;
	b635.2sd = 2*raw[[s]]$B635.SD;
	thresh635 = b635 + b635.2sd;
	check635 = as.logical(f635 < thresh635);
	f635bgsub[check635] = NA;
	dat635[, s] = f635bgsub;
	names(dat635)[s] = paste(names(raw)[s], '_635', sep='');
	
	f532 = raw[[s]]$F532.Mean;
	b532 = raw[[s]]$B532.Mean;
	f532bgsub = f532 - b532;
	b532.2sd = 2*raw[[s]]$B532.SD;
	thresh532 = b532 + b532.2sd;
	check532 = as.logical(f532 < thresh532);
	f532bgsub[check532] = NA;
	dat532[, s] = f532bgsub;
	names(dat532)[s] = paste(names(raw)[s], '_532', sep='');
	
}; rm(s, f635, b635, f635bgsub, b635.2sd, thresh635, check635, f532, b532, f532bgsub, b532.2sd, thresh532, check532);

dat = cbind(dat635, dat532);
dat = dat[, order(names(dat))];
save(dat, file='dat_FandC_2color.RData');

##################

design0 = read.csv('GSE10624_series_matrix_metadataClean.csv', sep='\t', row.names=1);
design = as.data.frame(t(design0));
names(design) = gsub('!', '', names(design));

traits0 = as.data.frame(matrix(nrow=ncol(dat), ncol=3));
rownames(traits0) = colnames(dat);
for (s in 1:nrow(traits0)) {
	tmp = strsplit(rownames(traits0)[s], '_')[[1]];
	s_id = tmp[1];
	cy = tmp[2];
	d_row = match(s_id, rownames(design));
	if (cy==635) {
		traits0[s, 1] = design$Sample_source_name_ch1[d_row];
		traits0[s, 2] = design$Sample_characteristics_ch1[d_row];
		traits0[s, 3] = 1;
	} else if (cy==532) {
		traits0[s, 1] = design$Sample_source_name_ch2[d_row];
		traits0[s, 2] = design$Sample_characteristics_ch2[d_row];
		traits0[s, 3] = 2;
	} else {
		stop('');	
	}
}; rm(s, tmp, s_id, cy, d_row);

traits = traits0;
colnames(traits) = c('sex', 'status', 'ch');
for (row in 1:nrow(traits)) {
	if (traits[row, 2]=='brooding female') {
		traits[row, 1:2] = c('female', NA);
	} else if (traits[row, 2]=='T-male') {
		traits[row, 1:2] = c('male', 'T');
	} else if (traits[row, 2]=='NT-male') {
		traits[row, 1:2] = c('male', 'NT');
	} else {
		stop('');
	}
}; rm(row);
traits$sex = as.numeric(as.factor(traits$sex));
traits$status = as.numeric(as.factor(traits$status));

traitCor = exn.computeAndPlotTraitCors(traits, net$MEs);

##################
dat = as.data.frame(t(dat));
datSplit = split(dat, as.factor(traits0[,1]));

lapply(datSplit, function(df){mean(cor(t(df), use='p')[upper.tri(cor(t(df), use='p'))])});

datAvg0 = as.data.frame(matrix(nrow=length(datSplit), ncol=ncol(datSplit[[1]])))
dimnames(datAvg0) = list(names(datSplit), names(datSplit[[1]]));
for (f in 1:length(datSplit)) {
	this = datSplit[[f]];
	avg = apply(this, 2, mean, na.rm=T);
	this = as.data.frame(rbind(this, avg));
	datAvg0[match(names(datSplit)[f], rownames(datAvg0)), ] = avg;
	datSplit[[f]] = this;
}; rm(f, this, avg);

outP = removeOutlierProbesIterate(t(datAvg));
outNA = removeTooManyNAs(outP$dataClean, probe_thresh=floor(ncol(outP$dataClean)/3), sample_thresh=floor(nrow(outP$dataClean)/3));
outS = outlierSamplesIterate(outNA$dataClean);
datAvg = as.data.frame(normalize.quantiles(as.matrix(outS$dataClean)));
names(datAvg) = names(outS$dataClean); rownames(datAvg) = rownames(outS$dataClean);
datAvg = as.data.frame(t(datAvg));
#datAvg[which(datAvg=='NaN', arr.ind=T)] = NA;
sft = pickSoftThreshold(datAvg, networkType='signed', verbose=3, blockSize=5000);


#####################
datM0 = dat[, traits$sex==2];
datM_T0 = dat[, traits$sex==2 & traits$status==2];
datM_NT0 = dat[, traits$sex==2 & traits$status==1];
datF0 = dat[, traits$sex==1];