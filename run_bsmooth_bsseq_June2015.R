###########################################################################################
###### setup
###########################################################################################
rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');
SUBJECTS = list.files()[grep('^3', list.files())];

###########################################################################################
###### load data for each subject and convert to BSseq objects
###########################################################################################

### columns where colClasses is NULL won't be read
#FILE = 'aligned.adapters.q30.m0.methratio_CpGcombined.CG.clean';
#COLCLASSES = c('character','numeric',rep('NULL',3),'numeric','numeric',rep('NULL',5));
FILE = 'aligned.adapters.q30.m0.methratio.clean'; # takes ~30min total to read in files and save RData object
COLCLASSES = c('character','numeric','character','character','NULL','numeric','numeric',rep('NULL',5));

### load data, ignore some columns
d0 = list();
for (s in SUBJECTS) {
	print(Sys.time())
	cat(s,'\n');
	d0[[s]] = read.table(paste(s, '/', FILE, sep=''), 
						 header=F, sep='\t',
						 colClasses=COLCLASSES,
						# nrows=1000
						 );
}; rm(s); gc();
print(Sys.time())
save(d0, file=paste(FILE, '_d0.RData', sep=''));
print(Sys.time())

### OR save each subject out separately
for (d in 1:length(d0)) {
	n = paste('d0_',names(d0)[d],sep='');
	cat(n,'\n');
	assign(n, d0[[d]]);
	save(list=n, file=paste(FILE, '_', n, '.RData', sep=''));
}; rm(d);



### make new list with BSseq objects, no strand info, no nucleotide context
# bsd = list();
# for (s in 1:length(d0)) {
	# cat(names(d0)[s], '\n');
	# #names(d0[[s]]) = c('chr','pos','Cov','M');
	# names(d0[[s]]) = c('chr','pos','strand','context','Cov','M');
	# if (names(d0)[s] %in% c('3157_TENNISON','3677_MONK')) {
		# group = 'ND';
	# } else if (names(d0)[s] %in% c('3165_BRISCOE','3581_LYNLEY')) {
		# group = 'D'
	# } else {
		# stop('Check samples');
	# }
	# bsd[[names(d0)[s]]] = BSseq(M=as.matrix(d0[[s]]$M, ncol=1), 
							    # Cov=as.matrix(d0[[s]]$Cov, ncol=1), 
							    # pos=d0[[s]]$pos, 
							    # chr=d0[[s]]$chr,
							    # sampleNames=names(d0)[s],
							    # pData=data.frame(group=as.character(group), row.names=names(d0)[s])
							    # );
# }; rm(s,group); gc();
# save(bsd, file=paste(FILE, '_bsd.RData', sep=''));
##########
### OR ###
##########
### make new list with BSseq objects, including strand info and nucleotide context
bsd = list();
for (s in 1:length(d0)) {
	cat(names(d0)[s], '\n');
	print(Sys.time())
	names(d0[[s]]) = c('chr','pos','strand','context','Cov','M');
	tmp = as.data.frame(cbind(d0[[s]][,1:2],d0[[s]][,2],d0[[s]][,3]));
	names(tmp)[2:4] = c('start', 'end', 'strand');
	if (names(d0)[s] %in% c('3157_TENNISON','3677_MONK')) {
		group = 'ND';
	} else if (names(d0)[s] %in% c('3165_BRISCOE','3581_LYNLEY')) {
		group = 'D'
	} else {
		stop('Check samples');
	}
	bsd[[names(d0)[s]]] = BSseq(M=as.matrix(d0[[s]]$M, ncol=1),
								Cov=as.matrix(d0[[s]]$Cov, ncol=1),
								gr=data.frame2GRanges(tmp),
								sampleNames=names(d0)[s],
								pData=data.frame(group=as.character(group), row.names=names(d0)[s]),
								rmZeroCov=TRUE
								);
	mcols(bsd[[names(d0)[s]]]) = d0[[s]]$context;
}; rm(tmp,s,group); gc();
print(Sys.time())
save(bsd, file=paste(FILE, '_bsd.RData', sep=''));
print(Sys.time())
###########################################################################################
###### combine subjects and filter scaffolds with too few methylation loci
###########################################################################################
print(Sys.time())
### combine into single BSseq object
bsdall = combine(bsd[[1]], bsd[[2]]);gc();print(Sys.time())
bsdall = combine(bsdall, bsd[[3]]);gc();print(Sys.time())
bsdall = combine(bsdall, bsd[[4]]);gc();print(Sys.time())
save(bsdall, file=paste(FILE, '_bsdCombined.RData', sep=''));

### set number of loci required to keep scaffold
### should be value of ns parameter to BSmooth
NUMLOCI = 70;
### get number of methylation loci per scaffold
scaffoldCounts = data.frame(scaffold=attributes(seqnames(bsdall))$values, 
							num=attributes(seqnames(bsdall))$lengths
							);
							
### get scaffolds to keep
scaffoldsCountsKeep = scaffoldCounts[scaffoldCounts$num > NUMLOCI,];
NAME = paste('bsdalln', NUMLOCI, sep='');
assign(NAME, bsdall[which(seqnames(bsdall) %in% scaffoldsCountsKeep$scaffold)]);
#scaffoldsCountsKeepNums = na.omit(as.numeric(unlist(strsplit(as.character(scaffoldsCountsKeep$scaffold), '_'))));

# can't get get(NAME) to work in save() so be careful of object names
save(bsdalln70, file=paste(FILE, '_bsdCombined_n', NUMLOCI, '.RData', sep=''));

###########################################################################################
###### check how much of genome is covered by remaining scaffolds
###########################################################################################

### read in scaffold lengths
sl = read.table('/Volumes/fishstudies/_Burtoni_genome_files/scaffold_lengths');

### filter to scaffolds from above
slFilt = sl[match(scaffoldsCountsKeep$scaffold, sl$V1), ];
slFilt = slFilt[order(as.numeric(rownames(slFilt))), ];

### compute percentage of genome covered by remaining scaffolds
sum(slFilt[,2]) / sum(sl[,2]);

### if kept scaffolds with >70 loci, this should be 0.9693055

###########################################################################################
###### run BSmooth, try different maxGap settings
###########################################################################################
Sys.time()
bsdalln70.fit.ns70h1000mg10e8 = BSmooth(bsdalln70, ns=70, h=1000, maxGap=10^8, mc.cores=4, parallelBy='sample');
Sys.time()
gc();
save(bsdalln70.fit.ns70h1000mg10e8, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8.RData'
	 );
Sys.time()	 
bsdalln70.fit.ns70h1000mg10e4 = BSmooth(bsdalln70, ns=70, h=1000, maxGap=10^4, mc.cores=4, parallelBy='sample');
Sys.time()					    				
gc();
save(bsdalln70.fit.ns70h1000mg10e4, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e4.RData'
	 );
Sys.time()
bsdalln70.fit.ns70h1000mg10e3 = BSmooth(bsdalln70, ns=70, h=1000, maxGap=10^3, mc.cores=4, parallelBy='sample');
Sys.time()					    				
gc();
save(bsdalln70.fit.ns70h1000mg10e3, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e3.RData'
	 );

### took ~10-30min for each, shorter maxGap meant longer runtime


######
Sys.time()
bsdalln70.fit.ns20h200mg500 = BSmooth(bsdalln70, ns=20, h=200, maxGap=500, mc.cores=4, parallelBy='sample');
Sys.time()
save(bsdalln70.fit.ns20h200mg500, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns20h200mg500.RData'
	 );

###########################################################################################
###### filter based on coverage and remove scaffolds again if needed
###########################################################################################

FIT = bsdalln70.fit.ns20h200mg500;

FIT.Cov = getCoverage(FIT);
# keep only sites that had 5x or higher coverage in all samples
reqCov = 5;
reqSamples = 4;
keep = which(apply(FIT.Cov, 1, function(f) sum(f>=reqCov)==reqSamples));
FIT.FiltX = FIT[keep, ];

NUMLOCI = 20;
### get number of methylation loci per scaffold
scaffoldCountsPostFit = data.frame(scaffold=attributes(seqnames(FIT.FiltX))$values, 
								   num=attributes(seqnames(FIT.FiltX))$lengths
							       );
### get scaffolds to keep
scaffoldsCountsPostFitKeep = scaffoldCountsPostFit[scaffoldCountsPostFit$num > NUMLOCI,];
FIT.FiltXn = FIT.FiltX[which(seqnames(FIT.FiltX) %in% scaffoldsCountsPostFitKeep$scaffold)];

### compute percentage of genome covered by remaining scaffolds
sl = read.table('/Volumes/fishstudies/_Burtoni_genome_files/scaffold_lengths');
slFilt = sl[match(scaffoldsCountsPostFitKeep$scaffold, sl$V1), ];
slFilt = slFilt[order(as.numeric(rownames(slFilt))), ];
sum(slFilt[,2]) / sum(sl[,2]);

##########
bsdalln70.fit.ns20h200mg500.5x.n20 = FIT.FiltXn;
save(bsdalln70.fit.ns20h200mg500.5x.n20, 
	 file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e3_5x_n70.RData'
	 );

### for 5x coverage and ns=70, was left with 2012 scaffolds covering 0.954575 of genome

### go forward with bsdalln70.fit.ns70h1000mg10e8.5x.n70


###########################################################################################
###### compute t-statistics
###########################################################################################
FIT.FiltXn = bsdalln70.fit.ns70h1000mg10e8.5x.n70;

groupD = c('3165_BRISCOE','3581_LYNLEY');
groupND = c('3157_TENNISON','3677_MONK');
Sys.time()
FIT.FiltXn.T = BSmooth.tstat(FIT.FiltXn, 
							 group1=groupND,
						     group2=groupD,
			       	 		 estimate.var='group2',
							 local.correct=T,
							 mc.cores=4,
							 verbose=T
							 );gc();
Sys.time()
### took ~15min

bsdalln70.fit.ns70h1000mg10e8.5x.n70.tstat_Dvar = FIT.FiltXn.T;
save(bsdalln70.fit.ns70h1000mg10e8.5x.n70.tstat_Dvar, 
     file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_Dvar.RData'
     );
     
###########################################################################################
###### find DMRs
###########################################################################################

FIT.FiltXn.T = bsdalln70.fit.ns70h1000mg10e8.5x.n70.tstat_Dvar;

plot(FIT.FiltXn.T);
apply(attributes(FIT.FiltXn.T)$stats, 2, summary);

for (pct in c(.01, .001, .0001)) {
	cat(pct, '\n');
	int = c(pct/2, 1-(pct/2));
	tvals = quantile(attributes(FIT.FiltXn.T)$stats[,5], int, na.rm=T);
	assign(paste('dmrs0', pct, sep=''), dmrFinder(FIT.FiltXn.T, cutoff=tvals));gc();
}; rm(pct,int,tvals);

dmrs01e04 = get('dmrs01e-04');
save(dmrs00.001, dmrs00.01, dmrs01e04, 
     file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_Dvar_dmrs0.RData'
     );
     
###########################################################################################
###### 
###########################################################################################
     
nThresh = 8;
meanDiffThresh = 0.1;
dmrs.01 = subset(dmrs00.01, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.001 = subset(dmrs00.001, n>=nThresh & abs(meanDiff)>=meanDiffThresh);
dmrs.1e04 = subset(dmrs01e04, n>=nThresh & abs(meanDiff)>=meanDiffThresh);

par(mfrow=c(3,3)); 
hist(abs(dmrs.01$meanDiff), breaks=nrow(dmrs.01)); 
hist(abs(dmrs.001$meanDiff), breaks=nrow(dmrs.001)); 
hist(abs(dmrs.1e04$meanDiff), breaks=nrow(dmrs.1e04));

hist(dmrs.01$width, breaks=nrow(dmrs.01)); 
hist(dmrs.001$width, breaks=nrow(dmrs.001)); 
hist(dmrs.1e04$width, breaks=nrow(dmrs.1e04));

hist(dmrs.01$invdensity, breaks=nrow(dmrs.01)); 
hist(dmrs.001$invdensity, breaks=nrow(dmrs.001)); 
hist(dmrs.1e04$invdensity, breaks=nrow(dmrs.1e04));


FIT.FiltXn = bsdalln70.fit.ns70h1000mg10e8.5x.n70;
FILE = 'aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_NDvar_dmrs.01_n10_md.2.pdf'

pData = pData(FIT.FiltXn);
pData$col = c('blue','red','red','blue');
pData(FIT.FiltXn) = pData;				 

pdf(file=FILE, width=10, height=5);
plotManyRegions(FIT.FiltXn, dmrs.001, extend=5000, addRegions=dmrs.001);
dev.off();

###########################################################################################
###### compare dmrs from different settings
###########################################################################################
rm(list=ls());

load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_dmrs0.RData');

dmrs00.01same = dmrs00.01;
dmrs00.001same = dmrs00.001;
dmrs01e04same = dmrs01e04;
rm(dmrs00.01, dmrs00.001, dmrs01e04);

load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_NDvar_dmrs0.RData');

dmrs00.01NDvar = dmrs00.01;
dmrs00.001NDvar = dmrs00.001;
dmrs01e04NDvar = dmrs01e04;
rm(dmrs00.01, dmrs00.001, dmrs01e04);

load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n70Smoothed.ns70h1000mg10e8_5x_n70_tstat_Dvar_dmrs0.RData');

dmrs00.01Dvar = dmrs00.01;
dmrs00.001Dvar = dmrs00.001;
dmrs01e04Dvar = dmrs01e04;
rm(dmrs00.01, dmrs00.001, dmrs01e04);

dmrs0 = list(S01=dmrs00.01same, S001=dmrs00.001same, S0001=dmrs01e04same,
			 ND01=dmrs00.01NDvar, ND001=dmrs00.001NDvar, ND0001=dmrs01e04NDvar, 
			 D01=dmrs00.01Dvar, D001=dmrs00.001Dvar, D0001=dmrs01e04Dvar
			 );
for (i in 1:length(dmrs0)) {
	rownames(dmrs0[[i]]) = paste(dmrs0[[i]]$chr, ':', dmrs0[[i]]$start, '-', dmrs0[[i]]$end, sep='');
}; rm(i);

dmrs01 = dmrs0[c(1,4,7)];

par(mfrow=c(1,3));
COL = which(names(dmrs01$S01)=='n');
WGCNA::verboseBoxplot(c(dmrs01$S01[,COL], 
						dmrs01$ND01[,COL], 
						dmrs01$D01[,COL]), 
				      c(rep('S',length(dmrs01$S01[,COL])),
				      	rep('ND',length(dmrs01$ND01[,COL])),
				      	rep('D',length(dmrs01$D01[,COL]))),
				     # ylim=c(0,15),
				      ylab=names(dmrs01$S01)[COL],
				      xlab='',
				      frame.plot=F,
				      col='grey'
				      );
COL = which(names(dmrs01$S01)=='width');
WGCNA::verboseBoxplot(c(dmrs01$S01[,COL], 
						dmrs01$ND01[,COL], 
						dmrs01$D01[,COL]), 
				      c(rep('S',length(dmrs01$S01[,COL])),
				      	rep('ND',length(dmrs01$ND01[,COL])),
				      	rep('D',length(dmrs01$D01[,COL]))),
				      #ylim=c(0,600),
				      ylab=names(dmrs01$S01)[COL],
				      xlab='',
				      frame.plot=F,
				      col='grey'
				      );
COL = which(names(dmrs01$S01)=='invdensity');
WGCNA::verboseBoxplot(c(dmrs01$S01[,COL], 
						dmrs01$ND01[,COL], 
						dmrs01$D01[,COL]), 
				      c(rep('S',length(dmrs01$S01[,COL])),
				      	rep('ND',length(dmrs01$ND01[,COL])),
				      	rep('D',length(dmrs01$D01[,COL]))),
				     # ylim=c(0,100),
				      ylab=names(dmrs01$S01)[COL],
				      xlab='',
				      frame.plot=F,
				      col='grey'
				      );



nThresh = 10;
meanDiffThresh = 0.2;

for (i in 1:length(dmrs01)) {
	dmrs01[[i]] = subset(dmrs01[[i]], n>=nThresh & abs(meanDiff)>=meanDiffThresh);
}; rm(i);



