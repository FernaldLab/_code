

#########################################################
###### load data and filter down to common hits
##############################################################
setwd('/Volumes/fishstudies/_methylation');
rm(list=ls());

# load differential methylation files
files = list.files()[grep('significant$', list.files())];
# files = list.files()[grep('untested',list.files())];


dat0 = list();
for (f in files) {
	dat0[[f]] = read.table(f, colClasses='character');
}; rm(f);

##############################################################
# working with "untested" files
# 3157 (nd1), 3677 (nd2), 3165 (d1), 3581 (d2)

dat = dat0;
for (d in 1:length(dat0)) {
	rownames(dat[[d]]) = apply(dat[[d]][,1:2], 1, function(f) paste(f, collapse=':'));
}; rm(d);

alldat0 = list();
# 3157 - nd1
tmp1 = dat[[1]][,c(3,4,6,7)];
tmp2 = dat[[2]][,c(3,4,6,7)];
alldat0[['nd1']] = as.data.frame(rbind(tmp1, tmp2[!(rownames(tmp2) %in% rownames(tmp1)), ]));
names(alldat0[['nd1']]) = c('strand','nucs','C','U');
x = unlist(strsplit(rownames(alldat0[['nd1']]),':'));
alldat0[['nd1']] = alldat0[['nd1']][order(x[seq(1,length(x),2)], as.numeric(x[seq(2,length(x),2)])),];

# 3677 - nd2
tmp1 = dat[[3]][,c(3,4,6,7)];
tmp2 = dat[[4]][,c(3,4,6,7)];
alldat0[['nd2']] = as.data.frame(rbind(tmp1, tmp2[!(rownames(tmp2) %in% rownames(tmp1)), ]));
names(alldat0[['nd2']]) = c('strand','nucs','C','U');
x = unlist(strsplit(rownames(alldat0[['nd2']]),':'));
alldat0[['nd2']] = alldat0[['nd2']][order(x[seq(1,length(x),2)], as.numeric(x[seq(2,length(x),2)])),];

# 3165 - d1
tmp1 = dat[[1]][,c(3,4,9,10)];
tmp2 = dat[[3]][,c(3,4,9,10)];
alldat0[['d1']] = as.data.frame(rbind(tmp1, tmp2[!(rownames(tmp2) %in% rownames(tmp1)), ]));
names(alldat0[['d1']]) = c('strand','nucs','C','U');
x = unlist(strsplit(rownames(alldat0[['d1']]),':'));
alldat0[['d1']] = alldat0[['d1']][order(x[seq(1,length(x),2)], as.numeric(x[seq(2,length(x),2)])),];

# 3581 - d2
tmp1 = dat[[2]][,c(3,4,9,10)];
tmp2 = dat[[4]][,c(3,4,9,10)];
alldat0[['d2']] = as.data.frame(rbind(tmp1, tmp2[!(rownames(tmp2) %in% rownames(tmp1)), ]));
names(alldat0[['d2']]) = c('strand','nucs','C','U');
x = unlist(strsplit(rownames(alldat0[['d2']]),':'));
alldat0[['d2']] = alldat0[['d2']][order(x[seq(1,length(x),2)], as.numeric(x[seq(2,length(x),2)])),];
rm(tmp1,tmp2,x);

alldat = alldat0;
all_count_thresh = 10;
alldat = lapply(alldat, function(f) f[.getRowsAboveCountThreshSingle(dat=f,thresh=all_count_thresh,col=c(3,4)), ]);
for (i in 1:length(alldat)) {
	print(names(alldat)[i]);
	alldat[[i]] = as.data.frame(cbind(alldat[[i]], val=as.numeric(alldat[[i]]$C)/as.numeric(alldat[[i]]$U)));
	# clear out sites with invalid nucleotides since too close to start of scaffold
	alldat[[i]] = alldat[[i]][alldat[[i]]$nucs!='0.000', ];
}; rm(i);


pos_counts = table(c(rownames(alldat$nd1), rownames(alldat$nd2), rownames(alldat$d1), rownames(alldat$d2)));
pos_comm = names(pos_counts)[pos_counts==4];

alldat_comm = alldat;
for (i in 1:length(alldat_comm)) {
	alldat_comm[[i]] = alldat_comm[[i]][rownames(alldat_comm[[i]]) %in% pos_comm, ];
}; rm(i);

# combine into data frame
adc = as.data.frame(cbind(alldat_comm$nd1, alldat_comm$nd2[,3:5], alldat_comm$d1[,3:5], alldat_comm$d2[,3:5]));
names(adc)[3:14] = c('nd1C','nd1U','nd1val','nd2C','nd2U','nd2val','d1C','d1U','d1val','d2C','d2U','d2val');

adc = as.data.frame(cbind(adc[,1:2], nucsRC=.addReverseComplementColumn(adc,2,1)$nucsRC, adc[,3:14]));

adc = .addFisherPvals(adc,c(4,5,7,8));			# nd1 vs nd2
adc = .addFisherPvals(adc,c(10,11,13,14));		# d1 vs d2
adc = .addFisherPvals(adc,c(4,5,10,11));		# nd1 vs d1
adc = .addFisherPvals(adc,c(4,5,13,14));		# nd1 vs d2
adc = .addFisherPvals(adc,c(7,8,10,11));		# nd2 vs d1
adc = .addFisherPvals(adc,c(7,8,13,14));		# nd2 vs d2
names(adc)[16:21] = c('nd1nd2','d1d2','nd1d1','nd1d2','nd2d1','nd2d2');

adc = as.data.frame(cbind(adc, ndval=apply(adc[,c(6,9)], 1, mean), dval=apply(adc[,c(12,15)], 1, mean)));
sigrows0 = apply(adc[,18:21], 1, function(f) all(f < .05));

# remove hits where nd1-nd2 or d1-d2 is significantly different
siginternal = apply(cbind(adc$nd1nd2<.05, adc$d1d2<.05)[,1:2], 1, function(f) any(f));
sigrows = sigrows & !siginternal;

# check that all hits are consistent 
tmp = cbind(adc$nd1val>adc$d1val, adc$nd1val>adc$d2val, adc$nd2val>adc$d1val, adc$nd2val>adc$d2val);
consistent = apply(tmp, 1, function(f) length(unique(f))==1);
sigrows = sigrows & consistent;


# add columns with methylation dinucleotide to use as factor
adc = .addMerCol(adc, strandCol=1, nucsCol=2, rcCol=which(names(adc)=='nucsRC'), l=2, colname='di');
adc = as.data.frame(cbind(adc, diH=adc$di));
adc$diH = as.character(adc$diH);

# add columns with methylation trinucleotide to use as factor
adc = .addMerCol(adc, strandCol=1, nucsCol=2, rcCol=which(names(adc)=='nucsRC'), l=3, colname='tri');
adc = as.data.frame(cbind(adc, triH=adc$tri));
adc$triH = as.character(adc$triH);

# edit diH and triH columns
for (r in 1:nrow(adc)) {
	adc[r, which(names(adc)=='diH')] = .flipNonGtoH(adc[r, which(names(adc)=='diH')], 2);
	adc[r, which(names(adc)=='triH')] = .flipNonGtoH(adc[r, which(names(adc)=='triH')], 2);
	#adc[r, which(names(adc)=='triHH')] = .flipNonGtoH(adc[r, which(names(adc)=='triHH')], 3);
}; rm(r);

adc = as.data.frame(cbind(adc, triHH=adc$triH));
adc$triHH = as.character(adc$triHH);
for (r in 1:nrow(adc)) {
	adc[r, which(names(adc)=='triHH')] = .flipNonGtoH(adc[r, which(names(adc)=='triHH')], 3);
}; rm(r);

# add column to denote DM positions
adc = as.data.frame(cbind(adc, sig=sigrows));

# add column to denote whether D or ND had higher avg methylation
adc = cbind(adc, higher=adc$ndval>adc$dval);
adc$higher[adc$higher] = 'ND';
adc$higher[adc$higher=='FALSE'] = 'D';
adc$higher = as.character(adc$higher);

# add column to denote whether all pairwise comparisons go in same direction
tmp = cbind(adc$nd1val>adc$d1val, adc$nd1val>adc$d2val, adc$nd2val>adc$d1val, adc$nd2val>adc$d2val, sigrows);
adc = as.data.frame(cbind(adc, cons=rep('D', nrow(adc))));


# check consistency



# # adc = cbind(adc, higher=rep('D', nrow(adc)));
# # adc$higher = as.character(adc$higher);
# # tmp = cbind(adc$nd1val>adc$d1val, adc$nd1val>adc$d2val, adc$nd2val>adc$d1val, adc$nd2val>adc$d2val);
# # # check consistency
# # if (all(apply(tmp,1,function(f) length(unique(f))) == 1)) {
	# # adc$higher[tmp[,1]] = 'ND';
# # }; rm(tmp);

##############################################################

# filter based on coverage, then get positions in at least 'thresh' files
dat = dat0;
count_thresh = 10;
file_thresh = 4;

dat = lapply(dat, function(f) f[.getRowsAboveCountThresh(dat=f,thresh=count_thresh,col=c(6,7,9,10)), ]);
pos = lapply(dat, function(f) f = paste(f[,1], f[,2], sep=':'));
hit_counts = table(unlist(pos));
common = names(hit_counts)[hit_counts >= file_thresh];
for (d in 1:length(dat)) {
	dat[[d]] = dat[[d]][paste(dat[[d]][,1], dat[[d]][,2], sep=':') %in% common, ];
}; rm(d);

# add fisher test pvals
dat = lapply(dat, function(f) .addFisherPvals(f, cols=c(6,7,9,10)));

###############################
### if require hit to be in all files
# combine into one matrix
datCommon = as.data.frame(dat[[1]][,1:4]);
names(datCommon) = c('scaffold','pos','strand','nucs');

# add reads in order 3157 (nd1), 3677 (nd2), 3165 (d1), 3581 (d2)
datCommon = cbind(datCommon, dat[[1]][,6:7], 
							 as.numeric(dat[[1]][,6]) / as.numeric(dat[[1]][,7]),
							 dat[[3]][,6:7],
							 as.numeric(dat[[3]][,6]) / as.numeric(dat[[3]][,7]),
							 dat[[1]][,9:10],
							 as.numeric(dat[[1]][,9]) / as.numeric(dat[[1]][,10]),
							 dat[[2]][,9:10],
							 as.numeric(dat[[2]][,9]) / as.numeric(dat[[2]][,10])
							 );
names(datCommon)[5:16] = c('nd1C','nd1U','nd1','nd2C','nd2U','nd2','d1C','d1U','d1','d2C','d2U','d2');

# add pvals
datCommon = cbind(datCommon, nd1d1=dat[[1]]$p, nd1d2=dat[[2]]$p, nd2d1=dat[[3]]$p, nd2d2=dat[[4]]$p);

# remove hits where nd1-nd2 or d1-d2 is significantly different
remove = c(which(.addFisherPvals(datCommon, cols=c(5,6,8,9))$p<.05), 
		   which(.addFisherPvals(datCommon, cols=c(11,12,14,15))$p<.05)
		   );
datCommon = datCommon[-remove, ];

################################
#########################################################
###### add information
##############################################################


###############################
### move forward with matrix where hits were required to be in all files

# add flipped nucs column for rev strand hits
datC = datCommon;
datC = .addReverseComplementColumn(datC, nucsCol=4, strandCol=3);

# add columns with methylation dinucleotide to use as factor
datC = .addMerCol(datC, strandCol=3, nucsCol=4, rcCol=which(names(datC)=='nucsRC'), l=2, colname='di');
datC = as.data.frame(cbind(datC, diH=datC$di));
datC$diH = as.character(datC$diH);

# add columns with methylation trinucleotide to use as factor
datC = .addMerCol(datC, strandCol=3, nucsCol=4, rcCol=which(names(datC)=='nucsRC'), l=3, colname='tri');
datC = as.data.frame(cbind(datC, triH=datC$tri));
datC$triH = as.character(datC$triH);

# edit diH and triH columns
for (r in 1:nrow(datC)) {
	datC[r, which(names(datC)=='diH')] = .flipNonGtoH(datC[r, which(names(datC)=='diH')], 2);
	datC[r, which(names(datC)=='triH')] = .flipNonGtoH(datC[r, which(names(datC)=='triH')], 2);
	#datC[r, which(names(datC)=='triHH')] = .flipNonGtoH(datC[r, which(names(datC)=='triHH')], 3);
}; rm(r);

datC = as.data.frame(cbind(datC, triHH=datC$triH));
datC$triHH = as.character(datC$triHH);
for (r in 1:nrow(datC)) {
	datC[r, which(names(datC)=='triHH')] = .flipNonGtoH(datC[r, which(names(datC)=='triHH')], 3);
}; rm(r);

# add column to denote who had higher methylation
datC = cbind(datC, higher=rep('D', nrow(datC)));
datC$higher = as.character(datC$higher);
tmp = cbind(datC$nd1>datC$d1, datC$nd1>datC$d2, datC$nd2>datC$d1, datC$nd2>datC$d2);
# check consistency
if (all(apply(tmp,1,function(f) length(unique(f))) == 1)) {
	datC$higher[tmp[,1]] = 'ND';
}; rm(tmp);

# add columns with nd and d averages
datC = cbind(datC, ND=rep(NA,nrow(datC)));
datC = cbind(datC, D=rep(NA,nrow(datC)));
datC$ND = apply(data.frame(datC$nd1,datC$nd2), 1, mean);
datC$D = apply(data.frame(datC$d1,datC$d2), 1, mean);

# add columns with higher val
datC = cbind(datC, highval=rep(NA,nrow(datC)));
datC$highval[datC$higher=='D'] = datC$D[datC$higher=='D'];
datC$highval[datC$higher=='ND'] = datC$ND[datC$higher=='ND'];

datC.bed = cbind(scaffold=datC[,1], pos0=as.numeric(datC$pos)-1, datC[, 2:ncol(datC)]); 
write.table(datC.bed, file='dm_all4_20X_fromR.bed', sep='\t', row.names=F, col.names=F, quote=F, na='.');
####################

###############################################
##### investigate overlaps with genomic features

# # # tmp=read.table('dm_all4_20X/Astatotilapia_burtoni.BROADAB2fix_CDS.gtf-closest_Dref',header=F,sep='\t')
# # # hits = paste(tmp[,7],tmp[,8],tmp[,9],tmp[,10],tmp[,11],sep=':');
# # # tmp = tmp[!duplicated(hits), ];

### files generated with "dm_closest_and_overlapsNEW_Feb.sh"
dir = 'dm_all4_10X'
files = list.files(dir);
files = files[grep('closest', files)];

ov0 = list();
for (f in files) {
	print(f)
	ov0[[f]] = read.table(paste(dir,'/',f,sep=''), header=F, sep='\t', fill=T, colClasses='character');
}; rm(f);

########################################
############# require same strand
########################################
#### TEs
# require same strand
TE = ov0[['Abur_final_TE.bed-closest_Dref_stranded']][,-c(2,4,5)];
TE = TE[TE$V7!='.', ];
TE = TE[abs(as.numeric(TE$V13)) <= 5000, ];
TE = .addStreamColumn(TE,3,10);
dTE = cbind(datC,TE);
dTEhit = dTE[dTE$stream=='hit', ];

#### miRNAs
# none within 5kb

#### UTRs
UTR = ov0[['Astatotilapia_burtoni.BROADAB1.UTRs.gff3-closest_Dref_stranded']][,-c(2,4,5)];
UTR = UTR[UTR$V7!='.', ];
UTR = UTR[abs(as.numeric(UTR$V16)) <= 5000, ];
UTR = .addStreamColumn(UTR,3,13);
#UTR = UTR[!duplicated(apply(UTR, 1, function(f) paste(f[1:2],collapse=' '))), ];
dUTR = cbind(datC, UTR[match(apply(datC, 1, function(f) paste(f[1:2],collapse=' ')), apply(UTR, 1, function(f) paste(f[1:2],collapse=' '))), ]);

#### lncs
LNC = ov0[['abur.lnc.final.gtf-closest_Dref_stranded']][,-c(2,4,5)];
LNC = LNC[LNC$V7!='.', ];
LNC = LNC[abs(as.numeric(LNC$V16)) <= 5000, ];
LNC = .addStreamColumn(LNC,3,13);

#### known SNPs
# use unstranded since vcf had no strand info
kSNP = ov0[['Assembly_SNPs.noHeader.gff3-closest_Dref']][,-c(2,4,5)];
kSNP = kSNP[abs(as.numeric(kSNP$V16)) <= 5000, ];
kSNP = .addStreamColumn(kSNP,3,13);
dkSNP = cbind(datC, kSNP[match(apply(datC, 1, function(f) paste(f[1:2],collapse=' ')), apply(kSNP, 1, function(f) paste(f[1:2],collapse=' '))), ]);
dkSNPhit = dkSNP[which(dkSNP$stream=='hit'),];

#### exons
CDS = ov0[['Astatotilapia_burtoni.BROADAB2fix_CDS.gtf-closest_Dref_stranded']][,-c(2,4,5)];
CDS = CDS[CDS$V7!='.', ];
CDS = CDS[abs(as.numeric(CDS$V16)) <= 5000, ];
CDS = .addStreamColumn(CDS,3,13);

#### introns
INT = ov0[['Astatotilapia_burtoni.BROADAB2fix.intron.gtf-closest_Dref_stranded']][,-c(2,4,5)];
INT = INT[INT$V7!='.', ];
INT = INT[abs(as.numeric(INT$V16)) <= 5000, ];
INT = .addStreamColumn(INT,3,13);

#### individual SNPs
# use unstranded since vcf had no strand info
nd1SNP = ov0[['ATCACG-accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf-closest_Dref']][,-c(2,4,5)];
nd1SNP = nd1SNP[nd1SNP$V7!=-1, ];
nd1SNP = nd1SNP[abs(as.numeric(nd1SNP$V17)) <= 5000, ];
nd1SNP = .addStreamColumn(nd1SNP,3,14);

nd2SNP = ov0[['TGACCA-accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf-closest_Dref']][,-c(2,4,5)];
nd2SNP = nd2SNP[nd2SNP$V7!=-1, ];
nd2SNP = nd2SNP[abs(as.numeric(nd2SNP$V17)) <= 5000, ];
nd2SNP = .addStreamColumn(nd2SNP,3,14);

ndSNP = intersect(apply(nd1SNP,1,function(f) paste(f[4:5],collapse=':')),   
				  apply(nd2SNP,1,function(f) paste(f[4:5],collapse=':'))
				  );
nd1SNP = nd1SNP[apply(nd1SNP,1,function(f) paste(f[4:5],collapse=':')) %in% ndSNP, ];
nd2SNP = nd2SNP[apply(nd2SNP,1,function(f) paste(f[4:5],collapse=':')) %in% ndSNP, ];

d1SNP = ov0[['CGATGT-accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf-closest_Dref']][,-c(2,4,5)];
d1SNP = d1SNP[d1SNP$V7!=-1, ];
d1SNP = d1SNP[abs(as.numeric(d1SNP$V17)) <= 5000, ];
d1SNP = .addStreamColumn(d1SNP,3,14);

d2SNP = ov0[['TTAGGC-accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf-closest_Dref']][,-c(2,4,5)];
d2SNP = d2SNP[d2SNP$V7!=-1, ];
d2SNP = d2SNP[abs(as.numeric(d2SNP$V17)) <= 5000, ];
d2SNP = .addStreamColumn(d2SNP,3,14);

dSNP = intersect(apply(d1SNP,1,function(f) paste(f[4:5],collapse=':')),   
				  apply(d2SNP,1,function(f) paste(f[4:5],collapse=':'))
				  );
d1SNP = d1SNP[apply(d1SNP,1,function(f) paste(f[4:5],collapse=':')) %in% dSNP, ];
d2SNP = d2SNP[apply(d2SNP,1,function(f) paste(f[4:5],collapse=':')) %in% dSNP, ];

#### CNVs
# none within 5kb

#### microsatellites
MIC = ov0[['microsatellites.noHeader.gff3-closest_Dref_stranded']][,-c(2,4,5)];
MIC = MIC[MIC$V7!='.', ];
MIC = MIC[abs(as.numeric(MIC$V16)) <= 5000, ];
MIC = .addStreamColumn(MIC,3,13);


ov = list(TE=TE, UTR=UTR, CDS=CDS, INT=INT, LNC=LNC, MIC=MIC, kSNP=kSNP, nd1SNP=nd1SNP, nd2SNP=nd2SNP, d1SNP=d1SNP, d2SNP=d2SNP);
dmsites = lapply(ov, function(ff) apply(ff,1,function(f) paste(f[1:2],collapse=':')));

for (f in 1:length(ov)) {
	ov[[f]] = cbind(ov[[f]], site=rep(NA,nrow(ov[[f]])));
	ov[[f]]$site = dmsites[[f]];
}; rm(f);

sitesCDSINT = intersect(dmsites$CDS,dmsites$INT);    # careful includes downstream hits
sitesCDSINTTE = intersect(dmsites$TE, sitesCDSINT);





d = datC;
# re-arrange d
rownames(d) = apply(d[,1:2], 1, function(f) paste(f, collapse=':'));
d = d[, -c(1,2)];
d = d[, c(1:2,19,3:4,6:7,9:10,12:13,5,8,11,14,26:27,25,28,15:18,20:24)];			


########################################################################
# investigate hits for CDS,INT,TE together

tmpTE = ov$TE;
tmpTE = cbind(tmpTE[,1:4], 
			  rep(NA,nrow(tmpTE)), 
			  tmpTE[,c(7,5,6)], 
			  rep(NA,nrow(tmpTE)), 
			  tmpTE[,9], 
			  rep(NA,nrow(tmpTE)),
			  tmpTE[,c(8,10:12)]
			  );
names(tmpTE) = names(ov$CDS);

CDSINTTE = list();
for (s in 1:length(sitesCDSINTTE)) {
	this = sitesCDSINTTE[s];
	CDSINTTE[[this]] = rbind(ov$CDS[ov$CDS$site == this, ], 
						     ov$INT[ov$INT$site == this, ], 
						     tmpTE[tmpTE$site == this, ]
						     );
}; rm(s,this);

		
########################################################################

######GET RID OF DOWNSTREAM HITS BEFORE COMBINING WITH d ########

tmp = TE[TE$stream != 'down', ];
#tmp = TE[match(rownames(d), apply(TE,1,function(f) paste(f[1:2],collapse=':'))), ];
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, TEtype=tmp$V10, TEdist=tmp$V13, TEstream=tmp$stream);

tmp = UTR[UTR$stream != 'down', ];
#tmp = UTR[match(rownames(d), apply(UTR,1,function(f) paste(f[1:2],collapse=':'))), ]; ## losing duplicates
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, UTRid=tmp$V15, UTRtype=tmp$V9, UTRdist=tmp$V16, UTRstream=tmp$stream);

tmp = CDS[CDS$stream != 'down', ];
#tmp = CDS[match(rownames(d), apply(CDS,1,function(f) paste(f[1:2],collapse=':'))), ]; ## losing duplicates and potential exon1 info
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, CDSid=apply(tmp, 1, function(f) paste(unlist(strsplit(f[12],'; '))[2:3], collapse='; ')), 		
			 CDSframe=tmp$V14, CDSdist=tmp$V16, CDSstream=tmp$stream
			 );

tmp = INT[INT$stream != 'down', ];			 
#tmp = INT[match(rownames(d), apply(INT,1,function(f) paste(f[1:2],collapse=':'))), ]; ## losing duplicates and intron1 info
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, INTid=apply(tmp, 1, function(f) paste(unlist(strsplit(f[12],'; '))[2:3], collapse='; ')), 
			 INTdist=tmp$V16, INTstream=tmp$stream
			 );

tmp = LNC[LNC$stream != 'down', ];		 
#tmp = LNC[match(rownames(d), apply(LNC,1,function(f) paste(f[1:2],collapse=':'))), ];
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, LNCid=apply(tmp, 1, function(f) paste(unlist(strsplit(f[12],'; '))[2:3], collapse='; ')), 
			 LNCdist=tmp$V16, LNCstream=tmp$stream);

tmp = MIC[MIC$stream != 'down', ];	 
#tmp = MIC[match(rownames(d), apply(MIC,1,function(f) paste(f[1:2],collapse=':'))), ];
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, MICinfo=tmp$V15, MICdist=tmp$V16, MICstream=tmp$stream);	

tmp = kSNP[kSNP$stream != 'down', ];
#tmp = kSNP[match(rownames(d), apply(kSNP,1,function(f) paste(f[1:2],collapse=':'))), ];
tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
d = cbind(d, kSNPinfo=tmp$V15, kSNPdist=tmp$V16, kSNPstream=tmp$stream);

for (snp in c('nd1SNP','nd2SNP','d1SNP','d2SNP')) {
	tmp = get(snp)[get(snp)$stream != 'down', ];
	#tmp = get(snp)[match(rownames(d), apply(get(snp),1,function(f) paste(f[1:2],collapse=':'))), ];
	tmp = tmp[match(rownames(d), apply(tmp,1,function(f) paste(f[1:2],collapse=':'))), ];
	d = cbind(d, paste(tmp$V10,',',tmp$V11,sep=''), tmp$V12, tmp$V16, tmp$V17, tmp$stream);
	names(d)[(ncol(d)-4):ncol(d)] = c(paste(snp,'nt',sep=''), paste(snp,'score',sep=''), paste(snp,'info',sep=''), paste(snp,'dist',sep=''), paste(snp,'stream',sep=''));
}; rm(snp)
		 
########################################
#####################################################
########################################

# get total number of 2mers and 3mers
mer2 = read.table('/Volumes/fishstudies/_Burtoni_genome_files/2mer_counts.jf.dump');
mer3 = read.table('/Volumes/fishstudies/_Burtoni_genome_files/3mer_counts.jf.dump');

numCG = mer2[mer2[,1] == 'CG', 2];
numCH = sum(mer2[mer2[,1] %in% c('CA','CC','CT'), 2]);
numDi = numCG + numCH;

fisher.test(matrix(c(sum(d$diH=='CG'), sum(d$diH=='CH'), numCG-sum(d$diH=='CG'), numCH-sum(d$diH=='CH')), ncol=2));
fisher.test(matrix(c(sum(d$diH=='CH'), sum(d$diH=='CG'), numCH-sum(d$diH=='CH'), numCG-sum(d$diH=='CG')), ncol=2));

numCGG = mer3[mer3[,1] == 'CGG', 2];
numCGH = sum(mer3[mer3[,1] %in% c('CGA','CGC','CGT'), 2]);
numCHG = sum(mer3[mer3[,1] %in% c('CAG','CCG','CTG'), 2]);
numCHH = sum(mer3[mer3[,1] != 'CGG', 2]);

chisq.test(matrix(c(sum(d$triHH=='CGG'), sum(d$triHH=='CGH'), sum(d$triHH=='CHG'), sum(d$triHH=='CHH'), numCGG-sum(d$triHH=='CGG'), numCGH-sum(d$triHH=='CGH'), numCHG-sum(d$triHH=='CHG'), numCHH-sum(d$triHH=='CHH')), ncol=2));

########################################
#####################################################
########################################

st = d[, grep('stream',names(d))];
hits = apply(st, 1, function(f) f[!is.na(f)]);
hitsnodown = lapply(hits, function(f) f[f!='down']);
hitsdown = names(sapply(hitsnodown,length))[sapply(hitsnodown,length)==0];

# see summary of how many features each site hit 
table(sapply(hitsnodown, length))

# find most common multihits
hittypesum = sapply(hitsnodown, function(f) paste(names(f),collapse=':'));
sort(table(hittypesum));

dnd = d[!(rownames(d) %in% hitsdown), ];
stnodown = dnd[, grep('stream',names(dnd))];
dtnodown = dnd[, grep('[A-Z]dist',names(dnd))];

dnd = .addClosestHitColumn(dnd, grep('dist',names(dnd)));
alldist = .getDistanceToUpstreamAndHits(dnd, grep('stream',names(dnd)), grep('[A-Z]dist',names(dnd)));
# get rid of downstream hit distances
alldist = lapply(alldist, function(f) f[f$type!='down',]);
.barplotStackedFromVecs(list(sort(sapply(alldist, function(f) nrow(f)/137)[c(1:4,6:11)])))
pie(sort(sapply(alldist, function(f) nrow(f)/137)[c(1:4,6:11)]))

# compute percentage of sites with upstream/hits to each feature type
sapply(alldist, nrow) / nrow(d);			# use nrow(d) to get percentage based on total number of DM

# get summary of all hits/upstream hit distances
tmpdists = as.numeric(unlist(alldist));
summary(tmpdists[!is.na(tmpdists)]);

lapply(alldist, function(f) summary(f$dist[f$type != 'down']));


# get all distances
tmp = as.matrix(dtnodown);
dlist = list();
for (j in 1:ncol(tmp)) {
	ttt = abs(as.numeric(tmp[,j]))
	dlist[[gsub('dist','',colnames(tmp)[j])]] = ttt[!is.na(ttt)];
}; rm(j,ttt);
dlist = dlist[sapply(dlist,length) > 0];


dndCDSINT = dnd[rownames(dnd) %in% sitesCDSINT, ];
# since sitesCDSINT includes downstream hit info
dndCDSINT = dndCDSINT[apply(dndCDSINT, 1, function(f) any(grepl('CDS|INT', f))), ];


########################################
#####################################################
########################################


verboseBoxplot(dnd$highval[dnd$strand=='+'], dnd$TEstream[dnd$strand=='+'] %in% c('hit','up'))


# transposons
dndTE = dnd[dnd$TEstream %in% c('up','hit'), ];
TEnd = TE[TE$stream %in% c('hit','up'), ];
TEids = apply(TEnd,1,function(f) paste(f[4:6],collapse=':'));

summary(abs(as.numeric(dndTE$TEdist)));
summary(abs(as.numeric(dndTE$TEdist[dndTE$closest=='TE'])));


cnames = c('+','-','CG','CH','D','ND','closest','DNA','LINE','LTR','SINE');
TEsites = as.data.frame(matrix(nrow=nrow(dndTE), ncol=length(cnames), dimnames=list(rownames(dndTE), cnames)));	
TEsites$'+' = as.numeric(dndTE$strand=='+');
TEsites$'-' = as.numeric(dndTE$strand=='-');
TEsites$CG = as.numeric(dndTE$diH=='CG');
TEsites$CH = as.numeric(dndTE$diH=='CH');
TEsites$D = as.numeric(dndTE$higher=='D');
TEsites$ND = as.numeric(dndTE$higher=='ND');
TEsites$closest = as.numeric(grepl('TE', dndTE$closest));
TEsites$DNA = as.numeric(dndTE$TEtype=='DNA');
TEsites$LINE = as.numeric(dndTE$TEtype=='LINE');
TEsites$SINE = as.numeric(dndTE$TEtype=='SINE');
TEsites$LTR = as.numeric(dndTE$TEtype=='LTR');

heatmap(as.matrix(TEsites[order(apply(TEsites, 1, sum)), order(apply(TEsites, 2, sum))]),scale='none',Rowv=NA,Colv=NA)
tmp = TEsites[TEsites$closest==1,]
tmp= tmp[,-which(names(tmp)=='closest')];
heatmap(as.matrix(tmp[order(apply(tmp, 1, sum)), order(apply(tmp, 2, sum))]),scale='none',Rowv=NA,Colv=NA);


dndTEcl = dndTE[grepl('TE',dndTE$closest), ];
summary(as.numeric(dndTEcl$closestdist)[as.numeric(dndTEcl$closestdist) != 0]);



tmp = dndTE$nucs;
tmp[dndTE$strand=='-'] = dndTE$nucsRC[dndTE$strand=='-'];
.seqLogoFromSeqVec(tmp);



#####################################################

# UTRs
dndUTR = dnd[!is.na(dnd$UTRstream),];
summary(abs(as.numeric(dndUTR$UTRdist)));

dndUTR[grep('three',dndUTR$UTRtype),];
dndUTR[grep('five',dndUTR$UTRtype),]

dndUTR3 = dndUTR[dndUTR$UTRtype=='three_prime_utr', ];
dndUTR5 = dndUTR[dndUTR$UTRtype=='five_prime_utr', ];


table(apply(dndUTR3[,c(72,74,76,78,80,82,84)], 1, function(f) paste(f[!is.na(f)],collapse=':')));
dndUTR3[dndUTR3$closest=='CDS',];

# get CDS ids
tmp = unlist(strsplit(unlist(strsplit(dndUTR3[dndUTR3$closest=='CDS',]$CDSid, '; ')), ' '));
tmp = tmp[grep('^ab', tmp)];
tmp = unlist(lapply(strsplit(tmp,'\\.'), function(f) paste(f[1:4],collapse='.')));
for(i in tmp){ print(CDS[grep(i,CDS$V15), ])  }; rm(i);

cnames = c('+','-','CG','CH','D','ND','closest','three_prime','five_prime');
UTRsites = as.data.frame(matrix(nrow=nrow(dndUTR), ncol=length(cnames), dimnames=list(rownames(dndUTR), cnames)));	
UTRsites$'+' = as.numeric(dndUTR$strand=='+');
UTRsites$'-' = as.numeric(dndUTR$strand=='-');
UTRsites$CG = as.numeric(dndUTR$diH=='CG');
UTRsites$CH = as.numeric(dndUTR$diH=='CH');
UTRsites$D = as.numeric(dndUTR$higher=='D');
UTRsites$ND = as.numeric(dndUTR$higher=='ND');
UTRsites$closest = as.numeric(grepl('UTR', dndUTR$closest));
UTRsites$three_prime = as.numeric(dndUTR$UTRtype=='three_prime_utr');
UTRsites$five_prime = as.numeric(dndUTR$UTRtype=='five_prime_utr');
heatmap(as.matrix(UTRsites[order(apply(UTRsites, 1, sum)), order(apply(UTRsites, 2, sum))]),scale='none',Rowv=NA,Colv=NA)
tmp = UTRsites[UTRsites$closest==1,]
tmp= tmp[,-which(names(tmp)=='closest')];
heatmap(as.matrix(tmp[order(apply(tmp, 1, sum)), order(apply(tmp, 2, sum))]),scale='none',Rowv=NA,Colv=NA);

#####################################################

# CDSs
dndCDS = dnd[!is.na(dnd$CDSstream),];
summary(abs(as.numeric(dndCDS$CDSdist)));
dndCDS[grepl('exon_number 1$', dndCDS$CDSid), ];

tmp = unlist(strsplit(unlist(strsplit(dndCDS[grepl('exon_number 1$', dndCDS$CDSid),]$CDSid, '; ')), ' '));
tmp = tmp[grep('^ab', tmp)];
tmp = unlist(lapply(strsplit(tmp,'\\.'), function(f) paste(f[1:4],collapse='.')));
for(i in tmp){ print(CDS[grep(i,CDS$V15), ])  }; rm(i);


dndCDSe1 = dndCDS[grepl('exon_number 1$', dndCDS$CDSid), ];

#####################################################
# kSNPs

tmp = kSNP[match(rownames(d), apply(kSNP,1,function(f) paste(f[1:2],collapse=':'))), ];
tmp2 = tmp[tmp$V16 %in% c(-2,-1,0,1,2),];

closekSNPs = apply(tmp2[,1:2],1,function(f) paste(f,collapse=':'));
dndSNP = dnd[rownames(dnd) %in% closekSNPs,];