setwd('/Volumes/fishstudies/_gatkbp_backup/tophat/');
files = paste(list.files(), '/accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf', sep='');

snps = list();
for (f in 1:length(files)) {
	print(files[f]);
	snps[[f]] = read.table(files[f], sep='\t', colClasses='character');
}; rm(f);

names(snps) = gsub('/accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf', '', files);


### find status specific snps

groupLookup = data.frame(num=c('ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA'),
						  group=c('ND', 'D', 'D', 'ND')
						  );
						 
# shared ND snps
ND1 = paste(snps$ATCACG$V1, snps$ATCACG$V2, sep=':');
ND2 = paste(snps$TGACCA$V1, snps$TGACCA$V2, sep=':');
ND = ND1[ND1 %in% ND2];

# shared D snps
D1 = paste(snps$CGATGT$V1, snps$CGATGT$V2, sep=':');
D2 = paste(snps$TTAGGC$V1, snps$TTAGGC$V2, sep=':');
D = D1[D1 %in% D2];

# snps common to both
comm = D[D %in% ND];
# D unique
Du = D[!(D %in% ND)];
# ND unique
NDu = ND[!(ND %in% D)];

# filter original data
snpsFilt = snps;
for (s in 1:length(snpsFilt)) {
	gp = groupLookup[match(names(snps)[s], groupLookup[,1]), 2];
	if (gp=='ND') {uniSnps = NDu};
	if (gp=='D') {uniSnps = Du};
	keep = paste(snps[[s]]$V1, snps[[s]]$V2, sep=':') %in% uniSnps;
	snpsFilt[[s]] = snpsFilt[[s]][keep, ];
}; rm(s,gp,uniSnps,keep);


bsnps0 = read.table('../../_Burtoni_annotations/Assembly_SNPs.noHeader.gff3', header=F, sep='\t', colClasses='character');
bsnps = paste(bsnps0$V1, bsnps0$V5, sep=':');

snpsFiltAssembly = snpsFilt;
snpsFiltNovel = snpsFilt;
for (s in 1:length(snpsFilt)) {
	assembly = paste(snpsFilt[[s]]$V1, snpsFilt[[s]]$V2, sep=':') %in% bsnps;
	snpsFiltAssembly[[s]] = snpsFilt[[s]][assembly, ];
	snpsFiltNovel[[s]] = snpsFilt[[s]][!assembly, ];
}; rm(s);

for (s in 1:length(snpsFilt)) {
	gp = groupLookup[match(names(snpsFilt)[s], groupLookup[,1]), 2];
	fileN = paste(files[s], '_', gp, 'novel', sep='');print(fileN)
	fileA = paste(files[s], '_', gp, 'assembly', sep='');print(fileA)
	write.table(snpsFiltNovel[[s]], file=fileN, quote=F, row.names=F, col.names=F, sep='\t');
	write.table(snpsFiltAssembly[[s]], file=fileA, quote=F, row.names=F, col.names=F, sep='\t');
}; rm(s,gp,fileN,fileA);

########### continued in /Volumes/fishstudies/_scripts/snps_bqsr.sh


######### load some result from /Volumes/fishstudies/_scripts/snps_bqsr.sh
rm(list=ls());
setwd('/Volumes/fishstudies/_gatkbp_backup/tophat/');
dirs = list.files();
snps = vector(mode='list', length=length(dirs));
names(snps) = dirs;

for (d in 1:length(dirs)) {
	if (dirs[d] %in% c('ATCACG','TGACCA')) {
		files = c(paste(dirs[d], 'accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf_NDassemblyQ50.vcf-GTF.hitFilt', sep='/'),
				  paste(dirs[d], 'accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf_NDnovelQ50.vcf-GTF.hitFilt', sep='/')
				  );print(files)
	} else if (dirs[d] %in% c('CGATGT','TTAGGC')) {
		files = c(paste(dirs[d], 'accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf_DassemblyQ50.vcf-GTF.hitFilt', sep='/'),
				  paste(dirs[d], 'accepted_hits.RG.DD.reorderSam.realign.BQSR.bam.SNPs.vcf_DnovelQ50.vcf-GTF.hitFilt', sep='/')
				  );print(files)
	} else {
		stop('WRONG DIR');
	}
	print(dirs[d])
	for (f in 1:length(files)) {
		print(files[f])
		snps[[d]][[f]] = read.table(files[f], header=F, sep='\t');
		names(snps[[d]])[f] = files[f];
	}
}; rm(d,f);

anno = read.table('../../_Burtoni_annotations/geneNamesTree_AB_noNONE',header=T,sep='\t', quote="");

snps2 = snps;
for (s in 1:length(snps)) {
	this_s = snps[[s]];
	for (f in 1:length(this_s)) {
		x = gsub('gene_id ', '', apply(this_s[[f]], 1, function(fxn) strsplit(fxn[19],';')[[1]][[1]]));
		this_s[[f]] = cbind(this_s[[f]], anno[match(x, anno$gene_id), ]);
		snps2[[s]][[f]] = this_s[[f]];
	}
}; rm(s,f,this_s);


snps2n = list(ATCACG=snps2$ATCACG[[2]], CGATGT=snps2$CGATGT[[2]], TGACCA=snps2$TGACCA[[2]], TTAGGC=snps2$TTAGGC[[2]]);
for (s in 1:length(snps2n)) {
	snps2n[[s]] = split(snps2n[[s]], as.factor(snps2n[[s]]$symbol));
}; rm(s);

snps2nND = names(snps2n$ATCACG)[names(snps2n$ATCACG) %in% names(snps2n$TGACCA)];
snps2nD = names(snps2n$CGATGT)[names(snps2n$CGATGT) %in% names(snps2n$TTAGGC)];

snps2nNDonly = snps2nND[!(snps2nND %in% snps2nD)];
snps2nDonly = snps2nD[!(snps2nD %in% snps2nND)];

snps2nNDonly.mat = matrix(nrow=length(snps2nNDonly), ncol=2);
rownames(snps2nNDonly.mat) = snps2nNDonly;
colnames(snps2nNDonly.mat) = c('ATCACG','TGACCA');
for (row in 1:nrow(snps2nNDonly.mat)) {
	sym = rownames(snps2nNDonly.mat)[row];
	for (col in 1:ncol(snps2nNDonly.mat)) {
		ind = match(colnames(snps2nNDonly.mat)[col], names(snps2n));
		x = snps2n[[ind]][names(snps2n[[ind]])==sym][[1]];#print(x)
		snps2nNDonly.mat[row,col] = length(unique(paste(x[,1],x[,2])));
	}
}; rm(row,sym,col,ind);

snps2nDonly.mat = matrix(nrow=length(snps2nDonly), ncol=2);
rownames(snps2nDonly.mat) = snps2nDonly;
colnames(snps2nDonly.mat) = c('CGATGT','TTAGGC');
for (row in 1:nrow(snps2nDonly.mat)) {
	sym = rownames(snps2nDonly.mat)[row];
	for (col in 1:ncol(snps2nDonly.mat)) {
		ind = match(colnames(snps2nDonly.mat)[col], names(snps2n));
		x = snps2n[[ind]][names(snps2n[[ind]])==sym][[1]];#print(x)
		snps2nDonly.mat[row,col] = length(unique(paste(x[,1],x[,2])));
	}
}; rm(row,sym,col,ind);

de = read.table('/Volumes/fishstudies/_gatkbp_backup/cuffdiff_postRealign_bqsr_fa/gene_exp.diff.sig.GeneNames',header=F,sep='\t');