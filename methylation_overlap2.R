setwd('/Volumes/handsfiler$/FishStudies/_methylation/dm_multiple-overlaps');
rm(list=ls());

files = list.files()[grepl('filt$', list.files(), ignore.case=T)];
ol0 = list();
for (f in 1:length(files)) {
	ol0[[f]] =  read.table(files[f], header=F, sep='\t');
}; rm(f);
names(ol0) = gsub('dm_multiple.bed-', '', files, fixed=T);

ol0$comboGTFFILT = ol0$comboGTFFILT[ol0$comboGTFFILT$V6!='CDS', ];
fpaste = paste(ol0$comboGTFFILT$V5,ol0$comboGTFFILT$V6,sep=':');

# no miRNA hits
ol0.cnv = ol0$comboGTFFILT[fpaste=='BROAD:cnv', ];
ol0.TE = ol0$comboGTFFILT[fpaste %in% c('BROAD:DNA','BROAD:LINE','BROAD:LTR','BROAD:RC','BROAD:SINE','BROAD:Unknown'), ];
ol0.utr5 = ol0$comboGTFFILT[fpaste=='BROAD:five_prime_utr', ];
ol0.utr3 = ol0$comboGTFFILT[fpaste=='BROAD:three_prime_utr', ];
ol0.SNPassembly = ol0$comboGTFFILT[fpaste=='A_burtoni_broad_scaffolds_v1:SNP', ];
ol0.lnc = ol0$comboGTFFILT[fpaste=='cuffLNC:exon', ];
ol0.cds = ol0$comboGTFFILT[fpaste=='protein_coding:exon', ];

ol = list(cnv=ol0.cnv, TE=ol0.TE, utr5=ol0.utr5, utr3=ol0.utr3, SNP.assembly=ol0.SNPassembly, lncRNA=ol0.lnc, CDS=ol0.cds, microsat=ol0$microsatFILT, DErealign=ol0$DE_postRealignFILT, DErealign.bqsr.fa=ol0$DE_postRealign_bqsr_faFILT, SNP.ATCACG=ol0$SNP_ATCACG_bqsrFilt, SNP.CGATGT=ol0$SNP_CGATGT_bqsrFilt, SNP.TGACCA=ol0$SNP_TGACCA_bqsrFilt, SNP.TTAGGC=ol0$SNP_TTAGGC_bqsrFilt);

sapply(ol, function(f) table(f[,ncol(f)]));

olPos = list();
for (f in 1:length(ol)) {
	#olPos[[f]] = unique(paste(ol[[f]][,1], ol[[f]][,2], ol[[f]][,3], sep=':'));
	olPos[[f]] = unique(paste(ol[[f]][,1], ol[[f]][,3], sep=':'));
}; rm(f);
names(olPos) = names(ol);

dm = read.table('../dm_multiple.bed');
dm = as.data.frame(cbind(dm, mat.or.vec(nr=nrow(dm), nc=length(olPos))));
names(dm)[4:ncol(dm)] = names(olPos);
for (fea in 1:length(olPos)) {
	if (length(olPos[[fea]]) > 0) {
		pos = match(olPos[[fea]], paste(dm[,1], dm[,3], sep=':'));
		dm[pos, fea+3] = 1;
	} else {
		next;
	}
}; rm(fea, pos);
dm = as.data.frame(cbind(dm, total=apply(dm[,4:ncol(dm)],1,sum)));


snp.ND = apply(dm, 1, function(f) as.numeric(f[14])+as.numeric(f[16]));
snp.D = apply(dm, 1, function(f) as.numeric(f[15])+as.numeric(f[17]));