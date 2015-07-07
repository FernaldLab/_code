#setwd('/Users/ath/Documents/_Fernald_lab/_methylation');
setwd('/Volumes/fishstudies/_methylation');
rm(list=ls());
diff = read.table('3157_vs_3165.filtered.50percent.significant.gtfintersect_with_genes.bed', sep='\t', header=F);

diff = diff[as.character(diff$V14)=='exon', ];
pvals = c();
syms = c();
for (r in 1:nrow(diff)) {
	p = fisher.test(matrix(as.numeric(c(diff[r,7:8], diff[r,10:11])), ncol=2))$p.value;
	pvals = c(pvals, p);
	s = toupper(gsub(' ','',unlist(lapply(strsplit(diff$V22[r],"\\("), function(xx) xx[1]))));
	syms = c(syms, s);
}; rm(r, p, s);
diff = cbind(diff, pvals, syms);
rm(pvals, syms);

diff.split = split(diff, as.character(diff[, 22]));

diff.mat = matrix(rep(0, ncol(diff)), ncol=ncol(diff));
colnames(diff.mat) = colnames(diff);
for (g in 1:length(diff.split)) {
	this_gene = diff.split[[g]];#print(this_gene)
	if (nrow(this_gene)==1) {
		diff.mat = rbind(diff.mat, this_gene);
	} else {
		tmp = split(this_gene, paste(this_gene$V1, this_gene$V2, sep=':'));
		for (s in 1:length(tmp)) {
			if (length(unique(paste(tmp[[s]]$V1, tmp[[s]]$V2, sep=':')))==1) {
				diff.mat = rbind(diff.mat, tmp[[s]][1, ]);
			} else {
				stop('ERROR');
			}
		}
	}
}; rm(g, this_gene, tmp, s);
rm(diff.split);
diff.mat = diff.mat[-1, ];
diff.mat = cbind(diff.mat, matrix(nrow=nrow(diff.mat), ncol=4));
names(diff.mat)[26:29] = c('p1', 'p2', 'higher', 'gene_id');
for (row in 1:nrow(diff.mat)) {
	diff.mat$p1[row] = diff.mat[row,8] / sum(diff.mat[row,7:8]);
	diff.mat$p2[row] = diff.mat[row,11] / sum(diff.mat[row,10:11]);
	if (diff.mat$p1[row] > diff.mat$p2[row]) {
		diff.mat$higher[row] = '3157-ND';
	} else {
		diff.mat$higher[row] = '3165-D';
	}
	diff.mat$gene_id[row] = gsub('gene_id ', '', strsplit(diff.mat[row, 20],';')[[1]][1]);
}; rm(row);

diff.matNA = diff.mat[which(diff.mat$syms=='NOANNOTATION'), ];
diff.mat = diff.mat[-which(diff.mat$syms=='NOANNOTATION'), ];

diff.mat.split = split(diff.mat, as.character(diff.mat$V22));
diff.matNA.split = split(diff.matNA, apply(diff.matNA,1,function(f) strsplit(f[20],';')[[1]][1]));

nrows = unlist(lapply(diff.mat.split,nrow));
nrowsNA = unlist(lapply(diff.matNA.split,nrow));

dmg = unique(toupper(gsub(' ','',unlist(lapply(strsplit(names(diff.mat.split),"\\("), function(xx) xx[1])))));
dmg = unlist(IDs[names(IDs) %in% dmg]);

############

# get cufflinks output
path = '../_LYNLEY_RNAseq/cufflinksOUT/GTF/orig/';
dirs = list.files(path);
cuff = list();
for (f in 1:length(dirs)) {
	cuff[[f]] = read.table(paste(path, dirs[f], '/genes.fpkm_tracking', sep=''),header=T,sep='\t');
	names(cuff)[f] = gsub('cufflinks', '', dirs[f]);
}; rm(f);

cuff = lapply(cuff, function(f) f = f[f$FPKM_status=='OK', ]);
cuff2 = lapply(cuff, function(f) f = f[, c(1, 7, 10)]);
cuff = cuff2; rm(cuff2);






#############
############
library(org.Hs.eg.db); IDs=as.list(org.Hs.egSYMBOL2EG);library(RDAVIDWebService)
aa=read.table('~/Documents/_Fernald_lab/geneNamesTree_AB_noNONE',header=T,sep='\t', quote="");
g=aa$symbol;
bg = unique(toupper(gsub(' ','',unlist(lapply(strsplit(g,"\\("), function(xx) xx[1])))));
bg2=unlist(IDs[names(IDs) %in% bg]);
#write.table(bg2,file='dmgBG.txt',quote=F,row.names=F)

d=DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d, inputIds=bg2, idType='ENTREZ_GENE_ID', listName='genome_bg', listType='Background');
addList(d, inputIds=dmg, idType='ENTREZ_GENE_ID', listName='dmg', listType='Gene');
setCurrentBackgroundPosition(d, 2);
chart = getFunctionalAnnotationChart(d);

#############


tfup=read.table('dmg_all_promoters_hs_2000_1.00_comp-trecount.list2ref_fdradjusted_pup.txt', header=T, sep='\t', row.names=1);
upnames0=rownames(which(tfup <.05, arr.ind=T));
upnames=unlist(lapply(strsplit(gsub('^[^$]*.', '', upnames0), '_'), function(f) f[[1]]));

tfdn=read.table('~/Downloads/dmg_all_promoters_hs_2000_1.00_comp-trecount.list2ref_fdradjusted_pdown.txt', header=T, sep='\t', row.names=1);
downnames0=rownames(which(tfdn <.05, arr.ind=T));
downnames=unlist(lapply(strsplit(gsub('^[^$]*.', '', downnames0), '_'), function(f) f[[1]]));
downnamesE = unlist(IDs[names(IDs) %in% downnames]);




###############

library(biomaRt); library(RDAVIDWebService)
mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');

aa=read.table('/Volumes/fishstudies/_Burtoni_annotations/geneNamesTree_AB_noNONE',header=T,sep='\t', quote="");
g=aa$symbol;
bg = unique(toupper(gsub(' ','',unlist(lapply(strsplit(g,"\\("), function(xx) xx[1])))));
bg_ens = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=bg, mart=mart);
bg_ens = bg_ens[,1];

g = unique(diff.mat$V22);
g = toupper(gsub(' ','',unlist(lapply(strsplit(g,"\\("), function(xx) xx[1]))));
g_ens = getBM(attributes='ensembl_gene_id', filters='hgnc_symbol', values=g, mart=mart)[,1];

d = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(d, inputIds=bg_ens, idType='ENSEMBL_GENE_ID', listName='BG', listType='Background');
addList(d, inputIds=g_ens, idType='ENSEMBL_GENE_ID', listName='dm_genes', listType='Gene');
setCurrentBackgroundPosition(d, which(getBackgroundListNames(d) == 'BG'));
fc = getFunctionalAnnotationChart(d)