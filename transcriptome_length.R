rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_Burtoni_annotations/');

gff = read.table('ref_AstBur1.0_scaffolds.clean.translate.final.gff3',sep='\t',header=F);

gffCDS = split(gff, gff$V3)$CDS;

tmp = unlist(strsplit(gffCDS$V9, 'gene='));
tmp = unlist(strsplit(tmp[grep('^LOC',tmp)], ';'));
tmp = tmp[grep('^LOC',tmp)];

gffCDS2 = split(gffCDS, tmp);

gffCDS2ls = c();
for (g in 1:length(gffCDS2)) {#print(names(gffCDS2)[g])
	if (i %% 1000 == 0) {print(i)}
	this = gffCDS2[[g]];
	tmp = split(this, grep('^ID', unlist(strsplit(this$V9, ';')), value=T));#print(tmp)
	tmpls = c();
	for (i in 1:length(tmp)) {
		tmpls = c(tmpls, sum(abs(tmp[[i]]$V5 - tmp[[i]]$V4)))
	}
	gffCDS2ls = c(gffCDS2ls, max(tmpls));
}; rm(g,this,tmp,i,tmpls)