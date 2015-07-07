# convert fish esembl IDs to human

library(biomaRt);

##### pulling ids from gtf annotated with gene syms and ensembl ids
dat = nd2;					#nd2 from _scripts/snp_gene_hits.py
ens = c();
for (row in 1:nrow(dat)) {
	this_str = dat[row, 20];
	tmp = strsplit(strsplit(this_str, ';')[[1]][2], '_')[[1]];
	if (length(tmp) == 4) {
		ens = c(ens, tmp[4]);
	}
}; rm(row,this_str,tmp);


##### for single list
ens = read.table('gene_exp.diff.sig.GeneNamesENS.txt')[,1];
#ens = read.table('../../_Burtoni_annotations/geneNamesTree_AB_noNONE',header=T,sep='\t',quote='')[,2];
ensSplit = split(ens,substring(ens,4,6));
sp = names(table(substring(ens,4,6)));
lookup = matrix(c(sp, 'drerio', 'gaculeatus', 'olatipes', 'tnigroviridis'),ncol=2);

marts = list();
for (s in 1:length(sp)) {
	print(paste(lookup[s,2], '_gene_ensembl', sep=''))
	marts[[s]] = useMart('ensembl', paste(lookup[s,2], '_gene_ensembl', sep=''));
	names(marts)[s] = lookup[s,2];
}; rm(s);

hits = list();
for (s in 1:length(marts)) {
	print(names(marts)[s]);
	vals = ensSplit[names(ensSplit)==lookup[lookup[,2]==names(marts)[s],1]][[1]];
	hits[[s]] = getBM(attributes='hsapiens_homolog_ensembl_gene', filters='ensembl_gene_id', values=vals, mart=marts[[s]]);
	names(hits)[s] = names(marts)[s];
}; rm(s);

write.table(unlist(hits), 'gene_exp.diff.sig.GeneNamesENS_hsapiens.txt',quote=F,row.names=F,col.names=F)
#write.table(unlist(hits), 'allGeneNamesENS_hsapiens.txt',quote=F,row.names=F,col.names=F)

##### DAVID GO
d = DAVIDWebService$new(email='ahilliar@stanford.edu');
#geneList = unlist(hits);
addList(d, inputIds=geneList, idType='ENSEMBL_GENE_ID', listName='background', listType='Background');
# for (s in 1:length(hits)) {
	# print(names(hits)[s]);
	# addList(d, inputIds=hits[[s]][,1], idType='ENSEMBL_GENE_ID', listName=names(hits)[s], listType='Gene');
# }; rm(s);
# setCurrentBackgroundPosition(d, 2);
# setCurrentGeneListPosition(d, 1);

# charts = list();
# for (s in 1:length(hits)) {
	# print(names(hits)[s]);
	# setCurrentGeneListPosition(d, s);
	# charts[[s]] = getFunctionalAnnotationChart(d);
	# names(charts)[s] = names(hits)[s];
# }; rm(s);