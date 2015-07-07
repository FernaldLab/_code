
library(biomaRt);

fish=c('gaculeatus_gene_ensembl', 'drerio_gene_ensembl', 'olatipes_gene_ensembl', 'tnigroviridis_gene_ensembl');
marts = vector(mode='list', length=length(fish));
for (f in 1:length(fish)) {
	marts[[f]] = useMart('ensembl', dataset=fish[f])
}; rm(f);
names(marts) = c('gac', 'dar', 'orl', 'tni');

# attributes = c('ensembl_gene_id','entrezgene','hgnc_symbol', 'external_gene_id');
# filters = 'ensembl_gene_id';
# values = mod$geneAnnos$ivory[,2];
# hits = matrix(ncol=length(attributes));
# colnames(hits) = attributes;
# for (f in 1:length(marts)) {
	# result = getBM(attributes=attributes, filters=filters, values=values, mart=marts[[f]]);
	# if (nrow(result) > 0) {
		# hits = rbind(hits, result);
	# }
# }; rm(f, result);
# if (sum(is.na(hits[1, ])) == ncol(hits)) {
	# hits = hits[-1, ];
# }

hsa = useMart('ensembl','hsapiens_gene_ensembl');
getBM(attributes=c('go_id'),filters='entrezgene',values=mod$geneIDs$ivory,hsa);
library(GO.db);
goterms = as.list(GOTERM);


netIDs = getAllGeneSymbolsAndHumanEntrezIDs(DATA,annos,IDs,F);
allgo = getBM(attributes=c('go_id'), filters='entrezgene', values=unlist(netIDs$allGenesIDs), hsa);

#entrezList=mod$geneIDs
modGO = getModGOfromEntrez(mod$geneIDs);