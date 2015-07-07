mart = useMart('ensembl', 'hsapiens_gene_ensembl');
a = listAttributes(mart);
species = 'ggorilla|ggallus|pabelii|mmusculus|mmulatta|ptroglodytes|mdomestica|oanatinus|ppaniscus';
#attr = a[grepl(species, a[,1]) & grepl('homolog', a[,1]) & (grepl('dN|dS', a[,2]) | grepl('gene', a[,1])), ];
#attr = a[grepl(species, a[,1]) & grepl('homolog', a[,1]) & grepl('gene', a[,1]), ];

attr = a[grepl(species, a[,1]) & grepl('homolog', a[,1]) & (grepl('\\%', a[,2]) | grepl('gene', a[,1])), ];


ids = c('ENSG00000174469', 'ENSG00000128573', 'ENSG00000227232');
getBM(attributes=attr[1:4, 1], filters='ensembl_gene_id', values=ids, mart);
getBM(attributes=attr[5:8, 1], filters='ensembl_gene_id', values=ids, mart);