###########################################################################################
###### setup
###########################################################################################
rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/');
library('bsseq');

load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20.RData');
load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20_tstat.RData');
load('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20_tstat_dmrs0.RData');
gc();

FIT.FiltXn = bsdalln20.fit.ns20h200mg500.5x.n20;
DMR0 = dmrs00.01;
###########################################################################################
###### 
###########################################################################################

# sampleNames(FIT.FiltXn);
# head(getCoverage(FIT.FiltXn));
# head(getMeth(FIT.FiltXn));

nTHRESH = 5;
diffTHRESH = 0.1;
DMR = subset(DMR0, n>=nTHRESH & abs(meanDiff)>=diffTHRESH);
table(DMR$direction);
summary(DMR$invdensity);

#write.table(DMR, file='aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20_tstat_dmrs_n5_md.1.txt',quote=F, row.names=F, col.names=T, sep='\t')



DMR = read.table('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20_tstat_dmrs_n3_md.1.txt',sep='\t',header=T)
par(mfrow=c(2,3))
for (i in 7:12) {
	WGCNA::verboseBoxplot(abs(DMR[,i]), DMR$direction, xlab='', ylab=names(DMR)[i], frame.plot=F, col='grey');
}; rm(i);



###########################################################################################
###### use bedtools to get overlaps
###########################################################################################

setwd('aligned.adapters.q30.m0.methratio.CG.clean.Combined_n20Smoothed.ns20h200mg500_5x_n20_tstat_dmrs_n3_md.1');

g0 = read.table('ref_AstBur1.0_scaffolds.clean.translate.final.genesOnly.gff3-closest', sep='\t');

MAXDIST = 5000;
g = subset(g0, grepl('^[a-zA-Z0-9]', V25) & abs(V26)<=MAXDIST);

#gS = split(g, paste(g[,1],":",g[,2],"-",g[,3],sep=""));

# add column for overlapping gene hits
.addGenesToDMRs = function (DMR, overlapResults) {
	gS = split(overlapResults, paste(overlapResults[,1],":",overlapResults[,2],"-",overlapResults[,3],sep=""));
	DMRgenes = as.data.frame(cbind(DMR, gene.overlap=rep(NA,nrow(DMR)), gene.dist=rep(NA,nrow(DMR))));
	for (row in 1:nrow(DMRgenes)) {
		thisDMR = paste(DMRgenes[row,1],":",DMRgenes[row,2],"-",DMRgenes[row,3],sep="");
		rownames(DMRgenes)[row] = thisDMR;
		hit = gS[match(thisDMR, names(gS))][[1]];
		if (is.null(hit)) {
			next;
		} else if (nrow(hit)==1) {
			tmp = strsplit(hit[,25], 'gene=')[[1]][2];
			tmp = strsplit(tmp, ';')[[1]][1];
			DMRgenes$gene.overlap[row] = tmp;
			DMRgenes$gene.dist[row] = hit[, ncol(hit)];
		} else {
			tmp = strsplit(hit[,25], 'gene=');
			tmp = sapply(tmp, function(f) strsplit(f[2],';')[[1]][1]);#print(thisDMR);print(tmp)
			DMRgenes$gene.overlap[row] = paste(tmp, collapse=',');
			DMRgenes$gene.dist[row] = paste(hit[, ncol(hit)], collapse=',');
		}
	}
	return(DMRgenes);
}

DMRgenes = .addGenesToDMRs(DMR, g);

### may want to curate DMRs with >1 overlap to choose one
DMRgenes[grepl(',',DMRgenes$gene.overlap), ];

# e.g. if DMR is over 5'utr of only one
# DMRgenes$gene.overlap[rownames(DMRgenes) == 'scaffold_124:214106-214361'] = 'LOC102297721';




# get homologs and human versions of gene hits
lookup = read.table('../../_Burtoni_annotations/H_burtoni_rna_blastx_FISH_ENS_top1_reciprocalBackFrom_Drer_Olat_Onil_Trub_ENS_pep_noComments_passRecipHsENS',sep='\t',header=T);

### very overwritten, could be factored, don't care at the moment
.addHumanGenes = function (DMRgenes, lookup) {
	DMRgenes = as.data.frame(cbind(DMRgenes, gene.homolog=rep(NA,nrow(DMRgenes)), gene.hs=rep(NA,nrow(DMRgenes))));
	for (row in 1:nrow(DMRgenes)) {
		thisLOC = DMRgenes$gene.overlap[row];
		# check if there no burtoni gene overlap
		if (is.na(thisLOC)) {
			next;
		# check if overlap covers two burtoni genes
		} else if (grepl(',', thisLOC, fixed=T)) {
			theseLOCs = unlist(strsplit(thisLOC, ',', fixed=T));#print(thisLOC)
			lrows = which(lookup$gene %in% theseLOCs);#print(lookup[lrows,])
			# check if neither is in lookup
			if (length(lrows) == 0) {
				next;
			} else {
				tmpLook = lookup[lrows, ];
				tmpHs = unlist(strsplit(tmpLook$hsENS, ','));
				DMRgenes$gene.homolog[row] = paste(unique(rownames(tmpLook)), collapse=',');
				if (length(tmpHs) == 0) {
					DMRgenes$gene.hs[row] = NA;
				} else if (length(tmpHs) == 1) {
					DMRgenes$gene.hs[row] = tmpHs;
				} else {
					DMRgenes$gene.hs[row] = paste(tmpHs, collapse=',');
				}
			}
		} else {
			lrow = which(lookup$gene == thisLOC);
			if (length(lrow) == 0) {
				next;
			} else {
				DMRgenes$gene.homolog[row] = paste(unique(rownames(lookup)[lrow]), collapse=',');
				tmpHs = unique(unlist(strsplit(lookup$hsENS[lrow], ',')));
				if (length(tmpHs) == 0) {
					DMRgenes$gene.hs[row] = NA;
				} else if (length(tmpHs) == 1) {
					DMRgenes$gene.hs[row] = tmpHs;
				} else {
					DMRgenes$gene.hs[row] = paste(as.vector(na.omit(tmpHs)), collapse=',');
				}
			} 
		}
	}
	return(DMRgenes);
}

DMRgenesHs = .addHumanGenes(DMRgenes, lookup);

DMRgenesHs[grepl(',',DMRgenesHs$gene.hs), ];


.checkMultiHsHits = function (DMRgenesHs, attributes, mart=NULL) {
	if (is.null(mart)) {
		mart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
	} else {
		mart = mart;
	}
	mrows = DMRgenesHs[grep(',',DMRgenesHs$gene.hs), ];
	for (r in 1:nrow(mrows)) {
		print(mrows[r, ]);
		tmp = getBM(attributes=attributes, filters='ensembl_gene_id', values=unlist(strsplit(mrows$gene.hs[r],',')), mart=mart);
		print(tmp);
		cat('\n--------------------------------------------------------------------------------------------------------------\n');
	}
	cat('----------------------------------------------------------------------\n',nrow(mrows), 'instances\n');
}

#hsMart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
attributes = c('ensembl_gene_id','description','external_gene_name','go_id','name_1006');
.checkMultiHsHits(DMRgenesHs, attributes=attributes, mart=hsMart);


####
library(RDAVIDWebService);

bg = as.vector(na.omit(unlist(strsplit(lookup$hsENS,',')))); 
genesForDavid = as.vector(na.omit(unlist(strsplit(DMRgenesHs$gene.hs,','))));
genesForDavid = genesForDavid[grepl('^ENS', genesForDavid)];
david = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(david, inputIds=bg, idType='ENSEMBL_GENE_ID', listName='bg', listType='Background');
addList(david, inputIds=genesForDavid, idType='ENSEMBL_GENE_ID', listName='DMRgenesHs', listType='Gene');
setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'bg'));
DMRgenesHsDavid0 = getFunctionalAnnotationChart(david);
DMRgenesHsDavid = subset(DMRgenesHsDavid0, FDR<=10);


genesForDavid_D = as.vector(na.omit(unlist(strsplit(DMRgenesHs$gene.hs[DMRgenesHs$direction=='hyper'],','))));
genesForDavid_D = genesForDavid_D[grepl('^ENS', genesForDavid_D)];
addList(david, inputIds=genesForDavid_D, idType='ENSEMBL_GENE_ID', listName='DMRgenesHs_D', listType='Gene');
DMRgenesHsDavid_D0 = getFunctionalAnnotationChart(david);
DMRgenesHsDavid_D = subset(DMRgenesHsDavid_D0, FDR<=10);

genesForDavid_ND = as.vector(na.omit(unlist(strsplit(DMRgenesHs$gene.hs[DMRgenesHs$direction=='hypo'],','))));
genesForDavid_ND = genesForDavid_ND[grepl('^ENS', genesForDavid_ND)];
addList(david, inputIds=genesForDavid_ND, idType='ENSEMBL_GENE_ID', listName='DMRgenesHs_ND', listType='Gene');
DMRgenesHsDavid_ND0 = getFunctionalAnnotationChart(david);
DMRgenesHsDavid_ND = subset(DMRgenesHsDavid_ND0, FDR<=10);


### remove genes that map to multiple hs ids
multiRows = which(grepl(',', DMRgenesHs$gene.hs, fixed=T));
#moreThanTwo = which(nchar(DMRgenesHs$gene.hs) > 31);
DMRgenesHsNoMulti = DMRgenesHs[-multiRows, ];
#DMRgenesHsNoMulti = DMRgenesHs[-moreThanTwo, ];


bgMultiRows = which(grepl(',',lookup$hsENS,fixed=T));
#moreThanTwoBG = which(nchar(lookup$hsENS) > 31);
bgNoMulti = as.vector(na.omit(unlist(strsplit(lookup$hsENS[-bgMultiRows],','))));
#bgNoMulti = as.vector(na.omit(unlist(strsplit(lookup$hsENS[-moreThanTwoBG],','))));

genesForDavid_D = as.vector(na.omit(unlist(strsplit(DMRgenesHsNoMulti$gene.hs[DMRgenesHsNoMulti$direction=='hyper'],','))));
genesForDavid_ND = as.vector(na.omit(unlist(strsplit(DMRgenesHsNoMulti$gene.hs[DMRgenesHsNoMulti$direction=='hypo'],','))));

david = DAVIDWebService$new(email='ahilliar@stanford.edu');
addList(david, inputIds=bgNoMulti, idType='ENSEMBL_GENE_ID', listName='bgNoMulti', listType='Background');
addList(david, inputIds=genesForDavid_D, idType='ENSEMBL_GENE_ID', listName='D', listType='Gene');
addList(david, inputIds=genesForDavid_ND, idType='ENSEMBL_GENE_ID', listName='ND', listType='Gene');

setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'bgNoMulti'));
setCurrentGeneListPosition(david, which(getGeneListNames(david) == 'D'));
DMRgenesHsDavid_D0 = getFunctionalAnnotationChart(david);
DMRgenesHsDavid_D = subset(DMRgenesHsDavid_D0, FDR<=10);

setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'bgNoMulti'));
setCurrentGeneListPosition(david, which(getGeneListNames(david) == 'ND'));
DMRgenesHsDavid_ND0 = getFunctionalAnnotationChart(david);
DMRgenesHsDavid_ND = subset(DMRgenesHsDavid_ND0, FDR<=10);

.getTermClustersFromDAVID = function (david, overlap=5L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=100L, EASEthresh=NULL) {
	cl = getClusterReport(david, overlap=overlap, initialSeed=initialSeed, finalSeed=finalSeed, linkage=linkage, kappa=kappa);
	cl = attributes(cl)$cluster;
	if (is.numeric(EASEthresh)) {
		cl = cl[sapply(cl, function(f) f$EnrichmentScore) >= EASEthresh];
	}
	return(cl);
}

BG = 'bg';
setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == BG));
setCurrentGeneListPosition(david, which(getGeneListNames(david) == 'D') );
DMRgenesHsDavid_Dclust = .getTermClustersFromDAVID(david, EASEthresh=1.2, overlap=5L, linkage=0.5, kappa=100L);

setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == BG));
setCurrentGeneListPosition(david, which(getGeneListNames(david) == 'ND') );
DMRgenesHsDavid_NDclust = .getTermClustersFromDAVID(david, EASEthresh=1.2, overlap=5L, linkage=0.5, kappa=100L);



.getTermGenes = function(chart,term) {
	if (is.numeric(term)) {
		return(unlist(strsplit(chart$Genes[term], ', ')));
	} else if (is.character(term)) {
		return(unlist(strsplit(chart$Genes[chart$Term==term], ', ')));
	} else {
		stop()
	}
}

.getTermKeywordGenes = function(chart, keyword, ignoreCase=T) {
	return(unique(unlist(strsplit(subset(chart, grepl(keyword,Term,ignore.case=ignoreCase))$Genes, ', '))));
}


hsMart = useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl');
#attributes = c('ensembl_gene_id','description','external_gene_name','go_id','name_1006');
attributes = c('ensembl_gene_id','description','external_gene_name');
getBM(attributes=attributes, 
	  filters='ensembl_gene_id', 
	  values=unique(.getDMRsWithGenes(DMRgenesHs, .getTermKeywordGenes(DMRgenesHsDavid_D, 'fibronectin'))$gene.hs), 
	  mart=hsMart
	  );




.getDMRsWithGenes = function(DMRs, genes, geneColName='gene.hs') {
	hits = c();
	colNum = which(names(DMRs) == geneColName);
	for (row in 1:nrow(DMRs)) {
		dmrgenes = unlist(strsplit(DMRs[,colNum][row], ','))
		if (length(intersect(genes, dmrgenes)) > 0) {
			hits = c(hits, TRUE);
		} else {
			hits = c(hits, FALSE);
		}
	}
	return(DMRs[hits, ]);
}


tmp = .getTermGenes(DMRgenesHsDavid_ND, 'GO:0030182~neuron differentiation');

rows = c();
for (row in 1:nrow(DMRgenesHs)) {
	sgenes = unlist(strsplit(DMRgenesHs$gene.hs[row], ',', fixed=T));
	if (any(tmp %in% sgenes)) {
		rows = c(rows, T)
	} else {
		rows = c(rows, F)
	}
	
}; rm(row,sgenes);



.convertIDsWithBiomaRt = function (biomart='ensembl', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', ids) {
	mart = useMart(biomart=biomart, dataset=dataset);
	return(getBM(attributes=c(fromID, toID), filters=fromID, values=ids, mart=mart));
}
bgNoMultiEntrez = .convertIDsWithBiomaRt(ids=bgNoMulti)$entrezgene;
bgEntrez = .convertIDsWithBiomaRt(ids=bg)$entrezgene;
genesForDavid_DEntrez = .convertIDsWithBiomaRt(ids=genesForDavid_D)$entrezgene;
genesForDavid_NDEntrez = .convertIDsWithBiomaRt(ids=genesForDavid_ND)$entrezgene;
library(GOFunction);
.GOFunctionAllOntologies = function (genes, refGenes, ...) {
		cat('\n==========  BP  ==========\n');
		tmpBP = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='BP', ...);print(tmpBP)
		if (!is.null(tmpBP)) {
			tmpBP = as.data.frame(cbind(ontology=rep('BP',nrow(tmpBP)), tmpBP));
		}
		cat('\n==========  CC  ==========\n');
		tmpCC = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='CC', ...);print(tmpCC)
		if (!is.null(tmpCC)) {
			tmpCC = as.data.frame(cbind(ontology=rep('CC',nrow(tmpCC)), tmpCC));
		}
		cat('\n==========  MF  ==========\n');
		tmpMF = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='MF', ...);print(tmpMF)
		if (!is.null(tmpMF)) {
			tmpMF = as.data.frame(cbind(ontology=rep('MF',nrow(tmpMF)), tmpMF));
		}
		tmpAll = as.data.frame(matrix(ncol=7, dimnames=list(1, c('ontology', 'goid', 'name', 'refnum', 'interestnum', 'pvalue', 'adjustp'))))
		for (res in c('tmpBP','tmpCC','tmpMF')) {
			if (!is.null(get(res))) {
				tmpAll = as.data.frame(rbind(tmpAll, get(res)))
			}
		}
		tmpAll = tmpAll[-1, ];
		return(tmpAll[order(tmpAll$adjustp, -tmpAll$interestnum), ]);
}

# wont get results if much stricter than this
.GOFunctionAllOntologies(genesForDavid_DEntrez, bgEntrez, fdrmethod='BH', fdrth=.5)

############################################################################################
############################################################################################
############################################################################################
### overlaps with TEs, miRNAs, lncRNAs...


te0 = read.table("Abur_final_TE.bed-closest", sep="\t", header=F);

# keep only hits within 5kb
te = te0[te0$V23<5000 & te0$V23>-5000, ];		# note different column

# get rid of lines where there was no hit
te = te[te$V20!=".",];				# note different column



# split into list where each element is a DMR
tes = split(te, paste(te[,1],":",te[,2],"-",te[,3],sep=""));

# get number of hits for each DMR
tesnrow = sapply(tes, nrow);

# some DMRs have >1 te hit, all are unique
tes[tesnrow>1]



# add column for overlapping TE hits
.addTEsToDMRs = function (DMR, overlapResults) {
	gS = split(overlapResults, paste(overlapResults[,1],":",overlapResults[,2],"-",overlapResults[,3],sep=""));
	DMRtes = as.data.frame(cbind(DMR, te.overlap=rep(NA,nrow(DMR)), te.dist=rep(NA,nrow(DMR))));
	for (row in 1:nrow(DMRtes)) {
		thisDMR = paste(DMRtes[row,1],":",DMRtes[row,2],"-",DMRtes[row,3],sep="");
		rownames(DMRtes)[row] = thisDMR;
		hit = gS[match(thisDMR, names(gS))][[1]];
		if (is.null(hit)) {
			next;
		} else if (nrow(hit)==1) {
			DMRtes$te.overlap[row] = hit[,20];
			DMRtes$te.dist[row] = hit[,ncol(hit)];
		} else {
			DMRtes$te.overlap[row] = paste(hit[,20], collapse=',');
			DMRtes$te.dist[row] = paste(hit[,ncol(hit)], collapse=',');
		}
	}
	return(DMRtes);
}

DMRgenesTEs = .addTEsToDMRs(DMRgenes, te);


############################################################################################
##### miRNAs

mir0 = read.table("abur_miRNAs-130326.fix.bed-closest", sep="\t", header=F);

# keep only hits within 5kb
mir = mir0[mir0$V26<5000 & mir0$V26>-5000, ];

# get rid of lines where there was no hit
mir = mir[mir$V17!=".",];	

#############################################################################################
##### lncRNAs

lnc0 = read.table("abur.lnc.final.gtf-closest", sep="\t", header=F);

# keep only hits within 5kb
lnc = lnc0[lnc0$V26<5000 & lnc0$V26>-5000, ];

# get rid of lines where there was no hit
lnc = lnc[lnc$V25!=".",];

# split into list where each element is a DMR
lncs = split(lnc, paste(lnc[,1],":",lnc[,2],"-",lnc[,3],sep=""));


###################################################################################################
#############################################################################################
#########################################################################################
.runFisherTest = function (DMR, col1, col2, testFunc='fisher.test') {
	cats1 = names(table(unlist(strsplit(DMR[,col1], ','))));#print(cats1)
	cats2 = names(table(unlist(strsplit(DMR[,col2], ','))));#print(cats2)
	mat = matrix(nrow=length(cats1), ncol=length(cats2), dimnames=list(cats1, cats2))
	for (i1 in 1:length(cats1)) {
		for (i2 in 1:length(cats2)) {
			#print(cats1[i1]); print(cats2[i2])
			mat[i1, i2] = sum(grepl(cats2[i2], DMR[DMR[,col1]==cats1[i1], col2]));
		}
	}
	return(list(mat, eval(call(testFunc, mat))));
}


.runFisherTest(DMRgenesTEs, 16, 19, 'chisq.test');




#########################################################################################
##### expression data from other animals
load('/Volumes/fishstudies-1/_Elim_RNAseq/4_expression/NCBI/htseq-count/rpkm_.25filt_na0_p20_mm100_mbs12000_ComBat3X/run3DATA.RData');

load('/Volumes/fishstudies-1/_Elim_RNAseq/4_expression/NCBI/htseq-count/rpkm_.25filt_na0_p20_mm100_mbs12000_ComBat3X/run3NET.RData');



tmp=unique(.getDMRsWithGenes(DMRgenesHs, .getTermKeywordGenes(DMRgenesHsDavid_D, 'fibronectin'))$gene.overlap)
tmp=tmp[tmp%in%names(DATA)]
.getkMERankAndQuantileForGenes(tmp,modkMEs)





