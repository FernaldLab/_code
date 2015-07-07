rm(list=ls());
setwd("~/Desktop/_mammalian_RNAseq/Alignments");
library(Rsubread);
basename='Bonobo';
firstreads=list.files()[grepl('BONOBO', list.files())];

rm(list=ls());
setwd("~/Desktop/_mammalian_RNAseq/Alignments");
library(Rsubread);
basename='Human';
firstreads=list.files()[grepl('HUMAN', list.files())];


rm(list=ls());
setwd("~/Desktop/_mammalian_RNAseq/Alignments");
library(Rsubread);
basename='Human';
firstreads=list.files()[grepl('HUMAN', list.files())];

species=c("Chicken", "Chimp", "Gorilla", "Macaque", "Mouse", "Opossum", "Orangutan", "Platypus")
for (s in species){
	firstreads=list.files()[grepl(toupper(s), list.files())];
	align(index=s,
	 readfile1=firstreads, readfile2=NULL, 
	 input_format='FASTQ', output_format='SAM', 
	 output_file=gsub('fastq', 'Rsubread.sam', firstreads),
	 nsubreads=10, TH1=3, TH2=1, nthreads=4, indels=5, phredOffset=33,
	 tieBreakQS=F, tieBreakHamming=T, unique=T, nBestLocations=1, 
	 minFragLength=50, maxFragLength=600, 
	 nTrim5=0, nTrim3=0, readGroupID=NULL, readGroup=NULL, color2base=F,
	 DP_GapOpenPenalty=-1, DP_GapExtPenalty=0, DP_MismatchPenalty=0, DP_MatchScore=2,
	 reportFusions=F
	 );
}
rm(s)


align(index=basename,
	 readfile1=firstreads, readfile2=NULL, 
	 input_format='FASTQ', output_format='SAM', 
	 output_file=gsub('fastq', 'Rsubread.sam', firstreads),
	 nsubreads=10, TH1=3, TH2=1, nthreads=4, indels=5, phredOffset=33,
	 tieBreakQS=F, tieBreakHamming=T, unique=T, nBestLocations=1, 
	 minFragLength=50, maxFragLength=600, 
	 nTrim5=0, nTrim3=0, readGroupID=NULL, readGroup=NULL, color2base=F,
	 DP_GapOpenPenalty=-1, DP_GapExtPenalty=0, DP_MismatchPenalty=0, DP_MatchScore=2,
	 reportFusions=F
	 );


#########
files = list.files()[grepl('sam$', list.files())];	 
for (f in 1:length(files)) {
	print(files[f])
	species = strsplit(files[f],'SRR')[[1]][1];
	if (species=='BONOBO') {
		gtffile = 'CHIMP.gtf'
	} else {
		gtffile = paste(species, '.gtf', sep='');
	}
	counts=featureCounts(files[f],
						annot.ext=gtffile, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id",
						useMetaFeatures=TRUE, allowMultiOverlap=FALSE, nthreads=4, strandSpecific=0, countMultiMappingReads=FALSE, 
						minMQS=0, isPairedEnd=FALSE, requireBothEndsMapped=FALSE, 
						checkFragLength=FALSE, minFragLength=50, maxFragLength=600,
						countChimericFragments=TRUE, chrAliases=NULL, reportReads=TRUE
						);
	out = paste(files[f], '.counts.RData', sep='');
	save(counts, file=out);
}; rm(f, species, gtffile, out);

############

rm(list=ls());
files = list.files()[grepl('RData$', list.files())];
all = list();
for (f in 1:length(files)) {
	print(files[f]);
	load(files[f]);
	r = apply(counts$counts, 2, function(x){x*(1000/counts$annotation$Length)*(1e6/sum(x))});
	all[[f]] = cbind(counts$counts, r);
	colnames(all[[f]]) = c('counts', 'rpkm');
	names(all)[f] = strsplit(files[f], '\\.')[[1]][1];
	rm(counts);
}; rm(f, r);

all = lapply(all, function(x) x[x[, 1]>0, ]);
all = all[unlist(lapply(all, nrow))>0];
all.split = split(all, as.factor(substr(names(all), 1, attributes(regexpr('^[A-Z]*', names(all)))[[1]]-3)));

ids = list();
for (sp in 1:length(all.split)) {
	x = all.split[[sp]];
	y = rownames(x[[1]]);
	if (length(x)==1) {
		ids[[sp]] = y;
	} else {
		for (s in 2:length(x)) {
			y = y[y %in% rownames(x[[s]])];
		}
		ids[[sp]] = y;
	}
	names(ids)[sp] = names(all.split)[sp];
}; rm(sp, x, y, s);

for (s in 1:length(all)) {
	ind = match(substr(names(all)[s], 1, attributes(regexpr('^[A-Z]*', names(all)[s]))[[1]]-3), names(ids));
	all[[s]] = all[[s]][rownames(all[[s]]) %in% ids[[ind]], ];
}; rm(s, ind);

library(biomaRt);
#species = 'ggorilla|ggallus|pabelii|mmusculus|mmulatta|hsapiens|mdomestica|oanatinus|ppaniscus';
mart = useMart('ensembl', 'ptroglodytes_gene_ensembl');
a = listAttributes(mart);

#attr = a[grepl(species, a[,1]) & grepl('homolog', a[,1]) & (grepl('\\%', a[,2]) | grepl('gene', a[,1])), ];
attr = a[grepl('hsapiens', a[,1]) & grepl('homolog', a[,1]) & (grepl('\\%', a[,2]) | grepl('gene', a[,1])), ];

# CHIMP has fewest genes 
###DONT DO THIS WAY JUST START WITH HUMAN TO AVOID CHIMP IDS MAPPING TO SAME HUMAN ID
genes = data.frame(CHIMP=ids$CHIMP, HUMAN=rep('', length(ids$CHIMP)), stringsAsFactors=F);
for (g in 1:nrow(genes)) {
	if (g %% 100 == 0) { print(g) }
	
	tmp = getBM(attributes=attr[1:2, 1], filters='ensembl_gene_id', values=genes[g, 1], mart);
	#print(is.na(tmp[1, 2]));
	if (nrow(tmp)==1) {
		if (tmp[1, 1]=='' | is.na(tmp[1, 2])) {
			genes[g, 2] = 'NA'; 
		} else {
			genes[g, 2] = as.character(tmp[1, 1]);
		}
	} else {
		if (length(unique(tmp[, 1]))==1) {
			genes[g, 2] = tmp[1, 1];
		} else {
			tmp = tmp[!(is.na(tmp[,2])),]
			maxid = max(tmp[, 2], na.rm=T);
			genes[g, 2] = tmp[which(tmp[, 2]==maxid)[1], 1];
		}
	}

}; rm(g, tmp, maxid);

genes = genes[genes$HUMAN!='NA', ]
human_genes = genes[, 2];

mart = useMart('ensembl', 'hsapiens_gene_ensembl');
a = listAttributes(mart);
#species = 'ggorilla|ggallus|pabelii|mmusculus|mmulatta|mdomestica|oanatinus';
species = 'ggorilla|pabelii|mmusculus|mmulatta|mdomestica|oanatinus';
attr = a[grepl(species, a[,1]) & grepl('homolog', a[,1]) & (grepl('\\%', a[,2]) | grepl('gene', a[,1])), ];

genes0 = genes;


genes = data.frame(genes0, GORILLA=rep('', nrow(genes0)), MACAQUE=rep('', nrow(genes0)), MOUSE=rep('', nrow(genes0)), OPOSSUM=rep('', nrow(genes0)), ORANGUTAN=rep('', nrow(genes0)), PLATYPUS=rep('', nrow(genes0)), stringsAsFactors=F);

#spVec = c('CHICKEN', 'GORILLA', 'MACAQUE', 'MOUSE', 'OPOSSUM', 'ORANGUTAN', 'PLATYPUS');
spVec = c('GORILLA', 'MACAQUE', 'MOUSE', 'OPOSSUM', 'ORANGUTAN', 'PLATYPUS');
for (g in 1:nrow(genes)) {
	if (g %% 100 == 0) {print(g)};
 	this_gene = genes$HUMAN[g];
	tmp = getBM(attributes=attr[-seq(3, nrow(attr), 3), 1], filters='ensembl_gene_id', values=this_gene, mart);
	tmpsym = tmp[, seq(1, ncol(tmp), 2)];
	tmpid = tmp[, seq(2, ncol(tmp), 2)];
	if (nrow(tmp)==1) {
		genes[g, 3:8] = tmpsym;
	} else {
		for (col in 1:ncol(tmpsym)) {
			if (length(unique(tmpsym[, col]))==1) {
				genes[g, col+2] = unique(tmpsym[, col]);
			} else {
				maxid = max(tmpid[, col], na.rm=T);
				goodid = tmpsym[which(tmpid[, col]==maxid)[1], col];
				genes[g, col+2] = goodid;
			}
		}
	}
}; rm(g, this_gene, tmp, col, maxid, goodid);

rowNAs = apply(genes, 1, function(g) sum(is.na(g)));
genesFilt = genes[rowNAs==0, ];

all0 = all;
for (s in 1:length(all)) {
	this_sp = substr(names(all)[s], 1, attributes(regexpr('^[A-Z]*', names(all)[s]))[[1]]-3);
	toKeep = genesFilt[, match(this_sp, names(genesFilt))];
	all[[s]] = all[[s]][rownames(all[[s]]) %in% toKeep, ];
}; rm(s, this_sp, toKeep);

# now MOUSE has fewest
mgenes = rownames(all[[18]]);
genesFilt.mgenes = genesFilt[genesFilt$MOUSE %in% mgenes, ];

# BECAUSE OF CHIMP HUMAN MAPPING ISSUE
dups = genesFilt.mgenes$HUMAN[duplicated(genesFilt.mgenes$HUMAN)];
genesFinal = genesFilt.mgenes[-(which(genesFilt.mgenes$HUMAN %in% dups)), ];

allFinal = all;
for (s in 1:length(allFinal)) {
	this_sp = substr(names(allFinal)[s], 1, attributes(regexpr('^[A-Z]*', names(allFinal)[s]))[[1]]-3);
	toKeep = genesFinal[, match(this_sp, names(genesFinal))];
	allFinal[[s]] = allFinal[[s]][rownames(allFinal[[s]]) %in% toKeep, ];
}; rm(s, this_sp, toKeep);

for (s in 1:length(allFinal)) {
	this_sp = substr(names(allFinal)[s], 1, attributes(regexpr('^[A-Z]*', names(allFinal)[s]))[[1]]-3);
	col = match(this_sp, names(genesFinal));
	allFinal[[s]] = allFinal[[s]][match(genesFinal[, col], rownames(allFinal[[s]])), ];
}; rm(s, this_sp, col);

dat0 = as.data.frame(matrix(nrow=nrow(allFinal[[1]]), ncol=length(allFinal)));
for (s in 1:length(allFinal)) {
	dat0[, s] = allFinal[[s]][, 2];
	names(dat0)[s] = names(allFinal)[s];
}; rm(s);

for (row in 1:nrow(dat0)) {
	if (is.na(rownames(allFinal$HUMANSRR306838)[row])) {
		rownames(dat0)[row] = rownames(allFinal$CHIMPSRR306811)[row];
	} else {
		rownames(dat0)[row] = rownames(allFinal$HUMANSRR306838)[row];
	}
}; rm(row);
save(dat0, file='mammal_dat0.RData');

x=apply(dat0, 1, function(g) sum(is.na(g)));
dat0 = dat0[-(which(x>10)), ];

####
# # [1] "BONOBOSRR306826.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o BONOBOSRR306826.Rsubread.sam                   ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file BONOBOSRR306826.Rsubread.sam...                           ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 13764                                        ||
# ||    Number of successfully assigned reads is : 0 (0.0%)                     ||
# ||    Running time : 0.00 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "BONOBOSRR306827.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o BONOBOSRR306827.Rsubread.sam                   ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file BONOBOSRR306827.Rsubread.sam...                           ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 24777783                                     ||
# ||    Number of successfully assigned reads is : 0 (0.0%)                     ||
# ||    Running time : 2.10 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "BONOBOSRR306828.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o BONOBOSRR306828.Rsubread.sam                   ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file BONOBOSRR306828.Rsubread.sam...                           ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 38196822                                     ||
# ||    Number of successfully assigned reads is : 0 (0.0%)                     ||
# ||    Running time : 3.19 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHICKENSRR306710.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHICKENSRR306710.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHICKEN.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHICKEN.gtf ...                                       ||
# ||    Number of features is 181261                                            ||
# ||    Number of meta-features is 17108                                        ||
# ||    Number of chromosomes is 934                                            ||
# ||                                                                            ||
# || Process SAM file CHICKENSRR306710.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 19476908                                     ||
# ||    Number of successfully assigned reads is : 0 (0.0%)                     ||
# ||    Running time : 1.64 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHICKENSRR306711.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHICKENSRR306711.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHICKEN.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHICKEN.gtf ...                                       ||
# ||    Number of features is 181261                                            ||
# ||    Number of meta-features is 17108                                        ||
# ||    Number of chromosomes is 934                                            ||
# ||                                                                            ||
# || Process SAM file CHICKENSRR306711.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 17557038                                     ||
# ||    Number of successfully assigned reads is : 0 (0.0%)                     ||
# ||    Running time : 1.46 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306811.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306811.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306811.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 20083064                                     ||
# ||    Number of successfully assigned reads is : 1763603 (8.8%)               ||
# ||    Running time : 1.71 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306812.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306812.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306812.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 13947644                                     ||
# ||    Number of successfully assigned reads is : 1393787 (10.0%)              ||
# ||    Running time : 1.28 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306813.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306813.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306813.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 20408261                                     ||
# ||    Number of successfully assigned reads is : 2704055 (13.2%)              ||
# ||    Running time : 1.85 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306814.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306814.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306814.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 17394854                                     ||
# ||    Number of successfully assigned reads is : 2005366 (11.5%)              ||
# ||    Running time : 1.56 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306815.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306815.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306815.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 22234086                                     ||
# ||    Number of successfully assigned reads is : 1922631 (8.6%)               ||
# ||    Running time : 2.01 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "CHIMPSRR306816.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o CHIMPSRR306816.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : CHIMP.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file CHIMP.gtf ...                                         ||
# ||    Number of features is 209271                                            ||
# ||    Number of meta-features is 28012                                        ||
# ||    Number of chromosomes is 666                                            ||
# ||                                                                            ||
# || Process SAM file CHIMPSRR306816.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 23317655                                     ||
# ||    Number of successfully assigned reads is : 2188029 (9.4%)               ||
# ||    Running time : 2.09 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "GORILLASRR306800.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o GORILLASRR306800.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : GORILLA.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file GORILLA.gtf ...                                       ||
# ||    Number of features is 296000                                            ||
# ||    Number of meta-features is 29216                                        ||
# ||    Number of chromosomes is 27                                             ||
# ||                                                                            ||
# || Process SAM file GORILLASRR306800.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 35257547                                     ||
# ||    Number of successfully assigned reads is : 5835930 (16.6%)              ||
# ||    Running time : 3.06 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "GORILLASRR306801.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o GORILLASRR306801.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : GORILLA.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file GORILLA.gtf ...                                       ||
# ||    Number of features is 296000                                            ||
# ||    Number of meta-features is 29216                                        ||
# ||    Number of chromosomes is 27                                             ||
# ||                                                                            ||
# || Process SAM file GORILLASRR306801.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 16254814                                     ||
# ||    Number of successfully assigned reads is : 2389418 (14.7%)              ||
# ||    Running time : 1.48 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306838.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306838.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306838.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 24513415                                     ||
# ||    Number of successfully assigned reads is : 15816493 (64.5%)             ||
# ||    Running time : 2.20 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306839.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306839.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306839.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 18850030                                     ||
# ||    Number of successfully assigned reads is : 10018848 (53.2%)             ||
# ||    Running time : 1.65 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306840.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306840.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306840.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 22576705                                     ||
# ||    Number of successfully assigned reads is : 12575897 (55.7%)             ||
# ||    Running time : 2.19 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306841.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306841.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306841.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 24325223                                     ||
# ||    Number of successfully assigned reads is : 12686771 (52.2%)             ||
# ||    Running time : 2.15 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306842.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306842.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306842.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 17422994                                     ||
# ||    Number of successfully assigned reads is : 11055088 (63.5%)             ||
# ||    Running time : 1.74 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "HUMANSRR306843.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o HUMANSRR306843.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : HUMAN.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file HUMAN.gtf ...                                         ||
# ||    Number of features is 1306656                                           ||
# ||    Number of meta-features is 63677                                        ||
# ||    Number of chromosomes is 265                                            ||
# ||                                                                            ||
# || Process SAM file HUMANSRR306843.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 7913181                                      ||
# ||    Number of successfully assigned reads is : 3518788 (44.5%)              ||
# ||    Running time : 0.70 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MACAQUESRR306777.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MACAQUESRR306777.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MACAQUE.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MACAQUE.gtf ...                                       ||
# ||    Number of features is 366946                                            ||
# ||    Number of meta-features is 30246                                        ||
# ||    Number of chromosomes is 1333                                           ||
# ||                                                                            ||
# || Process SAM file MACAQUESRR306777.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 19068947                                     ||
# ||    Number of successfully assigned reads is : 8060458 (42.3%)              ||
# ||    Running time : 1.65 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MACAQUESRR306778.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MACAQUESRR306778.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MACAQUE.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MACAQUE.gtf ...                                       ||
# ||    Number of features is 366946                                            ||
# ||    Number of meta-features is 30246                                        ||
# ||    Number of chromosomes is 1333                                           ||
# ||                                                                            ||
# || Process SAM file MACAQUESRR306778.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 22554234                                     ||
# ||    Number of successfully assigned reads is : 7706289 (34.2%)              ||
# ||    Running time : 1.96 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MACAQUESRR306779.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MACAQUESRR306779.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MACAQUE.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MACAQUE.gtf ...                                       ||
# ||    Number of features is 366946                                            ||
# ||    Number of meta-features is 30246                                        ||
# ||    Number of chromosomes is 1333                                           ||
# ||                                                                            ||
# || Process SAM file MACAQUESRR306779.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 21461283                                     ||
# ||    Number of successfully assigned reads is : 12082286 (56.3%)             ||
# ||    Running time : 2.05 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306757.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306757.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306757.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 25094445                                     ||
# ||    Number of successfully assigned reads is : 2675177 (10.7%)              ||
# ||    Running time : 2.18 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306758.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306758.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306758.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 18882745                                     ||
# ||    Number of successfully assigned reads is : 1950994 (10.3%)              ||
# ||    Running time : 1.61 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306759.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306759.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306759.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 20757817                                     ||
# ||    Number of successfully assigned reads is : 1792596 (8.6%)               ||
# ||    Running time : 1.78 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306760.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306760.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306760.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 17770683                                     ||
# ||    Number of successfully assigned reads is : 1903754 (10.7%)              ||
# ||    Running time : 1.53 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306761.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306761.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306761.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 18759557                                     ||
# ||    Number of successfully assigned reads is : 1894905 (10.1%)              ||
# ||    Running time : 1.62 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "MOUSESRR306762.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o MOUSESRR306762.Rsubread.sam                    ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : MOUSE.gtf (GTF)                                  ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file MOUSE.gtf ...                                         ||
# ||    Number of features is 628052                                            ||
# ||    Number of meta-features is 39179                                        ||
# ||    Number of chromosomes is 58                                             ||
# ||                                                                            ||
# || Process SAM file MOUSESRR306762.Rsubread.sam...                            ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 19726026                                     ||
# ||    Number of successfully assigned reads is : 1935587 (9.8%)               ||
# ||    Running time : 1.69 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "OPOSSUMSRR306742.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o OPOSSUMSRR306742.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : OPOSSUM.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file OPOSSUM.gtf ...                                       ||
# ||    Number of features is 206888                                            ||
# ||    Number of meta-features is 23899                                        ||
# ||    Number of chromosomes is 11                                             ||
# ||                                                                            ||
# || Process SAM file OPOSSUMSRR306742.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 11325682                                     ||
# ||    Number of successfully assigned reads is : 3982653 (35.2%)              ||
# ||    Running time : 0.98 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "OPOSSUMSRR306743.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o OPOSSUMSRR306743.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : OPOSSUM.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file OPOSSUM.gtf ...                                       ||
# ||    Number of features is 206888                                            ||
# ||    Number of meta-features is 23899                                        ||
# ||    Number of chromosomes is 11                                             ||
# ||                                                                            ||
# || Process SAM file OPOSSUMSRR306743.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 47574556                                     ||
# ||    Number of successfully assigned reads is : 15515441 (32.6%)             ||
# ||    Running time : 4.14 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "OPOSSUMSRR306744.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o OPOSSUMSRR306744.Rsubread.sam                  ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : OPOSSUM.gtf (GTF)                                ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file OPOSSUM.gtf ...                                       ||
# ||    Number of features is 206888                                            ||
# ||    Number of meta-features is 23899                                        ||
# ||    Number of chromosomes is 11                                             ||
# ||                                                                            ||
# || Process SAM file OPOSSUMSRR306744.Rsubread.sam...                          ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 22273667                                     ||
# ||    Number of successfully assigned reads is : 7873898 (35.4%)              ||
# ||    Running time : 1.96 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "ORANGUTANSRR306791.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o ORANGUTANSRR306791.Rsubread.sam                ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : ORANGUTAN.gtf (GTF)                              ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file ORANGUTAN.gtf ...                                     ||
# ||    Number of features is 223227                                            ||
# ||    Number of meta-features is 28443                                        ||
# ||    Number of chromosomes is 54                                             ||
# ||                                                                            ||
# || Process SAM file ORANGUTANSRR306791.Rsubread.sam...                        ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 36457958                                     ||
# ||    Number of successfully assigned reads is : 18935811 (51.9%)             ||
# ||    Running time : 3.17 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "ORANGUTANSRR306792.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o ORANGUTANSRR306792.Rsubread.sam                ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : ORANGUTAN.gtf (GTF)                              ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file ORANGUTAN.gtf ...                                     ||
# ||    Number of features is 223227                                            ||
# ||    Number of meta-features is 28443                                        ||
# ||    Number of chromosomes is 54                                             ||
# ||                                                                            ||
# || Process SAM file ORANGUTANSRR306792.Rsubread.sam...                        ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 17675725                                     ||
# ||    Number of successfully assigned reads is : 9766128 (55.3%)              ||
# ||    Running time : 1.70 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "PLATYPUSSRR306724.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o PLATYPUSSRR306724.Rsubread.sam                 ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : PLATYPUS.gtf (GTF)                               ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file PLATYPUS.gtf ...                                      ||
# ||    Number of features is 192553                                            ||
# ||    Number of meta-features is 26116                                        ||
# ||    Number of chromosomes is 13771                                          ||
# ||                                                                            ||
# || Process SAM file PLATYPUSSRR306724.Rsubread.sam...                         ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 10797205                                     ||
# ||    Number of successfully assigned reads is : 3436001 (31.8%)              ||
# ||    Running time : 0.97 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# [1] "PLATYPUSSRR306725.Rsubread.sam"

        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 1 SAM file                                       ||
# ||                           o PLATYPUSSRR306725.Rsubread.sam                 ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid36225               ||
# ||             Annotations : PLATYPUS.gtf (GTF)                               ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : no                                               ||
# ||         Strand specific : no                                               ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file PLATYPUS.gtf ...                                      ||
# ||    Number of features is 192553                                            ||
# ||    Number of meta-features is 26116                                        ||
# ||    Number of chromosomes is 13771                                          ||
# ||                                                                            ||
# || Process SAM file PLATYPUSSRR306725.Rsubread.sam...                         ||
# ||    Assign reads to features...                                             ||
# ||    Total number of reads is : 0                                            ||
# ||    Number of successfully assigned reads is : 0                            ||
# ||    Running time : 0.00 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

