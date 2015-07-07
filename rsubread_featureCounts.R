#####RUN 1

library(Rsubread);
library(limma);
setwd('~/Documents/Rsubread');
files = list.files('../_LYNLEY_RNAseq/')[grepl('sam$', list.files('../_LYNLEY_RNAseq/'))];
files = paste('../_LYNLEY_RNAseq/', files, sep='');
annot.ext = '../Astatotilapia_burtoni.BROADAB2.gtf';

counts = featureCounts(files = files,
					  annot.ext = annot.ext,
					  isGTFAnnotationFile = T,
					  GTF.featureType = 'exon',
					  GTF.attrType = 'transcript_id',
					  useMetaFeatures = T,
					  allowMultiOverlap = F,
					  nthreads = 4,
					  strandSpecific = 1,
					  countMultiMappingReads = F,
					  isPairedEnd = T,
					  reportReads = T
					  );
					  
save(counts, file='counts.RData');

colnames(counts$counts) = gsub('_Rsubread.sam', '', 
							   gsub('..._LYNLEY_RNAseq.130913_LYNLEY_0364_AC2HRPACXX_L6_', '', colnames(counts$counts))
							   );
rpkm = apply(counts$counts, 2, function(x){x*(1000/counts$annotation$Length)*(1e6/sum(x))});
expr = rowSums(rpkm >= .5) >= 4;
x = counts$counts[expr, ];
					  
status = factor(c('ND', 'D', 'ND', 'D'));		# assumes order is ATCACG, CGATGT, TGACCA, TTAGGC
design = model.matrix(~0+status);
colnames(design) = levels(status);

y = voom(x, design, plot=T);
plotMDS(y);

fit = lmFit(y, design);
contr = makeContrasts(DvsND=D-ND, levels=design);
fit.contr = eBayes(contrasts.fit(fit, contr));
dt = decideTests(fit.contr);

options(digits=3);
de = rownames(topTable(fit.contr,number=50));
#de.genename = substr(gsub('mrna','gene',de), 1, nchar(gsub('mrna','gene',de))-2);

annos = read.table('/Volumes/handsfiler$/FishStudies/_Burtoni_annotations/geneNamesTree_AB', sep='\t', header=T, na.strings='',fill=T, colClasses='character');


######################################################
splitOn = substr(gsub('mrna','gene',rownames(counts$counts)), 1, nchar(gsub('mrna','gene',rownames(counts$counts)))-2);
counts.gene = split(counts$counts, splitOn);
counts.gene = lapply(counts.gene, 
					 function(blah) matrix(blah, ncol=4, dimnames=list(rep('',length(blah)/4),colnames(counts$counts)))
					 );
counts.gene = lapply(counts.gene, colSums);
counts.geneMat = matrix(nrow=length(counts.gene), ncol=4, dimnames=list(names(counts.gene), names(counts.gene[[1]])));
for (g in 1:length(counts.gene)) {
	counts.geneMat[g, ] = counts.gene[[g]];
}; rm(g);

y.gene = voom(counts.geneMat, design, plot=T);
plotMDS(y.gene);
fit.gene = lmFit(y.gene, design);
fit.gene.contr = eBayes(contrasts.fit(fit.gene, contr));
dt.gene = decideTests(fit.gene.contr);
de.gene = rownames(topTable(fit.gene.contr,number=50));				  
					  
        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 4 SAM files                                      ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid16370               ||
# ||             Annotations : ../Astatotilapia_burtoni.BROADAB2.gtf (GTF)      ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : yes                                              ||
# ||         Strand specific : yes                                              ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# ||          Chimeric reads : counted                                          ||
# ||        Both ends mapped : not required                                     ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file ../Astatotilapia_burtoni.BROADAB2.gtf ...             ||
# ||    Number of features is 607912                                            ||
# ||    Number of meta-features is 52845                                        ||
# ||    Number of chromosomes is 2662                                           ||
# ||                                                                            ||
# || Process SAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_AT ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Total number of fragments is : 38738845                                 ||
# ||    Number of successfully assigned fragments is : 547180 (1.4%)            ||
# ||    Running time : 8.72 minutes                                             ||
# ||                                                                            ||
# || Process SAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_CG ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Total number of fragments is : 41726276                                 ||
# ||    Number of successfully assigned fragments is : 564813 (1.4%)            ||
# ||    Running time : 6.74 minutes                                             ||
# ||                                                                            ||
# || Process SAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TG ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Total number of fragments is : 38469019                                 ||
# ||    Number of successfully assigned fragments is : 882783 (2.3%)            ||
# ||    Running time : 4.60 minutes                                             ||
# ||                                                                            ||
# || Process SAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TT ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Total number of fragments is : 37594748                                 ||
# ||    Number of successfully assigned fragments is : 549609 (1.5%)            ||
# ||    Running time : 6.86 minutes                                             ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//


#####
##TRY WITH DUPLICATE REDAS REMOVED



rm(list=ls());
library(Rsubread);
library(limma);
setwd('~/Documents/Rsubread');
dir = '../_LYNLEY_RNAseq/'
files = list.files(dir)[grepl('markDuplicates.bam$', list.files(dir))];
files = paste(dir, files, sep='');
annot.ext = '../Astatotilapia_burtoni.BROADAB2.gtf';

counts = featureCounts(files = files,
					  annot.ext = annot.ext,
					  isGTFAnnotationFile = T,
					  GTF.featureType = 'exon',
					  GTF.attrType = 'transcript_id',
					  useMetaFeatures = T,
					  allowMultiOverlap = F,
					  nthreads = 4,
					  strandSpecific = 1,
					  countMultiMappingReads = F,
					  isPairedEnd = T,
					  reportReads = T
					  );
					  
					  
					  
					  
        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
# ||             Input files : 4 BAM files                                      ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                           o ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HR ... ||
# ||                                                                            ||
# ||             Output file : ./.Rsubread_featureCounts_pid16370               ||
# ||             Annotations : ../Astatotilapia_burtoni.BROADAB2.gtf (GTF)      ||
# ||      Assignment details : <input_file>.featureCounts                       ||
# ||                                                                            ||
# ||                 Threads : 4                                                ||
# ||                   Level : meta-feature level                               ||
# ||              Paired-end : yes                                              ||
# ||         Strand specific : yes                                              ||
# ||      Multimapping reads : not counted                                      ||
# || Multi-overlapping reads : not counted                                      ||
# ||                                                                            ||
# ||          Chimeric reads : counted                                          ||
# ||        Both ends mapped : not required                                     ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Load annotation file ../Astatotilapia_burtoni.BROADAB2.gtf ...             ||
# ||    Number of features is 607912                                            ||
# ||    Number of meta-features is 52845                                        ||
# ||    Number of chromosomes is 2662                                           ||
# ||                                                                            ||
# || Process BAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_AT ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Found reads that are not properly paired.                               ||
# ||    (missing mate or the mate is not the next read)                         ||
# ||    130948 reads have missing mates.                                        ||
# ||    Input was converted to a format accepted by featureCounts.              ||
# ||    Total number of fragments is : 25324816                                 ||
# ||    Number of successfully assigned fragments is : 370960 (1.5%)            ||
# ||    Running time : 15.93 minutes                                            ||
# ||                                                                            ||
# || Process BAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_CG ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Found reads that are not properly paired.                               ||
# ||    (missing mate or the mate is not the next read)                         ||
# ||    382610 reads have missing mates.                                        ||
# ||    Input was converted to a format accepted by featureCounts.              ||
# ||    Total number of fragments is : 29451971                                 ||
# ||    Number of successfully assigned fragments is : 412369 (1.4%)            ||
# ||    Running time : 19.31 minutes                                            ||
# ||                                                                            ||
# || Process BAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TG ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Found reads that are not properly paired.                               ||
# ||    (missing mate or the mate is not the next read)                         ||
# ||    229246 reads have missing mates.                                        ||
# ||    Input was converted to a format accepted by featureCounts.              ||
# ||    Total number of fragments is : 23689881                                 ||
# ||    Number of successfully assigned fragments is : 604466 (2.6%)            ||
# ||    Running time : 15.86 minutes                                            ||
# ||                                                                            ||
# || Process BAM file ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TT ... ||
# ||    Assign fragments (read pairs) to features...                            ||
# ||    Each fragment is counted once.                                          ||
# ||    Found reads that are not properly paired.                               ||
# ||    (missing mate or the mate is not the next read)                         ||
# ||    215661 reads have missing mates.                                        ||
# ||    Input was converted to a format accepted by featureCounts.              ||
# ||    Total number of fragments is : 26498951                                 ||
# ||    Number of successfully assigned fragments is : 414217 (1.6%)            ||
# ||    Running time : 16.63 minutes                                            ||
# ||                                                                            ||
# ||                         Read assignment finished.                          ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//
