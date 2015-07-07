setwd('/Volumes/fishstudies/_methylation/overlap_all');
rm(list=ls());

## number of lines in original methylation files
# 3157.wig.bed - 36868745
# 3165.wig.bed - 39926602
# 3581.wig.bed - 23365799
# 3677.wig.bed - 26630403

uniq0 = read.table('_uniq_ids_line_nums');
uniq0 = split(uniq0, substring(uniq0[, 2], 1, 4));
uniq = as.data.frame(matrix(nrow=nrow(uniq0[[1]]), ncol=length(uniq0)));
names(uniq) = names(uniq0);
rownames(uniq) = gsub('_uniq_ids|.hit|noCDS', '', gsub(paste(names(uniq0)[1], '.wig.bed-', sep=''), '', uniq0[[1]][, 2]));
for (s in 1:length(uniq0)) {
	uniq[, s] = uniq0[[s]][, 1];
}; rm(s);
rownames(uniq)[c(2,3,5,8)] = c('exon', 'intron', 'SNP_indiv', 'DE');

# > uniq
                    # 3157     3165    3581    3677
# CNVs               76020    75925   53578   58766
# exon             6239416  6533972 4338071 4685241
# intron          11883118 13049736 7217167 8396147
# SNPs_assembly      42581    46255   27443   31273
# SNP_indiv          13313    13608    8783    7641
# TEs              8387834  9034935 5441355 6154101
# UTRs             1991833  2155560 1284434 1451223
# DE                 43572    46048   27105   30523
# lncRNAs            61927    64708   38848   44104
# miRNAs              1094     1145     611     699
# microsatellites     9059    11171    7101    7499


total = matrix(c(3157, 3165, 3581, 3677, 36868745, 39926602, 23365799, 26630403), ncol=2);
uniq.pc = uniq;
for (col in 1:ncol(uniq)) {
	uniq.pc[, col] = uniq[, col] / total[total[,1]==names(uniq.pc)[col], 2] * 100;
}; rm(col);

nd = apply(uniq.pc[, c(1,4)], 1, mean);
d = apply(uniq.pc[, c(2,3)], 1, mean);

# > uniq.pc
                        # 3157         3165         3581         3677
# CNVs             0.206190908  0.190161437  0.229300954  0.220672590
# exon            16.923320824 16.364958881 18.565900528 17.593579038
# intron          32.230871976 32.684314082 30.887738956 31.528426363
# SNPs_assembly    0.115493489  0.115850079  0.117449440  0.117433446
# SNP_indiv        0.036109176  0.034082540  0.037589128  0.028692769
# TEs             22.750527581 22.628860327 23.287690697 23.109304805
# UTRs             5.402497427  5.398806540  5.497068600  5.449496953
# DE               0.118181403  0.115331628  0.116002881  0.114617116
# lncRNAs          0.167966119  0.162067386  0.166260097  0.165615218
# miRNAs           0.002967283  0.002867762  0.002614933  0.002624819
# microsatellites  0.024570948  0.027978840  0.030390572  0.028159544

uniq.pc2 = rbind(uniq.pc[c(2,3,6,7),], apply(uniq.pc[c(1,4,5,8:11),],2,sum));
rownames(uniq.pc2)[5] = 'other';

pie(sort(apply(uniq.pc2,1,mean)), main='% of methylated bases');  

dotchart(sort(apply(uniq.pc2,1,mean)), xlab='% of methylated bases', ylab='sequence feature',pch=19,cex=1.5,main='methylome(n=4)');
abline(v=sort(apply(uniq.pc2,1,mean)), col='red', lty='dashed');

dotchart(sort(apply(uniq.pc[c(1,4,5,8:11),],1,mean)),pch=19,cex=1.5,xlab='% of methylated bases',main='"other" features');
abline(v=sort(apply(uniq.pc[c(1,4,5,8:11),],1,mean)), col='red',lty='dashed');
##########

uniqOnly = list();
for (s in 1:length(names(uniq))) {
	print (names(uniq)[s])
	files = list.files()[grepl(paste(names(uniq)[s], '.*Only$', sep=''), list.files())];
	tmp = list();
	for (f in 1:length(files)) {
		print(files[f])
		tmp[[f]] = read.table(files[f])[, 1];
		names(tmp)[f] = gsub('.hit_uniq_ids_Only', '', gsub(paste(names(uniq)[s], '.wig.bed-', sep=''), '', files[f]));
	}
	uniqOnly[[s]] = tmp;
	names(uniqOnly)[s] = names(uniq)[s];
}; rm(s,tmp,f);


##########

files = list.files()[grepl('groupby',list.files())];

readOverlap_groupby = function(files, pattern, colClasses, col.names, len.name=4) {
	these_files = files[grepl(pattern, files)];
	out = list();
	for (f in 1:length(these_files)) {
		print(these_files[f]);
		out[[f]] = read.table(these_files[f], header=F, sep='\t', colClasses=colClasses, col.names=col.names);
	}
	names(out) = substring(these_files, 1, len.name);
	return(out);
}

exon_lvls = readOverlap_groupby(files,'noCDS',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);

meth_lvls = sapply(exon_lvls, function(f) mean(f[,1]));

intron_lvls = readOverlap_groupby(files,'intron',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);

meth_lvls = rbind(meth_lvls, sapply(intron_lvls, function(f) mean(f[,1])));
rownames(meth_lvls) = c('exon','intron');

TE_lvls = readOverlap_groupby(files,'TEs',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);

meth_lvls = rbind(meth_lvls, sapply(TE_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[3] = 'TEs';

CNV_lvls = readOverlap_groupby(files,'CNVs',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(CNV_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[4] = 'CNVs';

SNPs_assembly_lvls = readOverlap_groupby(files,'SNPs_assembly',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(SNPs_assembly_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[5] = 'SNPs_assembly';

UTR_lvls = readOverlap_groupby(files,'UTRs',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(UTR_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[6] = 'UTRs';

DE_lvls = readOverlap_groupby(files,'fa.hit',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(DE_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[7] = 'DE';

lncRNA_lvls = readOverlap_groupby(files,'lnc',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(lncRNA_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[8] = 'lncRNAs';

miRNA_lvls = readOverlap_groupby(files,'miRNA',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(miRNA_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[9] = 'miRNAs';

microsatellite_lvls = readOverlap_groupby(files,'microsatellite',
								colClasses=c(rep('NULL',7),'numeric',rep('NULL',3)),
								col.names=c(rep('NULL',7),'mean',rep('NULL',3)),
								len.name=4);
meth_lvls = rbind(meth_lvls, sapply(microsatellite_lvls, function(f) mean(f[,1])));
rownames(meth_lvls)[10] = 'microsatellites';

SNP_indiv_lvls = readOverlap_groupby(list.files('..'),'BQSR.hit$',
									 colClasses=c('character','numeric','numeric','character','numeric',rep('NULL',11)),
									 col.names=c('scaffold','start','end','id','mean',rep('NULL',11)),
									 len.name=4);
meth_lvls = rbind(meth_lvls, sapply(SNP_indiv_lvls, function(f) mean(f[,5])));
rownames(meth_lvls)[11] = 'SNP_indiv';


meth_lvls = meth_lvls[match(rownames(uniq.pc),rownames(meth_lvls)),];


# > meth_lvls
                      # 3157       3165       3581       3677
# CNVs            0.11741437 0.11514924 0.12080027 0.11883154
# exon            0.14069308 0.13749731 0.14722133 0.14567431
# intron          0.11176529 0.10728899 0.11831593 0.11558834
# SNPs_assembly   0.24472389 0.24429212 0.25736218 0.25246147
# SNP_indiv       0.18593871 0.19241885 0.20043773 0.19997751
# TEs             0.13402223 0.13050238 0.14268284 0.13983045
# UTRs            0.09604287 0.09241142 0.10056994 0.09911710
# DE              0.11172049 0.10634481 0.10925506 0.10814305
# lncRNAs         0.12057155 0.11942250 0.12753255 0.12797858
# miRNAs          0.09674384 0.07602358 0.08437568 0.08680473
# microsatellites 0.25465244 0.22327023 0.24065477 0.24453814


meth_lvls2 = t(meth_lvls);
meth_lvls2 = meth_lvls2[, order(colnames(meth_lvls2))];
##########################












##############################




files = list.files()[grepl('groupby',list.files())];

# these_files = files[grepl('noCDS',files)];
# these_classes = c('character',rep('numeric',9),'NULL');
# these_names = c('scaffold','start','end','count','count_distinct','min','max','mean','median','stdev','NULL');
# exon = list();
# for (f in 1:length(files[grepl('noCDS',files)])) {
	# print(these_files[f]);
	# exon[[f]] = read.table(these_files[f], header=F, sep='\t', colClasses=these_classes, col.names=these_names);
# }
# names(exon) = substring(these_files, 1, 4);
# rm(f, these_files, these_classes, these_names);



exon = readOverlap_groupby(files=files, pattern='noCDS', 	
						   colClasses=c('character',rep('numeric',9),'NULL'),
						   col.names=c('scaffold','start','end','count','count_distinct','min','max','mean','median','stdev','NULL'),
						   len.name=4
						   );

exon.hits = readOverlap_groupby(files=files, pattern='noCDS', 	
						        colClasses=c(rep('NULL',10),'character'),
						        col.names=c(rep('NULL',10), 'sites'),
						        len.name=4
						        );



























# can work on smaller bedtools intersect output files with no pre-filtering

groupLookup = data.frame(num=c('3157', '3165', '3581', '3677'),
						 group=c('ND', 'D', 'D', 'ND')
						 );

## DE from cuffdiff ##
de = list();
files = list.files()[grepl('fa.hit$', list.files())];
for (f in 1:length(files)) {
	print(files[f]);
	de[[f]] = read.table(files[f], header=F, sep='\t', stringsAsFactors=F);
	names(de)[f] = files[f];
}; rm(f);
names(de) = gsub('.wig.bed-cuffdiff_postRealign_bqsr_fa.hit', '', names(de));
names(de) = paste(names(de), groupLookup$group[match(names(de), groupLookup$num,)], sep='-');

numMethBases = sapply(de, nrow);

deSp = list();
for (f in 1:length(de)) {
	print(names(de)[f]);
	deSp[[f]] = split(de[[f]], as.factor(paste(de[[f]]$V6, de[[f]]$V7, de[[f]]$V8, sep=':')));
	names(deSp)[f] = names(de)[f];
}; rm(f);

numMethBasesByGene = lapply(deSp, function(f) sapply(f, nrow));
meanMethLevelsByGene = sapply(deSp, sapply, function(f) mean(f$V5));

fxn.pmb = function(x) {
	out = c();
	for (g in 1:length(x)) {
		tmp = strsplit(names(x)[g], ':')[[1]];
		l = as.numeric(tmp[3]) - as.numeric(tmp[2]);
		pmb = as.numeric(x[g]) / l;
		out = c(out, pmb);
	}
	names(out) = names(x);
	return(out);
}
pctMethBasesByGene = lapply(numMethBasesByGene, fxn.pmb);

library(WGCNA);

par(mfrow=c(2,2));
for (s in 1:length(de)) {
	verboseScatterplot(numMethBasesByGene[[s]], meanMethLevelsByGene[[s]],
					   abline=T, frame.plot=F, abline.color='red', abline.lty='dashed',
					   xlab='# of bases methylated', ylab='mean methylation level',
					   main=names(numMethBasesByGene)[s]
					   );
}; rm(s);

par(mfrow=c(2,2));
for (s in 1:length(de)) {
	verboseScatterplot(pctMethBasesByGene[[s]], meanMethLevelsByGene[[s]],
					   abline=T, frame.plot=F, abline.color='red', abline.lty='dashed',
					   xlab='% of bases methylated', ylab='mean methylation level',
					   main=names(pctMethBasesByGene)[s]
					   );
}; rm(s);

##   ##


##   ##

