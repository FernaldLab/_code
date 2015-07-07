setwd('/Volumes/fishstudies/_methylation/overlap_all');
rm(list=ls());

uniq0 = read.table('TEtype_uniq_ids_line_nums');
uniq0 = split(uniq0, substring(uniq0[, 2], 1, 4));
uniq = as.data.frame(matrix(nrow=nrow(uniq0[[1]]), ncol=length(uniq0)));
names(uniq) = names(uniq0);

rownames(uniq) = gsub('TEs.hit.|_uniq_ids_Only', '', gsub(paste(names(uniq0)[1], '.wig.bed-', sep=''), '', uniq0[[1]][, 2]));
for (s in 1:length(uniq0)) {
	uniq[, s] = uniq0[[s]][, 1];
}; rm(s);

allnums = read.table('_uniq_ids_line_nums');
allnums = allnums[grepl('TEs',allnums[,2]), ];

uniq.pc = uniq;
for (col in 1:ncol(uniq)) {
	uniq.pc[, col] = uniq.pc[, col] / allnums[match(names(uniq.pc)[col], substring(allnums[,2],1,4)), 1] * 100;
	#print(names(uniq)[col]);
	#print(allnums[match(names(uniq)[col], substring(allnums[,2],1,4)), 1])
}; rm(col);

uniq.pcmean = apply(uniq.pc,1,mean);

pie(sort(uniq.pcmean), main='% of methylated bases in TEs');

dotchart(sort(uniq.pcmean), xlab='% of methylated bases in TEs', ylab='',pch=19,cex=1.5,main='');
abline(v=sort(uniq.pcmean), col='red', lty='dashed');

# class I vs II vs unknown
WGCNA::verboseBoxplot(c(apply(uniq.pc[c(2,3,5),],2,sum), apply(uniq.pc[c(1,4),],2,sum), t(uniq.pc[6, ])), 
					  c(rep('Class I', 4), rep('Class II', 4), rep('Unknown', 4)),
					  xlab='', ylab='% of methylated cytosines'
					  );

###

uniqOnly = list();
for (s in 1:length(names(uniq))) {
	print (names(uniq)[s])
	files = list.files()[grepl(paste(names(uniq)[s], '.*TEs\\.hit\\..*Only$', sep=''), list.files())];
	tmp = list();
	for (f in 1:length(files)) {
		print(files[f])
		tmp[[f]] = read.table(files[f])[, 1];
		names(tmp)[f] = gsub('_uniq_ids_Only', '', gsub(paste(names(uniq)[s], '.wig.bed-TEs.hit.', sep=''), '', files[f]));
	}
	uniqOnly[[s]] = tmp;
	names(uniqOnly)[s] = names(uniq)[s];
}; rm(s,tmp,f);


# get levels, from output from _scripts/getMethLevelsTEtypes.py
lvls = read.table('TEtypes_methLevels', nrows=24);
tmp = c();
for (row in 1:nrow(lvls)) {
	tmp = c(tmp, strsplit(lvls[row,1], '\\.')[[1]][5]);
}; rm(row);
lvls = cbind(lvls, tmp);
rm(tmp);

lvlsmean = sapply(split(lvls,lvls$tmp), function(f) mean(f[,4]));
lvlsplit = split(lvls,lvls$tmp);
lvlmat = sapply(lvlsplit,function(f) f[,4]);

WGCNA::verboseBoxplot(c(apply(lvlmat[, c(2,3,5)],1,sum), apply(lvlmat[, c(1,4)],1,sum), t(lvlmat[, 6])), 
					  c(rep('Class I', 4), rep('Class II', 4), rep('Unknown', 4)),
					  xlab='', ylab='mean methylation level'
					  );



source('../../_code/bootstrapFunctions_6-16-13.R');
bout=bootstrap.ANOVA(lvlmat[,order(apply(lvlmat,2,mean))], Func='mean');

###

te = read.table('../../_Burtoni_annotations/Abur_final_TE.bed', header=F, sep='\t');

te_types = table(te$V4);
total_bp = c();
for (type in 1:length(te_types)) {
	print(names(te_types)[type]);
	tmp = te[te$V4==names(te_types)[type], ];
	total_bp = c(total_bp, sum( as.numeric(tmp$V3) - as.numeric(tmp$V2) ));
	names(total_bp)[type] = names(te_types)[type];
}; rm(type);

genome_bp_estimate = 841927056;
wigs = c(36868745, 39926602, 23365799, 26630403);





WGCNA::verboseScatterplot(apply(uniq.pc,1,mean), ((total_bp) / genome_bp_estimate*100),
						  abline=T,abline.lty='dashed',abline.col='red',
						  xlab='% of TEs',ylab='% of genome covered',
						  frame.plot=F,type='n'
						  );
text(apply(uniq.pc,1,mean), ((total_bp) / genome_bp_estimate*100),labels=names(apply(uniq.pc,1,mean)));

#########
######

# get info on TEs upstream of CDS

setwd('/Volumes/fishstudies/_methylation/overlap_all');
files = list.files()[grepl('^3157.*TEs.*upstream$', list.files())];
s3157 = list();
for (f in 1:length(files)) {
	print(files[f]);
	s3157[[f]] = read.table(files[f], header=F, sep='\t', colClasses='character');
}; rm(f);
names(s3157) = files;








x=read.table('TEupstream_methLevels',header=F,sep='\t');
m = matrix(nrow=nrow(x)/2, ncol=3, dimnames=list(c(x[seq(1,nrow(x),2),1]),c('all','noUp','up')))
for (row in seq(2,nrow(x),2)) {
	m[row/2, ] = as.numeric(unlist(strsplit(x[row,1],' ')));
}; rm(row);

for (i in 1:6) {
	no = m[seq(i,nrow(m),6), 2];
	up = m[seq(i,nrow(m),6), 3];
	print(rownames(m)[i])
	print(mean(no));
	print(mean(up))
	print(wilcox.test(no, up))
}; rm(i,no,up);


#ms = split(as.data.frame(m),substr(rownames(m),1,4));

#types = names(table(unlist(strsplit(rownames(m),'\\.'))))[table(unlist(strsplit(rownames(m),'\\.')))==4];









# get levels, from output from _scripts/getMethLevelsTEtypes.py
lvls = read.table('TEtypes_methLevels', nrows=24);
tmp = c();
for (row in 1:nrow(lvls)) {
	tmp = c(tmp, strsplit(lvls[row,1], '\\.')[[1]][5]);
}; rm(row);
lvls = cbind(lvls, tmp);
rm(tmp);

lvlsmean = sapply(split(lvls,lvls$tmp), function(f) mean(f[,4]));
lvlsplit = split(lvls,lvls$tmp);
lvlmat = sapply(lvlsplit,function(f) f[,4]);













