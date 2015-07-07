# study scaffold lengths

sl = read.table('/Volumes/fishstudies/_Burtoni_genome_files/scaffold_lengths');

total = sum(sl[,2]);
agg = c();
for (i in 1:nrow(sl)) {
	upto = sum(sl[1:i, 2]);
	pct = upto / total;
	agg = c(agg, pct);
}; rm(i,upto,pct);


an = read.table('/Volumes/fishstudies/_Burtoni_annotations/Astatotilapia_burtoni.BROADcombo.gtf_scaffold_counts');
tmp = unlist(strsplit(an$V2,'_'));
tmp = as.numeric(tmp[seq(2, length(tmp), 2)]);
an = an[order(tmp),];

total2 = sum(an[,1]);
agg2 = c();
for (i in 1:nrow(an)) {
	upto2 = sum(an[1:i, 1]);
	pct2 = upto2 / total2;
	agg2 = c(agg2, pct2);
}; rm(i,upto2,pct2);

par(mfrow=c(1,2)); plot(1:nrow(sl), agg, type='l'); plot(1:nrow(an), agg2, type='l');


sl[which(agg>.95)[1], ];
an[which(agg2>.95)[1], ];


par(mfrow=c(1,2)); 
plot(1:nrow(sl), agg, type='l'); 
abline(v=which(agg>.95)[1], col='red');
plot(1:nrow(an), agg2, type='l');
abline(v=which(agg2>.95)[1], col='red');

##########################################

load('aligned.adapters.q30.m0.methratio.CG.clean.upto1494_d0.RData');

common = intersect(d0[[1]]$chr, d0[[2]]$chr);
common = intersect(common, d0[[3]]$chr);
common = intersect(common, d0[[4]]$chr);

tmp = unlist(strsplit(common, '_'));
tmp = as.numeric(tmp[seq(2, length(tmp), 2)]);

##########

tNAs = which(is.na(attributes(bsdall.fit.tstat2x)$stats[,5]));
chrNAs = seqnames(granges(bsdall.fit2x)[tNAs]);