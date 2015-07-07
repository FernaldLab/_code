rm(list=ls());
#setwd('/Volumes/fishstudies/_methylation/newoverlap_all');
setwd('~/Documents/_methylation/');

subj = c('3157','3165','3581','3677');
#files = paste(subj, '_Astatotilapia_burtoni.BROADAB2fix_CDS.gtf-overlap_groupedV2_forR_CDScommon', sep='');
files = paste(subj, '_Astatotilapia_burtoni.BROADAB2fix_CDS.gtf-overlap_exon1s_groupedV2_forR_CDScommon', sep='');

for (f in 1:length(files)) {
	print(files[f]);
 	var = paste('f',substr(files[f],1,4),sep='');
 	assign(var, read.table(files[f],sep='\t',header=F,colClasses='character',col.names=c('cds', 'length', 'strand', 'gene', 'hit_count', 'pc_coverage', 'avg_lvl', 'dist_start', 'dist_start_pc', 'nucs', 'sense', 'levels')));
 	#print(head(get(var)));
}; rm(f, var);
#subj = paste('n',subj,sep='');


.getCNNinfo = function(subj, nucs) {
	out = list();
	lvls = list();
	sense = list();
	dist = list();
	dist_pc = list();
	for (s in 1:length(subj)) {
		pos = unlist(strsplit(get(subj[s])[[10]],',')) == nucs;
		lvl = as.numeric(unlist(strsplit(get(subj[s])[[12]],','))[pos]);
		sense = as.character(unlist(strsplit(get(subj[s])[[11]],','))[pos]);
		dist = as.numeric(unlist(strsplit(get(subj[s])[[8]],','))[pos]);
		dist_pc = as.numeric(unlist(strsplit(get(subj[s])[[9]],','))[pos]);
		out[[s]] = list(lvl=lvl,sense=sense,dist=dist,dist_pc=dist_pc);
		names(out)[s] = subj[s];
	}
	return(out);
}
CG = .getCNNinfo(subj,'CG');
CHG = .getCNNinfo(subj,'CHG');
CHH = .getCNNinfo(subj,'CHH');

# par(mfrow=c(4,3), oma=c(0,0,2,0));
# b = 1000;
# for (s in 1:4){
	# hist(CG[[s]]$dist_pc,breaks=b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	# hist(CHG[[s]]$dist_pc,breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	# hist(CHH[[s]]$dist_pc,breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
# }
# title(paste('both strands, breaks=',b,sep=''),outer=T);

par(mfrow=c(4,3), oma=c(0,0,2,0));
b = 1000; 
for (s in 1:4){
	cg = CG[[s]]$dist_pc;
	chg = CHG[[s]]$dist_pc;
	chh = CHH[[s]]$dist_pc;
	hist(cg[cg > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	hist(chg[chg > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	hist(chh[chh > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
}
title(paste('all exons, both strands, breaks=',b,sep=''),outer=T);

par(mfrow=c(4,3), oma=c(0,0,2,0));
b = 1000; 
for (s in 1:4){
	cg = CG[[s]]$dist_pc;
	chg = CHG[[s]]$dist_pc;
	chh = CHH[[s]]$dist_pc;
	hist(cg[cg > 0],breaks=length(cg)/b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	hist(chg[chg > 0],breaks=length(chg)/b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	hist(chh[chh > 0],breaks=length(chh)/b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
}
title(paste('all exons, both strands, breaks=length/',b,sep=''),outer=T);


par(mfrow=c(4,3), oma=c(0,0,2,0));
b = 1000; 
for (s in 1:4){
	cg = CGf[[s]]$dist_pc;
	chg = CHGf[[s]]$dist_pc;
	chh = CHHf[[s]]$dist_pc;
	hist(cg[cg > 0],breaks=length(cg)/b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	hist(chg[chg > 0],breaks=length(chg)/b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	hist(chh[chh > 0],breaks=length(chh)/b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
}
title(paste('first exons, both strands, breaks=length/',b,sep=''),outer=T);



par(mfrow=c(4,3), oma=c(0,0,2,0));
sense='as'; b=1000;
for (s in 1:4){
	cg = CG[[s]]$dist_pc[CG[[s]]$sense==sense]; 
	chg = CHG[[s]]$dist_pc[CHG[[s]]$sense==sense]; 
	chh = CHH[[s]]$dist_pc[CHH[[s]]$sense==sense]; 
	hist(cg[cg > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	hist(chg[chg > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	hist(chh[chh > 0],breaks=b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
}
title(paste('all exons, ', sense, ' strand, breaks=',b,sep=''),outer=T);

par(mfrow=c(4,3), oma=c(0,0,2,0));
sense='as'; b=1000;
for (s in 1:4){
	cg = CG[[s]]$dist_pc[CG[[s]]$sense==sense]; 
	chg = CHG[[s]]$dist_pc[CHG[[s]]$sense==sense]; 
	chh = CHH[[s]]$dist_pc[CHH[[s]]$sense==sense]; 
	hist(cg[cg > 0],breaks=length(cg)/b,col='grey',border='darkgrey',main=paste(s, ' CG'));
	hist(chg[chg > 0],breaks=length(chg)/b,col='grey',border='darkgrey',main=paste(s, ' CHG'));
	hist(chh[chh > 0],breaks=length(chh)/b,col='grey',border='darkgrey',main=paste(s, ' CHH'));
}
title(paste('all exons, ', sense, ' strand, breaks=length/',b,sep=''),outer=T);






.plotHistFromColumn = function(x, column, breakDiv=1000) {
	nums = as.numeric(x[,column]);
	n = nrow(x);
	m = signif(mean(nums),3);
	med = signif(median(nums),3);
	main = paste('n=', n, ', mean=', m, ', median=', med, sep='');
	hist(nums, breaks=n/breakDiv, col='grey', border='darkgrey', xlab=colnames(x)[column], main=main);
	abline(v=m, col='red');
	abline(v=med, col='blue');
}

.makeCDSdf = function(cds) {
	dists = as.numeric(unlist(strsplit(cds$dist_start, ',')));
	dists_pc = as.numeric(unlist(strsplit(cds$dist_start_pc, ',')));
	nucs = unlist(strsplit(cds$nucs, ','));
	sense = unlist(strsplit(cds$sense, ','));
	levels = unlist(strsplit(cds$levels, ','));
	df = data.frame(dists, dists_pc, nucs, sense, levels);
	return(df)
}
.makeCDSmat(x[1,])

.cdsStatsFromDataFrame = function(df) {
	pc.CG = sum(df$nucs=='CG') / nrow(df);
	pc.CHG = sum(df$nucs=='CHG') / nrow(df);
	pc.CHH = sum(df$nucs=='CHH') / nrow(df);
	#pc.CH = sum(df$nucs%in%c('CH','CHG','CHH')) / nrow(df);
	lvl.CG = mean(as.numeric(df$levels[df$nucs=='CG']));
	lvl.CHG = mean(as.numeric(df$levels[df$nucs=='CHG']));
	lvl.CHH = mean(as.numeric(df$levels[df$nucs=='CHH']));
	return(stats);
}

.collapseToGeneList = function(x) {
	genes = unique(x$gene);
	out = list()
	for (g in 1:length(genes)) {
		tmp = x[x$gene==genes[g], ];
		cds_ints = paste(tmp$cds, collapse=',');
		total_length = sum(as.numeric(tmp$length));
		total_hits = sum(as.numeric(tmp$hit_count));
		total_pc_coverage = total_hits / total_length;
		cds_avg_length = mean(as.numeric(tmp$length));
		cds_avg_hits = mean(as.numeric(tmp$hit_count));
		cds_avg_coverage = mean(as.numeric(tmp$pc_coverage));
		cds_avg_level = mean(as.numeric(tmp$avg_lvl));
		all_dists = unlist(strsplit(tmp$dist_start, ','));
		all_dists_pc = unlist(strsplit(tmp$dist_start_pc, ','));
		all_nucs = unlist(strsplit(tmp$nucs, ','));
		all_sense = unlist(strsplit(tmp$sense, ','));
		all_levels = unlist(strsplit(tmp$levels, ','));
		all = data.frame(all_dists, all_dists_pc, all_nucs, all_sense, all_levels);
		names(all) = gsub('all_', '', names(all));
		stats = data.frame(total_length, total_hits, total_pc_coverage, 
						   cds_avg_length, cds_avg_hits, cds_avg_coverage, cds_avg_level
						   );
		out[[g]] = list(cds=cds_ints, stats=stats, all=all);
		names(out)[g] = genes[g];
	}
	return(out);
}
.collapseToGeneList(head(x))



par(mfrow=c(1,2));
.plotHistFromColumn(x,7);
.plotHistFromColumn(x,6);




x=CG$n3157
#dists = c();
dists_pc = c();
for (row in 1:nrow(x)) {
	if (row %% 10000 == 0) {cat(signif(row/nrow(x)), '\n', sep='')}
	#dists = c(dists, unlist(strsplit(x$dist_start[row], ',')));
	dists_pc = c(dists_pc, unlist(strsplit(x$dist_start_pc[row], ',')));
}; rm(row);
#dists = as.numeric(dists);
dists_pc = as.numeric(dists_pc);
hist(dists_pc,breaks=length(dists_pc),col='grey',border='darkgrey');