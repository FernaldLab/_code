# analyze preseq output

setwd('/Volumes/OSX/_Fernald_lab/_LYNLEY_RNAseq/tophatOUT/GTF/padded');
lc_files = list.files()[grep('^lc', list.files())];
lc_files_verbose = lc_files[grepl('verbose', lc_files)];
lc_files = lc_files[!grepl('verbose', lc_files)];
lc = list();
for (file in 1:length(lc_files)) {
	lc[[file]] = read.table(lc_files[file], header=T); 
	names(lc)[file] = gsub('.txt', '', gsub('lc_extrap', '', lc_files[file]));
}; rm(file, lc_files);
lc_verbose = list();
for (file in 1:length(lc_files_verbose)) {
	lc_verbose[[file]] = read.table(lc_files_verbose[file], header=T, nrows=2); 
	names(lc_verbose)[file] = gsub('.txt', '', gsub('lc_extrap', '', lc_files_verbose[file]));
}; rm(file, lc_files_verbose);

upto = 1000;
ymax = max(unlist(lapply(lc, function(x){max(x$EXPECTED_DISTINCT)})));
col = colors(distinct=T)[grepl('blue|green|red|black|purple|orange', colors(distinct=T))];
col = col[!grepl('light', col)];
plot(lc[[1]]$TOTAL_READS[1:upto], lc[[1]]$EXPECTED_DISTINCT[1:upto], ylim=c(0,ymax), type='n', xlab='TOTAL_READS', ylab='DISTINCT_READS');
for (s in 1:length(lc)) {
	this_col = sample(col, 1);
 	lines(lc[[s]]$TOTAL_READS[1:upto], lc[[s]]$EXPECTED_DISTINCT[1:upto], col=this_col);
 	abline(v=lc_verbose[[s]][1, 3], col=this_col);
 	abline(h=lc_verbose[[s]][2, 3], col=this_col);
}; rm(s);
