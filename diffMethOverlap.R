setwd('~/Documents/_Fernald_lab/_methylation/');
files = list.files()[grep('significant$', list.files())];
raw = list();
allpos = c();
for (f in 1:length(files)) {
	raw[[f]] = read.table(files[f], sep='\t', header=F, colClasses='character');
	names(raw)[f] = files[f];
	allpos = c(allpos, paste(raw[[f]][,1], ':', raw[[f]][,2], sep=''));
}; rm(f);
allpostab = table(allpos);
overlap_all = names(allpostab)[allpostab==4];





for (f in 1:length(raw)) {
	
}; rm(f);