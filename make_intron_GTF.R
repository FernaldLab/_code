#x = read.table('~/Documents/_Fernald_lab/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.gtf', header=F, sep='\t');
#x = x[x$V3=='exon', ];
#save(x,file='~/Documents/_Fernald_lab/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.gtf_noCDS.RData')
#load('~/Documents/_Fernald_lab/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.gtf_noCDS.RData')
#x.int = as.data.frame(matrix(ncol=ncol(x)));
#xs = split(x, as.factor(apply(x, 1, function(f) strsplit(f[9],';')[[1]][1])));
#save(x.int, xs, file='~/Documents/_Fernald_lab/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.gtf_noCDS_split.RData');
rm(list=ls());
load(file='/Volumes/fishstudies/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.gtf_noCDS_split.RData');

for (g in 1:length(xs)) {
	print(g)
	gb = xs[[g]];
	gb = gb[order(gb$V4, gb$V5), ];
	gb.int = as.data.frame(matrix(ncol=ncol(gb)));
	gb = split(gb, as.factor(apply(gb, 1, function(f) strsplit(f[9],';')[[1]][2]))); 
	for (tr in gb) {#print(tr)
		ids = strsplit(tr[1,9], ';')[[1]];
		nrows = nrow(tr)-1;#print(nrows)
		if (nrows==0) {next}
		tr.int = as.data.frame(matrix(nrow=nrows, ncol=ncol(tr)));
		tr.int[, 1:3] = cbind(rep(tr[1,1], nrows), rep('ATH', nrows), rep('intron', nrows));
		tr.int[, 6:8] = tr[2:(nrows+1), 6:8];
		for (r in 2:nrow(tr)) {
			tr.int[r-1, 4:5] = c(tr[r-1, 5]+1, tr[r, 4]-1); 
			if (tr.int[r-1, 7] == '+') {
				tr.int[r-1, 9] = paste(ids[1], ids[2], paste(' intron_number ', r-1, ';', sep=''), sep=';');
			} else {
				tr.int[r-1, 9] = paste(ids[1], ids[2], paste(' intron_number ', nrows+2-r, ';', sep=''), sep=';');
			}
		}
		gb.int = rbind(gb.int, tr.int);
	}
	gb.int = gb.int[-1, ];
	gb.int = gb.int[order(gb.int$V4, gb.int$V5), ];#print(gb.int)
	x.int = rbind(x.int, gb.int);
}; rm(g, gb, ids, r, tr.int, nrows);

x.int = x.int[-1, ];
#print(x.int)
options(scipen=1)
x.int[70578,4] = 300000
x.int[grepl('e+',x.int$V4),]
save(x.int, file='~/Documents/_Fernald_lab/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.intron.gtf.RData')
write.table(x.int, file='/Volumes/fishstudies/_Burtoni_annotations/Astatotilapia_burtoni.BROADAB2fix.intron.gtf',row.names=F,col.names=F,quote=F,sep='\t');