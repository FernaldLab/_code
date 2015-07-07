setwd('/Volumes/fishstudies/_methylation/');
library(WGCNA);
.collapseFwdRevRows = function(mat) {
	mat2 = matrix(ncol=ncol(mat),
				  nrow=nrow(mat)/2,
				  dimnames=list(gsub('_fwd','',rownames(mat)[seq(1,nrow(mat),2)]), colnames(mat))
				  );
	for (r in 1:nrow(mat2)) {
		rr = 1 + (r-1)*2;
		mat2[r, ] = (mat[rr,] + mat[(rr+1),]) / 2;
	}
	return(mat2);
}
.verboseBoxplotFromMatrix = function(mat,...) {
	g = c();
	for (name in colnames(mat)) {
		g = c(g, rep(name, 4));
	}
	verboseBoxplot(as.vector(mat),as.factor(g),col='grey',xlab='',ylab='meth.level',...);
}


files_di = list.files()[grepl('di_lvls$',list.files())];
di = list();
for (f in 1:length(files_di)) {
	di[[f]] = t(read.table(files_di[f],row.names=1,colClasses=c('character',rep('numeric',4))));
	names(di)[f] = files_di[f];
}; rm(f);
di_mean = matrix(ncol=ncol(di[[1]]),nrow=length(di),dimnames=list(gsub('_Gflip|_di_lvls','',names(di)),colnames(di[[1]])));
for (f in 1:length(di)) {
	di_mean[f, ] = di[[f]][3,];
}; rm(f);
di_mean = di_mean[, -grep('N', colnames(di_mean))];
di_mean2 = .collapseFwdRevRows(di_mean);







files_tri = list.files()[grepl('tri_lvls$',list.files())];
tri = list();
for (f in 1:length(files_tri)) {
	tri[[f]] = t(read.table(files_tri[f],row.names=1,colClasses=c('character',rep('numeric',4))));
	names(tri)[f] = files_tri[f];
}; rm(f);
tri_mean = matrix(ncol=ncol(tri[[1]]),nrow=length(tri),dimnames=list(gsub('_Gflip|_tri_lvls','',names(tri)),colnames(tri[[1]])));
for (f in 1:length(tri)) {
	tri_mean[f, ] = tri[[f]][3,];
}; rm(f);
tri_mean = tri_mean[, -grep('N', colnames(tri_mean))];
tri_mean2 = .collapseFwdRevRows(tri_mean);






CHG = tri_mean2[, grep('C[ACT]G',colnames(tri_mean2))];
CHH = tri_mean2[, grep('C[ACT][ACT]',colnames(tri_mean2))];
CpG = tri_mean2[, grep('CG[ACGT]',colnames(tri_mean2))];