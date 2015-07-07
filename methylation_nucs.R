rm(list=ls());
#dir = '~/Documents/_LYNLEY_RNAseq/';
dir = '/Volumes/FishStudies/_methylation';
setwd(dir);
files = list.files()[grepl('wig.nucs$', list.files())];

for (f in 1:length(files)) {
	var = paste('n',substr(files[f],1,4),sep='');
	assign(var, read.table(files[f], sep='\t',header=F));
	print(head(get(var)));
}; rm(f, var);

all = list(n3157=n3157, n3165=n3165, n3581=n3581, n3677=n3677);
rm(n3157,n3165,n3581,n3677);
WGCNA::collectGarbage();
#save(all,file='all.wig.nucs.RData');
load(file='all.wig.nucs.RData');



flipNucsToPair = function(mat, nucCol) {
	newMat = mat;
	newMat[,nucCol] = as.character(newMat[,nucCol]);
	complement = matrix(c('A','T','C','G','T','A','G','C'), ncol=2, nrow=4);
	for (r in 1:nrow(mat)) {
		if(r %% 1000000 == 0){print(r)}
		n = as.character(mat[r, nucCol]);#print(n)
		if (substr(n,3,3)=='G'){
			tmp = complement[match(substr(n,2,2), complement[,1]), 2];#print(tmp)
			new_n = paste('C', tmp, sep='');#print(new_n)
			newMat[r, nucCol] = new_n;
		}
	}
	return(newMat);
}



cod = vector(mode='list', length=4);
names(cod) = names(all);
for (f in 1:length(all)) {
	all[[f]] = all[[f]][is.na(all[[f]]$V5), ];
	cod[[f]] = table(substr(all[[f]]$V4, 2, 4));
}; rm(f);


# for (s in 1:length(cod)) {
	# for (tr in 1:length(names(cod[[s]]))) {
		# b = substr(names(cod[[s]])[tr],2,2);
		# if (b=='G') {
			# t1 = complement[match(substr(names(cod[[s]])[tr],1,1), complement[,1]), 2];
			# t2 = 'C';
			# t3 = complement[match(substr(names(cod[[s]])[tr],3,3), complement[,1]), 2];
			# names(cod[[s]])[tr] = paste(t1,t2,t3,sep='');
		# }
	# }
# }; rm(s,tr,b,t1,t2,t3)


m = data.frame(matrix(nrow=length(cod[[1]]), ncol=4, dimnames=list(c(),names(cod))));
for (f in 1:length(cod)) {
	m[, f] = names(sort(cod[[f]], decreasing=T));
}; rm(f);


par(mfrow=c(2,2));for (f in 1:length(cod)) {
	b=barplot(sort(cod[[f]]),xaxt='n',xlab='',border='grey');
	text(x=b, y=-3, labels=names(sort(cod[[f]])), xpd=T, srt=45,cex=.8,adj=c(1,1));
}; rm(f);


cpg = list();
for (f in 1:length(cod)) {
	tmp = rev(names(sort(cod[[f]])));
	vec = c();
	for (j in 1:length(tmp)) {
		d1 = substr(tmp[j], 1, 2); d2 = substr(tmp[j], 2, 3);
		c1 = d1 %in% c('GC', 'CG'); c2 = d2 %in% c('GC', 'CG');
		if (any(c1, c2)) {
			vec = c(vec, j);
			names(vec)[length(vec)] = tmp[j];
		}
	}
	cpg[[f]] = vec;
	names(cpg)[f] = names(cod)[f];
}; rm(f, tmp, vec, j, d1, d2, c1, c2);

x = apply(data.frame(cod)[,c(2,4,6,8)],1,mean);
dotchart(sort(x),xlab='mean # (n=4)',main='triplets with methylation at mid base');