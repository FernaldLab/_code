complement = matrix(c('A','T','C','G','T','A','G','C'), ncol=2, nrow=4)

n = c();
for (r in 1:nrow(x)) {
	if (x[r,4] == '+') {
		tmp = substring(x[r,5], 3, 4);
		n = c(n, tmp);
	} else if (x[r,4] == '-') {
		tmp = substring(x[r,5], 2, 3);
		n1 = substring(tmp,1,1);
		n2 = substring(tmp,2,2);
		n1 = complement[match(n1, complement[,1]), 2];
		n2 = complement[match(n2, complement[,1]), 2];
		tmp = paste(n2, n1, sep='')
		n = c(n, tmp);
	} else {
		print('invalid strand info');
	}
}; rm(r,n1,n2,tmp);