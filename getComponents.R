# dat: rows are genes, cols are samples
getComponents = function (dat)
{
	tmp=impute.knn(as.matrix(dat));
	dat=as.matrix(tmp$data);
	dat = t(scale(t(dat)));
	svd1=svd(dat,nu=ncol(dat), nv=ncol(dat));
	return(svd1);
}