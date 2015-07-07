# # X is data

# X = scale(X);

# V <- t( eigen( t(X) %*% X , T )$vectors ) 
# U <- eigen( X %*% t( X ) , T )$vectors 
# D <- diag( sqrt( eigen( X %*% t( X ) , T )$values ) ) 

# X - U %*% D %*% t(V) 

# svd(X)$v - V 
# svd(X)$u - U 
# svd(X)$d - sqrt( eigen( X %*% t( X ) , T )$values ) 

# X-svd(X)$u %*% diag( svd(X)$d ) %*% t( svd(X)$v )

#########

runPCA = function(matAdj = 'centered and scaled matrix')
{
	return(eigen(cov(matAdj)))
}

varExplained = function(eigenList)
{
	par(mfrow=c(1,2));
	
	plot(
		eigenList$values / sum(eigenList$values),
		pch=21,
		col='black', bg = '#549cc4',
		ylim=c(0,1),
		xlab='PC', ylab='Variance Explained'
		) + abline(h=0.9)
		
	plot(
		cumsum(eigenList$values) / sum(eigenList$values),
		pch=21,
		col='black', bg = '#549cc4',
		ylim=c(0,1),
		xlab='PC', ylab='Cumulative Variance Explained'
		) + abline(h=0.9)
}

afterPCA = function(
				   matAdj = 'centered and scaled matrix',
				   meanList = 'list of col means of unadjusted matrix',
				   eigenList = 'List of eigenvals and eigenvec of matAdj cov matrix',
				   n = 'selected PCs',
				   specific_select = 'if T: n==1:n, else just n\'th cols'
				   )
{
	if (length(n) > ncol(matAdj))
	{
		stop('N is higher than number of PCs')
	}
	if (!specific_select & length(n) > 1)
	{
		stop('Use single number when selecting up to n\'th PC')
	}
	temp1 = t(eigenList$vectors[,n] %*% (t(eigenList$vectors[,n]) %*% t(matAdj))); #print(temp1)
	temp2 = t(matrix(meanList, nrow = nrow(matAdj), ncol = ncol(matAdj))); #print(temp2)
	temp1+temp2
}


# reconstMatrix = afterPCA(
 # matAdjust = apply(originalData, 2, function(i) i - mean(i)),
 # meanList = apply(originalData, 2, mean),
 # eigenList = pca,
 # n = 5,
 # specific_select = FALSE
# )
