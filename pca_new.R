# you probably want to use stats::prcomp for PCA on big matrices
load('~/face.rda')
runPCA <- function(mat = 'Unadjusted matrix') {
	eigen(cov(apply(mat, 2, function(i) i - mean(i))))
}
pca <- runPCA(faceData)
 
# > str(pca)
# List of 2
# $ values : num [1:32] 18915 10067 4724 2506 1366 ...
# $ vectors: num [1:32, 1:32] -0.236 -0.261 -0.192 -0.164 -0.204 ...



varExplained <- function(eigenList) {
 
par(mfrow = c(1,2))
 
plot(
 eigenList$value / sum(eigenList$value), pch = 21, col = 'black',
 bg = '#549cc4', ylim = c(0, 1), xlab = 'Principal Component',
 ylab = 'Variance Explained'
 ) + abline(h = 0.9)
 
plot(
 cumsum(eigenList$value) / sum(eigenList$value), pch = 21,
 col = 'black', bg = '#549cc4', ylim = c(0, 1), xlab = 'Principal Component',
 ylab = 'Cumulative Variance Explained'
 ) + abline(h = 0.9)
}
 
varExplained(pca)


afterPCA <- function(
 matAdjust = 'Centered matrix',
 meanList = 'List of column means of original (unadjusted) matrix',
 eigenList = 'List of eigenvalues and eigenvectors of adjust matrix covariance matrix',
 n = 'selected PC\'s',
 specific_select = 'If True: n == 1:n, if False: just n\'th columns') {
 
 if (length(n) > ncol(matAdjust)) stop('N is higher than the number of PC\'s')
 if (!specific_select & length(n) > 1) stop('Use a single number when selecting up to n\'th PC')
 if (!specific_select) n <- 1:n
 
 t(eigenList$vectors[,n] %*% (t(eigenList$vectors[,n]) %*% t(matAdjust))) + t(matrix(meanList, nrow = nrow(matAdjust), ncol = ncol(matAdjust)))
}
 
# ColorBrewer palette
library(RColorBrewer)
showMatrix <- function(x, ...) image(t(x[nrow(x):1,]), xaxt = 'none', yaxt = 'none', col = rev(colorRampPalette(brewer.pal(7, 'Blues'))(100)), ...)
 
reconstMatrix <- afterPCA(
 matAdjust = apply(faceData, 2, function(i) i - mean(i)),
 meanList = apply(faceData, 2, mean),
 eigenList = pca,
 n = 5,
 specific_select = FALSE
)
 
par(mfrow = c(1,2), mar = c(0, 0, 1, 0), bty = 'n')
showMatrix(faceData, main = 'Original Matrix')
showMatrix(reconstMatrix, main = 'First 5 PC\'s')
