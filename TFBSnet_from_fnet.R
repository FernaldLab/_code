fnet = read.table('~/Downloads/dmg_all_promoters_hs_2000_1.00_comp-trecount_list2ref_fdr_hypothesis.fnet.txt',header=T);
fnet[duplicated(fnet[,1]),1]=paste(fnet[duplicated(fnet[,1]),1],'.2',sep='');
fnet[duplicated(fnet[,1]),1]=paste(fnet[duplicated(fnet[,1]),1],'.2',sep='');
rownames(fnet)=fnet[,1];
fnet=fnet[,-1];

library(WGCNA);

ps=pickSoftThreshold(t(fnet),networkType='signed');
par( mfrow = c(1, 3) );
set = ps$fitIndices;
x = set$Power;
y = set$SFT.R.sq;
plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '');
text(x, y, labels = x);
abline( h = 0.8, col = 'red', lty = 'dashed' );

y = set$mean.k.;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '');
text(x, y, labels = x);
abline( h = 50, col = 'red', lty = 'dashed' );
y = set$slope;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = '');
text(x, y, labels = x);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

par(mfrow=c(1,2));hist(softConnectivity(t(fnet),type='signed',power=15));scaleFreePlot(softConnectivity(t(fnet),type='signed',power=15));

adj = adjacency(t(fnet), type='signed', power=15)
top = TOMdist(adj)
dendro = flashClust(as.dist(top), method = "average");
labels=cutreeDynamic(dendro,minClusterSize=10,distM=top);
MEs=moduleEigengenes(t(fnet),labels)$eigengenes;
#mergeCloseModules()
colors=labels2colors(labels);
plotDendroAndColors(dendro,colors,dendroLabels=F,rowText=labels);

MEs = moduleEigengenes(t(fnet), colors)$eigengenes;

kME = as.data.frame(cor(t(fnet), MEs, use='p'));

tkME=kME[colors=='turquoise',];
tkME=tkME[order(-tkME$MEturquoise),];
x=rownames(tkME[1:50,]);
fnet[rownames(fnet)%in%x,]
x=gsub('.2','',x,fixed=T);
myb=names(IDs)[IDs %in% x]
diff.mat[diff.mat$sym %in% myb,]