n=10000; mu=100; size=1; 
#xlim=c(0,250); ylim=c(0,450);
z=rnbinom(n, size=size, mu=mu); 
y=rpois(n, mu)

par(mfrow=c(1,2))
hist(z, breaks=n, main='negative binomial'
#,xlim=xlim, ylim=ylim
); 
abline(v=mu, col='red'); 
abline(v=median(z), col='grey');

hist(y, breaks=n, main='poisson'
#,xlim=xlim, ylim=ylim
); 
abline(v=mean(y), col='red'); 
abline(v=median(y), col='grey');

###############

rna=read.csv('~/Documents/_Fernald_lab/_RNAseq/mRNA_Diff_Exp_with_miRNA_regulators-1-FPKMonly.csv', header=T);
names(rna)[3:4]=c('dom','sub');
rna = rna[order(rna$Gene.model), ];

nToSim = 10;
rnaSIM = cbind(rna, matrix(nrow=nrow(rna), ncol=nToSim));
rownames(rnaSIM) = paste( rnaSIM[,1], '_', rnaSIM[,2], sep='');
rnaSIM = rnaSIM[, -c(1,2)];
names(rnaSIM)[3:length(names(rnaSIM))]=paste('sim', names(rnaSIM)[3:length(names(rnaSIM))], sep='');

for (transcript in 1:nrow(rna))
{
	
}