# rm(list=ls());
# setwd('/Volumes/fishstudies/_Burtoni_genome_files');

# files = list.files()[grepl('dump$',list.files())];
# for (f in 1:length(files)) {
	# tmp = read.table(files[f]);
	# vals = tmp[,2];
	# names(vals) = tmp[,1];
	# assign(paste('mer',f,sep=''), vals);
# }; rm(f,tmp,vals)

z1=read.table('../_Burtoni_genome_files/1mer_counts.jf.dump');
z2=read.table('../_Burtoni_genome_files/2mer_counts.jf.dump');
probs = matrix(nrow=4,ncol=4,dimnames=list(c('A','C','G','T'),c('A','C','G','T')))
for (row in 1:nrow(z2)) {
	leader = substr(z2[row,1],1,1);
	follower = substr(z2[row,1],2,2);
	probs[match(leader,rownames(probs)),match(follower,colnames(probs))] = z2[row,2] / z1[match(leader, z1[,1]),2];
}; rm(row,leader,follower)