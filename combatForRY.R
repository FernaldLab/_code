### remove batch effects
source('/Volumes/fishstudies/_code/ComBat.R');

# format expression data and write to file
datToWrite = as.data.frame(datImpute);
datToWrite = as.data.frame(cbind(gene=rownames(datToWrite), datToWrite));
write.table(datToWrite, file='datForComBat.txt', quote=F, sep='\t', row.names=F);

# format sample info table and write to file
# assumes samples in rows of rg match samples in cols of datToWrite
rg0 = read.table('../../../readGroupsLibs.txt',sep='\t',header=T);
rg = rg0; 
rg$BatchSeq = gsub('2015', '', rg$BatchSeq);
rg = cbind(names(datToWrite)[2:ncol(datToWrite)], rg);
rg = rg[, c(1,6,2,3,7)];
names(rg) = c('Array name', 'Sample name', 'Batch', 'Covariate 1', 'Covariate 2');
write.table(rg, file='rgForComBat.txt', quote=F, sep='\t', row.names=F);

# run ComBat to correct for library batch effect
combatout = ComBat(expression_xls='datForComBat.txt', 
				   sample_info_file='rgForComBat.txt', 
				   filter=F, write=F, skip=1
				   );