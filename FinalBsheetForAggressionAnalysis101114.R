setwd('/Volumes/fishstudies/_behaviorRA/');
dat0 = read.csv('FinalBsheetForAggressionAnalysis101114.csv');

dat.stable = dat0[dat0$timepoint=='stable', ];

cor.test(dat.stable$Days.in.dyad.without.kill, dat.stable$Total.Kills);

killcounts = as.numeric(names(table(dat.stable$Total.Kills)));
mat = as.data.frame(matrix(nrow=length(killcounts), ncol=3, 
						   dimnames=list(killcounts, 
						   				 c('Total.Kills','median-Days.in.dyad.without.kill', 'n')
						   				 )
						   )
					);
for (k in killcounts) {
	mat[k+1, ] = c(as.numeric(k),
				   median(dat.stable$Days.in.dyad.without.kill[dat.stable$Total.Kills==k]),
				   sum(dat.stable$Total.Kills==k));
}; rm(k);





datM = dat.stable[dat.stable$X.males.killed==2 & dat.stable$Total.Kills==2, ]