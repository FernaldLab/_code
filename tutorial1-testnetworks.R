#################################
###### build test networks ######
#################################

# start up
rm(list = ls());
setwd('~/Documents/_analysis_compare');
library('WGCNA');
options(stringsAsFactors=F);
load(file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed');

datList = DATAcollapsedList;
rm(DATAcollapsedList);
collectGarbage();

# compute scale-free diagnostics to choose power for adjacency 
pickSoft = list();
for( j in 1:3 ) {
	dat = as.data.frame( t( datList[[j]]$datETcollapsed ) );
	pickSoft[[j]] = pickSoftThreshold(dat, networkType = 'signed', blockSize = 6000, verbose = 2);
	collectGarbage();
	}
	rm(j, dat);
names(pickSoft) = names(datList);
#save(pickSoft, file = '_Rdata/VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed_pickSoft')

# plot scale-free fit, mean k, and slope of fit as function of power

#jpeg(file = 'test_networksVSP.jpg', width = 8, height = 10, units = 'in', quality = 100, type = 'quartz', res = 150);

par( mfrow = c(3, 3) );
for( j in 1:3 ) {
	set = pickSoft[[j]]$fitIndices;
	x = set$Power;
	y = set$SFT.R.sq;
	plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = 0.8, col = 'red', lty = 'dashed' );
	
	y = set$mean.k.;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = 50, col = 'red', lty = 'dashed' );
	
	y = set$slope;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = names(pickSoft)[j]);
	text(x, y, labels = x);
	abline( h = -2, col = 'red', lty = 'dashed' );
	abline( h = -1, col = 'red', lty = 'dashed' );
	}
	rm(j, set, x, y);
	
#dev.off();


# or to do one at a time
set0 = 'deviate2';

set0 = match(set0, names(pickSoft));
par( mfrow = c(1, 3) );
set = pickSoft[[set0]]$fitIndices;
x = set$Power;
y = set$SFT.R.sq;
plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = 0.8, col = 'red', lty = 'dashed' );

y = set$mean.k.;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = 50, col = 'red', lty = 'dashed' );
	
y = set$slope;
plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = names(pickSoft)[set0]);
text(x, y, labels = x);
abline( h = -2, col = 'red', lty = 'dashed' );
abline( h = -1, col = 'red', lty = 'dashed' );

rm(set0, set, x, y);




# use deviate3, power = 14
datVSP = datList$deviate3$datETcollapsed;
datVSP = as.data.frame( t(datVSP) );
k14 = softConnectivity(datVSP, type = 'signed', power = 14, blockSize = 6000, verbose = 3);
names(k14) = names(datVSP);
#save(k14, file = '_Rdata/VSPprobes_removed_deviate3_tooManyNAs_no1S3S_Qnorm_collapsed_k14');
hist(k14);
collectGarbage();
