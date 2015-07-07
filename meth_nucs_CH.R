rm(list=ls());
options(stringsAsFactors=F);
dir = '/Volumes/FishStudies/_methylation';
setwd(dir);
files = list.files()[grepl('_CH$', list.files())];

for (f in 1:length(files)) {
	print(files[f]);
 	var = paste('n',substr(files[f],1,4),sep='');
 	assign(var, read.table(files[f],sep='\t',header=F,
 						   colClasses=c('character','numeric','numeric','character',
 						   				'numeric','character','character','character',
 						   				'character','character')
 						   				)
 						   );
 	#print(head(get(var)));
}; rm(f, var);
WGCNA::collectGarbage();
save(n3157, file='meth_nucs_CH_n3157.RData');
save(n3165, file='meth_nucs_CH_n3165.RData');
save(n3581, file='meth_nucs_CH_n3581.RData');
save(n3677, file='meth_nucs_CH_n3677.RData');

.filteredHist = function(dat, histCol=5, filtCol=10, filt=NULL, main=NULL, ...) {
	if(is.null(main)) {
		if(is.null(filt)) {
			main=paste(gsub('()','',match.call()[2]), ', all', sep='');
		} else {
			main=paste(gsub('()','',match.call()[2]), ', ', 
					   paste(filt,collapse='+'), 
					   sep='');
		}
	} else {
		main=main;
	}
	if(is.null(filt)) {
		hist(as.numeric(dat[,histCol]), main=main, ...);
	} else if(length(filt)==1){
		hist(as.numeric(dat[dat[,filtCol]==filt, histCol]), main=main, ...);
	} else {
		hist(as.numeric(dat[dat[,filtCol]%in%filt, histCol]), main=main, ...);
	}
}


col='grey';border='darkgrey';xlab='meth.lvl';ylab='';breaks=100;xlim=c(0,1.1);ylim=c(0,4000000);
par(mfrow=c(2,3));
filts=list(NULL,'CG',c('CH','CHG','CHH'),'CHG','CHH');
for (filt in filts) {
	.filteredHist(n3157,filt=filt,histCol=5,filtCol=10,
				  col=col,border=border,breaks=breaks,
				 # xlim=xlim,ylim=ylim,
				  xlab=xlab,ylab=ylab
				  );
}; rm(filt);


hist(as.numeric(dat$V5),);
hist(as.numeric(dat$V5[dat$V10=='CG']));
hist(as.numeric(dat$V5[dat$V10%in%c('CH','CHG','CHH')]));
hist(as.numeric(dat$V5[dat$V10=='CHG']));
hist(as.numeric(dat$V5[dat$V10=='CHH']));