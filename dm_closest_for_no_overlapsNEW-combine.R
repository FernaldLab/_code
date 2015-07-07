rm(list=ls());
setwd('/Volumes/fishstudies/_methylation/new3157v3165overlaps/closest_to_noHits');

# read in data
files = list();
for (f in 1:length(list.files())) {
	this = read.table(list.files()[f],sep='\t',header=F,quote='',fill=T,colClasses='character');
	this = this[apply(this,1,function(f) !any(f==-1)), ];
	files[[f]] = this;
}; rm(f,this);
names(files) = list.files();

# get all methylation sites
sites = c();
for (f in files) {
	for (r in 1:nrow(f)) {
		sites = c(sites, paste(f[r,1:11], collapse=':'));
	}
}; rm(f,r);
sites = sort(unique(sites));

siteMat = data.frame(t(matrix(unlist(strsplit(sites, ':')), nrow=11)));
siteMat = as.data.frame(cbind(siteMat, closest=rep('NONE', nrow(siteMat))));
for (f in 1:length(files)) {
	this_f = files[[f]];
	print(names(files)[f]);
	for (r in 1:nrow(this_f)) {
		pos = paste(this_f[r, 1:11], collapse=':');
		hit = paste(this_f[r, 12:ncol(this_f)],collapse='\t');
		if (siteMat[sites==pos, 12] == 'NONE') {
			siteMat[sites==pos, 12] = hit;
		} else {
			siteMat[sites==pos, 12] = paste(siteMat[sites==pos, 12], hit, sep=';;');
		}
	}
}; rm(f, this_f, r, pos, hit);

thresh = 5000
siteMat2 = as.data.frame(cbind(siteMat, 
							   top=rep('', nrow(siteMat)),  
							   top1=rep('', nrow(siteMat)),
							   min_d=rep(0, nrow(siteMat))
							   )
						);
for (r in 1:nrow(siteMat2)) {
	h = strsplit(strsplit(as.character(siteMat2[r, 12]),';;')[[1]], '\t');
	dists = c();
	for (hh in h) {
		dists = c(dists, as.numeric(hh[length(hh)]));
	}
	siteMat2$min_d[r] = min(dists)
	ind = which(dists == min(dists));
	siteMat2$top1[r] = paste(h[ind][[1]],collapse='\t')
	if (any(dists<thresh)) {
		 for (d in h[dists<thresh]) {
		 	siteMat2$top[r] = paste(siteMat2$top[r], paste(d,collapse='\t'), sep=';;')
		 }
	}
}; rm(r, h, hh, dists, d, ind);

x = siteMat2[siteMat2$min_d==1,c(1:11,14,15)];
save.image(file='dm_closest_for_no_overlapsNEW-combineWORKSPACE.RData');


siteMat2cds_utr = siteMat2[grepl('CDS|five_prime_utr', siteMat2$top1, ignore.case=T), c(1:11,14:15)];
hstrands = c()
for(r in 1:nrow(siteMat2cds_utr)){
	hit = unlist(strsplit(siteMat2cds_utr$top1[r], '\t'));
	strand = hit[hit=='-' | hit=='+'];
	hstrands = c(hstrands, strand);
}; rm(r,hit,strand)
siteMat2cds_utr = as.data.frame(cbind(siteMat2cds_utr,hstrands));

siteMat2cds_utr_sense = siteMat2cds_utr[siteMat2cds_utr$X4==siteMat2cds_utr$hstrands,];
keep = c();
for(r in 1:nrow(siteMat2cds_utr_sense)){
	mpos = as.numeric(siteMat2cds_utr_sense[r,3]);
	h = unlist(strsplit(siteMat2cds_utr_sense$top1[r], '\t'));
	if (siteMat2cds_utr_sense$hstrands[r]=='+') {
		start = as.numeric(h[4]);
		if (start > mpos) {
			keep = c(keep, T);
			#print(siteMat2cds_utr_sense[r,])
			#print('up')
		} else {
			keep = c(keep, F);
		}
	} else {
		start = as.numeric(h[5]);
		if (start < mpos) {
			keep = c(keep, T);
			#print(siteMat2cds_utr_sense[r,])
			#print('up')
		} else {
			keep = c(keep, F);
		}
	}
}; rm(r,h,start,mpos);
siteMat2cds_utr_senseUp = siteMat2cds_utr_sense[keep,]
for (r in 1:nrow(siteMat2cds_utr_senseUp)) {
	h = unlist(strsplit(siteMat2cds_utr_senseUp$top1[r],'\t'));
	siteMat2cds_utr_senseUp$top1[r] = paste(h[1:(length(h)-1)],collapse='\t')
}; rm(r,h)








siteMat2cds_utr_antisense = siteMat2cds_utr[siteMat2cds_utr$X4!=siteMat2cds_utr$hstrands,];
keep = c();
for(r in 1:nrow(siteMat2cds_utr_antisense)){
	mpos = as.numeric(siteMat2cds_utr_antisense[r,3]);
	h = unlist(strsplit(siteMat2cds_utr_antisense$top1[r], '\t'));
	if (siteMat2cds_utr_antisense$hstrands[r]=='+') {
		start = as.numeric(h[4]);
		if (start > mpos) {
			keep = c(keep, T);
			#print(siteMat2cds_utr_antisense[r,])
			#print('up')
		} else {
			keep = c(keep, F);
		}
	} else {
		start = as.numeric(h[5]);
		if (start < mpos) {
			keep = c(keep, T);
			#print(siteMat2cds_utr_antisense[r,])
			#print('up')
		} else {
			keep = c(keep, F);
		}
	}
}; rm(r,h,start,mpos);
siteMat2cds_utr_antisenseUp = siteMat2cds_utr_antisense[keep,]
for (r in 1:nrow(siteMat2cds_utr_antisenseUp)) {
	h = unlist(strsplit(siteMat2cds_utr_antisenseUp$top1[r],'\t'));
	siteMat2cds_utr_antisenseUp$top1[r] = paste(h[1:(length(h)-1)],collapse='\t')
}; rm(r,h)









ugenes = c();
for (r in 1:nrow(siteMat2cds_utr_senseUp)) {
	h = unlist(strsplit(siteMat2cds_utr_senseUp$top1[r],'\t'));
	if (h[3]=='five_prime_utr') {
		hh = unlist(strsplit(h[9],'Parent='))[2];
		hhh = unlist(strsplit(hh,'.',fixed=T));
		hh = paste('ab.gene.',hhh[3],'.',hhh[4],sep='')
	} else if (h[3]=='CDS') {
		hh = unlist(strsplit(h[9],'\"'))[2]
	} else {
		stop('error')
	}
	ugenes = c(ugenes, hh)
}; rm(r,hh);

for (r in 1:nrow(siteMat2cds_utr_antisenseUp)) {
	h = unlist(strsplit(siteMat2cds_utr_antisenseUp$top1[r],'\t'));
	if (h[3]=='five_prime_utr') {
		hh = unlist(strsplit(h[9],'Parent='))[2];
		hhh = unlist(strsplit(hh,'.',fixed=T));
		hh = paste('ab.gene.',hhh[3],'.',hhh[4],sep='')
	} else if (h[3]=='CDS') {
		hh = unlist(strsplit(h[9],'\"'))[2]
	} else {
		stop('error')
	}
	ugenes = c(ugenes, hh)
}; rm(r,hh);