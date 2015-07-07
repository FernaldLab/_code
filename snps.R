setwd('~/Documents/_Fernald_lab/_LYNLEY_RNAseq/snps/');
files = list.files()[grepl('gtfintersect', list.files())];

snps = vector(mode='list', length=4);
names(snps) = substr(files, 1, 6);

for (f in 1:length(files)) {
	snps[[f]] = read.table(files[f], header=F, sep='\t');
	cat('\n', files[f], '\n', sep='');
	print(head(snps[[f]]));
}; rm(f);
#save(snps, file='gtfintersectsJan24.RData');

qscore = 2000;
snps2 = snps;
for (f in 1:length(snps2)) {
	snps2[[f]] = snps2[[f]][snps2[[f]][,6] > qscore, ];
}; rm(f, qscore);
#save(snps2, file='gtfintersectsJan24_QS2000.RData');

# snps2common = snps2;
# pos = paste(snps2[[1]][,1], ':', snps2[[1]][,2], sep='');
# for (f in 2:length(snps2)) {
	# tmp = paste(snps2[[f]][,1], ':', snps2[[f]][,2], sep='');
	# pos = intersect(pos, tmp);
# }; rm(f, tmp);
# for (f in 1:length(snps2common)) {
	# tmp = paste(snps2common[[f]][,1], ':', snps2[[f]][,2], sep='');
	# snps2common[[f]] = snps2common[[f]][tmp %in% pos, ];
# }; rm(f, tmp);

scaffoldPos = function (df) {
	return(paste(df[,1], ':', df[,2], sep=''));
}


D = c('CGATGT', 'TTAGGC');
ND = c('ATCACG', 'TGACCA');
snps2D = snps2[names(snps2) %in% D];
pos1 = unique(scaffoldPos(snps2D[[1]]));
pos2 = unique(scaffoldPos(snps2D[[2]]));
tmp = intersect(pos1, pos2);
snps2D[[1]] = snps2D[[1]][scaffoldPos(snps2D[[1]]) %in% tmp, ];
snps2D[[2]] = snps2D[[2]][scaffoldPos(snps2D[[2]]) %in% tmp, ];
rm(pos1, pos2, tmp);

snps2ND = snps2[names(snps2) %in% ND];
pos1 = unique(scaffoldPos(snps2ND[[1]]));
pos2 = unique(scaffoldPos(snps2ND[[2]]));
tmp = intersect(pos1, pos2);
snps2ND[[1]] = snps2ND[[1]][scaffoldPos(snps2ND[[1]]) %in% tmp, ];
snps2ND[[2]] = snps2ND[[2]][scaffoldPos(snps2ND[[2]]) %in% tmp, ];
rm(pos1, pos2, tmp, snps2);
#save(snps2D, snps2ND, file='gtfintersectsJan24_QS2000_D-ND.RData');

posD = unique(scaffoldPos(snps2D[[1]]));
posND = unique(scaffoldPos(snps2ND[[1]]));
posCom = intersect(posD, posND);
posDonly = posD[!(posD %in% posCom)];
posNDonly = posND[!(posND %in% posCom)];

snps2Donly = snps2D;
snps2Donly[[1]] = snps2Donly[[1]][scaffoldPos(snps2Donly[[1]]) %in% posDonly, ];
snps2Donly[[2]] = snps2Donly[[2]][scaffoldPos(snps2Donly[[2]]) %in% posDonly, ];
snps2NDonly = snps2ND;
snps2NDonly[[1]] = snps2NDonly[[1]][scaffoldPos(snps2NDonly[[1]]) %in% posNDonly, ];
snps2NDonly[[2]] = snps2NDonly[[2]][scaffoldPos(snps2NDonly[[2]]) %in% posNDonly, ];
rm(posD, posND, posCom, posDonly, posNDonly, snps2D, snps2ND);

check = c();
for (r in 1:nrow(snps2Donly[[1]])) {
	tmp1 = strsplit(as.character(snps2Donly[[1]][r, 10]), ':')[[1]][2];	#print(tmp1)
	tmp2 = strsplit(as.character(snps2Donly[[2]][r, 10]), ':')[[1]][2];
	check1 = as.numeric(strsplit(tmp1, ',')[[1]][1]) > as.numeric(strsplit(tmp1, ',')[[1]][2]);
	check2 = as.numeric(strsplit(tmp2, ',')[[1]][1]) > as.numeric(strsplit(tmp2, ',')[[1]][2]);
	check = c(check, check1==check2);
}; rm(r, tmp1, tmp2, check1, check2);
snps2Donly[[1]] = snps2Donly[[1]][check, ];
snps2Donly[[2]] = snps2Donly[[2]][check, ];
rm(check);

check = c();
for (r in 1:nrow(snps2NDonly[[1]])) {
	tmp1 = strsplit(as.character(snps2NDonly[[1]][r, 10]), ':')[[1]][2];	#print(tmp1)
	tmp2 = strsplit(as.character(snps2NDonly[[2]][r, 10]), ':')[[1]][2];
	check1 = as.numeric(strsplit(tmp1, ',')[[1]][1]) > as.numeric(strsplit(tmp1, ',')[[1]][2]);
	check2 = as.numeric(strsplit(tmp2, ',')[[1]][1]) > as.numeric(strsplit(tmp2, ',')[[1]][2]);
	check = c(check, check1==check2);
}; rm(r, tmp1, tmp2, check1, check2);
snps2NDonly[[1]] = snps2NDonly[[1]][check, ];
snps2NDonly[[2]] = snps2NDonly[[2]][check, ];
rm(check);
