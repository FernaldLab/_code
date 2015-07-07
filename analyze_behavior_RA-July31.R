rm(list=ls());
setwd('/Volumes/fishstudies/_behaviorRA');
#setwd('/Volumes/fishstudies/_behaviorRA/Nov2013Behavior');
source('/Volumes/fishstudies/_code/ethograms_from_scorevideo.R');

# load kill data and sheet to translate project id to fish id
dat0 = read.csv('masterDBfiltered072914_edited_in_R.csv');
ids0 = read.csv('Experience_Aggression_hormones061914_sheet-database.csv')[, 1:2];

### load establishment day data ###

### Nov 6 with April 23,nov 7 with April 24, nov 25 with May 12, nov 26 with May 13

# combine day1's (nov6, apr23) for establishment then day2's (nov7, apr24)
dir = 'Nov2013Behavior';
filesNov6 = paste(dir, '/', list.files(dir)[grepl('.*nov6.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
filesNov7 = paste(dir, '/', list.files(dir)[grepl('.*nov7.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
dir = 'EstablishmentAprilAllmales';
filesApr23 = paste(dir, '/', list.files(dir)[grepl('.*apr23.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
filesApr24 = paste(dir, '/', list.files(dir)[grepl('.*apr24.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');


files_establishment1 = c(filesNov6, filesApr23);
files_establishment2 = c(filesNov7, filesApr24);

rawdataEst1 = .getDataBatch(files_establishment1);
names(rawdataEst1) = gsub('Nov2013Behavior/|EstablishmentAprilAllmales/', '', names(rawdataEst1));
rawdataEst2 = .getDataBatch(files_establishment2);
names(rawdataEst2) = gsub('Nov2013Behavior/|EstablishmentAprilAllmales/', '', names(rawdataEst2));

rawdataEst1_ids = as.character(.getProjectIdFromFilenames(names(rawdataEst1)));
rawdataEst2_ids = as.character(.getProjectIdFromFilenames(names(rawdataEst2)));

# animals that were reused in multiple dyads
# # # toRemove = c('33','41','49','52');
# # # rawdataEst1 = rawdataEst1[-(which(rawdataEst1_ids %in% toRemove))];
# # # rawdataEst2 = rawdataEst2[-(which(rawdataEst2_ids %in% toRemove))];
# # # rawdataEst1_ids = as.character(.getProjectIdFromFilenames(names(rawdataEst1)));
# # # rawdataEst2_ids = as.character(.getProjectIdFromFilenames(names(rawdataEst2)));

rawdataEst1Codes = .combineBehaviorCodesFromFishList(rawdataEst1);
rawdataEst2Codes = .combineBehaviorCodesFromFishList(rawdataEst2);

ids = ids0[ids0[,1] %in% c(rawdataEst1_ids, rawdataEst2_ids), ];
dat = dat0[dat0$FISH_ID %in% ids[, 2], ];
dat = cbind(ids, dat[match(ids[, 2], dat$FISH_ID), ]);
#unlist(sapply(rawdata,lapply,length)[1,])
kill0 = dat$'Project.ID'[dat$'Total.Kills'==0];
kill01 = dat$'Project.ID'[dat$'Total.Kills'<2];
killmult = dat$'Project.ID'[dat$'Total.Kills'>=2];
kill3plus = dat$'Project.ID'[dat$'Total.Kills'>2];

male_directed = c('a','b','d','e','q','r','j');
female_directed = c('c','i','s','t');
territorial = c('h','k','l','z');

# work on day 1
countsEst1 = .getBehaviorCountMatrixFromFishList(rawdataEst1)$cMat;	#removes '12_log032014_nov6_A4C_D_NC.txt', '48_log052014_Apr23_A8A_D_NC.txt'
percentsEst1 = .getBehaviorCountMatrixFromFishList(rawdataEst1)$pMat;	#removes '12_log032014_nov6_A4C_D_NC.txt', '48_log052014_Apr23_A8A_D_NC.txt'
countsEst1 = cbind(countsEst1, 
				   '12_log032014_nov6_A4C_D_NC.txt'=rep(0,nrow(countsEst1)), 
				   '48_log052014_Apr23_A8A_D_NC.txt'=rep(0,nrow(countsEst1))
				   );
countsEst1 = countsEst1[, match(rawdataEst1_ids, .getProjectIdFromFilenames(colnames(countsEst1)))];

countsEst1C = matrix(nrow=3, ncol=ncol(countsEst1), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(countsEst1))));
countsEst1C[1, ] = apply(countsEst1[rownames(countsEst1) %in% male_directed, ], 2, sum);
countsEst1C[2, ] = apply(countsEst1[rownames(countsEst1) %in% female_directed, ], 2, sum);
#countsEst1C[3, ] = apply(countsEst1[rownames(countsEst1) %in% territorial, ], 2, sum);
countsEst1C[3, ] = countsEst1[rownames(countsEst1) %in% territorial, ];

percentsEst1C = matrix(nrow=3, ncol=ncol(percentsEst1), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(percentsEst1))));
percentsEst1C[1, ] = apply(percentsEst1[rownames(percentsEst1) %in% male_directed, ], 2, sum);
percentsEst1C[2, ] = apply(percentsEst1[rownames(percentsEst1) %in% female_directed, ], 2, sum);
percentsEst1C[3, ] = percentsEst1[rownames(percentsEst1) %in% territorial, ];


TOPLOT = countsEst1;
#TOPLOT = TOPLOT / 15;
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
						  main=rawdataEst1Codes[match(rownames(TOPLOT)[i], rawdataEst1Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Establishment day 1', outer=T);

TOPLOT = apply(countsEst1, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors, total', names=xnames, col='grey', notch=F,
					  main='Establishment day 1: '
					  );
rm(TOPLOT, tofactor, xnames);

heatmap(cor(countsEst1[, apply(countsEst1, 2, sum)>0]), symm=T);

YLIM = c(0,135)
TOPLOT = countsEst1C;
#TOPLOT = TOPLOT / 15;
par(mfrow=c(1,4), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F, ylim=YLIM,
						  main=paste(rownames(TOPLOT)[i], ': ', sep='')
						  );
}; rm(i, TOPLOT, tofactor, xnames);
TOPLOT = apply(countsEst1C, 2, sum);
#TOPLOT = TOPLOT / 15;
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
					  main='total: ', ylim=YLIM
					  );
rm(TOPLOT, tofactor, xnames);
title(main='Establishment day 1', outer=T);

#######
# work on day 2

countsEst2 = .getBehaviorCountMatrixFromFishList(rawdataEst2)$cMat;

countsEst2C = matrix(nrow=3, ncol=ncol(countsEst2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(countsEst2))));
countsEst2C[1, ] = apply(countsEst2[rownames(countsEst2) %in% male_directed, ], 2, sum);
countsEst2C[2, ] = apply(countsEst2[rownames(countsEst2) %in% female_directed, ], 2, sum);
countsEst2C[3, ] = apply(countsEst2[rownames(countsEst2) %in% territorial, ], 2, sum);
#countsEst2C[3, ] = countsEst2[rownames(countsEst2) %in% territorial, ];

TOPLOT = countsEst2;
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
						  main=rawdataEst2Codes[match(rownames(TOPLOT)[i], rawdataEst2Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Establishment day 2', outer=T);

TOPLOT = apply(countsEst2, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors, total', names=xnames, col='grey', notch=F,
					  main='Establishment day 2: '
					  );
rm(TOPLOT, tofactor, xnames);

heatmap(cor(countsEst2[, apply(countsEst2, 2, sum)>0]), symm=T);

YLIM = c(0,175)
TOPLOT = countsEst2C;
par(mfrow=c(1,4), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F, ylim=YLIM,
						  main=paste(rownames(TOPLOT)[i], ': ', sep='')
						  );
}; rm(i, TOPLOT, tofactor, xnames);
TOPLOT = apply(countsEst2C, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
					  main='total: ', ylim=YLIM
					  );
rm(TOPLOT, tofactor, xnames);
title(main='Establishment day 2', outer=T);


# # # .sumRow = function (mat, rowname) {
	# # # return(sum(mat[rownames(mat)==rowname, ]));
# # # }
# # # C1 = 'kill01';
# # # C2 = 'killmult';
# # # BEH = 'c';
# # # X1 = sapply(rawdataEst1Combined[names(rawdataEst1Combined) %in% get(C1)], .sumRow, BEH);
# # # X2 = sapply(rawdataEst1Combined[names(rawdataEst1Combined) %in% get(C2)], .sumRow, BEH);
# # # WGCNA::verboseBoxplot(c(X1, X2), 
					  # # # c(rep(C1, length(X1)), rep(C2, length(X2))),
					  # # # xlab='', ylab=paste(rawdataEst1Codes[match(BEH, rawdataEst1Codes[, 1]), 2], ' (', BEH, ')', sep=''), 
					  # # # col='grey', frame.plot=F
					  # # # );

# # # par(mfrow=c(2,7))
# # # for (BEH in rawdataEst1Codes[, 1]) {
	# # # X1 = sapply(rawdataEst1Combined[names(rawdataEst1Combined) %in% get(C1)], .sumRow, BEH);
	# # # X2 = sapply(rawdataEst1Combined[names(rawdataEst1Combined) %in% get(C2)], .sumRow, BEH);
	# # # WGCNA::verboseBoxplot(c(X1, X2), 
					  # # # c(rep(C1, length(X1)), rep(C2, length(X2))),
					  # # # xlab='', ylab=paste(rawdataEst1Codes[match(BEH, rawdataEst1Codes[, 1]), 2], ' (', BEH, ')', sep=''), 
					  # # # col='grey', frame.plot=F
					  # # # );
# # # }


###################


######################################
### load intruder day data ###

### nov 25 with May 12, nov 26 with May 13

# combine day1's (nov25, may12) for intruder then day2's (nov26, may13)
dir = 'Nov2013Behavior';
filesNov25 = paste(dir, '/', list.files(dir)[grepl('.*nov25.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
filesNov26 = paste(dir, '/', list.files(dir)[grepl('.*nov26.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
dir = 'maintenanceIntruderMay12_13---remember to edit 54_log052814_May12_A8B_S_C';
filesMay12 = paste(dir, '/', list.files(dir)[grepl('.*may12.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');
filesMay13 = paste(dir, '/', list.files(dir)[grepl('.*may13.*_D_.*txt$', list.files(dir), ignore.case=T)], sep='');

files_intruder1 = c(filesNov25, filesMay12);
files_intruder2 = c(filesNov26, filesMay13);

rawdataInt1 = .getDataBatch(files_intruder1);
names(rawdataInt1) = gsub('Nov2013Behavior/|maintenanceIntruderMay12_13---remember to edit 54_log052814_May12_A8B_S_C/', '', names(rawdataInt1));
rawdataInt2 = .getDataBatch(files_intruder2);
names(rawdataInt2) = gsub('Nov2013Behavior/|maintenanceIntruderMay12_13---remember to edit 54_log052814_May12_A8B_S_C/', '', names(rawdataInt2));

rawdataInt1_ids = as.character(.getProjectIdFromFilenames(names(rawdataInt1)));
rawdataInt2_ids = as.character(.getProjectIdFromFilenames(names(rawdataInt2)));

# animals that were reused in multiple dyads
# # # toRemove = c('33','41','49','52');
# # # rawdataInt1 = rawdataInt1[-(which(rawdataInt1_ids %in% toRemove))];
# # # rawdataInt2 = rawdataInt2[-(which(rawdataInt2_ids %in% toRemove))];
# # # rawdataInt1_ids = as.character(.getProjectIdFromFilenames(names(rawdataInt1)));
# # # rawdataInt2_ids = as.character(.getProjectIdFromFilenames(names(rawdataInt2)));


rawdataInt1Codes = .combineBehaviorCodesFromFishList(rawdataInt1);
rawdataInt2Codes = .combineBehaviorCodesFromFishList(rawdataInt2);

ids = ids0[ids0[,1] %in% c(rawdataInt1_ids, rawdataInt2_ids), ];
dat = dat0[dat0$FISH_ID %in% ids[, 2], ];
dat = cbind(ids, dat[match(ids[, 2], dat$FISH_ID), ]);

kill0 = dat$'Project.ID'[dat$'Total.Kills'==0];
kill01 = dat$'Project.ID'[dat$'Total.Kills'<2];
killmult = dat$'Project.ID'[dat$'Total.Kills'>=2];
kill3plus = dat$'Project.ID'[dat$'Total.Kills'>2];

male_directed = c('a','b','d','e','q','r');
female_directed = c('c','i','s','t');
territorial = c('h','k','l','z');

# work on day 1
countsInt1 = .getBehaviorCountMatrixFromFishList(rawdataInt1)$cMat;	

countsInt1C = matrix(nrow=3, ncol=ncol(countsInt1), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(countsInt1))));
countsInt1C[1, ] = apply(countsInt1[rownames(countsInt1) %in% male_directed, ], 2, sum);
countsInt1C[2, ] = apply(countsInt1[rownames(countsInt1) %in% female_directed, ], 2, sum);
#countsInt1C[3, ] = apply(countsInt1[rownames(countsInt1) %in% territorial, ], 2, sum);
countsInt1C[3, ] = countsInt1[rownames(countsInt1) %in% territorial, ];


TOPLOT = countsInt1;
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
						  main=rawdataInt1Codes[match(rownames(TOPLOT)[i], rawdataInt1Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Intruder day 1', outer=T);

TOPLOT = apply(countsInt1, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors, total', names=xnames, col='grey', notch=F,
					  main='Intruder day 1: '
					  );
rm(TOPLOT, tofactor, xnames);

heatmap(cor(countsInt1[, apply(countsInt1, 2, sum)>0]), symm=T);

YLIM = c(0,350)
TOPLOT = countsInt1C;
par(mfrow=c(1,4), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F, ylim=YLIM,
						  main=paste(rownames(TOPLOT)[i], ': ', sep='')
						  );
}; rm(i, TOPLOT, tofactor, xnames);
TOPLOT = apply(countsInt1C, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
					  main='total: ', ylim=YLIM
					  );
rm(TOPLOT, tofactor, xnames);
title(main='Intruder day 1', outer=T);

# work on day 2
countsInt2 = .getBehaviorCountMatrixFromFishList(rawdataInt2)$cMat;	

countsInt2C = matrix(nrow=3, ncol=ncol(countsInt2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(countsInt2))));
countsInt2C[1, ] = apply(countsInt2[rownames(countsInt2) %in% male_directed, ], 2, sum);
countsInt2C[2, ] = apply(countsInt2[rownames(countsInt2) %in% female_directed, ], 2, sum);
countsInt2C[3, ] = apply(countsInt2[rownames(countsInt2) %in% territorial, ], 2, sum);
#countsInt2C[3, ] = countsInt2[rownames(countsInt2) %in% territorial, ];

TOPLOT = countsInt2;
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
						  main=rawdataInt2Codes[match(rownames(TOPLOT)[i], rawdataInt2Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Intruder day 2', outer=T);

YLIM = c(0,175)
TOPLOT = countsInt2C;
par(mfrow=c(1,4), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F, ylim=YLIM,
						  main=paste(rownames(TOPLOT)[i], ': ', sep='')
						  );
}; rm(i, TOPLOT, tofactor, xnames);
TOPLOT = apply(countsInt2C, 2, sum);
tofactor = .getProjectIdFromFilenames(names(TOPLOT)) %in% killmult;
xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
WGCNA::verboseBoxplot(TOPLOT, as.factor(tofactor),
					  xlab='', ylab='num behaviors', names=xnames, col='grey', notch=F,
					  main='total: ', ylim=YLIM
					  );
rm(TOPLOT, tofactor, xnames);
title(main='Intruder day 2', outer=T);


###################

# generate separate rawdata list and matrix of behavior codes for each day
# # days = c('nov6', 'nov7', 'nov25', 'nov26');
# # for (day in days) {
	# # name = paste(day, 'rawdata', sep='');
	# # assign(name, rawdata[grepl(day, names(rawdata), ignore.case=T)]);
	# # assign(paste(name, 'Codes', sep=''), .combineBehaviorCodesFromFishList(get(name)));
	# # assign(gsub('rawdata', 'behCounts', name), .getBehaviorCountMatrixFromFishList(get(name)));
	# # print(day);
	# # print(unlist(sapply(get(name), lapply, length)[1, ]));
# # }; rm(day, name);

days = c('rawdataEst1', 'rawdataEst2', 'rawdataInt1', 'rawdataInt2');
# loop through each day and build dataSummary lists
for (day in 1:length(days)) {
	print(days[day])
	# assign current day to variable
	#raw = get(paste(days[day], 'rawdata', sep=''));
	raw = get(days[day])
	tmp = raw;
	# loop through fish in current day
	for (f in 1:length(raw)) {
		print(names(raw)[f])
		if (length(raw[[f]]$data) < 3) {
			cat('skipping ', names(raw)[f], ' because <3 behaviors\n', sep='');
			next;
		} else {
			tmp[[f]] = .dataSummary(raw[[f]], plot=F);
		}
	}
	assign(paste(gsub('raw', '', days[day]), sep=''), tmp);
}; rm(day,raw,tmp,f);


###################

source('/Volumes/fishstudies-1/_code/bootstrapFunctions_6-16-13.R');


BEH = 'c';
YLIM = c(0,150)
allcodes = rbind(rawdataEst1Codes, rawdataEst2Codes, rawdataInt1Codes, rawdataInt2Codes)
YLAB = allcodes[allcodes[,1] == BEH, 2][[1]]

par(mfrow=c(2,2));
GP = kill01;
out = bootstrap.2paired(countsEst1[BEH, .getProjectIdFromFilenames(colnames(countsEst1)) %in% GP], 
 						countsEst2[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(countsEst1[BEH, .getProjectIdFromFilenames(colnames(countsEst1)) %in% GP], 
 						countsEst2[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM)
GP = kill01;
tmp = countsInt1[, .getProjectIdFromFilenames(colnames(countsInt1)) %in% .getProjectIdFromFilenames(colnames(countsInt2))];
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						countsInt2[BEH, .getProjectIdFromFilenames(colnames(countsInt2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						countsInt2[BEH, .getProjectIdFromFilenames(colnames(countsInt2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM);
 						
 						
 						
 						
 						
 						
 						
 						
 						
 						
BEH = 1;
YLIM = c(0,225)
YLAB = 'num male_directed'

par(mfrow=c(2,2));
GP = kill01;
out = bootstrap.2paired(countsEst1C[BEH, .getProjectIdFromFilenames(colnames(countsEst1)) %in% GP], 
 						countsEst2C[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(countsEst1C[BEH, .getProjectIdFromFilenames(colnames(countsEst1)) %in% GP], 
 						countsEst2C[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM)
GP = kill01;
tmp = countsInt1C[, .getProjectIdFromFilenames(colnames(countsInt1C)) %in% .getProjectIdFromFilenames(colnames(countsInt2C))];
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						countsInt2C[BEH, .getProjectIdFromFilenames(colnames(countsInt2C)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						countsInt2C[BEH, .getProjectIdFromFilenames(colnames(countsInt2C)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM)
 						
 						
 						
 						
 						
par(mfrow=c(2,3))						
for (GP in c('killmult','kill01')) {	
	for (BEH in 1:3) {				
		#BEH = 1;									 						
		#GP = 'killmult';
		if (BEH==1) {YLIM=c(0,250)}
		if (BEH==2) {YLIM=c(0,150)}
		if (BEH==3) {YLIM=c(0,30)}
		COL = 'lightgrey';
		TOPLOT=cbind(countsEst1C[BEH, .getProjectIdFromFilenames(colnames(countsEst1)) %in% get(GP)],
			 		 countsEst2C[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% get(GP)],
			 		 countsInt1C[BEH, .getProjectIdFromFilenames(colnames(countsInt1)) %in% get(GP)]
			 		 );
		boxplot(TOPLOT, boxcol=COL, col=COL, 
			    ylab=paste('# ', rownames(countsEst1C)[BEH], ' behaviors', sep=''), 
			    names=c('Establishment day 1','Establishment day 2','Intruder day 1'), 
			    frame.plot=F, 
			    main=paste(GP, ' (n=', nrow(TOPLOT), ')', sep=''), cex.axis=1.1, cex.lab=1.3, ylim=YLIM
			    );
		stripchart(data.frame(TOPLOT), vertical = T, add = T, pch = 21, cex = 1.5, bg='black');
		col = 'darkgrey'
		for (row in 1:nrow(TOPLOT)) {
			if (TOPLOT[row, 1] < TOPLOT[row, 2] & TOPLOT[row, 2] < TOPLOT[row, 3]) {
				col = 'red';
			} else {
				col = COL;
			}
			segments(1, TOPLOT[row, 1], 2, TOPLOT[row, 2], col = col)
		}
		for (row in 1:nrow(TOPLOT)) {
			if (TOPLOT[row, 1] < TOPLOT[row, 2] & TOPLOT[row, 2] < TOPLOT[row, 3]) {
				col = 'red';
			} else {
				col = COL;
			}
			segments(2, TOPLOT[row, 2], 3, TOPLOT[row, 3], col = col)
		}
	}
}




# establishment day2 vs intruder day1
BEH = 1;
YLIM = c(0,250)
YLAB = 'num male_directed'
par(mfrow=c(1,2))
GP = kill01;
out = bootstrap.2paired(countsEst2C[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP], 
 						countsInt1C[BEH, .getProjectIdFromFilenames(colnames(countsInt1)) %in% GP],
 						noDist=T,
 						conditionNames=c('Establishment day 2','Intruder day 1'),
 						main='kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM,col='grey', border='darkgrey');
 						
 GP = killmult;			
 tmp = countsInt1C[BEH, .getProjectIdFromFilenames(colnames(countsInt1)) %in% GP];	
 out = bootstrap.2paired(countsEst2C[BEH, .getProjectIdFromFilenames(colnames(countsEst2)) %in% GP], 
 						tmp[.getProjectIdFromFilenames(names(tmp)) %in% .getProjectIdFromFilenames(colnames(countsEst2))],
 						noDist=T,
 						conditionNames=c('Establishment day 2','Intruder day 1'),
 						main='kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM,col='grey', border='darkgrey');
##############

# duration

STATSTAT = 'median'

STAT=STATSTAT;
TOPLOT = .getBehaviorDurationMatrixFromFishList(rawdataEst1,STAT);
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab=paste(STAT, ' duration (frames)', sep=''), names=xnames, col='grey', notch=F,
						  main=rawdataEst1Codes[match(rownames(TOPLOT)[i], rawdataEst1Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Establishment day 1', outer=T);

STAT=STATSTAT;
TOPLOT = .getBehaviorDurationMatrixFromFishList(rawdataEst2,STAT);
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab=paste(STAT, ' duration (frames)', sep=''), names=xnames, col='grey', notch=F,
						  main=rawdataEst2Codes[match(rownames(TOPLOT)[i], rawdataEst2Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Establishment day 2', outer=T);

STAT=STATSTAT;
TOPLOT = .getBehaviorDurationMatrixFromFishList(rawdataInt1,STAT);
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab=paste(STAT, ' duration (frames)', sep=''), names=xnames, col='grey', notch=F,
						  main=rawdataInt1Codes[match(rownames(TOPLOT)[i], rawdataInt1Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Intruder day 1', outer=T);

STAT=STATSTAT;
TOPLOT = .getBehaviorDurationMatrixFromFishList(rawdataInt2,STAT);
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in 1:nrow(TOPLOT)) {
	tofactor = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
	xnames = c(paste('kill0/1 (n=', sum(tofactor==F), ')', sep=''), paste('kill>=2 (n=', sum(tofactor), ')', sep=''));
	WGCNA::verboseBoxplot(TOPLOT[i,], as.factor(tofactor),
						  xlab='', ylab=paste(STAT, ' duration (frames)', sep=''), names=xnames, col='grey', notch=F,
						  main=rawdataInt2Codes[match(rownames(TOPLOT)[i], rawdataInt2Codes[,1]),2]
						  );
}; rm(i, TOPLOT, tofactor, xnames);
title(main='Intruder day 2', outer=T);


#######
STAT = 'total';
durationsEst1 = .getBehaviorDurationMatrixFromFishList(rawdataEst1, stat=STAT);	

durationsEst1C = matrix(nrow=3, ncol=ncol(durationsEst1), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(durationsEst1))));
durationsEst1C[1, ] = apply(durationsEst1[rownames(durationsEst1) %in% male_directed, ], 2, sum);
durationsEst1C[2, ] = apply(durationsEst1[rownames(durationsEst1) %in% female_directed, ], 2, sum);
#durationsEst1C[3, ] = apply(durationsEst1[rownames(durationsEst1) %in% territorial, ], 2, sum);
durationsEst1C[3, ] = durationsEst1[rownames(durationsEst1) %in% territorial, ];


durationsEst2 = .getBehaviorDurationMatrixFromFishList(rawdataEst2, stat=STAT);

durationsEst2C = matrix(nrow=3, ncol=ncol(durationsEst2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(durationsEst2))));
durationsEst2C[1, ] = apply(durationsEst2[rownames(durationsEst2) %in% male_directed, ], 2, sum);
durationsEst2C[2, ] = apply(durationsEst2[rownames(durationsEst2) %in% female_directed, ], 2, sum);
durationsEst2C[3, ] = apply(durationsEst2[rownames(durationsEst2) %in% territorial, ], 2, sum);
#durationsEst2C[3, ] = durationsEst2[rownames(durationsEst2) %in% territorial, ];

durationsInt1 = .getBehaviorDurationMatrixFromFishList(rawdataInt1, stat=STAT);	

durationsInt1C = matrix(nrow=3, ncol=ncol(durationsInt1), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(durationsInt1))));
durationsInt1C[1, ] = apply(durationsInt1[rownames(durationsInt1) %in% male_directed, ], 2, sum);
durationsInt1C[2, ] = apply(durationsInt1[rownames(durationsInt1) %in% female_directed, ], 2, sum);
#durationsInt1C[3, ] = apply(durationsInt1[rownames(durationsInt1) %in% territorial, ], 2, sum);
durationsInt1C[3, ] = durationsInt1[rownames(durationsInt1) %in% territorial, ];

durationsInt2 = .getBehaviorDurationMatrixFromFishList(rawdataInt2, stat=STAT);	

durationsInt2C = matrix(nrow=3, ncol=ncol(durationsInt2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(durationsInt2))));
durationsInt2C[1, ] = apply(durationsInt2[rownames(durationsInt2) %in% male_directed, ], 2, sum);
durationsInt2C[2, ] = apply(durationsInt2[rownames(durationsInt2) %in% female_directed, ], 2, sum);
durationsInt2C[3, ] = apply(durationsInt2[rownames(durationsInt2) %in% territorial, ], 2, sum);
#durationsInt2C[3, ] = durationsInt2[rownames(durationsInt2) %in% territorial, ];

BEH = 1;
YLIM = c(0,7500)
YLAB = paste('duration (frames) ', rownames(durationsEst1C)[BEH], sep='');

par(mfrow=c(2,2));
GP = kill01;
out = bootstrap.2paired(durationsEst1C[BEH, .getProjectIdFromFilenames(colnames(durationsEst1)) %in% GP], 
 						durationsEst2C[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(durationsEst1C[BEH, .getProjectIdFromFilenames(colnames(durationsEst1)) %in% GP], 
 						durationsEst2C[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM)
GP = kill01;
tmp = durationsInt1C[, .getProjectIdFromFilenames(colnames(durationsInt1C)) %in% .getProjectIdFromFilenames(colnames(durationsInt2C))];
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						durationsInt2C[BEH, .getProjectIdFromFilenames(colnames(durationsInt2C)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM); 
GP = killmult;
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						durationsInt2C[BEH, .getProjectIdFromFilenames(colnames(durationsInt2C)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM)
 						
 						
par(mfrow=c(2,3))						
for (GP in c('killmult','kill01')) {	
	for (BEH in 1:3) {				
		#BEH = 1;									 						
		#GP = 'killmult';
		if (BEH==1 & GP=='killmult') {YLIM=c(0,25000)} else {YLIM=c(0,5000)}
		if (BEH==2) {YLIM=c(0,2500)}
		if (BEH==3) {YLIM=c(0,10000)}
		COL = 'lightgrey';
		TOPLOT=cbind(durationsEst1C[BEH, .getProjectIdFromFilenames(colnames(durationsEst1)) %in% get(GP)],
			 		 durationsEst2C[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% get(GP)],
			 		 durationsInt1C[BEH, .getProjectIdFromFilenames(colnames(durationsInt1)) %in% get(GP)]
			 		 );
		boxplot(TOPLOT, boxcol=COL, col=COL, 
			    ylab=paste('duration (frames) ', rownames(durationsEst1C)[BEH], ' behaviors', sep=''), 
			    names=c('Establishment day 1','Establishment day 2','Intruder day 1'), 
			    frame.plot=F, 
			    main=paste(GP, ' (n=', nrow(TOPLOT), ')', sep=''), cex.axis=1.1, cex.lab=1.3, ylim=YLIM
			    );
		stripchart(data.frame(TOPLOT), vertical = T, add = T, pch = 21, cex = 1.5, bg='black');
		col = 'darkgrey'
		for (row in 1:nrow(TOPLOT)) {
			if (TOPLOT[row, 1] < TOPLOT[row, 2] & TOPLOT[row, 2] < TOPLOT[row, 3]) {
				col = 'red';
			} else {
				col = COL;
			}
			segments(1, TOPLOT[row, 1], 2, TOPLOT[row, 2], col = col)
		}
		for (row in 1:nrow(TOPLOT)) {
			if (TOPLOT[row, 1] < TOPLOT[row, 2] & TOPLOT[row, 2] < TOPLOT[row, 3]) {
				col = 'red';
			} else {
				col = COL;
			}
			segments(2, TOPLOT[row, 2], 3, TOPLOT[row, 3], col = col)
		}
	}
}



# establishment day2 vs intruder day1
BEH = 1;
YLIM = c(0,25000)
YLAB = 'duration (frames) male_directed'
par(mfrow=c(1,2))
GP = kill01;
out = bootstrap.2paired(durationsEst2C[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP], 
 						durationsInt1C[BEH, .getProjectIdFromFilenames(colnames(durationsInt1)) %in% GP],
 						noDist=T,
 						conditionNames=c('Establishment day 2','Intruder day 1'),
 						main='kill0/1',
 						dataDescriptor=YLAB, ylim=YLIM,col='grey', border='darkgrey');
 						
 GP = killmult;			
 tmp = durationsInt1C[BEH, .getProjectIdFromFilenames(colnames(durationsInt1)) %in% GP];	
 out = bootstrap.2paired(durationsEst2C[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP], 
 						tmp[.getProjectIdFromFilenames(names(tmp)) %in% .getProjectIdFromFilenames(colnames(durationsEst2))],
 						noDist=T,
 						conditionNames=c('Establishment day 2','Intruder day 1'),
 						main='kill>=2',
 						dataDescriptor=YLAB, ylim=YLIM,col='grey', border='darkgrey');
 						
 						
 						
 						
 						
BEH = 'c';
YLIM = c(0,1500)
allcodes = rbind(rawdataEst1Codes, rawdataEst2Codes, rawdataInt1Codes, rawdataInt2Codes)
YLAB = allcodes[allcodes[,1] == BEH, 2][[1]]

par(mfrow=c(2,2));
GP = kill01;
out = bootstrap.2paired(durationsEst1[BEH, .getProjectIdFromFilenames(colnames(durationsEst1)) %in% GP], 
 						durationsEst2[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill0/1',
 						ylim=YLIM,
 						dataDescriptor=YLAB); 
GP = killmult;
out = bootstrap.2paired(durationsEst1[BEH, .getProjectIdFromFilenames(colnames(durationsEst1)) %in% GP], 
 						durationsEst2[BEH, .getProjectIdFromFilenames(colnames(durationsEst2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Establishment, kill>=2',
 						ylim=YLIM,
 						dataDescriptor=YLAB)
GP = kill01;
tmp = durationsInt1[, .getProjectIdFromFilenames(colnames(durationsInt1)) %in% .getProjectIdFromFilenames(colnames(durationsInt2))];
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						durationsInt2[BEH, .getProjectIdFromFilenames(colnames(durationsInt2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill0/1',
 						ylim=YLIM,
 						dataDescriptor=YLAB); 
GP = killmult;
out = bootstrap.2paired(tmp[BEH, .getProjectIdFromFilenames(colnames(tmp)) %in% GP], 
 						durationsInt2[BEH, .getProjectIdFromFilenames(colnames(durationsInt2)) %in% GP],
 						noDist=T,
 						conditionNames=c('day1','day2'),
 						main='Intruder, kill>=2',
 						ylim=YLIM,
 						dataDescriptor=YLAB);
 						
 						
###### get behavior count differences from Est.2 to Int.1

# limit to common animals and behaviors

keeprows = intersect(rownames(countsEst2), rownames(countsInt1));
keepcols = intersect(.getProjectIdFromFilenames(colnames(countsEst2)), .getProjectIdFromFilenames(colnames(countsInt1)));

tmpest2 = countsEst2[match(keeprows,rownames(countsEst2)), match(keepcols, .getProjectIdFromFilenames(colnames(countsEst2)))];
tmpint1 = countsInt1[match(keeprows,rownames(countsInt1)), match(keepcols, .getProjectIdFromFilenames(colnames(countsInt1)))];
diffcountsInt1_Est2 = tmpint1 - tmpest2;

par(mfrow=c(2,6), oma=c(0,0,2,0));
for(i in 1:nrow(diffcountsInt1_Est2)) {  
	beh = rawdataEst2Codes[match(rownames(diffcountsInt1_Est2)[i], rawdataEst2Codes[,1]), 2];
	WGCNA::verboseBoxplot(diffcountsInt1_Est2[i, ], 
					      .getProjectIdFromFilenames(colnames(diffcountsInt1_Est2)) %in% killmult, 
					      xlab='', 
					      ylab=paste(beh, sep=''), 
					      names=c('kill01 (n=10)','kill>=2 (n=9)'), col='grey', frame.plot=F, cex.axis=1.2		##### CAREFUL OF n's in names
					      );  
	abline(h=0, col='red')
}; rm(i)
title('Change in number of behaviors from Est.2 to Int.1', outer=T);


tmpest2 = durationsEst2[match(keeprows,rownames(durationsEst2)), match(keepcols, .getProjectIdFromFilenames(colnames(durationsEst2)))];
tmpint1 = durationsInt1[match(keeprows,rownames(durationsInt1)), match(keepcols, .getProjectIdFromFilenames(colnames(durationsInt1)))];
diffdurationsInt1_Est2 = tmpint1 - tmpest2;

#TOPLOT = diffdurationsInt1_Est2;
# if doing percent change
TOPLOT = diffcountsInt1_Est2 / tmpest2 * 100;
TOPLOT[!is.finite(TOPLOT)] = NaN;
TOPLOT = TOPLOT[, !(apply(TOPLOT,2,function(f) sum(is.nan(f))) == nrow(TOPLOT))];
# 
FACTOR = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
par(mfrow=c(2,6), oma=c(0,0,2,0));
for(i in 1:nrow(TOPLOT)) {  
	beh = rawdataEst2Codes[match(rownames(TOPLOT)[i], rawdataEst2Codes[,1]), 2];
	# check when doing %change
	check1 = length(unique(TOPLOT[i, FACTOR]));print(check1)
	check2 = length(unique(TOPLOT[i, !FACTOR]));print(check2)
	names = c(paste('kill01 (n=', check2, ')', sep=''), paste(paste('kill>=2 (n=', check1, ')', sep='')));
	if (check1==1 | check2==1) {next}
	#
	#print(TOPLOT[i,])
	WGCNA::verboseBoxplot(TOPLOT[i, ], FACTOR, 
					      xlab='', ylab=paste(beh, sep=''), 
					      names=names, col='grey', frame.plot=F, cex.axis=1	, notch=T	
					      );  
	abline(h=0, col='red')
}; rm(i)
title('% Change in number of behaviors from Est.2 to Int.1', outer=T);

######################

diffcountsInt1_Est2C = matrix(nrow=3, ncol=ncol(diffcountsInt1_Est2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(diffcountsInt1_Est2))));
diffcountsInt1_Est2C[1, ] = apply(diffcountsInt1_Est2[rownames(diffcountsInt1_Est2) %in% male_directed, ], 2, sum);
diffcountsInt1_Est2C[2, ] = apply(diffcountsInt1_Est2[rownames(diffcountsInt1_Est2) %in% female_directed, ], 2, sum);
#diffcountsInt1_Est2C[3, ] = apply(diffcountsInt1_Est2[rownames(diffcountsInt1_Est2) %in% territorial, ], 2, sum);
diffcountsInt1_Est2C[3, ] = diffcountsInt1_Est2[rownames(diffcountsInt1_Est2) %in% territorial, ];


diffdurationsInt1_Est2C = matrix(nrow=3, ncol=ncol(diffdurationsInt1_Est2), dimnames=list(c('male_directed','female_directed','territorial_neutral'), c(colnames(diffdurationsInt1_Est2))));
diffdurationsInt1_Est2C[1, ] = apply(diffdurationsInt1_Est2[rownames(diffdurationsInt1_Est2) %in% male_directed, ], 2, sum);
diffdurationsInt1_Est2C[2, ] = apply(diffdurationsInt1_Est2[rownames(diffdurationsInt1_Est2) %in% female_directed, ], 2, sum);
#diffdurationsInt1_Est2C[3, ] = apply(diffdurationsInt1_Est2[rownames(diffdurationsInt1_Est2) %in% territorial, ], 2, sum);
diffdurationsInt1_Est2C[3, ] = diffdurationsInt1_Est2[rownames(diffdurationsInt1_Est2) %in% territorial, ];

par(mfrow=c(1,3), oma=c(0,0,2,0));
for(i in 1:nrow(diffcountsInt1_Est2C)) {  
	#beh = rawdataEst2Codes[match(rownames(diffcountsInt1_Est2C)[i], rawdataEst2Codes[,1]), 2];
	WGCNA::verboseBoxplot(diffcountsInt1_Est2C[i, ], 
					      .getProjectIdFromFilenames(colnames(diffcountsInt1_Est2C)) %in% killmult, 
					      xlab='', 
					      ylab=rownames(diffcountsInt1_Est2C)[i], 
					      names=c('kill01 (n=10)','kill>=2 (n=9)'), col='grey', frame.plot=F, cex.axis=1.2		##### CAREFUL OF n's in names
					      );  
	abline(h=0, col='red')
}; rm(i)
title('Change in number of behaviors from Est.2 to Int.1', outer=T);



# if doing percent change
TOPLOT = diffdurationsInt1_Est2C / durationsEst2C * 100;
TOPLOT[!is.finite(TOPLOT)] = NaN;
TOPLOT = TOPLOT[, !(apply(TOPLOT,2,function(f) sum(is.nan(f))) == nrow(TOPLOT))];
# 
FACTOR = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
par(mfrow=c(1,3), oma=c(0,0,2,0));
for(i in 1:nrow(TOPLOT)) {  
	#beh = rawdataEst2Codes[match(rownames(TOPLOT)[i], rawdataEst2Codes[,1]), 2];
	# check when doing %change
	check1 = length(unique(TOPLOT[i, FACTOR]));print(check1)
	check2 = length(unique(TOPLOT[i, !FACTOR]));print(check2)
	names = c(paste('kill01 (n=', check2, ')', sep=''), paste(paste('kill>=2 (n=', check1, ')', sep='')));
	if (check1==1 | check2==1) {next}
	#
	#print(TOPLOT[i,])
	WGCNA::verboseBoxplot(TOPLOT[i, ], FACTOR, 
					      xlab='', ylab=rownames(TOPLOT)[i], 
					      names=names, col='grey', frame.plot=F, cex.axis=1	, notch=T	
					      );  
	abline(h=0, col='red')
}; rm(i)
title('% Change in total duration (frames) of behaviors from Est.2 to Int.1', outer=T);

#V2
TOPLOT = diffdurationsInt1_Est2C / durationsEst2C * 100;
TOPLOT[!is.finite(TOPLOT)] = NaN;
TOPLOT = TOPLOT[, !(apply(TOPLOT,2,function(f) sum(is.nan(f))) == nrow(TOPLOT))];
# 
FACTOR = .getProjectIdFromFilenames(colnames(TOPLOT)) %in% killmult;
par(mfrow=c(1,3), oma=c(0,2,0,0));
for(i in 1:nrow(TOPLOT)) {  
	#beh = rawdataEst2Codes[match(rownames(TOPLOT)[i], rawdataEst2Codes[,1]), 2];
	# check when doing %change
	check1 = length(unique(TOPLOT[i, FACTOR]));print(check1)
	check2 = length(unique(TOPLOT[i, !FACTOR]));print(check2)
	#names = c(paste('Level 1-2 (n=', check2, ')', sep=''), paste(paste('Level 3 (n=', check1, ')', sep='')));
	names = c('Level 1-2', 'Level 3');
	if (check1==1 | check2==1) {next}
	#
	#print(TOPLOT[i,])
	if (rownames(TOPLOT)[i] =='male_directed') {MAIN = 'Male-directed: '}
	if (rownames(TOPLOT)[i] =='female_directed') {MAIN = 'Female-directed: '}
	if (rownames(TOPLOT)[i] =='territorial_neutral') {MAIN = 'Territorial/Neutral: '}
	if (i==1) {
		YLAB = '% change total duration of behavior\nfrom Est.2 to Int.1';
		YLIM = c(0,12500);
		par(mai=c(.7, .8, .4, .3));
	} else {
		YLAB='';
		YLIM = c(-200,800);
		par(mai=c(.7, .25, .4, .3));
	}
	WGCNA::verboseBoxplot(TOPLOT[i, ], FACTOR, 
					      xlab='', ylab=YLAB, ylim=YLIM,
					      names=names, col='grey', frame.plot=F, cex.axis=1.3, cex.lab=1.4, notch=F, main=MAIN	
					      );  
	abline(h=0, col='red')
}; rm(i)
#title('% Change in total duration (frames) of behaviors from Est.2 to Int.1', outer=T);


##################
FACTOR = .getProjectIdFromFilenames(colnames(durationsEst1C)) %in% killmult;
COLS = FACTOR;
COLS[COLS] = 'blue';
COLS[COLS=='FALSE'] = 'darkgrey';

par(mfrow=c(3,3))
for (i in 1:nrow(durationsEst1C)) {
	WGCNA::verboseScatterplot(durationsEst1C[i, ], diffdurationsInt1_Est2C[i, ], 
		   abline=T, abline.col='red',
		   xlab='total duration of behaviors, Est.1', ylab='change from Est2 to Int1', 
		   main=paste(rownames(durationsEst1C)[i], '\n', sep=''), 
		   frame.plot=F, 
		   col='black',bg=COLS,
		   pch=21,cex=2, cex.axis=1
		   )
}
for (i in 1:nrow(durationsEst1C)) {
	WGCNA::verboseScatterplot(durationsEst1C[i, FACTOR], diffdurationsInt1_Est2C[i, FACTOR], 
		   abline=T, abline.col='red',
		   xlab='total duration of behaviors, Est.1', ylab='change from Est2 to Int1', 
		   main='kill>=2', 
		   frame.plot=F, 
		   col='black',bg=COLS[FACTOR],
		   pch=21,cex=2, cex.axis=1
		   )
}
for (i in 1:nrow(durationsEst1C)) {
	WGCNA::verboseScatterplot(durationsEst1C[i, !FACTOR], diffdurationsInt1_Est2C[i, !FACTOR], 
		   abline=T, abline.col='red',
		   xlab='total duration of behaviors, Est.1', ylab='change from Est2 to Int1', 
		   main='kill01', 
		   frame.plot=F, 
		   col='black',bg=COLS[!FACTOR],
		   pch=21,cex=2, cex.axis=1
		   )
}