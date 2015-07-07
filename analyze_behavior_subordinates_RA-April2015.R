### see also:
##### analyze_killbehavior_hormones_RA-July31.R
##### analyze_behavior_RA-July31.R
rm(list=ls());
source('/Volumes/fishstudies/_code/ethograms_from_scorevideo.R');
load("/Volumes/fishstudies/_behaviorRA/WORKSPACE_analyze_behavior_RA-Aug28.RData");
dir = '/Volumes/fishstudies/_behaviorRA/Behavior of subordinates';
setwd(dir);
allfiles = list.files(dir);

# get subordinate files for specific day, e.g. establishment day 2
doms = colnames(countsEst2);																	#doms = colnames(countsInt1);
date_tank = unlist(lapply(strsplit(doms,'_'), function(f) paste(f[3:4],collapse='_')));
# be careful of filenames with different formatting for date or tank, e.g. 'Apr_24' vs Apr24
inds = c();
for (i in date_tank) {
	inds = c(inds, grep(i, allfiles, ignore.case=T));
}; rm(i)
files = allfiles[inds];
subs = files[grep('_S_', files)];

date_tank_subs = unlist(lapply(strsplit(subs,'_'), function(f) paste(f[3:4],collapse='_')));
inds = c();
for (i in date_tank_subs) {
	inds = c(inds, grep(i, doms, ignore.case=T));
}; rm(i)
dyads = cbind(doms[inds],subs);

dyads = cbind(doms,subs);

rawdataEst2sub = .getDataBatch(subs);
rawdataEst2sub_ids = as.character(.getProjectIdFromFilenames(names(rawdataEst2sub)));
rawdataEst2subCodes = .combineBehaviorCodesFromFishList(rawdataEst2sub);

# rawdataInt1sub = .getDataBatch(subs);
# rawdataInt1sub_ids = as.character(.getProjectIdFromFilenames(names(rawdataInt1sub)));
# rawdataInt1subCodes = .combineBehaviorCodesFromFishList(rawdataInt1sub);


countsEst2sub = .getBehaviorCountMatrixFromFishList(rawdataEst2sub)$cMat;
# countsInt1sub = .getBehaviorCountMatrixFromFishList(rawdataInt1sub)$cMat;
# countsInt1 = countsInt1[, colnames(countsInt1) %in% dyads[,1]]

countsEst2id = as.numeric(unlist(lapply(strsplit(colnames(countsEst2), '_'), function(f) f[1])));
# countsInt1id = as.numeric(unlist(lapply(strsplit(colnames(countsInt1), '_'), function(f) f[1])));

par(mfrow=c(2,4),oma=c(0,0,2,0));
for (row in 1:nrow(countsEst2sub)) {
	WGCNA::verboseBoxplot(c(countsEst2sub[row, which(countsEst2id %in% kill01)], 
							countsEst2sub[row, which(countsEst2id %in% killmult)]), 
						  c(rep('kill 0/1',10),rep('kill > 1',9)),
						  ylab=rawdataEst2subCodes[row,2],
						  xlab='',frame.plot=F, col='lightgrey'
						  );
}; rm(row);
title('subordinate males, separated by group of their dominant dyad partner',outer=T);



par(mfrow=c(2,4),oma=c(0,0,2,0));
for (row in 1:nrow(countsInt1sub)) {
	WGCNA::verboseBoxplot(c(countsInt1sub[row, which(countsInt1id %in% kill01)], 
							countsInt1sub[row, which(countsInt1id %in% killmult)]), 
						  c(rep('kill 0/1',7),rep('kill > 1',8)),
						  ylab=rawdataInt1subCodes[row,2],
						  xlab='',frame.plot=F, col='lightgrey'
						  );
}; rm(row);

##########################################################

#### get latency to first attack for doms
.getFirstInstancesOfAllBehaviors = function (rawdataListEntry) {
	sData = split(rawdataListEntry$data, rawdataListEntry$data$behavior);
	firsts = as.data.frame(matrix(nrow=length(sData), ncol=2));
	for (b in 1:length(sData)) {
		tmp = sData[[b]];
		tmp = tmp[order(as.numeric(tmp$start_frame)), ];
		firsts[b, ] = c(tmp$behavior[1], tmp$start_frame[1]);
	}
	names(firsts) = c('behavior', 'start_frame');
	return(list(sData=sData, firsts=firsts));
}

.getFirstInstancesOfAllBehaviorsList = function (rawdataList) {
	firsts = list();
	for (s in 1:length(rawdataList)) {
		subname = names(rawdataList)[s];
		tmp = rawdataList[[s]];
		if (nrow(tmp$data) > 0) {
			firsts[[subname]] = .getFirstInstancesOfAllBehaviors(tmp)$firsts;
		} else {
			warning(paste('Skipping ', subname, ' because no data', sep=''));
			next;
		}
	}
	return(firsts)
}

.combineFirstsList = function (getFirstInstancesOfAllBehaviorsListOutput) {
	x = getFirstInstancesOfAllBehaviorsListOutput;
	behs = unique(unlist(lapply(x, function(f) f$behavior)));
	mat = as.data.frame(matrix(nrow=length(behs), ncol=length(x), dimnames=list(behs, names(x))));
	for (s in 1:length(x)) {
		sfirsts = x[[s]][match(rownames(mat), x[[s]]$behavior), ];
		mat[, s] = sfirsts$start_frame;
	}
	return(mat)
}

firstsEst1 = .combineFirstsList(.getFirstInstancesOfAllBehaviorsList(rawdataEst1));
firstsEst2 = .combineFirstsList(.getFirstInstancesOfAllBehaviorsList(rawdataEst2));

names(firstsEst1) = .getProjectIdFromFilenames(names(firstsEst1));
names(firstsEst2) = .getProjectIdFromFilenames(names(firstsEst2));

x = dat[dat$Project.ID %in% killmult, ]$Total.Kills;
names(x) = rownames(dat[dat$Project.ID %in% killmult, ]);

y = 





firstsEst1_01 = firstsEst1[, .getProjectIdFromFilenames(names(firstsEst1)) %in% kill01];
firstsEst1_mult = firstsEst1[, .getProjectIdFromFilenames(names(firstsEst1)) %in% killmult];




#FIRSTS = firstsEst1[match(aggEst1[,1], rownames(firstsEst1)), ];
FIRSTS = firstsEst1;
CODES = rawdataEst1Codes;
CEX = 1.5;
par(mfrow=c(3,5));
for (i in 1:nrow(FIRSTS)) {
	main = paste(CODES[match(rownames(FIRSTS)[i], CODES[,1]), 2],'\n',sep='');
	f01 = .getProjectIdFromFilenames(names(FIRSTS)) %in% kill01;#print(f01)
	na1 = sum(!is.na(as.numeric(FIRSTS[i, f01])));#print(na1)
	na2 = sum(!is.na(as.numeric(FIRSTS[i, !f01])));#print(na2)
	if (na1<2) {
		boxplot(as.numeric(FIRSTS[i, !f01]), main=main, xlab=paste('kill>1 (n=',na2,')',sep=''), cex.lab=CEX, cex.main=CEX, frame.plot=F, col='lightgrey');
	} else if (na2<2) {
		boxplot(as.numeric(FIRSTS[i, f01]), main=main, xlab=paste('kill<=1 (n=',na1,')',sep=''), cex.lab=CEX, cex.main=CEX, frame.plot=F, col='lightgrey');
	} else {
		verboseBoxplot(as.numeric(FIRSTS[i,]), f01, 
					   frame.plot=F, 
					   xlab='', ylab='frame #',
					   names=c(paste('kill>1 (n=',na2,')',sep=''),paste('kill<=1 (n=',na1,')',sep='')), 
					   col='lightgrey',
					   main=main, notch=F
					   );
	}
}; rm(FIRSTS, i);









FIRSTS = firstsEst2;
CODES = rawdataEst2Codes;
CEX = 1.5;
par(mfrow=c(2,3),oma=c(0,0,2,0));
for (i in 1:nrow(FIRSTS)) {
	main = paste(CODES[match(rownames(FIRSTS)[i], CODES[,1]), 2],'\n',sep='');
	f01 = .getProjectIdFromFilenames(names(FIRSTS)) %in% kill01;#print(f01)
	na1 = sum(!is.na(as.numeric(FIRSTS[i, f01])));#print(na1)
	na2 = sum(!is.na(as.numeric(FIRSTS[i, !f01])));#print(na2)
	if (na1<2 || na2<2){
		next
	} else {
		verboseBoxplot(as.numeric(FIRSTS[i,]), f01, 
					   frame.plot=F, 
					   xlab='', ylab='frame #',
					   names=c(paste('kill>1 (n=',na2,')',sep=''),paste('kill<=1 (n=',na1,')',sep='')), 
					   col='lightgrey',
					   main=main, notch=F
					   );
	}

}; rm(FIRSTS, i);
#title(main='Est.2', outer=T)












################################

dat2 = dat[,c(5,6,10,11,12,15,23:25,31:35)];
####################################

tmp = firstsEst2[, colnames(firstsEst2) %in% colnames(firstsEst1)];

tmpc1 = countsEst1[, .getProjectIdFromFilenames(colnames(countsEst1)) %in% colnames(firstsEst1)];
colnames(tmpc1) = .getProjectIdFromFilenames(colnames(tmpc1))
tmpc2 = countsEst2[, .getProjectIdFromFilenames(colnames(countsEst2)) %in% colnames(firstsEst1)];
colnames(tmpc2) = .getProjectIdFromFilenames(colnames(tmpc2))

forlm = as.data.frame(cbind(as.numeric(firstsEst1[rownames(firstsEst1)=='d', ]),
							as.numeric(tmp[rownames(tmp)=='s', ]),
							as.numeric(tmpc1[rownames(tmpc1)=='d', ]),
							as.numeric(tmpc2[rownames(tmpc2)=='s', ]),
							as.numeric(colnames(tmp) %in% kill01)
							)
					  );


names(forlm) = c('est1_chase_l', 'est2_quiver_l', 'est1_chase_c', 'est2_quiver_c', 'kill01');
forlm$est1_chase = as.numeric(forlm$est1_chase);
forlm$est2_quiver = as.numeric(forlm$est2_quiver);
forlm$kill01 = as.numeric(forlm$kill01);

#lm_kill01 = lm(kill01 ~ est1_chase_l + est2_quiver_l + est1_chase_c + est2_quiver_c, data=forlm);
lm_kill01 = lm(kill01 ~ est1_chase_c + est2_quiver_c, data=forlm);

# lm_kill01.est = lm_kill01[[1]][[1]] + lm_kill01[[1]][[2]]*forlm$est1_chase_l + lm_kill01[[1]][[3]]*forlm$est2_quiver_l + lm_kill01[[1]][[4]]*forlm$est1_chase_c + lm_kill01[[1]][[5]]*forlm$est2_quiver_c;
lm_kill01.est = lm_kill01[[1]][[1]] + lm_kill01[[1]][[2]]*forlm$est1_chase_c + lm_kill01[[1]][[3]]*forlm$est2_quiver_c;


