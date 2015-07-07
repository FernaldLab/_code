rm(list=ls());
setwd('/Volumes/fishstudies/_behaviorRA');
library(WGCNA);
source('/Volumes/fishstudies/_code/bootstrapFunctions_6-16-13.R');

dat0 = read.csv('Database 06-5-14-1_sheet-AggressionDmales.csv');
dat = dat0;
rownames(dat) = dat$'Fish..';
dat = dat[, -c(3,4,6,8,11,14,16)];
dat$'X.females.killed'[is.na(dat$'X.females.killed')] = 0;
dat$'X.males.killed'[is.na(dat$'X.males.killed')] = 0;
dat$'Total.Kills'[is.na(dat$'Total.Kills')] = 0;

# all killers vs non-killers
#datKill = dat[!(is.na(dat$'Total.Kills')), ];

par(mfrow=c(1,2));
verboseBoxplot(dat$'Age.at.dyad.creation', as.factor(dat$'Total.Kills'), col='grey', xlab='',ylab='',main='Age.at.dyad.creation',ylim=c(120,260));
verboseBoxplot(dat$'Age.at.dissection', as.factor(dat$'Total.Kills'), col='grey', xlab='',ylab='',main='Age.at.dyad.dissection',ylim=c(120,260));



datH0 = read.csv('AggressionSummarySteroids061814.csv');
datH0[,10] = as.numeric(datH0[,10])
par(mfrow=c(2,3),oma=c(0,0,2,0))
for( col in 7:10 ) {
	verboseBoxplot(datH0[,col],as.factor(datH0$Aggression),xlab='',ylab='',main=names(datH0)[col],col='grey')
}
for( col in 9:10 ) {
	verboseBoxplot(datH0[,col],as.factor(datH0$timepoint),xlab='',ylab='',main=names(datH0)[col],col='grey')
}
title(main='AggressionSummarySteroids061814', outer=T)

datH = datH0[-60, ];		# duplicate FISH.ID?
rownames(datH) = datH$'Fish.ID';
datH = datH[, -c(1,2)];

tmp = dat[rownames(dat) %in% rownames(datH), ];
tmp = tmp[match(rownames(datH), rownames(tmp)), ];
DAT = cbind(tmp, datH); rm(tmp);

DAT.PKT = DAT[DAT$Status=='PKT', ];
DAT.KT = DAT[DAT$Status=='KT', ];
DAT.AT = DAT[DAT$Status=='AT', ];
DAT.T = DAT[DAT$Status=='T', ];

par(mfrow=c(1,3))
verboseBoxplot(c(DAT.PKT$'Days.in.dyad.until.1st.kill',DAT.PKT$'latency.2nd.Kill',DAT.PKT$'latency.3rd.kill'),
			   as.factor(c(rep('1st',nrow(DAT.PKT)),rep('2nd',nrow(DAT.PKT)),rep('3rd',nrow(DAT.PKT)))),
			   #numberStandardErrors=2,
			   ylim=c(0,50),
			   ylab='latency to kill (days)', xlab='kill #', 
			   col='grey',
			   notch=F,
			   main='PKT: ');
verboseBoxplot(c(DAT.KT$'Days.in.dyad.until.1st.kill',DAT.KT$'latency.2nd.Kill',DAT.KT$'latency.3rd.kill'),
			   as.factor(c(rep('1st',nrow(DAT.KT)),rep('2nd',nrow(DAT.KT)),rep('3rd',nrow(DAT.KT)))),
			   #numberStandardErrors=2,
			   ylim=c(0,50),
			   #ylab='latency to kill (days)', xlab='kill #', 
			   ylab='', xlab='',
			   col='grey',
			   notch=F,
			   main='KT: ');
verboseBoxplot(c(DAT.AT$'Days.in.dyad.until.1st.kill',DAT.AT$'latency.2nd.Kill',DAT.AT$'latency.3rd.kill'),
			   as.factor(c(rep('1st',nrow(DAT.AT)),rep('2nd',nrow(DAT.AT)),rep('3rd',nrow(DAT.AT)))),
			   #numberStandardErrors=2,
			   ylim=c(0,50),
			   #ylab='latency to kill (days)', xlab='kill #', 
			   ylab='', xlab='',
			   col='grey',
			   notch=F,
			   main='AT: ');


COL=7;
tmp = cbind(DAT.PKT[, COL], c(DAT.KT[,COL],rep(NA,4)), c(DAT.AT[,COL],rep(NA,13)));
x=bootstrap.ANOVA(tmp, groupNames=c('PKT','KT','AT'), dataDescriptor='latency to 1st kill (days)');

COL=9;
tmp = cbind(DAT.PKT[, COL], c(DAT.KT[,COL],rep(NA,4)), c(DAT.AT[,COL],rep(NA,13)));
x=bootstrap.ANOVA(tmp, groupNames=c('PKT','KT','AT'), dataDescriptor='latency to 2nd kill (days)');

COL=10;
tmp = cbind(DAT.PKT[, COL], c(DAT.KT[,COL],rep(NA,4)), c(DAT.AT[,COL],rep(NA,13)));
x=bootstrap.ANOVA(tmp, groupNames=c('PKT','KT','AT'), dataDescriptor='latency to 3rd kill (days)');


# OR
par(mfrow=c(1,3))
COL=7;
tmp = c(DAT.PKT[, COL], DAT.KT[, COL], DAT.AT[, COL]);
grp = c(rep('PKT',nrow(DAT.PKT)), rep('KT',nrow(DAT.KT)), rep('AT',nrow(DAT.AT)));
verboseBoxplot(tmp, grp, col='grey', notch=F, xlab='', ylab='latency to 1st kill (days)',ylim=c(0,50));
COL=9;
tmp = c(DAT.PKT[, COL], DAT.KT[, COL], DAT.AT[, COL]);
grp = c(rep('PKT',nrow(DAT.PKT)), rep('KT',nrow(DAT.KT)), rep('AT',nrow(DAT.AT)));
verboseBoxplot(tmp, grp, col='grey', notch=F, xlab='', ylab='latency to 2nd kill (days)',ylim=c(0,50));
COL=10;
tmp = c(DAT.PKT[, COL], DAT.KT[, COL], DAT.AT[, COL]);
grp = c(rep('PKT',nrow(DAT.PKT)), rep('KT',nrow(DAT.KT)), rep('AT',nrow(DAT.AT)));
verboseBoxplot(tmp, grp, col='grey', notch=F, xlab='', ylab='latency to 3rd kill (days)',ylim=c(0,50));





X = DAT
cors = corAndPvalue(X[,c(3,4,6:15,18:21),],use='p')
par(oma=c(7,0,0,7)); heatmap(cors$cor ,symm=T);		

# estradiol looks to have some significant correlations

# plot across all fish
COL = 20;
TOPLOT = DAT;
par(mfrow=c(3,5), oma=c(0,0,2,0));
for (i in c(3,4,6:15,18,19,21)) {
	if (sum(!is.na(TOPLOT[, i]),na.rm=T)<3 | sum(TOPLOT[, i]!=0,na.rm=T)<3) {
		next;
	}
	verboseScatterplot(TOPLOT[, COL], TOPLOT[, i], 
					   abline=T, abline.lty='dashed', abline.col='red', frame.plot=F, 
					   ylab=names(TOPLOT)[i], xlab=names(TOPLOT)[COL], 
					   pch=21, bg='grey', cex=1.2
					   );
}; rm(i);
title(paste('(n=', sum(!is.na(TOPLOT[,COL])), ')', sep=''), outer=T);

# plot across different statuses
COL = 20;
STATUS = 'PKT'
TOPLOT = DAT[DAT$Status==STATUS, ];
par(mfrow=c(4,5), oma=c(0,0,2,0));
for (i in c(3,4,6:15,18,19,21:23)) {
	if (sum(!is.na(TOPLOT[, i]))<3 | length(unique(TOPLOT[, i]))==1) {
		next;
	}
	verboseScatterplot(TOPLOT[, COL], TOPLOT[, i], 
					   abline=T, abline.lty='dashed', abline.col='red', frame.plot=F, 
					   ylab=names(TOPLOT)[i], xlab=names(TOPLOT)[COL], 
					   pch=21, bg='grey', cex=1.2,
					   corOptions="use='p'"
					   );
}; rm(i);
title(paste(STATUS, ' (n=', sum(!is.na(TOPLOT[,COL])), ')', sep=''), outer=T);

COL1 = 20; COL2 = 21;
par(mfrow=c(1,4))
for (STATUS in c('PKT','KT','AT','T')) {
	TOPLOT = get(paste('DAT.',STATUS,sep=''));
	if (STATUS=='PKT') {
		YLAB = names(TOPLOT)[COL2];
		XLAB = names(TOPLOT)[COL1];
	} else {
		YLAB = ''; XLAB = '';
	}
	verboseScatterplot(TOPLOT[, COL1], TOPLOT[, COL2], 
				   abline=T, abline.lty='dashed', abline.col='red', frame.plot=F, 
				   ylab=YLAB, xlab=XLAB, 
				   pch=21, bg='grey', cex=1.2,
				   corOptions="use='p'",
				   ylim=c(0,250000), xlim=c(0,400000),
				   main=paste(STATUS,': ',sep='')
				   );
}


dat0 = read.csv('Experience_Aggression_hormones061914_sheet-database.csv');

#dat = dat0[!is.na(dat0$'Testo.pg.mL'), ];
# or
dat = dat0[!is.na(dat0$'Estradiolpg.mL'), ];
dat$'Cortisol.pg.mL'[dat$'Cortisol.pg.mL'=='TOO LOW'] = NA;

dat = dat[match(rownames(DAT),dat[,2]),]
DAT = cbind(DAT, length=dat$'body.length.mm.', mass=dat$'body.mass.g.')
DAT.PKT = DAT[DAT$Status=='PKT', ];
DAT.KT = DAT[DAT$Status=='KT', ];
DAT.AT = DAT[DAT$Status=='AT', ];
DAT.T = DAT[DAT$Status=='T', ];










.corHeatmapWithPvalLabels = function(data, mar=c(8,9,3,3), col=blueWhiteRed(50), main='', cex.text=.5, thresh=.05/(ncol(data)^2/2-(ncol(data)/2)), ...) {
	tmp = corAndPvalue(data, ...);
	cors = tmp$cor; 
	pvals = tmp$p;
	textMat = paste(signif(cors, 2), '\n(', signif(pvals, 1), ')', sep='');
	dim(textMat) = dim(cors);
	textMat[pvals>thresh] = '';
	par(mar=mar);
	labeledHeatmap(Matrix=cors, xLabels=names(data), yLabels=names(data), ySymbols=names(data), colorLabels=F, colors=col, textMatrix=textMat, setStdMargins=F, cex.text=cex.text, zlim=c(-1, 1), main=main, ...)
}

.corHeatmapWithPvalLabels(DAT[,c(3,4,6:15,18:21),],main='all fish');
.corHeatmapWithPvalLabels(DAT[DAT$Status=='PKT',c(3,4,6:15,18:21),],main='PKT');
.corHeatmapWithPvalLabels(DAT[DAT$Status=='KT',c(3,4,6:15,18:21),],main='KT');
.corHeatmapWithPvalLabels(DAT[DAT$Status=='T',c(3,4,6:15,18:21),],main='T');

# # dat = read.csv('Database 06-5-14-1_sheet-AggressionDmales.csv');
# # for(i in c(5,7,10,12,13,15,17,18,19,20)){dat[,i] = as.numeric(dat[,i])}; rm(i);
# # dat2 = read.csv('Experience_Aggression_hormones061914_sheet-database.csv');
# # dat2$'Fish.ID' = as.numeric(dat2$'Fish.ID');

# # dat2f = dat2[!is.na(dat2$'Fish.ID'), ];
# # datf = dat[dat$'Fish..' %in% dat2f$'Fish.ID', ];
# # datf = datf[match(dat2f$'Fish.ID', datf$'Fish..'), ];
# # datf = datf[!is.na(datf$Run), ];


dat0 = read.csv('Experience_Aggression_hormones061914_sheet-database.csv');

#dat = dat0[!is.na(dat0$'Testo.pg.mL'), ];
# or
dat = dat0[!is.na(dat0$'Estradiolpg.mL'), ];
dat$'Cortisol.pg.mL'[dat$'Cortisol.pg.mL'=='TOO LOW'] = NA;

par(mfrow=c(1,4),oma=c(0,0,2,0));
for(m in 4:7) {
	verboseBoxplot(as.numeric(dat[,m]), as.factor(dat$Status), xlab='', ylab='', main=names(dat)[m], col='grey')
}; rm(m);
title(main='Experience_Aggression_hormones061914_sheet-database',outer=T)

par(mfrow=c(2,1),oma=c(0,0,2,0));
for(m in 6:7) {
	verboseBoxplot(as.numeric(dat[,m]), as.factor(paste(dat$Status, '-', dat$'plasma.collected')), xlab='', ylab='', main=names(dat)[m], col='grey',notch=F)
}; rm(m);
title(main='Experience_Aggression_hormones061914_sheet-database',outer=T)


par(mfrow=c(2,2),oma=c(0,0,2,0))
notch=F;
m=6;
verboseBoxplot(as.numeric(dat[dat$Status=='T',m]), as.factor(dat$'plasma.collected'[dat$Status=='T']), xlab='', ylab='', main=paste('T-',names(dat)[m]), col='grey',notch=notch);
verboseBoxplot(as.numeric(dat[dat$Status=='PKT',m]), as.factor(dat$'plasma.collected'[dat$Status=='PKT']), xlab='', ylab='', main=paste('PKT-',names(dat)[m]), col='grey',notch=notch);
m=7;
verboseBoxplot(as.numeric(dat[dat$Status=='T',m]), as.factor(dat$'plasma.collected'[dat$Status=='T']), xlab='', ylab='', main=paste('T-',names(dat)[m]), col='grey',notch=notch);
verboseBoxplot(as.numeric(dat[dat$Status=='PKT',m]), as.factor(dat$'plasma.collected'[dat$Status=='PKT']), xlab='', ylab='', main=paste('PKT-',names(dat)[m]), col='grey',notch=notch);
title(main='Experience_Aggression_hormones061914_sheet-database',outer=T)

verboseBoxplot(as.numeric(dat[,m]), as.factor(dat$'plasma.collected'), xlab='', ylab='', main=paste('PKT-',names(dat)[m]), col='grey');


par(mfrow=c(2,2),oma=c(0,0,2,0))
notch=F;
m=6;
verboseBoxplot(as.numeric(dat[dat$'plasma.collected'=='NOSTRESS',m]), as.factor(dat$Status[dat$'plasma.collected'=='NOSTRESS']), xlab='', ylab='', main=paste('NO STRESS-',names(dat)[m]), col='grey',notch=notch);
verboseBoxplot(as.numeric(dat[dat$'plasma.collected'=='DURING_ATTACK',m]), as.factor(dat$Status[dat$'plasma.collected'=='DURING_ATTACK']), xlab='', ylab='', main=paste('DURING_ATTACK-',names(dat)[m]), col='grey',notch=notch);
m=7;
verboseBoxplot(as.numeric(dat[dat$'plasma.collected'=='NOSTRESS',m]), as.factor(dat$Status[dat$'plasma.collected'=='NOSTRESS']), xlab='', ylab='', main=paste('NO STRESS-',names(dat)[m]), col='grey',notch=notch);
verboseBoxplot(as.numeric(dat[dat$'plasma.collected'=='DURING_ATTACK',m]), as.factor(dat$Status[dat$'plasma.collected'=='DURING_ATTACK']), xlab='', ylab='', main=paste('DURING_ATTACK-',names(dat)[m]), col='grey',notch=notch);
title(main='Experience_Aggression_hormones061914_sheet-database',outer=T)

# # par(mfrow=c(3,4))
# # for (m in 4:7) {
	# # for (n in 4:7) {
		# # if(m==n){next}
		# # verboseBoxplot(as.numeric(c(dat[,m],dat[,n])),
						# # as.factor(c( rep(names(dat)[m], nrow(dat)), rep(names(dat)[n], nrow(dat)) )),
						# # xlab='',ylab='',col='grey')
	# # }	
# # }; rm(m,n);



######################

bigDB = read.csv('Database 06-5-14-1_NewDB.csv',header=T);
bigDB = bigDB[bigDB$'FISH_ID.1' %in% rownames(DAT), ];
bigDB = bigDB[bigDB$GSI!='', ];
bigDB = bigDB[match(rownames(DAT), bigDB$'FISH_ID.1'), ];
bigDAT = cbind(DAT, bigDB);