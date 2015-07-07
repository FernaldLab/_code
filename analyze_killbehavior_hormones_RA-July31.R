# # rm(list=ls());
# # setwd('/Volumes/fishstudies/_behaviorRA');

# # dbNew = read.csv('Database 072314update_NewDB.csv');
# # dbAgg = read.csv('Database 072314update_AggressionDmales.csv');
# # datH = read.csv('AggressionSummarySteroids072114UPDATE.csv');

# # dbNew$SEX = gsub(' ', '', dbNew$SEX);
# # dbNew = dbNew[dbNew$SEX=='male', ];
# # dbNew$TISSUE_TYPE = gsub(' ', '', dbNew$TISSUE_TYPE);
# # dbNew = dbNew[dbNew$TISSUE_TYPE=='testes', ];

# # DAT = cbind(dbNew, 
			# # dbAgg[match(dbNew$FISH_ID, dbAgg$'Fish..'), ], 
			# # datH[match(dbNew$FISH_ID, datH$'Fish.ID'), ]
			# # );
			
# # write.csv(DAT, file='masterDB_updated_July23.csv', row.names=F);

# # #################################################################

# # DATnum = DAT[, c(4,5,18,25,27,30,32,33,35,37:40,48:51)];
# # for (i in 1:ncol(DATnum)) {
	# # DATnum[, i] = as.numeric(DATnum[, i]);
# # }; rm(i);


rm(list=ls());
setwd('/Volumes/fishstudies/_behaviorRA');
source('/Volumes/fishstudies/_code/bootstrapFunctions_6-16-13.R');
dat = read.csv('masterDBfiltered072914.csv');
dat = dat[, -c(1,2,6:8,11:15,17,18,20:25,27,29,30,44,45,55:57)];
dat = cbind(dat, 'time.in.dyad'=dat$'Age.at.dissection'-dat$'Age.at.dyad.creation', gp=rep(NA,nrow(dat)));
dat$gp[dat$'Total.Kills'==0] = 'kill0';
dat$gp[dat$'Total.Kills' %in% c(1,2)] = 'kill1/2';
dat$gp[dat$'Total.Kills'>=3] = 'kill>=3';
dat = cbind(dat, gp2=dat$gp);
dat$gp2[dat$'Total.Kills'==1] = 'kill1';
dat$gp2[dat$'Total.Kills'==2] = 'kill2';
#write.csv(dat, file='masterDBfiltered072914_edited_in_R.csv');

# check for relationship between total kills and total time in dyad
WGCNA::verboseScatterplot(dat$'time.in.dyad', dat$'Total.Kills', abline=T, abline.col='red', frame.plot=F, xlab='Total time in dyad (days)', ylab='Total kills', cex=1.2, main=paste('(n=', length(dat$'time.in.dyad'), ')', sep=''));

##########

longterm = c(which(dat$'paradigm..NA.longterm.'=='longterm'), which(is.na(dat$'paradigm..NA.longterm.')));
stable = c(which(dat$'timepoint..NA.stable.'=='stable'), which(is.na(dat$'timepoint..NA.stable.')));
postattack = which(grepl('POST_ATTACK',dat$'timepoint..NA.stable.'))
kill0 = which(dat$'Total.Kills'==0);

# compare latency to 1st kill in 1-2kill fish vs >=3kill fish
DAT = dat[longterm, ];
DAT = dat[c(longterm,postattack), ];

#DAT = DAT[DAT$'Total.Kills'>0, ];
#DAT = DAT[DAT$'Total.Kills'>1, ];
#DAT = dat[grepl('ATTACK',dat$'timepoint..NA.stable.'),];
#DAT = dat[grepl('DURING_ATTACK',dat$'timepoint..NA.stable.'),];

# within kill0, compare during vs post attack
par(mfrow=c(3,2));
YLIM=c(0,3e+5);
DAT = dat[grepl('ATTACK',dat$'timepoint..NA.stable.'),];
DAT = DAT[DAT$gp2=='kill0', ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Estradiolpg.mL', main='kill0: ',ylim=YLIM);
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Cortisol.pg.mL', main='kill0: ',ylim=YLIM);
DAT = dat[grepl('ATTACK',dat$'timepoint..NA.stable.'),];
DAT = DAT[DAT$gp2 %in% c('kill2','kill>=3'), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Estradiolpg.mL', main='kill>1: ',ylim=YLIM);
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Cortisol.pg.mL', main='kill>1: ',ylim=YLIM);
DAT = dat[grepl('ATTACK',dat$'timepoint..NA.stable.'),];
DAT = DAT[DAT$gp2=='kill>=3', ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Estradiolpg.mL', main='kill>=3: ',ylim=YLIM);
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$'timepoint..NA.stable.'), col='grey', xlab='', ylab='Cortisol.pg.mL', main='kill>=3: ',ylim=YLIM);


# 1/2 kills vs >=3 kills
par(mfrow=c(1,3));
YLIM = c(0,60);
XLIM = c(0,60);
BREAKS = c(0,10,20,30,40,50,60);
CEX=1
TOPLOT = DAT[DAT$'Total.Kills' >= 3,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('>=3 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);
TOPLOT = DAT[DAT$'Total.Kills' < 3,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('1/2 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);
WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill', as.factor(DAT$'Total.Kills' >= 3),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', names=c('1/2 kills','>=3 kills'), cex=CEX);


# 1 kill vs 2 kills vs >=3 kills
par(mfrow=c(2,3));
YLIM = c(0,60);
XLIM = c(0,60);
BREAKS = c(0,10,20,30,40,50,60);
CEX=1
TOPLOT = DAT[DAT$'Total.Kills' >= 3,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('>=3 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);

TOPLOT = DAT[DAT$'Total.Kills' == 2,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('2 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);

TOPLOT = DAT[DAT$'Total.Kills' == 1,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('1 kill (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);

WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill', as.factor(DAT$gp2),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', cex=CEX, frame.plot=F,ylim=YLIM);
WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill'[!DAT$gp2=='kill1'], as.factor(DAT$gp2[!DAT$gp2=='kill1']),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', cex=CEX, frame.plot=F,ylim=YLIM);
WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill'[!DAT$gp2=='kill2'], as.factor(DAT$gp2[!DAT$gp2=='kill2']),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', cex=CEX, frame.plot=F,ylim=YLIM);

#
par(mfrow=c(1,2))
WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill',as.factor(DAT$'Total.Kills'),notch=T,col='grey', xlab='Total kills', ylab='Days in dyad until first kill', main='', cex=CEX, frame.plot=F,ylim=YLIM)
WGCNA::verboseScatterplot(DAT$'Total.Kills', DAT$'Days.in.dyad.until.1st.kill', xlab='Total kills', ylab='Days in dyad until first kill', main='', cex=CEX, frame.plot=F,ylim=YLIM, abline=T)


#### figure
par(oma=c(0,1,0,0))
WGCNA::verboseScatterplot(DAT$'Total.Kills'[DAT$'Total.Kills'>0], DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'>0], xlab='Fish damaged by dominant', ylab='Days in dyad until first damage', main='', cex=CEX, frame.plot=F,ylim=c(0,60), abline=T,xlim=c(.5,5.5));
segments(.75,median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==1]), 1.25, median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==1]),lwd=2)
segments(1.75,median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==2]), 2.25, median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==2]),lwd=2)
segments(2.75,median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==3]), 3.25, median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==3]),lwd=2)
segments(3.75,median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==4]), 4.25, median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==4]),lwd=2)
segments(4.75,median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==5]), 5.25, median(DAT$'Days.in.dyad.until.1st.kill'[DAT$'Total.Kills'==5]),lwd=2)

# 1/2 kills vs >=3 kills
par(mfrow=c(2,3));
#par(mfrow=c(1,3));
YLIM = c(0,60);
XLIM = c(0,60);
BREAKS = c(0,10,20,30,40,50,60);
CEX=1
TOPLOT = DAT[DAT$'Total.Kills' >= 3,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('>=3 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);
TOPLOT = DAT[DAT$'Total.Kills' < 3,]$'Days.in.dyad.until.1st.kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='Days in dyad until first kill', main=paste('1-2 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);
WGCNA::verboseBoxplot(DAT$'Days.in.dyad.until.1st.kill', as.factor(DAT$'Total.Kills' >= 3),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', names=c('1-2 kills','>=3 kills'), cex=CEX);

#par(mfrow=c(1,3));
YLIM = c(0,60);
XLIM = c(0,60);
BREAKS = c(0,10,20,30,40);
#BREAKS=10
CEX=1
TOPLOT = DAT[DAT$'Total.Kills' >= 3,]$'latency.2nd.Kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='latency to 2nd kill (days)', main=paste('>=3 kills (n=', length(TOPLOT), ')', sep=''), breaks=BREAKS, cex=CEX);
TOPLOT = DAT[DAT$'Total.Kills' < 3,]$'latency.2nd.Kill';
hist(TOPLOT, col='grey', border='darkgrey', xlim=XLIM, ylim=YLIM, xlab='latency to 2nd kill (days)', main=paste('1-2 kills (n=', sum(!is.na(TOPLOT)), ')', sep=''), breaks=BREAKS, cex=CEX);
WGCNA::verboseBoxplot(DAT$'latency.2nd.Kill', as.factor(DAT$'Total.Kills' >= 3),notch=T,col='grey', xlab='', ylab='latency to 2nd kill (days)', main='', names=c('1-2 kills','>=3 kills'), cex=CEX, ylim=YLIM);

##############

par(mfrow=c(1,2));
gp = 'kill1/2';
WGCNA::verboseBoxplot(as.numeric(c(DAT$'Days.in.dyad.until.1st.kill'[DAT$gp==gp], DAT$'latency.2nd.Kill'[DAT$gp==gp])), as.factor(c(rep('1st',sum(DAT$gp==gp)), rep('2nd',sum(DAT$gp==gp)))),xlab='',ylab='latency to kill (days)', col='grey', main=paste(gp, ': ', sep=''));
gp = 'kill>=3';

WGCNA::verboseBoxplot(as.numeric(c(DAT$'Days.in.dyad.until.1st.kill'[DAT$gp==gp], DAT$'latency.2nd.Kill'[DAT$gp==gp], DAT$'latency.3rd.kill'[DAT$gp==gp])), 
					  as.factor(c(rep('1st',sum(DAT$gp==gp)), rep('2nd',sum(DAT$gp==gp)), rep('3rd',sum(DAT$gp==gp)))),
					  xlab='',ylab='latency to kill (days)', 
					  col='grey', main=paste(gp, ': ', sep=''), 
					  ylim=c(0,60));
					  
					  
par(mfrow=c(1,3));
WGCNA::verboseBoxplot(as.numeric(c(DAT$'Days.in.dyad.until.1st.kill'[DAT$gp==gp], DAT$'latency.2nd.Kill'[DAT$gp==gp])), 
					  as.factor(c(rep('1st',sum(DAT$gp==gp)), rep('2nd',sum(DAT$gp==gp)))),
					  xlab='',ylab='latency to kill (days)', 
					  col='grey', main=paste(gp, ': ', sep=''), 
					  ylim=c(0,30));

WGCNA::verboseBoxplot(as.numeric(c(DAT$'latency.2nd.Kill'[DAT$gp==gp], DAT$'latency.3rd.kill'[DAT$gp==gp])), 
					  as.factor(c(rep('2nd',sum(DAT$gp==gp)), rep('3rd',sum(DAT$gp==gp)))),
					  xlab='',ylab='latency to kill (days)', 
					  col='grey', main=paste(gp, ': ', sep=''), 
					  ylim=c(0,30));

WGCNA::verboseBoxplot(as.numeric(c(DAT$'Days.in.dyad.until.1st.kill'[DAT$gp==gp], DAT$'latency.3rd.kill'[DAT$gp==gp])), 
					  as.factor(c(rep('1st',sum(DAT$gp==gp)), rep('3rd',sum(DAT$gp==gp)))),
					  xlab='',ylab='latency to kill (days)', 
					  col='grey', main=paste(gp, ': ', sep=''), 
					  ylim=c(0,30));

DAT = DAT[-which(DAT$'latency.3rd.kill'=='no info due to vacation'),]					# "no info due to vacation"
par(mfrow=c(1,2));
out=bootstrap.2paired(DAT[DAT$gp==gp, ]$'Days.in.dyad.until.1st.kill', 
					  DAT[DAT$gp==gp, ]$'latency.2nd.Kill',
					  noDist=T,
					  dataDescriptor='latency to kill (days)', 
					  conditionNames=c('1st','2nd'),xlab='',
					  main=gp,
					  ylim=c(0,25));
out=bootstrap.2paired(DAT[DAT$gp==gp, ]$'latency.2nd.Kill', 
					  as.numeric(DAT[DAT$gp==gp, ]$'latency.3rd.kill'),
					  noDist=T,
					  dataDescriptor='latency to kill (days)', 
					  conditionNames=c('2nd','3rd'),xlab='',
					  main=gp,
					  ylim=c(0,25));					  

LINECOL = 'darkgrey';				  
TOPLOT = DAT[, colnames(DAT) %in% c('Days.in.dyad.until.1st.kill', 'latency.2nd.Kill', 'latency.3rd.kill')];
TOPLOT = apply(TOPLOT,2,as.numeric);
TOPLOT = TOPLOT[DAT$gp=='kill>=3',]
boxplot(TOPLOT, boxcol='lightgrey', col='lightgrey', ylab='latency to kill (days)', names=c('1st','2nd','3rd'), frame.plot=F, ylim=c(0,30), main='Aggression level 3');
stripchart(data.frame(TOPLOT), vertical = T, add = T, pch = 21, cex = 1.5, bg='black');
for (row in 1:nrow(TOPLOT)) {
	if (TOPLOT[row, 1] > TOPLOT[row, 2] & TOPLOT[row, 2] > TOPLOT[row, 3]) {
		col = 'red';
	} else {
		col = LINECOL;
	}
	segments(1, TOPLOT[row, 1], 2, TOPLOT[row, 2], col = col)
}
for (row in 1:nrow(TOPLOT)) {
	if (TOPLOT[row, 1] > TOPLOT[row, 2] & TOPLOT[row, 2] > TOPLOT[row, 3]) {
		col = 'red';
	} else {
		col = LINECOL;
	}
	segments(2, TOPLOT[row, 2], 3, TOPLOT[row, 3], col = col)
}
segments(1,29,2,29);
text(1.5,30,'****');
segments(2,22,3,22);
text(2.5,23,'*');



		# # data.lab = dataDescriptor;
			# # # boxplot of conditions with individual data points
			# # toPlot = data[, 1:2];
			
			# # if (ptemp$p == 0) {toPaste = paste('p < 1e-5 (n=', length(diffs), ')', sep = '')}
			# # else {toPaste = paste('p = ', signif(ptemp$p, 2), ' (n=', length(diffs), ')', sep = '')}
			# # toPaste=paste(main, ': ', toPaste, sep='')
			# # # draw boxes
			# # boxplot(toPlot, ylab = data.lab,
					# # #medlty = medlty,
					# # main = toPaste, frame.plot = F,
					# # #cex.lab = cex.lab,
					# # #cex.axis = cex.axis,
					# # ...);
					
			# # # add data points
			# # stripchart(toPlot, vertical = T, add = T, pch = pch, bg = col, cex = 1.5);
					   
			# # # connect paired points across conditions
			# # for (row in 1:nrow(data)) {segments(1, data[row, 1], 2, data[row, 2], col = border)}
###################

DAT = dat[longterm, ];

par(mfrow=c(2,2));
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp=='kill>=3'), col='grey', xlab='', ylab='Estradiol(pg/mL)',names=c('kill<3','kill>=3'),ylim=c(0,4e+5))

WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp=='kill>=3'), col='grey', xlab='', ylab='Cortisol(pg/mL)',ylim=c(0,4e+5),names=c('kill<3','kill>=3'));

WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp), col='grey', xlab='', ylab='Estradiol(pg/mL)',ylim=c(0,4e+5))

WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp), col='grey', xlab='', ylab='Cortisol(pg/mL)',ylim=c(0,4e+5))

DAT = dat[c(longterm,postattack), ];
DAT = DAT[DAT$'Total.Kills'>0, ];
par(mfrow=c(2,3));
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp), col='grey', xlab='', ylab='Estradiol(pg/mL)',ylim=c(0,4e+5))
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp), col='grey', xlab='', ylab='Cortisol(pg/mL)',ylim=c(0,4e+5));
DAT = dat[c(longterm,postattack), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp2), col='grey', xlab='', ylab='Estradiol(pg/mL)',ylim=c(0,4e+5))
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp2), col='grey', xlab='', ylab='Cortisol(pg/mL)',ylim=c(0,4e+5))
DAT = dat[c(longterm,postattack), ];
DAT = DAT[DAT$gp %in% c('kill0','kill>=3'), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp2), col='grey', xlab='', ylab='Estradiol(pg/mL)',ylim=c(0,4e+5))
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp2), col='grey', xlab='', ylab='Cortisol(pg/mL)',ylim=c(0,4e+5))

par(mfrow=c(1,2));
WGCNA::verboseScatterplot(DAT$'Total.Kills', DAT$'Estradiolpg.mL',abline=T,frame.plot=F,xlab='Total kills', ylab='Estradiol(pg/mL)');
WGCNA::verboseScatterplot(DAT$'Total.Kills', as.numeric(DAT$'Cortisol.pg.mL'),abline=T,frame.plot=F,xlab='Total kills',ylab='Cortisol(pg/mL)');

par(mfrow=c(1,2));
WGCNA::verboseScatterplot(DAT$'time.in.dyad', DAT$'Estradiolpg.mL',abline=T,frame.plot=F,xlab='Total time in dyad (days)', ylab='Estradiol(pg/mL)',ylim=c(0,4e+5));
WGCNA::verboseScatterplot(DAT$'time.in.dyad', as.numeric(DAT$'Cortisol.pg.mL'),abline=T,frame.plot=F,xlab='Total time in dyad (days)',ylab='Cortisol(pg/mL)',ylim=c(0,4e+5));



par(mfrow=c(1,3))
WGCNA::verboseScatterplot(DAT$'Estradiolpg.mL'[DAT$gp=='kill0'], as.numeric(DAT$'Cortisol.pg.mL'[DAT$gp=='kill0']),abline=T,frame.plot=F,xlab='Estradiol(pg/mL)', ylab='Cortisol(pg/mL)',ylim=c(0,2.5e+5), xlim=c(0,4e+5),main='kill0: ');
WGCNA::verboseScatterplot(DAT$'Estradiolpg.mL'[DAT$gp=='kill1/2'], as.numeric(DAT$'Cortisol.pg.mL'[DAT$gp=='kill1/2']),abline=T,frame.plot=F,xlab='Estradiol(pg/mL)', ylab='Cortisol(pg/mL)',ylim=c(0,2.5e+5), xlim=c(0,4e+5),main='kill1/2: ');
WGCNA::verboseScatterplot(DAT$'Estradiolpg.mL'[DAT$gp=='kill>=3'], as.numeric(DAT$'Cortisol.pg.mL'[DAT$gp=='kill>=3']),abline=T,frame.plot=F,xlab='Estradiol(pg/mL)', ylab='Cortisol(pg/mL)',ylim=c(0,2.5e+5), xlim=c(0,4e+5),main='kill>=3: ');

###

DAT = dat[grepl('ATTACK',dat$'timepoint..NA.stable.'),];

par(mfrow=c(1,2))
#WGCNA::verboseBoxplot(DAT$'Total.Kills', as.factor(DAT$'timepoint..NA.stable.'), xlab='', ylab='Total kills', col='grey', frame.plot=F);
#WGCNA::verboseBoxplot(DAT$'Testo.pg.mL', as.factor(DAT$'timepoint..NA.stable.'), xlab='', ylab='Testo.pg.mL', col='grey', frame.plot=F);
#WGCNA::verboseBoxplot(DAT$'X11keto.pg.mL', as.factor(DAT$'timepoint..NA.stable.'), xlab='', ylab='X11keto.pg.mL', col='grey', frame.plot=F);
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$'timepoint..NA.stable.'), xlab='', ylab='Estradiolpg.mL', col='grey', frame.plot=F);
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$'timepoint..NA.stable.'), xlab='', ylab='Cortisol.pg.mL', col='grey', frame.plot=F);

# basal hormone levels in multiple kill vs 0 kill fish
par(mfrow=c(2,3))
DAT = dat[union(stable,postattack), ];
DAT = DAT[DAT$gp2 %in% c('kill0','kill2','kill>=3'), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp2=='kill0'), xlab='', ylab='Estradiolpg.mL', col='grey', frame.plot=F, names=c('kill>1','kill0'));
DAT = dat[union(stable,postattack), ];
DAT = DAT[DAT$gp2 %in% c('kill0','kill>=3'), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp2), xlab='', ylab='Estradiolpg.mL', col='grey', frame.plot=F);
DAT = dat[union(stable,postattack), ];
WGCNA::verboseBoxplot(DAT$'Estradiolpg.mL', as.factor(DAT$gp2), xlab='', ylab='Estradiolpg.mL', col='grey', frame.plot=F);

DAT = dat[union(stable,postattack), ];
DAT = DAT[DAT$gp2 %in% c('kill0','kill2','kill>=3'), ];
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp2=='kill0'), xlab='', ylab='Cortisol.pg.mL', col='grey', frame.plot=F, names=c('kill>1','kill0'));
DAT = dat[union(stable,postattack), ];
DAT = DAT[DAT$gp2 %in% c('kill0','kill>=3'), ];
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp2), xlab='', ylab='Cortisol.pg.mL', col='grey', frame.plot=F);
DAT = dat[union(stable,postattack), ];
WGCNA::verboseBoxplot(as.numeric(DAT$'Cortisol.pg.mL'), as.factor(DAT$gp2), xlab='', ylab='Cortisol.pg.mL', col='grey', frame.plot=F);


###############
#############
###############
par(mfrow=c(1,2))

TMPdat = dat[intersect(longterm, kill12), ];
WGCNA::verboseBoxplot(as.numeric(c(TMPdat$'Days.in.dyad.until.1st.kill',TMPdat$'latency.2nd.Kill')),
			   as.factor(c(rep('1st',nrow(TMPdat)),rep('2nd',nrow(TMPdat)))),
			   #numberStandardErrors=2,
			   ylim=c(0,50),
			   ylab='latency to kill (days)', xlab='kill #', 
			   col='grey',
			   notch=T,
			   main='1-2 kills: ');
			   
			   
TMPdat = dat[intersect(longterm, kill3more), ];
WGCNA::verboseBoxplot(as.numeric(c(TMPdat$'Days.in.dyad.until.1st.kill',TMPdat$'latency.2nd.Kill',TMPdat$'latency.3rd.kill')),
			   as.factor(c(rep('1st',nrow(TMPdat)),rep('2nd',nrow(TMPdat)),rep('3rd',nrow(TMPdat)))),
			   #numberStandardErrors=2,
			   ylim=c(0,50),
			   ylab='latency to kill (days)', xlab='kill #', 
			   col='grey',
			   notch=T,
			   main='>=3 kills: ');
			   
			   
			   
#####################
####################
######################
newdat = read.csv('AggressionUpdatedSheet090414_fixformulas.csv');
# has sex of killed fish
DAT = newdat[newdat$'sex.of.1st.kill'=='M' & newdat$'sex.of.2nd.kill'=='M', ];
gp = DAT$'Total.Kills' >= 2		
		
par(mfrow=c(1,2));
out=bootstrap.2paired(DAT[gp, ]$'Days.in.dyad.until.1st.kill', 
					  DAT[gp, ]$'latency.2nd.Kill',
					  noDist=T,
					  dataDescriptor='latency to kill (days)', 
					  conditionNames=c('1st','2nd'),xlab='',
					  main='',
					  ylim=c(0,50));
out=bootstrap.2paired(DAT[gp, ]$'latency.2nd.Kill', 
					  as.numeric(DAT[gp, ]$'latency.3rd.kill'),
					  noDist=T,
					  dataDescriptor='latency to kill (days)', 
					  conditionNames=c('2nd','3rd'),xlab='',
					  main='',
					  ylim=c(0,50));		


C1 = match('Days.in.dyad.until.1st.kill', colnames(DAT))
C2 = match('latency.2nd.Kill', colnames(DAT));
TOPLOT = DAT[gp,c(C1,C2)]	  
boxplot(TOPLOT,notch=T, col='lightgrey', border='darkgrey', frame.plot=F,ylim=c(0,50),
	    ylab='latency to damage other male (days)', 
	    names=c('1st','2nd'), 
	    main=paste('p=', signif(wilcox.test(TOPLOT[,1], TOPLOT[,2], paired=T)$p.value, 3), ' (n=', nrow(TOPLOT), ')',sep='')
	    );
stripchart(TOPLOT, vertical = T, add = T, pch=21, bg='black')
for (i in 1:nrow(TOPLOT)) {
	if (TOPLOT[i,1] > TOPLOT[i,2]) {
		COL='red'
	} else {
		COL='darkgrey'
	}
	segments(1,TOPLOT[i,1],2,TOPLOT[i,2],col=COL);
}; rm(i);




		  

LINECOL = 'darkgrey';				  
TOPLOT = DAT[, colnames(DAT) %in% c('Days.in.dyad.until.1st.kill', 'latency.2nd.Kill', 'latency.3rd.kill')];
TOPLOT = apply(TOPLOT,2,as.numeric);
TOPLOT = TOPLOT[gp,]
boxplot(TOPLOT, boxcol='lightgrey', col='lightgrey', ylab='latency to damage (days)', names=c('1st','2nd','3rd'), frame.plot=F, ylim=c(0,30), main='Aggression level 3');
stripchart(data.frame(TOPLOT), vertical = T, add = T, pch = 21, cex = 1.5, bg='black');
for (row in 1:nrow(TOPLOT)) {
	if (TOPLOT[row, 1] > TOPLOT[row, 2] & TOPLOT[row, 2] > TOPLOT[row, 3]) {
		col = 'red';
	} else {
		col = LINECOL;
	}
	segments(1, TOPLOT[row, 1], 2, TOPLOT[row, 2], col = col)
}
for (row in 1:nrow(TOPLOT)) {
	if (TOPLOT[row, 1] > TOPLOT[row, 2] & TOPLOT[row, 2] > TOPLOT[row, 3]) {
		col = 'red';
	} else {
		col = LINECOL;
	}
	segments(2, TOPLOT[row, 2], 3, TOPLOT[row, 3], col = col)
}
# segments(1,29,2,29);
# text(1.5,30,'****');
# segments(2,22,3,22);
# text(2.5,23,'*');

##################
# redo old plot but split into MvsF first kills
# # # DAT = dat[longterm, ];
# # # > TOPLOT = DAT[DAT$gp2 %in% c('kill1','kill>=3'), ]
# # # > WGCNA::verboseBoxplot(TOPLOT$'Days.in.dyad.until.1st.kill', as.factor(TOPLOT$gp2),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main='', cex=CEX);


DAT=cbind(DAT,gp3=rep(NA,nrow(DAT)));
DAT$gp3[DAT$'Total.Kills'<2] = 'kill01';
DAT$gp3[DAT$'Total.Kills'>=2] = 'kill>=2';

Mfirst = datnew$FISH_ID[datnew$sex.of.1st.kill=='M'];
Ffirst = datnew$FISH_ID[datnew$sex.of.1st.kill=='F'];
DATm = DAT[DAT$FISH_ID %in% Mfirst, ];
DATf = DAT[DAT$FISH_ID %in% Ffirst, ];

par(mfrow=c(2,2));
YLIM=c(0,50)
TOPLOT = DATm[DATm$gp2 %in% c('kill1','kill>=3'), ]
WGCNA::verboseBoxplot(TOPLOT$'Days.in.dyad.until.1st.kill', as.factor(TOPLOT$gp2),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('M first (n=',nrow(TOPLOT),'): ',sep=''), cex=CEX, ylim=YLIM);
TOPLOT = DATf[DATf$gp2 %in% c('kill1','kill>=3'), ]
WGCNA::verboseBoxplot(TOPLOT$'Days.in.dyad.until.1st.kill', as.factor(TOPLOT$gp2),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('F first (n=',nrow(TOPLOT),'): ',sep=''), cex=CEX, ylim=YLIM);
tmp1 = DATm[DATm$gp2=='kill>=3',]$'Days.in.dyad.until.1st.kill';
tmp2 = DATf[DATf$gp2=='kill>=3',]$'Days.in.dyad.until.1st.kill';
WGCNA::verboseBoxplot(c(tmp1,tmp2), as.factor(c(rep('M first',length(tmp1)), rep('F first',length(tmp2)))),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('kill>=3 (n=',length(c(tmp1,tmp2)),'):',sep=''), cex=CEX, ylim=YLIM);
tmp1 = DATm[DATm$gp2=='kill1',]$'Days.in.dyad.until.1st.kill';
tmp2 = DATf[DATf$gp2=='kill1',]$'Days.in.dyad.until.1st.kill';
WGCNA::verboseBoxplot(c(tmp1,tmp2), as.factor(c(rep('M first',length(tmp1)), rep('F first',length(tmp2)))),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('kill1 (n=',length(c(tmp1,tmp2)),'):',sep=''), cex=CEX, ylim=YLIM);


par(mfrow=c(2,2));
YLIM=c(0,50)
TOPLOT = DATm
WGCNA::verboseBoxplot(TOPLOT$'Days.in.dyad.until.1st.kill', as.factor(TOPLOT$gp3),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('M first (n=',nrow(TOPLOT),'): ',sep=''), cex=CEX, ylim=YLIM);
TOPLOT = DATf
WGCNA::verboseBoxplot(TOPLOT$'Days.in.dyad.until.1st.kill', as.factor(TOPLOT$gp3),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('F first (n=',nrow(TOPLOT),'): ',sep=''), cex=CEX, ylim=YLIM);
tmp1 = DATm[DATm$gp3=='kill>=2',]$'Days.in.dyad.until.1st.kill';
tmp2 = DATf[DATf$gp3=='kill>=2',]$'Days.in.dyad.until.1st.kill';
WGCNA::verboseBoxplot(c(tmp1,tmp2), as.factor(c(rep('M first',length(tmp1)), rep('F first',length(tmp2)))),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('kill>=2 (n=',length(c(tmp1,tmp2)),'):',sep=''), cex=CEX, ylim=YLIM);
tmp1 = DATm[DATm$gp3=='kill01',]$'Days.in.dyad.until.1st.kill';
tmp2 = DATf[DATf$gp3=='kill01',]$'Days.in.dyad.until.1st.kill';
WGCNA::verboseBoxplot(c(tmp1,tmp2), as.factor(c(rep('M first',length(tmp1)), rep('F first',length(tmp2)))),notch=T,col='grey', xlab='', ylab='Days in dyad until first kill', main=paste('kill01 (n=',length(c(tmp1,tmp2)),'):',sep=''), cex=CEX, ylim=YLIM);




##################
##############

# final all hormone data compared to number of days between last kill and dissection

newdat = read.csv('AggressionUpdatedSheet090414_fixformulas.csv');
skint = c();
for ( r in 1:nrow(newdat) ) {
	kcols = newdat[r, grepl('kill.date', names(newdat))];
	if (all(kcols=='')) { next };
	check = which(kcols=='');
	if (length(check)==0) {
		lastdate = as.Date(as.character(kcols[length(kcols)]), '%m/%d/%Y');
	} else {
		lastdate = as.Date(as.character(kcols[min(check) - 1]), '%m/%d/%Y');
	}
	sacdate = as.Date(newdat[r, names(newdat)=='DISSECTION_DATE'], '%m/%d/%Y');
	days = as.numeric(sacdate - lastdate);
	skint = c(skint, days);
	names(skint)[r] = newdat[r, names(newdat)=='FISH_ID'];
}; rm(r);


# cortisol
cort0 = read.csv('combined_summaryelisaruns7-13Austin_cortisol.csv');
cort = cort0[ as.character(cort0$Fish.ID) %in% names(skint),];
cort = cbind(cort, skint=skint[match(as.character(cort$Fish.ID), names(skint))]);
WGCNA::verboseBoxplot(cort$skint, cort$timepoint, notch=F, col='grey', xlab='', ylab='days', frame.plot=F);
WGCNA::verboseScatterplot(cort$skint, cort$Cortisol.pg.mL, abline=T, frame.plot=F, xlab='days', ylab='cortisol');

par(mfrow=c(2,2));
YMAX=200000
XMAX=30
for (GP in c('10 days POST_ATTACK', 'DURING_ATTACK', 'stable')) {
	ROWS = cort$timepoint==GP;
	WGCNA::verboseScatterplot(cort$skint[ROWS], cort$Cortisol.pg.mL[ROWS], abline=T, frame.plot=F, xlab='days', ylab='cortisol', main=paste(GP,'\n',sep=''), cex=1.5, pch=21, ylim=c(0,YMAX), xlim=c(0,XMAX), abline.col='red');
}
ROWS = cort$timepoint %in% c('10 days POST_ATTACK', 'stable')
WGCNA::verboseScatterplot(cort$skint[ROWS], cort$Cortisol.pg.mL[ROWS], abline=T, frame.plot=F, xlab='days', ylab='cortisol', main='10 days POST_ATTACK + stable\n', cex=1.5, pch=21, ylim=c(0,YMAX), xlim=c(0,XMAX), abline.col='red');

par(mfrow=c(2,2));
WGCNA::verboseBoxplot(cort$skint, cort$timepoint, notch=F, col='grey', xlab='', ylab='days', frame.plot=F, cex.axis=1.1);
WGCNA::verboseBoxplot(cort$skint, cort$timepoint=='DURING_ATTACK', notch=F, col='grey', xlab='', ylab='days', frame.plot=F, cex.axis=1.1, names=c('10 days POST_ATTACK + stable', 'DURING_ATTACK'));
WGCNA::verboseBoxplot(cort$Cortisol.pg.mL, cort$timepoint, notch=F, col='grey', xlab='', ylab='cortisol', frame.plot=F, cex.axis=1.1);
WGCNA::verboseBoxplot(cort$Cortisol.pg.mL, cort$timepoint=='DURING_ATTACK', notch=F, col='grey', xlab='', ylab='cortisol', frame.plot=F, cex.axis=1.1, names=c('10 days POST_ATTACK + stable', 'DURING_ATTACK'));

ckills = newdat[newdat$FISH_ID %in% cort$Fish.ID, c(3,65)];
cort = cbind(cort, kills=ckills[match(cort$Fish.ID, ckills$FISH_ID),2]);

##############

# hand-curated by RA, number of days alone data (number of days between last male kill and dissection)
rm(list=ls());
setwd('/Volumes/fishstudies/_behaviorRA');
ddat0 = read.csv('dayswithoutsubordinate.csv');
ddat = ddat0[, c(63,64,73:88)];

# tmp = read.csv('dayswithoutsubordinate100114update.csv')[,c(3,79)];
# tmp = tmp[tmp$FISH_ID %in% ddat$Fish.ID.1,];
# ddat = ddat[ddat$Fish.ID.1 %in% tmp$FISH_ID, ];
# tmp = tmp[match(ddat$Fish.ID.1, tmp$FISH_ID), ];
# ddat$days.without.subordinate = tmp$days.between.last.kill.and.dissection;

par(mfrow=c(1,3), oma=c(0,0,3,0))
C1 = 9;
YLIM = c(0, 300000); XLIM = c(0, 10);
CEX = 1.1;
for (C2 in c(5,6,12)) {
	if(C2==12){YLIM=c(0,6)}
	WGCNA::verboseScatterplot(as.numeric(ddat[,C1]), as.numeric(ddat[,C2]), xlim=XLIM, ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab=names(ddat)[C2],
						 type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
	text(ddat[,C1], ddat[,C2], labels=ddat$Total.Kills.1);
}
title(main=paste('stable (n=', nrow(ddat), ')\n', sep=''), outer=T);

TOPLOT = ddat[ddat0$plate=='2', ]
par(mfrow=c(1,3), oma=c(0,0,3,0))
C1 = 9;
YLIM = c(0, 300000); XLIM = c(0, 10);
CEX = 1.1;
for (C2 in c(5,6,12)) {
	if(C2==12){YLIM=c(0,6)}
	WGCNA::verboseScatterplot(as.numeric(TOPLOT[,C1]), as.numeric(TOPLOT[,C2]), xlim=XLIM, ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(TOPLOT)[C1], ylab=names(TOPLOT)[C2],
						 type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
	text(TOPLOT[,C1], TOPLOT[,C2], labels=TOPLOT$Total.Kills.1);
}
title(main=paste('stable (n=', nrow(TOPLOT), ')\n', sep=''), outer=T);


kt = read.csv('combined_summaryelisaruns7-13Austin_11keto_run_10_old.csv');
tt = read.csv('combined_summaryelisaruns7-13Austin_testo_plate5_run11.csv');

kt = kt[match(ddat$Fish.ID.1,kt$Fish.ID),];
tt = tt[match(ddat$Fish.ID.1,tt$FishID),];

ktcombo0 = cbind(ddat$X11keto.pg.mL, kt$X11keto.pg.mL);
ktcombo = c();
for (r in 1:nrow(ktcombo0)) {
	if (is.na(ktcombo0[r,1])) {
		if (is.na(ktcombo0[r,2])) {
			ktcombo=c(ktcombo, NA);
			next;
		} else {
			ktcombo=c(ktcombo, ktcombo0[r,2]);
			next;
		}
	}
	if (is.numeric(ktcombo0[r,1])) {
		if (is.na(ktcombo0[r,2])) {
			ktcombo=c(ktcombo, ktcombo0[r,1]);
			next;
		} else {
			ktcombo=c(ktcombo, ktcombo0[r,2]);
			next;
		}
	}
}; rm(r);
names(ktcombo) = ddat$Fish.ID.1;

ttcombo0 = cbind(ddat$Testo.pg.mL, tt$tesstosterone.pg.mL);
ttcombo = c();
for (r in 1:nrow(ttcombo0)) {
	if (is.na(ttcombo0[r,1])) {
		if (is.na(ttcombo0[r,2])) {
			ttcombo=c(ttcombo, NA);
			next;
		} else {
			ttcombo=c(ttcombo, ttcombo0[r,2]);
			next;
		}
	}
	if (is.numeric(ttcombo0[r,1])) {
		if (is.na(ttcombo0[r,2])) {
			ttcombo=c(ttcombo, ttcombo0[r,1]);
			next;
		} else {
			ttcombo=c(ttcombo, ttcombo0[r,2]);
			next;
		}
	}
}; rm(r);
names(ttcombo) = ddat$Fish.ID.1;

par(mfrow=c(2,3), oma=c(0,0,3,0))
C2=5;
WGCNA::verboseScatterplot(ddat[,C1], kt[,C2], xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab=names(kt)[C2],
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], kt[,C2], labels=ddat$Total.Kills.1);
C2=4;
WGCNA::verboseScatterplot(ddat[,C1], ddat[,C2], xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab=names(ddat)[C2],
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], ddat[,C2], labels=ddat$Total.Kills.1);
WGCNA::verboseScatterplot(ddat[,C1], ktcombo, xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab='11kt-both plates',
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], ktcombo, labels=ddat$Total.Kills.1);

#title(main=paste('stable (n=', nrow(tt)-sum(is.na(kt$X11keto.pg.mL)), ')\n', sep=''), outer=T);
C2=5;
WGCNA::verboseScatterplot(ddat[,C1], tt[,C2], xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab=names(tt)[C2],
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], tt[,C2], labels=ddat$Total.Kills.1);
C2=3;
WGCNA::verboseScatterplot(ddat[,C1], ddat[,C2], xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab=names(ddat)[C2],
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], ddat[,C2], labels=ddat$Total.Kills.1);
WGCNA::verboseScatterplot(ddat[,C1], ttcombo, xlim=XLIM, #ylim=YLIM,
						  abline=T, abline.col='red', frame.plot=F, 
						  xlab=names(ddat)[C1], ylab='T-both plates',
						  type='n',
						  cex.lab=CEX, cex.axis=CEX, cex.main=CEX);
text(ddat[,C1], ttcombo, labels=ddat$Total.Kills.1);



#####
# 0 day mult kill vs nonkillers

kt = read.csv('combined_summaryelisaruns7-13Austin_11keto_run_10_old.csv');
tt = read.csv('combined_summaryelisaruns7-13Austin_testo_plate5_run11.csv');
newdat = read.csv('AggressionUpdatedSheet090414_fixformulas.csv');
kills = newdat[, c(19,65)];
zeroday = ddat$Fish.ID.1[ddat$days.without.subordinate==0];

kt = cbind(kt, kills[match(kt$Fish.ID, kills$FISH_ID.1),]);
kt = cbind(kt, days.without.subordinate=ddat[match(kt$Fish.ID, ddat$Fish.ID.1),9]);
kt = kt[kt$timepoint=='stable',];
kt$days.without.subordinate[is.na(kt$days.without.subordinate)] = 0;
kt = kt[kt$days.without.subordinate==0, ];

tt = cbind(tt, kills[match(tt$FishID, kills$FISH_ID.1),]);
tt = cbind(tt, days.without.subordinate=ddat[match(tt$FishID, ddat$Fish.ID.1),9]);
tt = tt[tt$timepoint=='stable',];
tt$days.without.subordinate[is.na(tt$days.without.subordinate)] = 0;
tt = tt[tt$days.without.subordinate==0, ];


par(mfrow=c(1,2), oma=c(0,0,2,0))
WGCNA::verboseBoxplot(kt$X11keto.ng.ml, kt$Total.Kills>=2, 
	   				  notch=F, col='grey',xlab='',ylab='11keto.ng.ml', frame.plot=F,
	   				  names=c(paste('kill01 (n=', sum(kt$Total.Kills<2), ')', sep=''),
	   				  		  paste('kill>=2 (n=', sum(kt$Total.Kills>=2), ')', sep=''))
	   				  )
WGCNA::verboseBoxplot(tt$testosterone.ng.ml, tt$Total.Kills>=2, 
					  notch=F, col='grey',xlab='', ylab='testosterone.ng.ml', frame.plot=F,
					  names=c(paste('kill01 (n=', sum(tt$Total.Kills<2), ')', sep=''),
	   				  		  paste('kill>=2 (n=', sum(tt$Total.Kills>=2), ')', sep=''))
	   				  )
title('stable, only fish with 0 days without subordinate', outer=T);

x=read.csv('androgenshighlowaggression.csv');
par(mfrow=c(1,2));
WGCNA::verboseBoxplot(x$X11keto.ng.ml, x$timepoint,notch=F,frame.plot=F,xlab='',ylab='11keto.ng.ml', col='grey', names=c('Level1-2 (n=3)','Level3 (n=6)'))
WGCNA::verboseBoxplot(x$testosterone.ng.ml, x$timepoint,notch=F,frame.plot=F,xlab='',ylab='testosterone.ng.ml', col='grey', names=c('Level1-2 (n=3)','Level3 (n=6)'))
#############################

cort = read.csv('combined_summaryelisaruns7-13Austin_cortisol.csv');
estra = read.csv('combined_summaryelisaruns7-13Austin_estradiol_plate4_run13.csv');
testo = read.csv('combined_summaryelisaruns7-13Austin_testo_plate6_run12old.csv');
keto = read.csv('combined_summaryelisaruns7-13Austin_11keto_run9.csv');
newdat = read.csv('AggressionUpdatedSheet090414_fixformulas.csv');

tmp = newdat[newdat$Fish_ID %in% cort$Fish.ID, ];
tmp = tmp[match(cort$Fish.ID, tmp$Fish_ID), ];
cort = cbind(cort, kills=tmp$Total.Kills);
gp = rep('x',nrow(cort));
gp[cort$kills>=2] = 'kill>=2';
gp[gp=='x'] = 'kill01';
cort = cbind(cort, gp);
cort = cort[!(cort$timepoint=='stable'), ];

tmp = newdat[newdat$Fish_ID %in% estra$Fish.ID, ];
tmp = tmp[match(estra$Fish.ID, tmp$Fish_ID), ];
estra = cbind(estra, kills=tmp$Total.Kills);
gp = rep('x',nrow(estra));
gp[estra$kills>=2] = 'kill>=2';
gp[gp=='x'] = 'kill01';
estra = cbind(estra, gp);
estra = estra[!(estra$timepoint=='stable'), ];


tmp = newdat[newdat$Fish_ID %in% testo$Fish.ID, ];
tmp = tmp[match(testo$Fish.ID, tmp$Fish_ID), ];
testo = cbind(testo, kills=tmp$Total.Kills);
gp = rep('x',nrow(testo));
gp[testo$kills>=2] = 'kill>=2';
gp[gp=='x'] = 'kill01';
testo = cbind(testo, gp);
testo = testo[!(testo$timepoint=='stable'), ];


tmp = newdat[newdat$Fish_ID %in% keto$Fish.ID, ];
tmp = tmp[match(keto$Fish.ID, tmp$Fish_ID), ];
keto = cbind(keto, kills=tmp$Total.Kills);
gp = rep('x',nrow(keto));
gp[keto$kills>=2] = 'kill>=2';
gp[gp=='x'] = 'kill01';
keto = cbind(keto, gp);
keto = keto[!(keto$timepoint=='stable'), ];

par(mfrow=c(2,2));
YLIM = c(0,300);
COL = 11;
FACT = 'kill>=2';
WGCNA::verboseBoxplot(cort[cort$gp==FACT,COL], cort$timepoint[cort$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(cort)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);
FACT = 'kill01';
WGCNA::verboseBoxplot(cort[cort$gp==FACT,COL], cort$timepoint[cort$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(cort)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = 'DURING_ATTACK';
WGCNA::verboseBoxplot(cort[cort$timepoint==FACT,COL], cort$gp[cort$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(cort)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = '10 days POST_ATTACK';
WGCNA::verboseBoxplot(cort[cort$timepoint==FACT,COL], cort$gp[cort$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(cort)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);			
					  
					  
					  
par(mfrow=c(2,2));
YLIM = c(0,500);
COL = 11;
names(estra)[COL] = 'Estradiol.ng.ml'
FACT = 'kill>=2';
WGCNA::verboseBoxplot(estra[estra$gp==FACT,COL], estra$timepoint[estra$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(estra)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);
FACT = 'kill01';
WGCNA::verboseBoxplot(estra[estra$gp==FACT,COL], estra$timepoint[estra$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(estra)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = 'DURING_ATTACK';
WGCNA::verboseBoxplot(estra[estra$timepoint==FACT,COL], estra$gp[estra$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(estra)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = '10 days POST_ATTACK';
WGCNA::verboseBoxplot(estra[estra$timepoint==FACT,COL], estra$gp[estra$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(estra)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);		  



par(mfrow=c(2,2));
YLIM = c(0,800);
COL = 11;
FACT = 'kill>=2';
WGCNA::verboseBoxplot(testo[testo$gp==FACT,COL], testo$timepoint[testo$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(testo)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);
FACT = 'kill01';
WGCNA::verboseBoxplot(testo[testo$gp==FACT,COL], testo$timepoint[testo$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(testo)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = 'DURING_ATTACK';
WGCNA::verboseBoxplot(testo[testo$timepoint==FACT,COL], testo$gp[testo$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(testo)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = '10 days POST_ATTACK';
WGCNA::verboseBoxplot(testo[testo$timepoint==FACT,COL], testo$gp[testo$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(testo)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);	
					  
par(mfrow=c(2,2));
YLIM = c(0,12);
COL = 11;
FACT = 'kill>=2';
WGCNA::verboseBoxplot(keto[keto$gp==FACT,COL], keto$timepoint[keto$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(keto)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);
FACT = 'kill01';
WGCNA::verboseBoxplot(keto[keto$gp==FACT,COL], keto$timepoint[keto$gp==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(keto)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = 'DURING_ATTACK';
WGCNA::verboseBoxplot(keto[keto$timepoint==FACT,COL], keto$gp[keto$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(keto)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);					  
FACT = '10 days POST_ATTACK';
WGCNA::verboseBoxplot(keto[keto$timepoint==FACT,COL], keto$gp[keto$timepoint==FACT],
					  col='grey', frame.plot=F, notch=F,
					  ylab=names(keto)[COL], xlab='',main=paste(FACT,'\n',sep=''), ylim=YLIM);	