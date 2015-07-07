rm(list=ls());
setwd('~/Documents/_Fernald_lab/_dyadStats');

#db = read.csv('RaiseFishAustin.csv');
#names(db) = toupper(c('date', 'run', 'created', 'tank', 'length.mm', 'mass.g', 'sex', 'density', 'age.days'));
#db0=db;
#db=db0[db0$RUN>12, ];
##OR##
db = read.csv('Growth_Rate_Data_11-14-13LB.csv');
names(db) = toupper(c('date', 'run', 'population', 'created', 'tank', 'length.mm', 'mass.g', 'sex', 'density', 'age.days'));
#db0=db;
#db=db[]


lm.len = lm(LENGTH.MM ~ AGE.DAYS + DENSITY, data=db);
lm.len.int = lm(LENGTH.MM ~ AGE.DAYS * DENSITY, data=db);

lm.len.est = lm.len[[1]][[1]] + lm.len[[1]][[2]]*db$AGE.DAYS + lm.len[[1]][[3]]*db$DENSITY;
lm.len.int.est = lm.len.int[[1]][[1]] + lm.len.int[[1]][[2]]*db$AGE.DAYS + lm.len.int[[1]][[3]]*db$DENSITY + lm.len.int[[1]][[4]]*db$DENSITY*db$AGE.DAYS;

par(mfrow=c(1,2));
WGCNA::verboseScatterplot(db$LENGTH.MM, lm.len.est, 
						  xlab='measured length', ylab='estimated length', main='lm.len:', 
						  abline=T, abline.col='red', abline.lty='dashed', 
						  xlim=c(0,50), ylim=c(0,50),
						  frame.plot=F);
WGCNA::verboseScatterplot(db$LENGTH.MM, lm.len.int.est, 
						  xlab='measured length', ylab='estimated length', main='lm.len.int:', 
						  abline=T, abline.col='red', abline.lty='dashed', 
						  xlim=c(0,50), ylim=c(0,50),
						  frame.plot=F);



par(mfrow=c(1,3), oma=c(0,1,2,0));

WGCNA::verboseScatterplot(db$LENGTH.MM, db$AGE.DAYS, 
						  xlab='Measured length (mm)', ylab='Age (days)', main='', cex.main=1.3,
						  abline=T, abline.col='red', abline.lty='dashed', 
						 xlim=c(20,60), ylim=c(40,140),
						  frame.plot=F, bg='grey', pch=21, col='black');
WGCNA::verboseScatterplot(db$LENGTH.MM, db$DENSITY, 
						  xlab='Measured length (mm)', ylab='Density (# of fish)', main='', cex.main=1.3,
						  abline=T, abline.col='red', abline.lty='dashed', 
						 xlim=c(20,60), ylim=c(10,50),
						  frame.plot=F, bg='grey', pch=21, col='black');
WGCNA::verboseScatterplot(db$LENGTH.MM, lm.len.est, 
						  xlab='Measured length (mm)', ylab='Estimated length (mm)', main='', cex.main=1.3,
						  abline=T, abline.col='red', abline.lty='dashed', 
						 xlim=c(20,60), ylim=c(20,50),
						  frame.plot=F, bg='grey', pch=21, col='black');
title('runs 20-21 & 24-33', outer=T)
						  
						  
						  
outk = kmeans(db[, names(db) %in% c('LENGTH.MM', 'AGE.DAYS', 'DENSITY')], 200);					  
outsurf = mba.surf(outk$centers, nrow(outk$centers), nrow(outk$centers));		
persp3d(outsurf$xyz.est, theta=135, phi=120, col="green3", scale=F,
  ltheta=-120, shade=1.5, expand=1, axes=T, box=F, xlab='LENGTH.MM', ylab='AGE.DAYS', zlab='DENSITY');			  
						  
###################
db.M = db[db$SEX=='M', ];

lm.len.M = lm(LENGTH.MM ~ AGE.DAYS + DENSITY, data=db.M);
lm.len.int.M = lm(LENGTH.MM ~ AGE.DAYS * DENSITY, data=db.M);

lm.len.est.M = lm.len.M[[1]][[1]] + lm.len.M[[1]][[2]]*db.M$AGE.DAYS + lm.len.M[[1]][[3]]*db.M$DENSITY;
lm.len.int.est.M = lm.len.int.M[[1]][[1]] + lm.len.int.M[[1]][[2]]*db.M$AGE.DAYS + lm.len.int.M[[1]][[3]]*db.M$DENSITY + lm.len.int.M[[1]][[4]]*db.M$DENSITY*db.M$AGE.DAYS;

par(mfrow=c(1,2));
WGCNA::verboseScatterplot(db.M$LENGTH.MM, lm.len.est.M, 
						  xlab='measured length', ylab='estimated length', main='lm.len.M:', 
						  abline=T, abline.col='red', abline.lty='dashed', 
						  xlim=c(0,50), ylim=c(0,50),
						  frame.plot=F);
WGCNA::verboseScatterplot(db.M$LENGTH.MM, lm.len.int.est.M, 
						  xlab='measured length', ylab='estimated length', main='lm.len.int.M:', 
						  abline=T, abline.col='red', abline.lty='dashed', 
						  xlim=c(0,50), ylim=c(0,50),
						  frame.plot=F);